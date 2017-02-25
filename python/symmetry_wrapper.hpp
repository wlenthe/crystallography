/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * Copyright (c) 2017, William C. Lenthe                                           *
 * All rights reserved.                                                            *
 *                                                                                 *
 * Redistribution and use in source and binary forms, with or without              *
 * modification, are permitted provided that the following conditions are met:     *
 *                                                                                 *
 * 1. Redistributions of source code must retain the above copyright notice, this  *
 *    list of conditions and the following disclaimer.                             *
 *                                                                                 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,    *
 *    this list of conditions and the following disclaimer in the documentation    *
 *    and/or other materials provided with the distribution.                       *
 *                                                                                 *
 * 3. Neither the name of the copyright holder nor the names of its                *
 *    contributors may be used to endorse or promote products derived from         *
 *    this software without specific prior written permission.                     *
 *                                                                                 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"     *
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE       *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE  *
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE    *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL      *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR      *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER      *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,   *
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE   *
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef _symmetry_wrapper_h_
#define _symmetry_wrapper_h_

#include "quaternion_wrapper.hpp"
#include <string>
#include <algorithm>
#include <cctype>

#include "symmetry.hpp"

typedef struct {
	PyObject_HEAD
	Symmetry<double> const * symPtr;
} PySymmetry;

static PySymmetry* Symmetry_create(Symmetry<double> const * sym);
static bool PySymmetry_Check(PyObject* obj);

static void Symmetry_dealloc(PySymmetry* self){
	Py_TYPE(self)->tp_free((PyObject*)self);
}

static int Symmetry_init(PySymmetry* self, PyObject* args, PyObject* kwds) {
	//parse arguments
	char* buff = NULL;
	const char *kwlist[] = {"name", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s", const_cast<char**>(kwlist), &buff)) return -1;
	self->symPtr = Symmetry<double>::FromName(buff);
	if(self->symPtr == NULL) {
		std::string err("unknown symmetry name '");
		err += buff;
		err += "'";
		PyErr_SetString(PyExc_ValueError, err.c_str());
		return -1;
	}
	return 0;
}

static PyObject* Symmetry_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
	PySymmetry* self = (PySymmetry*)type->tp_alloc(type, 0);
	Symmetry_init(self, args, kwds);
	return (PyObject*)self;
}

static PyObject* Symmetry_quOperators(PySymmetry* self) {
	npy_intp dims[1];
	std::vector<Quaternion<double> > const * const opPtr = self->symPtr->quOperators();
	dims[0] = opPtr->size();
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_OBJECT, 0);
	PyQuaternion** buff = (PyQuaternion**)PyArray_DATA(output);
	for(size_t i = 0; i < opPtr->size(); i++)	buff[i] = Quaternion_create(opPtr->at(i));
	return (PyObject*) output;
}

static PyObject* Symmetry_fzQu(PySymmetry* self, PyObject* array) {
	//get array object and its dimensions
	PyArrayObject* input = (PyArrayObject*)PyArray_FROM_OTF(array, NPY_QUAT, NPY_ARRAY_IN_ARRAY);
	if(input == NULL) {
		Py_XDECREF(input);
		return NULL;
	}
	int ndims = PyArray_NDIM(input);
	npy_intp* dims = PyArray_DIMS(input);
	Symmetry<double> const * sym = ((PySymmetry*)self)->symPtr;
	Quaternion<double>* from = (Quaternion<double>*)PyArray_DATA(input);

	//get number of points and compute fz quaternions
	int totalPoints = 1;
	for(npy_intp i = 0; i < ndims; i++) totalPoints *= (int)dims[i];
	if(1 == totalPoints && !PyArray_Check(array)) {
		return (PyObject*)Quaternion_create(sym->fzQu(from[0]));
	} else {
		PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(ndims, dims, NPY_QUAT, 0);
		Quaternion<double>* to = (Quaternion<double>*)PyArray_DATA(output);
		for(int i = 0; i < totalPoints; i++) to[i] = sym->fzQu(from[i]);
		return (PyObject*)output;
	}
}

//wrapper for members taking 2 quaternions and returning a quaternion
static PyObject* Symmetry_qu_wrapper(PySymmetry* self, PyObject* args, PyObject* kwds, Quaternion<double>(Symmetry<double>::*f)(const Quaternion<double>&,const Quaternion<double>&)const) {
	PyObject* q1;
	PyObject* q2;
	const char *kwlist[] = {"qu1", "qu2", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "OO", const_cast<char**>(kwlist), &q1, &q2)) return NULL;

	//attempt to convert both arguments to numpy arrays
	PyArrayObject* input1 = (PyArrayObject*)PyArray_FROM_OTF(q1, NPY_QUAT, NPY_ARRAY_IN_ARRAY);
	PyArrayObject* input2 = (PyArrayObject*)PyArray_FROM_OTF(q2, NPY_QUAT, NPY_ARRAY_IN_ARRAY);
	if(input1 == NULL) {
		PyErr_SetString(PyExc_ValueError, "this function requires 2 Quaternion arguments");
		Py_XDECREF(input1);
		return NULL;
	}
	if(input2 == NULL) {
		PyErr_SetString(PyExc_ValueError, "this function requires 2 Quaternion arguments");
		Py_XDECREF(input2);
		return NULL;
	}
	int ndims1 = PyArray_NDIM(input1);
	int ndims2 = PyArray_NDIM(input2);
	npy_intp* dims1 = PyArray_DIMS(input1);
	npy_intp* dims2 = PyArray_DIMS(input2);
	Symmetry<double> const * sym = ((PySymmetry*)self)->symPtr;
	Quaternion<double>* from1 = (Quaternion<double>*)PyArray_DATA(input1);
	Quaternion<double>* from2 = (Quaternion<double>*)PyArray_DATA(input2);

	//get number of points
	int totalPoints1 = 1;
	int totalPoints2 = 1;
	for(npy_intp i = 0; i < ndims1; i++) totalPoints1 *= (int)dims1[i];
	for(npy_intp i = 0; i < ndims2; i++) totalPoints2 *= (int)dims2[i];
	bool single1 = (1 == totalPoints1 && !PyArray_Check(q1));
	bool single2 = (1 == totalPoints2 && !PyArray_Check(q2));

	if(single1 && single2) {
		//2 single quaternions
		return (PyObject*)Quaternion_create((sym->*f)(from1[0], from2[0]));
	} else if(single1 && !single2) {
		//single quaternion and array
		PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(ndims2, dims2, NPY_QUAT, 0);
		Quaternion<double>* to = (Quaternion<double>*)PyArray_DATA(output);
		for(int i = 0; i < totalPoints2; i++)
			to[i] = (sym->*f)(from1[0], from2[i]);
		return (PyObject*)output;
 	} else if(!single1 && single2) {
		//array and single quaternion
		PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(ndims1, dims1, NPY_QUAT, 0);
		Quaternion<double>* to = (Quaternion<double>*)PyArray_DATA(output);
		for(int i = 0; i < totalPoints1; i++)
			to[i] = (sym->*f)(from1[i], from2[0]);
		return (PyObject*)output;
	} else {
		//2 arrays of quaternions
		if(ndims1 != ndims2) {
			PyErr_SetString(PyExc_ValueError, "inputs must have same dimensions");
			return NULL;
		}
		for(int i = 0; i < ndims1; i++) {
			if(dims1[i] != dims2[i]) {
				PyErr_SetString(PyExc_ValueError, "inputs must have same shape");
				return NULL;
			}
		}
		PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(ndims1, dims1, NPY_QUAT, 0);
		Quaternion<double>* to = (Quaternion<double>*)PyArray_DATA(output);
		for(int i = 0; i < totalPoints1; i++)
			to[i] = (sym->*f)(from1[i], from2[i]);
		return (PyObject*)output;
	}
}

static PyObject* Symmetry_Repr(PyObject* object) {
	PySymmetry* self = (PySymmetry*)object;
	std::string repr = self->symPtr->name() + " symmetry";
	return Py_BuildValue("s", repr.c_str());
}

static PyObject* Symmetry_disoQu  (PySymmetry* self, PyObject* args, PyObject* kwds) {return Symmetry_qu_wrapper(self, args, kwds, &Symmetry<double>::disoQu);}
static PyObject* Symmetry_nearbyQu(PySymmetry* self, PyObject* args, PyObject* kwds) {return Symmetry_qu_wrapper(self, args, kwds, &Symmetry<double>::nearbyQu);}

static PyObject* Symmetry_cubic        (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::Cubic        );}
static PyObject* Symmetry_cubiclow     (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::CubicLow     );}
static PyObject* Symmetry_hexagonal    (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::Hexagonal    );}
static PyObject* Symmetry_hexagonallow (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::HexagonalLow );}
static PyObject* Symmetry_tetragonal   (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::Tetragonal   );}
static PyObject* Symmetry_tetragonallow(PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::TetragonalLow);}
static PyObject* Symmetry_trigonal     (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::Trigonal     );}
static PyObject* Symmetry_trigonallow  (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::TrigonalLow  );}
static PyObject* Symmetry_orthorhombic (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::Orthorhombic );}
static PyObject* Symmetry_monoclinic   (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::Monoclinic   );}
static PyObject* Symmetry_triclinic    (PySymmetry* self) {return (PyObject*)Symmetry_create(&symmetry::Groups<double>::Triclinic    );}

//method table
static PyMethodDef Symmetry_methods[] = {
	{"quOperators"  , (PyCFunction)Symmetry_quOperators, METH_NOARGS                 , "sym.quOperators()\n\treturn sym's symmetry operators as Quaternion objects"                               },
	{"fzQu"         , (PyCFunction)Symmetry_fzQu       , METH_O                      , "sym.fzQu(qu)\n\tcompute the equivalent representation of qu in the fundamental zone"                      },
	{"disoQu"       , (PyCFunction)Symmetry_disoQu     , METH_VARARGS | METH_KEYWORDS, "sym.disoQu(qu1, qu2)\n\tcompute the disorientation between qu1 and qu2"                                   },
	{"nearbyQu"     , (PyCFunction)Symmetry_nearbyQu   , METH_VARARGS | METH_KEYWORDS, "sym.nearbyQu(qu1, qu2)\n\tcompute the equivalent representation of qu2 closest to qu1"                    },
	{NULL}  /* Sentinel */
};
                                 
//type definition
static PyTypeObject SymmetryType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"Symmetry",                    //tp_name
	sizeof(PySymmetry),            //tp_basicsize
	0,                             //tp_itemsize
	(destructor)Symmetry_dealloc,  //tp_dealloc
	0,                             //tp_print
	0,                             //tp_getattr
	0,                             //tp_setattr
	0,                             //tp_reserved
	Symmetry_Repr,                 //tp_repr
	0,                             //tp_as_number
	0,                             //tp_as_sequence
	0,                             //tp_as_mapping
	0,                             //tp_hash 
	0,                             //tp_call
	0,                             //tp_str
	0,                             //tp_getattro
	0,                             //tp_setattro
	0,                             //tp_as_buffer
	Py_TPFLAGS_DEFAULT,            //tp_flags
	//tp_doc
	"Symmetry(name)\n\
\t-name must be one of:\n\
\t\tcubic\n\
\t\tcubic low\n\
\t\thexagonal\n\
\t\thexagonal low\n\
\t\ttetragonal\n\
\t\ttetragonal low\n\
\t\ttrigonal\n\
\t\ttrigonal low\n\
\t\torthorhombic\n\
\t\tmonoclinic\n\
\t\ttriclinic\n",
	0,                             //tp_traverse
	0,                             //tp_clear
	0,                             //tp_richcompare
	0,                             //tp_weaklistoffset
	0,                             //tp_iter
	0,                             //tp_iternext
	Symmetry_methods,              //tp_methods
	0,                             //tp_members
	0,                             //tp_getset
	0,                             //tp_base
	0,                             //tp_dict
	0,                             //tp_descr_get
	0,                             //tp_descr_set
	0,                             //tp_dictoffset
	(initproc)Symmetry_init,       //tp_init
	0,                             //tp_alloc
	Symmetry_new,                  //tp_new
};
static PySymmetry* Symmetry_create(Symmetry<double> const * sym) {
	PySymmetry* obj = (PySymmetry*)SymmetryType.tp_alloc(&SymmetryType, 0);
	obj->symPtr = sym;
	return obj;
}
static bool PySymmetry_Check(PyObject* obj) {return (PyTypeObject*)(obj->ob_type) == &SymmetryType;}

#endif//_symmetry_wrapper_h_