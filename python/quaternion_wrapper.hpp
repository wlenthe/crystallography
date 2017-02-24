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
#ifndef _quaternion_wrapper_h_
#define _quaternion_wrapper_h_

//python says to include its headers first
#include "Python.h"
#include "structmember.h"
#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#include "numpy/arrayobject.h"
#include "numpy/npy_math.h"
#include "numpy/ufuncobject.h"

#include <sstream>
#include <iomanip>
#include <vector>

#include "quaternion.hpp"

typedef struct {
	PyObject_HEAD
	Quaternion<double> qu;
} PyQuaternion;

PyArray_Descr *quaternion_descr;//numpy type descriptor
PyQuaternion* Quaternion_create(Quaternion<double> qu);//create a PyQuaternion from a Quaternion<double>
static bool PyQuaternion_Check(PyObject* obj);//check if a PyObject is a PyQuaternion

static PyObject* Quaternion_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
	Quaternion<double> qu = Quaternion<double>::Identity();
	if(!PyArg_ParseTuple(args, "|dddd", &qu.w, &qu.x, &qu.y, &qu.z)) return NULL;
	return (PyObject*) Quaternion_create(qu);
}

//wrapper for members taking no arguments and returning another quaternion
static PyObject* Quaternion_member_wrapper(PyQuaternion* self, Quaternion<double>(Quaternion<double>::*f)()const) {
	return (PyObject*)Quaternion_create((self->qu.*f)());
}

//wrapper for members taking a scalar argument
static PyObject* Quaternion_scalar_wrapper(PyObject* self, PyObject* other, Quaternion<double>(Quaternion<double>::*f)(const double&)const) {
	if(!PyFloat_Check(other)) {
		PyErr_SetString(PyExc_ValueError, "this function requires a float argument");
		return NULL;
	}
	return (PyObject*)Quaternion_create( (((PyQuaternion*)self)->qu.*f)(PyFloat_AS_DOUBLE(other)) );
}

//wrapper for members taking a quaternion argument
static PyObject* Quaternion_quaternion_wrapper(PyObject* self, PyObject* other, double(Quaternion<double>::*f)(const Quaternion<double>&)const) {
	if(!PyQuaternion_Check(other)) {
		PyErr_SetString(PyExc_ValueError, "this function requires a Quaternion argument");
		return NULL;
	}
	return Py_BuildValue("d", (((PyQuaternion*)self)->qu.*f)(((PyQuaternion*)other)->qu));
}

//wrapper for binary operators
static PyObject* Quaternion_binary_wrapper(PyObject* self, PyObject* other, Quaternion<double>(Quaternion<double>::*f)(const Quaternion<double>&)const) {
	if(!PyQuaternion_Check(self) || !PyQuaternion_Check(other)) {
		PyErr_SetString(PyExc_ValueError, "Binary operations are defined for 2 quaternion objects only. Use the scalar* functions for binary operations with a scalar value.");
		return NULL;
	}
	return (PyObject*)Quaternion_create( (((PyQuaternion*)self)->qu.*f)(((PyQuaternion*)other)->qu) );
}

static PyObject* Quaternion_normalize(PyQuaternion* self) {return Quaternion_member_wrapper(self, &Quaternion<double>::normalize);}
static PyObject* Quaternion_conjugate(PyQuaternion* self) {return Quaternion_member_wrapper(self, &Quaternion<double>::conjugate);}
static PyObject* Quaternion_inverse(PyQuaternion* self) {return Quaternion_member_wrapper(self, &Quaternion<double>::inverse);}
static PyObject* Quaternion_explement(PyQuaternion* self) {return Quaternion_member_wrapper(self, &Quaternion<double>::explement);}
static PyObject* Quaternion_negate(PyQuaternion* self) {return Quaternion_member_wrapper(self, &Quaternion<double>::operator-);}
static PyObject* Quaternion_abs(PyQuaternion* self) {return Quaternion_member_wrapper(self, &Quaternion<double>::abs);}

static PyObject* Quaternion_scalar_add(PyObject* self, PyObject* other) {return Quaternion_scalar_wrapper(self, other, &Quaternion<double>::operator+);}
static PyObject* Quaternion_scalar_sub(PyObject* self, PyObject* other) {return Quaternion_scalar_wrapper(self, other, &Quaternion<double>::operator-);}
static PyObject* Quaternion_scalar_mul(PyObject* self, PyObject* other) {return Quaternion_scalar_wrapper(self, other, &Quaternion<double>::operator*);}
static PyObject* Quaternion_scalar_div(PyObject* self, PyObject* other) {return Quaternion_scalar_wrapper(self, other, &Quaternion<double>::operator/);}

static PyObject* Quaternion_dot(PyObject* self, PyObject* other) {return Quaternion_quaternion_wrapper(self, other, &Quaternion<double>::dot);}
static PyObject* Quaternion_mag(PyQuaternion* self) {return Py_BuildValue("d", ((PyQuaternion*)self)->qu.mag());}

static PyObject* Quaternion_binary_add(PyObject* self, PyObject* other) {return Quaternion_binary_wrapper(self, other, &Quaternion<double>::operator+);}
static PyObject* Quaternion_binary_sub(PyObject* self, PyObject* other) {return Quaternion_binary_wrapper(self, other, &Quaternion<double>::operator-);}
static PyObject* Quaternion_binary_mul(PyObject* self, PyObject* other) {return Quaternion_binary_wrapper(self, other, &Quaternion<double>::operator*);}
static PyObject* Quaternion_binary_div(PyObject* self, PyObject* other) {return Quaternion_binary_wrapper(self, other, &Quaternion<double>::operator/);}
static PyObject* Quaternion_urnary_negate(PyObject* self) {return Quaternion_negate((PyQuaternion*) self);}
static PyObject* Quaternion_urnary_abs(PyObject* self) {return Quaternion_abs((PyQuaternion*) self);}

static PyObject* Quaternion_Repr(PyObject* object) {
	PyQuaternion* self = (PyQuaternion*)object;
	std::stringstream ss;
	ss << "(" << std::fixed << std::setprecision(8);
	ss << (self->qu.w >= 0.0 ? " " : "") << self->qu.w << ", ";
	ss << (self->qu.x >= 0.0 ? " " : "") << self->qu.x << ", ";
	ss << (self->qu.y >= 0.0 ? " " : "") << self->qu.y << ", ";
	ss << (self->qu.z >= 0.0 ? " " : "") << self->qu.z << ")";
	return Py_BuildValue("s", ss.str().c_str());//"y" for bytes object
}

static PyObject* Quaternion_tuple(PyObject* self) {
	PyQuaternion* o = (PyQuaternion*)self;
	return PyTuple_Pack(4,
		PyFloat_FromDouble(o->qu.w),
		PyFloat_FromDouble(o->qu.x),
		PyFloat_FromDouble(o->qu.y),
		PyFloat_FromDouble(o->qu.z));
}

//method table
static PyMethodDef Quaternion_methods[] = {
	{"normalize", (PyCFunction)Quaternion_normalize , METH_NOARGS, "qu.normalize()\n\treturn a normalized copy qu"},
	{"conjugate", (PyCFunction)Quaternion_conjugate , METH_NOARGS, "qu.conjugate()\n\treturn a conjugated copy qu"},
	{"inverse"  , (PyCFunction)Quaternion_inverse   , METH_NOARGS, "qu.inverse()\n\treturn an inversed copy qu"},
	{"explement", (PyCFunction)Quaternion_explement , METH_NOARGS, "qu.explement()\n\treturn the smaller qu and its explement"},
	{"magnitude", (PyCFunction)Quaternion_mag       , METH_NOARGS, "qu.magnitude()\n\treturn magnitude of qu"},
	{"dot"      , (PyCFunction)Quaternion_dot       , METH_O     , "qu1.dot(qu2)\n\treturns dot product of qu1 and qu2"},
	{"scalarAdd", (PyCFunction)Quaternion_scalar_add, METH_O     , "qu.scalarAdd(s)\n\treturn a copy of qu s added to each element"},
	{"scalarSub", (PyCFunction)Quaternion_scalar_sub, METH_O     , "qu.scalarAdd(s)\n\treturn a copy of qu s subracted from each element"},
	{"scalarMul", (PyCFunction)Quaternion_scalar_mul, METH_O     , "qu.scalarAdd(s)\n\treturn a copy of qu with each element multiplied by s"},
	{"scalarDiv", (PyCFunction)Quaternion_scalar_div, METH_O     , "qu.scalarAdd(s)\n\treturn a copy of qu with each element divided by s"},
	{"toTuple"  , (PyCFunction)Quaternion_tuple     , METH_NOARGS, "qu.toTuple()\n\treturn (w,x,y,z)"},
	{NULL}//sentinel
};

//member table
PyMemberDef Quaternion_members[] = {
	{(char*)"w", T_DOUBLE, offsetof(PyQuaternion, qu) + 0*sizeof(double), 0, (char*)"cos(angle/2)"},
	{(char*)"x", T_DOUBLE, offsetof(PyQuaternion, qu) + 1*sizeof(double), 0, (char*)"nx*sin(angle/2)"},
	{(char*)"y", T_DOUBLE, offsetof(PyQuaternion, qu) + 2*sizeof(double), 0, (char*)"nx*sin(angle/2)"},
	{(char*)"z", T_DOUBLE, offsetof(PyQuaternion, qu) + 3*sizeof(double), 0, (char*)"nx*sin(angle/2)"},
	{NULL}//sentinel
};

//number method table
static PyNumberMethods QuaternionNumberMethods {
	Quaternion_binary_add,   //nb_add
	Quaternion_binary_sub,   //nb_subtract
	Quaternion_binary_mul,   //nb_multiply
	0,                       //nb_remainder
	0,                       //nb_divmod
	0,                       //nb_power
	Quaternion_urnary_negate,//nb_negative
	0,                       //nb_positive
	Quaternion_urnary_abs,   //nb_absolute
	0,                       //nb_bool
	0,                       //nb_invert
	0,                       //nb_lshift
	0,                       //nb_rshift
	0,                       //nb_and
	0,                       //nb_xor
	0,                       //nb_or
	0,                       //nb_int
	0,                       //nb_reserved
	0,                       //nb_float
	0,                       //nb_inplace_add
	0,                       //nb_inplace_subtract
	0,                       //nb_inplace_multiply
	0,                       //nb_inplace_remainder
	0,                       //nb_inplace_power
	0,                       //nb_inplace_lshift
	0,                       //nb_inplace_rshift
	0,                       //nb_inplace_and
	0,                       //nb_inplace_xor
	0,                       //nb_inplace_or
	0,                       //nb_floor_divide
	Quaternion_binary_div,   //nb_true_divide
	0,                       //nb_inplace_floor_divide
	0,                       //nb_inplace_true_divide
	0,                       //nb_index
	0,                       //nb_matrix_multiply
	0,                       //nb_inplace_matrix_multiply
};

PyTypeObject QuaternionType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"Quaternion",                 //tp_name
	sizeof(PyQuaternion),         //tp_basicsize
	0,                            //tp_itemsize
	0,                            //tp_dealloc
	0,                            //tp_print
	0,                            //tp_getattr
	0,                            //tp_setattr
	0,                            //tp_reserved
	Quaternion_Repr,              //tp_repr
	&QuaternionNumberMethods,     //tp_as_number
	0,                            //tp_as_sequence
	0,                            //tp_as_mapping
	0,                            //tp_hash 
	0,                            //tp_call
	Quaternion_Repr,                            //tp_str
	0,                            //tp_getattro
	0,                            //tp_setattro
	0,                            //tp_as_buffer
	Py_TPFLAGS_DEFAULT,           //tp_flags
	"Quaternion(w=1,x=0,y=0,z=0)",//tp_doc: this seems to cause issues runtime with numpy if len the string is shorter than 25 or longer than 32 bytes including null. extremely concerning...
	0,                            //tp_traverse
	0,                            //tp_clear
	0,                            //tp_richcompare
	0,                            //tp_weaklistoffset
	0,                            //tp_iter
	0,                            //tp_iternext
	Quaternion_methods,           //tp_methods
	Quaternion_members,           //tp_members
	0,                            //tp_getset
	0,                            //tp_base
	0,                            //tp_dict
	0,                            //tp_descr_get
	0,                            //tp_descr_set
	0,                            //tp_dictoffset
	0,                            //tp_init
	0,                            //tp_alloc
	Quaternion_new,               //tp_new
};
static bool PyQuaternion_Check(PyObject* obj) {return (PyTypeObject*)(obj->ob_type) == &QuaternionType;}

PyQuaternion* Quaternion_create(Quaternion<double> qu) {
	PyQuaternion* self = (PyQuaternion*)QuaternionType.tp_alloc(&QuaternionType, 0);
	self->qu = qu;
	return self;
}

//////////////////////////////////////////////////////////////////////
//                 numpy wrapping code begins here                  //
//////////////////////////////////////////////////////////////////////
static int NPY_QUAT;
static PyArray_ArrFuncs quaternionArrFuncs;

//numpy requires the following function to be defined for a new type: nonzero, copyswap, copyswapn, setitem, getitem, and cast
static bool Quaternion_nonzero(void* data, void* arr) {
	Quaternion<double> qu;
	if(PyArray_ISBEHAVED_RO((PyArrayObject*)arr)) {
		qu = *(Quaternion<double>*)data;
	} else {
		//must handle poorly behaved arrays
		char* cData = (char*) data;
		PyArray_Descr* descr = PyArray_DescrFromType(NPY_DOUBLE);
		descr->f->copyswap(&qu.w, cData + offsetof(Quaternion<double>, w), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(&qu.x, cData + offsetof(Quaternion<double>, x), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(&qu.y, cData + offsetof(Quaternion<double>, y), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(&qu.z, cData + offsetof(Quaternion<double>, z), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		Py_DECREF(descr);
	}
		 if(qu.w != 0.0) return true;
	else if(qu.x != 0.0) return true;
	else if(qu.y != 0.0) return true;
	else if(qu.z != 0.0) return true;
	else                return false;
}

static void Quaternion_copyswap(void* dest, void* src, int swap, void* arr) {
	PyArray_Descr* descr = PyArray_DescrFromType(NPY_DOUBLE);
	descr->f->copyswapn(dest, sizeof(double), src, sizeof(double), 4, swap, NULL);
	Py_DECREF(descr);
}

static void Quaternion_copyswapn(void* dest, npy_intp dstride, void* src, npy_intp sstride, npy_intp n, int swap, void* arr) {
	PyArray_Descr* descr = PyArray_DescrFromType(NPY_DOUBLE);
	descr->f->copyswapn(dest, sizeof(double), src, sizeof(double), 4, swap, NULL);
	descr->f->copyswapn((void*) &(((Quaternion<double>*)dest)->w), dstride, (void*) &(((Quaternion<double>*)src)->w), sstride, n, swap, NULL);
	descr->f->copyswapn((void*) &(((Quaternion<double>*)dest)->x), dstride, (void*) &(((Quaternion<double>*)src)->x), sstride, n, swap, NULL);
	descr->f->copyswapn((void*) &(((Quaternion<double>*)dest)->y), dstride, (void*) &(((Quaternion<double>*)src)->y), sstride, n, swap, NULL);
	descr->f->copyswapn((void*) &(((Quaternion<double>*)dest)->z), dstride, (void*) &(((Quaternion<double>*)src)->z), sstride, n, swap, NULL);
	Py_DECREF(descr);
}

static int Quaternion_setitem(PyObject* item, void* data, void* arr) {
	if(!PyQuaternion_Check(item)) {
		PyErr_SetString(PyExc_ValueError, "cannot set value with non Quaternion object");
		return 1;
	}
	Quaternion<double> qu = ((PyQuaternion*)item)->qu;
	if(arr == NULL || PyArray_ISBEHAVED((PyArrayObject*)arr)) {
		*((Quaternion<double>*)data) = qu;
	} else {
		//must handle poorly behaved arrays
		char* cData = (char*) data;
		PyArray_Descr* descr = PyArray_DescrFromType(NPY_DOUBLE);
		descr->f->copyswap(cData + offsetof(Quaternion<double>, w), (void*) &qu.w, !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(cData + offsetof(Quaternion<double>, x), (void*) &qu.x, !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(cData + offsetof(Quaternion<double>, y), (void*) &qu.y, !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(cData + offsetof(Quaternion<double>, z), (void*) &qu.z, !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		Py_DECREF(descr);
	}
	return PyErr_Occurred() ? -1 : 0;
}

static PyObject* Quaternion_getitem(void* data, void* arr) {
	Quaternion<double> qu;
	if((arr == NULL) || PyArray_ISBEHAVED_RO((PyArrayObject*)arr)) {
		qu = *((Quaternion<double>*)data);
	} else {
		//must handle poorly behaved arrays
		char* cData = (char*) data;
		PyArray_Descr* descr = PyArray_DescrFromType(NPY_DOUBLE);
		descr->f->copyswap(&qu.w, cData + offsetof(Quaternion<double>, w), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(&qu.x, cData + offsetof(Quaternion<double>, x), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(&qu.y, cData + offsetof(Quaternion<double>, y), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		descr->f->copyswap(&qu.z, cData + offsetof(Quaternion<double>, z), !PyArray_ISNOTSWAPPED((PyArrayObject*)arr), NULL);
		Py_DECREF(descr);
	}
	return (PyObject*) Quaternion_create(qu);
}

//helpers to cast from scalars
template <typename T>
void Quaternion_cast_wrapper(void* from, void* to, npy_intp n, void* fromarr, void* toarr) {
	static const double fillVal = 0.0;
	double* const destData = (double*) from;
	std::fill(destData, destData + 4 * n, fillVal);//scalar to quaternion doesn't make sense, fill with 0
}
template <typename T> void Quaternion_complex_cast_wrapper(void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper<T>(from, to, n, fromarr, toarr);}

//functions to cast from scalars
void Quaternion_from_BOOL       (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_bool      >(from, to, n, fromarr, toarr);}
void Quaternion_from_BYTE       (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_byte      >(from, to, n, fromarr, toarr);}
void Quaternion_from_UBYTE      (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_ubyte     >(from, to, n, fromarr, toarr);}
void Quaternion_from_SHORT      (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_short     >(from, to, n, fromarr, toarr);}
void Quaternion_from_USHORT     (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_ushort    >(from, to, n, fromarr, toarr);}
void Quaternion_from_INT        (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_int       >(from, to, n, fromarr, toarr);}
void Quaternion_from_UINT       (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_uint      >(from, to, n, fromarr, toarr);}
void Quaternion_from_LONG       (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_long      >(from, to, n, fromarr, toarr);}
void Quaternion_from_LONGLONG   (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_longlong  >(from, to, n, fromarr, toarr);}
void Quaternion_from_ULONGLONG  (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_ulonglong >(from, to, n, fromarr, toarr);}
void Quaternion_from_ULONG      (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_ulong     >(from, to, n, fromarr, toarr);}
void Quaternion_from_FLOAT      (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_uint32    >(from, to, n, fromarr, toarr);}
void Quaternion_from_DOUBLE     (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_uint64    >(from, to, n, fromarr, toarr);}
void Quaternion_from_LONGDOUBLE (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_cast_wrapper        <npy_longdouble>(from, to, n, fromarr, toarr);}
void Quaternion_from_CFLOAT     (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_complex_cast_wrapper<npy_uint32    >(from, to, n, fromarr, toarr);}
void Quaternion_from_CDOUBLE    (void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_complex_cast_wrapper<npy_uint64    >(from, to, n, fromarr, toarr);}
void Quaternion_from_CLONGDOUBLE(void* from, void* to, npy_intp n, void* fromarr, void* toarr) {Quaternion_complex_cast_wrapper<npy_longdouble>(from, to, n, fromarr, toarr);}

//ufunc wrappers (numpy versions of Quaternion_member_wrapper, Quaternion_scalar_wrapper, Quaternion_quaternion_wrapper, Quaternion_binary_wrapper)

//wrapper for members taking no arguments and returning another quaternion
static void Quaternion_member_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata, Quaternion<double>(Quaternion<double>::*f)()const) {
 char* const source = args[0];
 char* const dest = args[1];
 const npy_intp sourceStride = strides[0];
 const npy_intp destStride = strides[1];
 const npy_intp n = dimensions[0];
 for(npy_intp i = 0; i < n; i++)
	*((Quaternion<double>*)(dest+i*destStride)) = (((Quaternion<double>*)(source+sourceStride))->*f)();
}
static void Quaternion_conj_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_member_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::conjugate);}
static void Quaternion_inv_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_member_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::inverse);}
static void Quaternion_neg_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_member_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator-);}

//wrapper for binary operators
static void Quaternion_binary_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata, Quaternion<double>(Quaternion<double>::*f)(const Quaternion<double>&)const) {
 char* const source1 = args[0];
 char* const source2 = args[1];
 char* const dest = args[2];
 const npy_intp source1Stride = strides[0];
 const npy_intp source2Stride = strides[1];
 const npy_intp destStride = strides[2];
 const npy_intp n = dimensions[0];
 for(npy_intp i = 0; i < n; i++)
	*((Quaternion<double>*)(dest+i*destStride)) = (((Quaternion<double>*)(source1+source1Stride))->*f)(*((Quaternion<double>*)(source2+source2Stride)));
}
static void Quaternion_mul_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_binary_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator*);}

//wrapper for binary comparison operators
static void Quaternion_cmp_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata, bool(Quaternion<double>::*f)(const Quaternion<double>&)const) {
 char* const source1 = args[0];
 char* const source2 = args[1];
 char* const dest = args[2];
 const npy_intp source1Stride = strides[0];
 const npy_intp source2Stride = strides[1];
 const npy_intp destStride = strides[2];
 const npy_intp n = dimensions[0];
 for(npy_intp i = 0; i < n; i++)
	*((bool*)(dest+i*destStride)) = (((Quaternion<double>*)(source1+source1Stride))->*f)(*((Quaternion<double>*)(source2+source2Stride)));
}
static void Quaternion_eq_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_cmp_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator==);}
static void Quaternion_ne_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_cmp_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator!=);}
static void Quaternion_lt_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_cmp_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator<);}
static void Quaternion_le_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_cmp_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator<=);}
static void Quaternion_gt_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_cmp_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator>);}
static void Quaternion_ge_ufunc(char** args, npy_intp* dimensions, npy_intp* strides, void* innerloopdata) {return Quaternion_cmp_ufunc(args, dimensions, strides, innerloopdata, &Quaternion<double>::operator>=);}

static int Quaternion_registerNumpy() {
	//assign required array functions
	PyArray_InitArrFuncs(&quaternionArrFuncs);
	quaternionArrFuncs.getitem   = (PyArray_GetItemFunc*)  Quaternion_getitem;
	quaternionArrFuncs.setitem   = (PyArray_SetItemFunc*)  Quaternion_setitem;
	quaternionArrFuncs.copyswap  = (PyArray_CopySwapFunc*) Quaternion_copyswap;
	quaternionArrFuncs.copyswapn = (PyArray_CopySwapNFunc*)Quaternion_copyswapn;
	quaternionArrFuncs.nonzero   = (PyArray_NonzeroFunc*)  Quaternion_nonzero;

	//build array descriptor
	quaternion_descr = PyObject_New(PyArray_Descr, &PyArrayDescr_Type);
	quaternion_descr->typeobj   = &QuaternionType;
	quaternion_descr->kind      = 'q';
	quaternion_descr->type      = 'j';
	quaternion_descr->byteorder = '=';
	quaternion_descr->type_num  = 0;
	quaternion_descr->elsize    = 4*sizeof(double);
	quaternion_descr->alignment = sizeof(double);
	quaternion_descr->subarray  = NULL;
	quaternion_descr->fields    = NULL;
	quaternion_descr->names     = NULL;
	quaternion_descr->f         = &quaternionArrFuncs;

	//register data type with numpy and get type number
	NPY_QUAT = PyArray_RegisterDataType(quaternion_descr);
	if(NPY_QUAT < 0) return 1;
	
	//initialize cast array functions
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_BOOL       ), NPY_QUAT, Quaternion_from_BOOL       );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_BYTE       ), NPY_QUAT, Quaternion_from_BYTE       );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_UBYTE      ), NPY_QUAT, Quaternion_from_UBYTE      );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_SHORT      ), NPY_QUAT, Quaternion_from_SHORT      );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_USHORT     ), NPY_QUAT, Quaternion_from_USHORT     );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_INT        ), NPY_QUAT, Quaternion_from_INT        );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_UINT       ), NPY_QUAT, Quaternion_from_UINT       );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_LONG       ), NPY_QUAT, Quaternion_from_LONG       );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_ULONG      ), NPY_QUAT, Quaternion_from_LONGLONG   );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_LONGLONG   ), NPY_QUAT, Quaternion_from_ULONGLONG  );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_ULONGLONG  ), NPY_QUAT, Quaternion_from_ULONG      );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_FLOAT      ), NPY_QUAT, Quaternion_from_FLOAT      );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_DOUBLE     ), NPY_QUAT, Quaternion_from_DOUBLE     );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_LONGDOUBLE ), NPY_QUAT, Quaternion_from_LONGDOUBLE );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_CFLOAT     ), NPY_QUAT, Quaternion_from_CFLOAT     );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_CDOUBLE    ), NPY_QUAT, Quaternion_from_CDOUBLE    );
	PyArray_RegisterCastFunc(PyArray_DescrFromType(NPY_CLONGDOUBLE), NPY_QUAT, Quaternion_from_CLONGDOUBLE);

	//disable casting
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_BOOL       ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_BYTE       ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_UBYTE      ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_SHORT      ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_USHORT     ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_INT        ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_UINT       ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_LONG       ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_ULONG      ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_LONGLONG   ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_ULONGLONG  ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_FLOAT      ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_DOUBLE     ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_LONGDOUBLE ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_CFLOAT     ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_CDOUBLE    ), NPY_QUAT, NPY_NOSCALAR);
	PyArray_RegisterCanCast(PyArray_DescrFromType(NPY_CLONGDOUBLE), NPY_QUAT, NPY_NOSCALAR);

	//register ufuncs
	int urnaryArgs[3];
	urnaryArgs[0] = NPY_QUAT;
	urnaryArgs[1] = NPY_QUAT;
	PyObject* ufuncDict = PyModule_GetDict(PyImport_ImportModule("numpy"));
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "conjugate"), NPY_QUAT, Quaternion_conj_ufunc, urnaryArgs, NULL);
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "reciprocal"), NPY_QUAT, Quaternion_inv_ufunc, urnaryArgs, NULL);
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "negative"), NPY_QUAT, Quaternion_neg_ufunc, urnaryArgs, NULL);

	int binaryArgs[3];
	urnaryArgs[0] = NPY_QUAT;
	urnaryArgs[1] = NPY_QUAT;
	binaryArgs[2] = NPY_QUAT;
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "multiply"), NPY_QUAT, Quaternion_mul_ufunc, binaryArgs, NULL);

	binaryArgs[2] = NPY_BOOL;
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "equal"),         NPY_QUAT, Quaternion_eq_ufunc, binaryArgs, NULL);
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "not_equal"),     NPY_QUAT, Quaternion_ne_ufunc, binaryArgs, NULL);
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "less"),          NPY_QUAT, Quaternion_lt_ufunc, binaryArgs, NULL);
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "less_equal"),    NPY_QUAT, Quaternion_le_ufunc, binaryArgs, NULL);
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "greater"),       NPY_QUAT, Quaternion_gt_ufunc, binaryArgs, NULL);
	PyUFunc_RegisterLoopForType((PyUFuncObject *)PyDict_GetItemString(ufuncDict, "greater_equal"), NPY_QUAT, Quaternion_ge_ufunc, binaryArgs, NULL);
	return 0;
}

//function to convert a numpy array of doubles to a numpy array of Quaternion objects
static PyObject* double2object(PyObject* self, PyObject* args, PyObject* kwds) {
	//parse arguments
	PyObject* array = NULL;
	int axis = -1;
	static char const* kwlist[] = {"array", "axis", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|i", const_cast<char**>(kwlist), &array, &axis)) return NULL;

	//get array object and its dimensions
	PyArrayObject* input = (PyArrayObject*)PyArray_FROM_OTF(array, NPY_FLOAT64, NPY_ARRAY_IN_ARRAY);
	if(input == NULL) {
		Py_XDECREF(input);
		return NULL;
	}
	int ndims = PyArray_NDIM(input);
	npy_intp* dims = PyArray_DIMS(input);

	//check shape
	if(axis == -1) axis += ndims;
	if(axis >= ndims || axis < 0) {
		std::stringstream ss;
		ss << "axis must be >= 0 and < " << ndims;
		PyErr_SetString(PyExc_ValueError, ss.str().c_str());
		Py_DECREF(input);
		return NULL;
	}

	if(dims[axis] < 4) {
		std::stringstream ss;
		ss << "axis " << axis << " is too small (expected 4, got " << dims[axis] << ")";
		PyErr_SetString(PyExc_ValueError, ss.str().c_str());
		Py_XDECREF(input);
		return NULL;
	} else if(dims[axis] > 4) {
		std::stringstream ss;
		ss << "axis " << axis << " is too large (expected 4, got " << dims[axis] << ")";
		PyErr_WarnEx(NULL, ss.str().c_str(), 1);
	}

	// get number of points and stride in axis
	int totalPoints = 1;
	for(int i = 0; i < ndims; i++) totalPoints *= dims[i];
	totalPoints /= dims[axis];
	npy_intp stride = PyArray_STRIDE(input, axis) / sizeof(double);

	//create output array of objects
	std::vector<npy_intp> newDims(dims, dims+ndims);
	newDims[axis] = 1;
	if(ndims > 1) newDims.erase(newDims.begin()+axis);
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(newDims.size(), newDims.data(), NPY_QUAT, 0);
	
	//get data
	double* to = (double*)PyArray_DATA(output);
	if(stride == 1) {
		//if the data is contigous pass pointers directly
		double* from = (double*)PyArray_DATA(input);
		const size_t delta = dims[axis];
		for(int i = 0; i < totalPoints; i++) {
			to[0] = from[0];
			to[1] = from[1];
			to[2] = from[2];
			to[3] = from[3];
			from += delta;
			to += 4;
		}
	} else {
		//if the data isn't contigous, copy to buffer
		std::vector<npy_intp> index(ndims, 0);
		for(int i = 0; i < totalPoints; i++) {
			double* from = (double*)PyArray_GetPtr(input, index.data());
			to[0] = from[0       ];
			to[1] = from[  stride];
			to[2] = from[2*stride];
			to[3] = from[3*stride];
			to += 4;
			index.back() += 1;
			for(int j = ndims-1; j > 0; j--) {
				if((index[j] > 0 && j == axis) || index[j] == dims[j]) {
					index[j] = 0;
					index[j-1] += 1;
				}
			}
		}
	}
	return (PyObject*)output;
}

//function to convert a numpy array of Quaternion objects to a numpy array of double
static PyObject* object2double(PyObject* self, PyObject* array) {
	//get array object and its dimensions
	PyArrayObject* input = (PyArrayObject*)PyArray_FROM_OTF(array, NPY_QUAT, NPY_ARRAY_IN_ARRAY);
	if(input == NULL) {
		Py_XDECREF(input);
		return NULL;
	}
	int ndims = PyArray_NDIM(input);
	npy_intp* dims = PyArray_DIMS(input);

	//get number of points
	int totalPoints = 4;
	for(int i = 0; i < ndims; i++) totalPoints *= dims[i];

	//create output array of doubles
	std::vector<npy_intp> newDims(dims, dims+ndims);
	newDims.push_back(4);
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(newDims.size(), newDims.data(), NPY_FLOAT64, 0);
	
	//get pointers and copy data
	double* from = (double*)PyArray_DATA(input);
	double* to = (double*)PyArray_DATA(output);
	for(int i = 0; i < totalPoints; i++) to[i] = from[i];
	return (PyObject*)output;
}

#endif//_quaternion_wrapper_h_