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
#ifndef _rotations_wrapper_h_
#define _rotations_wrapper_h_

#include "Python.h"
#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#include "numpy/arrayobject.h"

#include "quaternion_wrapper.hpp"

#include "rotations.hpp"

#include <vector>
#include <sstream>
#include <random>
#include <limits>

//templated wrapper function for all possible rotation conversions except qu2*
template <size_t N_FROM, size_t N_TO, bool QuatOut = false>
static PyObject* rotation_wrapper(PyObject* self, PyObject* args, PyObject* kwds, void(*f)(double const * const, double * const)) {
	//parse arguments
	PyObject* array = NULL;
	int axis = std::numeric_limits<int>::min();
	static char const* kwlist[] = {"array", "axis", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|i", const_cast<char**>(kwlist), &array, &axis)) return NULL;
	bool axisSet = axis != std::numeric_limits<int>::min();
	if(!axisSet) axis = -1;

	//get array object as doubles and its dimensions
	PyArrayObject* input = (PyArrayObject*)PyArray_FROM_OTF(array, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
	if(input == NULL) {
		Py_XDECREF(input);
		return NULL;
	}
	int ndims = PyArray_NDIM(input);
	npy_intp* dims = PyArray_DIMS(input);

	//check shape
	if(axis < 0) axis += ndims;
	if(axis >= ndims || axis < 0) {
		std::stringstream ss;
		ss << "axis must be >= 0 and < " << ndims;
		PyErr_SetString(PyExc_ValueError, ss.str().c_str());
		Py_DECREF(input);
		return NULL;
	}

	//handle orientation matrix
	bool asMatrix = false;
	if(9 == N_FROM) {
		if(axisSet) {
			//axis specified
			asMatrix = dims[axis] < 9 && axis+1 < ndims;
		} else {
			//default axis
			if(dims[axis] < 9 && ndims > 1) {
				asMatrix = true;
				--axis;
			}
		}
	}

	std::vector<npy_intp> vdims;
	if(asMatrix) {
		//check first axis of matrix
		if(dims[axis] < 3) {
			std::stringstream ss;
			ss << "axis " << axis << " is too small (expected 3, got " << dims[axis] << ")";
			PyErr_SetString(PyExc_ValueError, ss.str().c_str());
			Py_XDECREF(input);
			return NULL;
		} else if(dims[axis] > 3) {
			std::stringstream ss;
			ss << "axis " << axis << " is too large (expected 3, got " << dims[axis] << ")";
			PyErr_WarnEx(NULL, ss.str().c_str(), 1);
		}

		//check second axis of matrix (must be exactly 3 or strides will be long)
		if(dims[axis+1] != 3) {
			std::stringstream ss;
			ss << "axis " << axis+1 << " is wrong (expected 3, got " << dims[axis+1] << ")";
			PyErr_SetString(PyExc_ValueError, ss.str().c_str());
			Py_XDECREF(input);
			return NULL;
		}

		//reshape from 3x3 to 9 component vector
		for(int i = 0; i <= axis; i++) vdims.push_back(dims[i]);
		vdims.back() *= 3;
		for(int i = axis+2; i < ndims; i++) vdims.push_back(dims[i]);
		PyArray_Dims newDims;
		newDims.ptr = vdims.data();
		newDims.len = (int)vdims.size();
		input = (PyArrayObject*)PyArray_Newshape(input, &newDims, NPY_CORDER);

		ndims -= 1;
		dims = vdims.data();
		npy_intp* dims =  PyArray_DIMS(input);
	} else {
		if(dims[axis] < N_FROM) {
			std::stringstream ss;
			ss << "axis " << axis << " is too small (expected " << N_FROM << ", got " << dims[axis] << ")";
			PyErr_SetString(PyExc_ValueError, ss.str().c_str());
			Py_XDECREF(input);
			return NULL;
		} else if(dims[axis] > N_FROM) {
			std::stringstream ss;
			ss << "axis " << axis << " is too large (expected " << N_FROM << ", got " << dims[axis] << ")";
			PyErr_WarnEx(NULL, ss.str().c_str(), 1);
		}
	}

	//get number of points and stride in axis
	int totalPoints = 1;
	for(npy_intp i = 0; i < ndims; i++) totalPoints *= (int)dims[i];
	totalPoints /= (int)dims[axis];
	npy_intp stride = PyArray_STRIDE(input, axis) / sizeof(double);

	//create output array
	std::vector<npy_intp> newDims(dims, dims+ndims);
	if(QuatOut) {
		newDims.erase(newDims.begin() + axis);
		if(newDims.empty()) newDims.push_back(1);
	} else newDims[axis] = N_TO;
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY((int)newDims.size(), newDims.data(), QuatOut ? NPY_QUAT : NPY_DOUBLE, 0);
	
	//get data
	if(stride == 1) {
		//if the data is contigous pass pointers directly
		double* from = (double*)PyArray_DATA(input);
		double* to = (double*)PyArray_DATA(output);
		for(int i = 0; i < totalPoints; i++) {
			f(from, to);//execute conversion
			from += dims[axis];
			to += N_TO;
		}
	} else {
		//if the data isn't contigous, copy to buffer
		double fromBuff[N_FROM], toBuff[N_TO];
		std::vector<npy_intp> index(ndims, 0);
		if(QuatOut) {
			double* to = (double*)PyArray_DATA(output);
			for(int i = 0; i < totalPoints; i++) {
				double* from = (double*)PyArray_GetPtr(input, index.data());
				for(int j = 0; j < N_FROM; j++) fromBuff[j] = from[j*stride];
				f(fromBuff, to);//execute conversion
				to += N_TO;

				index.back() += 1;
				for(int j = ndims-1; j > 0; j--) {
					if((index[j] > 0 && j == axis) || index[j] == dims[j]) {
						index[j] = 0;
						index[j-1] += 1;
					}
				}
			}

		} else {
			for(int i = 0; i < totalPoints; i++) {
				double* from = (double*)PyArray_GetPtr(input, index.data());
				double* to = (double*)PyArray_GetPtr(output, index.data());
				for(int j = 0; j < N_FROM; j++) fromBuff[j] = from[j*stride];
				f(fromBuff, toBuff);//execute conversion
				for(int j = 0; j < N_TO; j++) to[j*stride] = toBuff[j];

				index.back() += 1;
				for(int j = ndims-1; j > 0; j--) {
					if((index[j] > 0 && j == axis) || index[j] == dims[j]) {
						index[j] = 0;
						index[j-1] += 1;
					}
				}
			}
		}
	}

	if(9 == N_TO) {
		//reshape orientation matricies to 3x3
		for(int i = 0; i < axis; i++) vdims.push_back(dims[i]);
		vdims.push_back(3);
		vdims.push_back(3);
		for(int i = axis+1; i < ndims; i++) vdims.push_back(dims[i]);
		PyArray_Dims newDims;
		newDims.ptr = vdims.data();
		newDims.len = (int)vdims.size();
		return PyArray_Newshape(output, &newDims, NPY_CORDER);
	} else {
		return (PyObject*)output;
	}
}

//templated wrapper function for qu2* rotation conversions
template <size_t N_TO>
static PyObject* rotation_wrapper_qu(PyObject* self, PyObject* other, void(*f)(double const * const, double * const)) {
	//get array object as doubles and its dimensions
	PyArrayObject* input = (PyArrayObject*)PyArray_FROM_OTF(other, NPY_QUAT, NPY_ARRAY_IN_ARRAY);
	if(input == NULL) {
		Py_XDECREF(input);
		return NULL;
	}
	int ndims = PyArray_NDIM(input);
	npy_intp* dims =  PyArray_DIMS(input);

	//create output array
	std::vector<npy_intp> newDims(dims, dims+ndims);
	newDims.push_back(N_TO);
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY((int)newDims.size(), newDims.data(), NPY_DOUBLE, 0);

	//get number of points
	int totalPoints = 1;
	for(npy_intp i = 0; i < ndims; i++) totalPoints *= (int)dims[i];
	
	//convert
	double* from = (double*)PyArray_DATA(input);
	double* to = (double*)PyArray_DATA(output);
	for(int i = 0; i < totalPoints; i++) {
		f(from, to);//execute conversion
		from += 4;
		to += N_TO;
	}

	if(9 == N_TO) {
		//reshape orientation matricies to 3x3
		newDims.back() = 3;
		newDims.push_back(3);
		PyArray_Dims reshapedDims;
		reshapedDims.ptr = newDims.data();
		reshapedDims.len = (int)newDims.size();
		return PyArray_Newshape(output, &reshapedDims, NPY_CORDER);
	} else {
		return (PyObject*)output;
	}
}

static PyObject* setConvention(PyObject* self, PyObject* args) {
	//parse arguments
	char* convention = NULL;
	if(!PyArg_ParseTuple(args, "s", &convention)) return NULL;

	//set convention
	if(0 == strcmp(convention, "active")) rotations::Constants<double>::convention = rotations::Constants<double>::active;
	else if(0 == strcmp(convention, "passive")) rotations::Constants<double>::convention = rotations::Constants<double>::passive;
	else {
		std::stringstream ss;
		ss << "invalid convention '" << convention << "'. Convention must be 'active' or 'passive'";
		PyErr_SetString(PyExc_ValueError, ss.str().c_str());
	}
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject* getConvention(PyObject* self, PyObject* arg) {
	return Py_BuildValue("s", rotations::Constants<double>::passive == rotations::Constants<double>::convention ? "passive" : "active");
}

static PyObject* random_orientations(PyObject* self, PyObject* args, PyObject* kwds) {
	//parse arguments
	unsigned int count = 0;
	char* representation = NULL;
	static char const* kwlist[] = {"count", "representation", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "I|s", const_cast<char**>(kwlist), &count, &representation)) return NULL;
	size_t repr = 6;
	if(NULL == representation) repr = 6;
	else if(0 == strcmp("eu", representation)) repr = 0;
	else if(0 == strcmp("om", representation)) repr = 1;
	else if(0 == strcmp("ax", representation)) repr = 2;
	else if(0 == strcmp("ro", representation)) repr = 3;
	else if(0 == strcmp("qu", representation)) repr = 4;
	else if(0 == strcmp("ho", representation)) repr = 5;
	else if(0 == strcmp("cu", representation)) repr = 6;
	else {
		std::stringstream ss;
		ss << "invalid representation '" << representation << "'. Convention must be 'eu', 'om', 'ax', 'ro', 'qu', 'ho', 'cu', or None ('cu')";
		PyErr_SetString(PyExc_ValueError, ss.str().c_str());
		return NULL;
	}

	//create output array
	npy_intp dims[2];
	dims[0] = count;
	dims[1] = 3;
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(2, dims, NPY_DOUBLE, 0);
	
	//fill with random cubochoric vectors
	double* buff = (double*)PyArray_DATA(output);
	std::random_device seed;
	std::mt19937 generator(std::random_device{}());//seed mersenne twister
	double a2 = rotations::Constants<double>::cuA / 2.0;
	std::uniform_real_distribution<double> distribution(-a2, a2);
	for(size_t i = 0; i < 3*count; i++) buff[i] = distribution(generator);

	//convert to appropriate representation
	if(6 == repr) {
		return (PyObject*)output;//return unchanged cu's
	} else if(0 == repr) {
		for(size_t i = 0; i < count; i++) rotations::cu2eu(buff+3*i, buff+3*i);//convert to eu in place
		return (PyObject*)output;
	}else if(5 == repr) {
		for(size_t i = 0; i < count; i++) rotations::cu2ho(buff+3*i, buff+3*i);//convert to ho in place
		return (PyObject*)output;
	} else {
		//create new array of appropriate size
		if(4 == repr) {
			PyArrayObject* converted = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_QUAT, 0);
			Quaternion<double>* buffConv = (Quaternion<double>*)PyArray_DATA(converted);
			for(size_t i = 0; i < count; i++) rotations::cu2qu(buff+3*i, buffConv[i].data());
			return (PyObject*) converted;
		} else {
			if(1 == repr) dims[1] = 9;
			else dims[1] = 4;
			PyArrayObject* converted = (PyArrayObject*)PyArray_EMPTY(2, dims, NPY_DOUBLE, 0);
			double* buffConv = (double*)PyArray_DATA(converted);

			if(1 == repr) {
				for(size_t i = 0; i < count; i++) rotations::cu2om(buff+3*i, buffConv+9*i);
				npy_intp d[3];
				d[0] = count;
				d[1] = d[2] = 3;	
				PyArray_Dims newDims;
				newDims.len = 3;
				newDims.ptr = d;
				return PyArray_Newshape(converted, &newDims, NPY_CORDER);
			} else if(2 == repr) {
				for(size_t i = 0; i < count; i++) rotations::cu2ax(buff+3*i, buffConv+4*i);
			} else if(3 == repr) {
				for(size_t i = 0; i < count; i++) rotations::cu2ro(buff+3*i, buffConv+4*i);
			}
			return (PyObject*)converted;
		}
	}
}

//wrapper functions for each conversion routine
//eu -->
static PyObject* eu2om_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,9     >(self, args, kwds, rotations::eu2om<double>);}
static PyObject* eu2ax_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4     >(self, args, kwds, rotations::eu2ax<double>);}
static PyObject* eu2ro_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4     >(self, args, kwds, rotations::eu2ro<double>);}
static PyObject* eu2qu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4,true>(self, args, kwds, rotations::eu2qu<double>);}
static PyObject* eu2ho_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,3     >(self, args, kwds, rotations::eu2ho<double>);}
static PyObject* eu2cu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,3     >(self, args, kwds, rotations::eu2cu<double>);}

//om -->
static PyObject* om2eu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<9,3     >(self, args, kwds, rotations::om2eu<double>);}
static PyObject* om2ax_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<9,4     >(self, args, kwds, rotations::om2ax<double>);}
static PyObject* om2ro_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<9,4     >(self, args, kwds, rotations::om2ro<double>);}
static PyObject* om2qu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<9,4,true>(self, args, kwds, rotations::om2qu<double>);}
static PyObject* om2ho_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<9,3     >(self, args, kwds, rotations::om2ho<double>);}
static PyObject* om2cu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<9,3     >(self, args, kwds, rotations::om2cu<double>);}

//ax -->
static PyObject* ax2eu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,3     >(self, args, kwds, rotations::ax2eu<double>);}
static PyObject* ax2om_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,9     >(self, args, kwds, rotations::ax2om<double>);}
static PyObject* ax2ro_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,4     >(self, args, kwds, rotations::ax2ro<double>);}
static PyObject* ax2qu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,4,true>(self, args, kwds, rotations::ax2qu<double>);}
static PyObject* ax2ho_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,3     >(self, args, kwds, rotations::ax2ho<double>);}
static PyObject* ax2cu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,3     >(self, args, kwds, rotations::ax2cu<double>);}

//ro -->
static PyObject* ro2eu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,3     >(self, args, kwds, rotations::ro2eu<double>);}
static PyObject* ro2om_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,9     >(self, args, kwds, rotations::ro2om<double>);}
static PyObject* ro2ax_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,4     >(self, args, kwds, rotations::ro2ax<double>);}
static PyObject* ro2qu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,4,true>(self, args, kwds, rotations::ro2qu<double>);}
static PyObject* ro2ho_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,3     >(self, args, kwds, rotations::ro2ho<double>);}
static PyObject* ro2cu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<4,3     >(self, args, kwds, rotations::ro2cu<double>);}

//qu -->
static PyObject* qu2eu_wrapper(PyObject* self, PyObject* other) {return rotation_wrapper_qu<3>(self, other, rotations::qu2eu<double>);}
static PyObject* qu2om_wrapper(PyObject* self, PyObject* other) {return rotation_wrapper_qu<9>(self, other, rotations::qu2om<double>);}
static PyObject* qu2ax_wrapper(PyObject* self, PyObject* other) {return rotation_wrapper_qu<4>(self, other, rotations::qu2ax<double>);}
static PyObject* qu2ro_wrapper(PyObject* self, PyObject* other) {return rotation_wrapper_qu<4>(self, other, rotations::qu2ro<double>);}
static PyObject* qu2ho_wrapper(PyObject* self, PyObject* other) {return rotation_wrapper_qu<3>(self, other, rotations::qu2ho<double>);}
static PyObject* qu2cu_wrapper(PyObject* self, PyObject* other) {return rotation_wrapper_qu<3>(self, other, rotations::qu2cu<double>);}

//ho -->
static PyObject* ho2eu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,3     >(self, args, kwds, rotations::ho2eu<double>);}
static PyObject* ho2om_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,9     >(self, args, kwds, rotations::ho2om<double>);}
static PyObject* ho2ax_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4     >(self, args, kwds, rotations::ho2ax<double>);}
static PyObject* ho2ro_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4     >(self, args, kwds, rotations::ho2ro<double>);}
static PyObject* ho2qu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4,true>(self, args, kwds, rotations::ho2qu<double>);}
static PyObject* ho2cu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,3     >(self, args, kwds, rotations::ho2cu<double>);}

//cu -->
static PyObject* cu2eu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,3     >(self, args, kwds, rotations::cu2eu<double>);}
static PyObject* cu2om_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,9     >(self, args, kwds, rotations::cu2om<double>);}
static PyObject* cu2ax_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4     >(self, args, kwds, rotations::cu2ax<double>);}
static PyObject* cu2ro_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4     >(self, args, kwds, rotations::cu2ro<double>);}
static PyObject* cu2qu_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,4,true>(self, args, kwds, rotations::cu2qu<double>);}
static PyObject* cu2ho_wrapper(PyObject* self, PyObject* args, PyObject* kwds) {return rotation_wrapper<3,3     >(self, args, kwds, rotations::cu2ho<double>);}

#endif//_rotations_wrapper_h_