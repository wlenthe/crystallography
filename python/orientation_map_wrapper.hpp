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
#ifndef _orientation_map_wrapper_h_
#define _orientation_map_wrapper_h_

#include "Python.h"
#include "structmember.h"//PyMemberDef full definition
#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#include "numpy/arrayobject.h"

#include "quaternion_wrapper.hpp"
#include "symmetry_wrapper.hpp"

#include "orientation_map.hpp"

#include <vector>

typedef struct {
	PyObject_HEAD
	OrientationMap<double> om;
	PyArrayObject* quats;
	PyArrayObject* phases;
	PyArrayObject* quality;
	PyArrayObject* syms;
	std::uint32_t rows, cols;
	double xRes, yRes;
} PyOrientationMap;
static bool PyOrientationMap_Check(PyObject* obj);

static void OrientationMap_dealloc(PyOrientationMap* self){
	//PyExc_AttributeError
	if(NULL != self->quats  ) Py_XDECREF(self->quats);
	if(NULL != self->phases ) Py_XDECREF(self->phases);
	if(NULL != self->quality) Py_XDECREF(self->quality);
	if(NULL != self->syms   ) Py_XDECREF(self->syms);
	Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* OrientationMap_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
	PyOrientationMap* self = (PyOrientationMap*)type->tp_alloc(type, 0);
	if(NULL != self) {
		self->quats   = (PyArrayObject*)NULL;
		self->phases  = (PyArrayObject*)NULL;
		self->quality = (PyArrayObject*)NULL;
		self->syms    = (PyArrayObject*)NULL;
	}
	return (PyObject*)self;
}

int PyOrientationMap_readFile(PyOrientationMap* map, std::string fileName) {
	try {
		map->om = OrientationMap<double>(fileName);
	} catch (std::exception& e) {
		std::string error = "error reading " + fileName + ": " + e.what();
		PyErr_SetString(PyExc_IOError, error.c_str());
		return 1;
	}
	//copy meta data
	map->rows = map->om.rows;
	map->cols = map->om.cols;
	map->xRes = map->om.xRes;
	map->yRes = map->om.yRes;

	//copy symmetry
	std::vector<npy_intp> dims(1, map->om.syms.size());
	map->syms = (PyArrayObject*)PyArray_EMPTY(dims.size(), dims.data(), NPY_OBJECT, 0);
	PySymmetry** pSyms = (PySymmetry**)PyArray_DATA(map->syms);
	for(size_t i = 0; i < map->om.syms.size(); i++) pSyms[i] = Symmetry_create(map->om.syms[i]);

	//wrap data
	dims[0] = map->om.rows;
	dims.push_back(map->om.cols);
	int typenum = NPY_NOTYPE;
	if(std::is_same<size_t, std::uint32_t>::value) typenum = NPY_UINT32;
	else if(std::is_same<size_t, std::uint64_t>::value) typenum = NPY_UINT64;
	static_assert(std::is_same<size_t, std::uint32_t>::value || std::is_same<size_t, std::uint64_t>::value, "could not determine numpy type for size_t");
	map->phases  = (PyArrayObject*)PyArray_SimpleNewFromData(dims.size(), dims.data(), typenum   , (void*)map->om.phases->data());
	map->quality = (PyArrayObject*)PyArray_SimpleNewFromData(dims.size(), dims.data(), NPY_DOUBLE, (void*)map->om.quality->data());
	map->quats   = (PyArrayObject*)PyArray_SimpleNewFromData(dims.size(), dims.data(), NPY_QUAT  , (void*)map->om.quats->data());
	return 0;
}

static int OrientationMap_init(PyOrientationMap* self, PyObject* args, PyObject* kwds) {
	char* fileName;
	const char *kwlist[] = {"fileName", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "s", const_cast<char**>(kwlist), &fileName)) return -1;
	return PyOrientationMap_readFile(self, fileName);
}

//method table
static PyMethodDef OrientationMap_methods[] = {
	{NULL}//sentinel
};

//member table
static PyMemberDef OrientationMap_members[] = {
	{(char*)"quats"     , T_OBJECT, offsetof(PyOrientationMap, quats)  , READONLY, (char*)"Quaternoins"                  },
	{(char*)"phases"    , T_OBJECT, offsetof(PyOrientationMap, phases) , READONLY, (char*)"phase (index into symmetries)"},
	{(char*)"quality"   , T_OBJECT, offsetof(PyOrientationMap, quality), READONLY, (char*)"index quality"                },
	{(char*)"symmetries", T_OBJECT, offsetof(PyOrientationMap, syms)   , READONLY, (char*)"phase symmetries"             },
	{(char*)"rows"      , T_UINT  , offsetof(PyOrientationMap, rows)   , READONLY, (char*)"map rows"                     },
	{(char*)"cols"      , T_UINT  , offsetof(PyOrientationMap, cols)   , READONLY, (char*)"map columns"                  },
	{(char*)"xRes"      , T_DOUBLE, offsetof(PyOrientationMap, xRes)   , READONLY, (char*)"map x resolution"             },
	{(char*)"yRes"      , T_DOUBLE, offsetof(PyOrientationMap, yRes)   , READONLY, (char*)"map y resolution"             },
	{NULL}//sentinel
};

//type definition
static PyTypeObject OrientationMapType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"OrientationMap",                      //tp_name
	sizeof(PyOrientationMap),              //tp_basicsize
	0,                                     //tp_itemsize
	(destructor)OrientationMap_dealloc,    //tp_dealloc
	0,                                     //tp_print
	0,                                     //tp_getattr
	0,                                     //tp_setattr
	0,                                     //tp_reserved
	0,                                     //tp_repr
	0,                                     //tp_as_number
	0,                                     //tp_as_sequence
	0,                                     //tp_as_mapping
	0,                                     //tp_hash 
	0,                                     //tp_call
	0,                                     //tp_str
	0,                                     //tp_getattro
	0,                                     //tp_setattro
	0,                                     //tp_as_buffer
	Py_TPFLAGS_DEFAULT|Py_TPFLAGS_BASETYPE,//tp_flags
	"OrientationMap(filename)",            //tp_doc
	0,                                     //tp_traverse
	0,                                     //tp_clear
	0,                                     //tp_richcompare
	0,                                     //tp_weaklistoffset
	0,                                     //tp_iter
	0,                                     //tp_iternext
	OrientationMap_methods,                //tp_methods
	OrientationMap_members,                //tp_members
	0,                                     //tp_getset
	0,                                     //tp_base
	0,                                     //tp_dict
	0,                                     //tp_descr_get
	0,                                     //tp_descr_set
	0,                                     //tp_dictoffset
	(initproc)OrientationMap_init,         //tp_init
	0,                                     //tp_alloc
	OrientationMap_new,                    //tp_new
};
static PyOrientationMap* OrientationMap_create() {return (PyOrientationMap*)OrientationMapType.tp_alloc(&OrientationMapType, 0);}
static bool PyOrientationMap_Check(PyObject* obj) {return (PyTypeObject*)(obj->ob_type) == &OrientationMapType;}

typedef struct {
	// PyObject_HEAD
	PyOrientationMap super;
	//orientationmap members
	PyArrayObject* quats;
	PyArrayObject* phases;
	PyArrayObject* quality;
	PyArrayObject* syms;
	std::uint32_t rows, cols;
	double xRes, yRes;

	//segmentation members
	std::shared_ptr< Segmentation<double > > seg;//use shared_ptr instead of instance so the struct stays POD type and offsetof can be used
	PyArrayObject* ids;
	PyArrayObject* avgQuats;
	PyArrayObject* avgPhases;
	PyArrayObject* numPixels;
	std::uint32_t numGrains;
} PySegmentation;
static bool PySegmentation_Check(PyObject* obj);

static void Segmentation_dealloc(PySegmentation* self){
	if(NULL != self->ids      ) Py_XDECREF(self->ids);
	if(NULL != self->avgQuats ) Py_XDECREF(self->avgQuats);
	if(NULL != self->avgPhases) Py_XDECREF(self->avgPhases);
	if(NULL != self->numPixels) Py_XDECREF(self->numPixels);
	Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* Segmentation_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
	PySegmentation* self = (PySegmentation*)type->tp_alloc(type, 0);
	if(NULL != self) {
		self->ids       = (PyArrayObject*)NULL;
		self->avgQuats  = (PyArrayObject*)NULL;
		self->avgPhases = (PyArrayObject*)NULL;
		self->numPixels = (PyArrayObject*)NULL;
	}
	return (PyObject*)self;
}

static int Segmentation_init(PySegmentation* self, PyObject* args, PyObject* kwds) {
	PyObject* omObj;
	double tolerance = 5.0;
	const char *kwlist[] = {"orientation_map", "tolerance", NULL};
	if(!PyArg_ParseTupleAndKeywords(args, kwds, "O|d", const_cast<char**>(kwlist), &omObj, &tolerance)) return -1;
	if(!PyOrientationMap_Check(omObj)) {
		PyErr_SetString(PyExc_ValueError, "orientation_map must be OrientationMap class");
		return -1;
	}
	//copy parent class
	PyOrientationMap* om = (PyOrientationMap*) omObj;
	self->super.om      = om->om;
	self->super.quats   = om->quats;
	self->super.phases  = om->phases;
	self->super.quality = om->quality;
	self->super.syms    = om->syms;
	self->super.rows    = om->rows;
	self->super.cols    = om->cols;
	self->super.xRes    = om->xRes;
	self->super.yRes    = om->yRes;
	Py_INCREF(om->quats);
	Py_INCREF(om->phases);
	Py_INCREF(om->quality);
	Py_INCREF(om->syms);

	//perform segmentation
	self->seg = std::make_shared< Segmentation<double> >(om->om, tolerance);
	self->numGrains = self->seg->numGrains;

	//wrap grain level data
	std::vector<npy_intp> dims(1, self->seg->numGrains);
	int typenum = NPY_NOTYPE;
	if(std::is_same<size_t, std::uint32_t>::value) typenum = NPY_UINT32;
	else if(std::is_same<size_t, std::uint64_t>::value) typenum = NPY_UINT64;
	static_assert(std::is_same<size_t, std::uint32_t>::value || std::is_same<size_t, std::uint64_t>::value, "could not determine numpy type for size_t");
	self->avgPhases = (PyArrayObject*)PyArray_SimpleNewFromData(dims.size(), dims.data(), typenum , (void*)self->seg->avgPhases.data());
	self->numPixels = (PyArrayObject*)PyArray_SimpleNewFromData(dims.size(), dims.data(), typenum , (void*)self->seg->numPixels.data());
	self->avgQuats  = (PyArrayObject*)PyArray_SimpleNewFromData(dims.size(), dims.data(), NPY_QUAT, (void*)self->seg->avgQuats.data());
	
	//wrap voxel level data
	dims[0] = self->seg->rows;
	dims.push_back(self->seg->cols);
	self->ids = (PyArrayObject*)PyArray_SimpleNewFromData(dims.size(), dims.data(),  typenum, (void*)self->seg->ids.data());
	return 0;
}

static PyObject* Symmetry_mapColor(PySegmentation* self) {
	npy_intp dims[2];
	dims[0] = self->seg->rows;
	dims[1] = self->seg->cols;
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(2, dims, NPY_UINT8, 0);
	std::uint8_t* buff = (std::uint8_t*)PyArray_DATA(output);
	std::vector<std::uint8_t> colors = self->seg->mapColor();
	std::copy(colors.begin(), colors.end(), buff);
	return (PyObject*) output;
}

static PyObject* Symmetry_neighbors(PySegmentation* self) {
	//get neighbors and determine type of size_t
	int typenum = NPY_NOTYPE;
	std::vector< std::set<size_t> > neighbors = self->seg->findNeighbors();
	if(std::is_same<size_t, std::uint32_t>::value) typenum = NPY_UINT32;
	else if(std::is_same<size_t, std::uint64_t>::value) typenum = NPY_UINT64;
	static_assert(std::is_same<size_t, std::uint32_t>::value || std::is_same<size_t, std::uint64_t>::value, "could not determine numpy type for size_t");

	//create array to hold array for each grain
	npy_intp dims[1];
	dims[0] = neighbors.size();
	PyArrayObject* output = (PyArrayObject*)PyArray_EMPTY(1, dims, NPY_OBJECT, 0);
	PyArrayObject** pArray = (PyArrayObject**)PyArray_DATA(output);
	for(size_t i = 0; i < neighbors.size(); i++) {
		size_t j = 0;
		dims[0] = neighbors[i].size();
		pArray[i] = (PyArrayObject*)PyArray_EMPTY(1, dims, typenum, 0);
		size_t* buff = (size_t*)PyArray_DATA(pArray[i]);
		for(std::set<size_t>::const_iterator iter = neighbors[i].cbegin(); iter != neighbors[i].cend(); ++iter) buff[j++] = *iter;
	}
	return (PyObject*) output;
}

//method table
static PyMethodDef Segmentation_methods[] = {
	{"mapColor", (PyCFunction)Symmetry_mapColor, METH_NOARGS, "om.mapColor()\n\tcolor grains such that no touching grains are the same"},
	{"neighbors", (PyCFunction)Symmetry_neighbors, METH_NOARGS, "om.mapGrainToPixel(grainArray)\n\tcopy values from a per grain array to a per pixel array"},
	{NULL}//sentinel
};

//member table
static PyMemberDef Segmentation_members[] = {
	{(char*)"ids"      , T_OBJECT, offsetof(PySegmentation, ids      ), READONLY, (char*)"grain ids"                           },
	{(char*)"avgQuats" , T_OBJECT, offsetof(PySegmentation, avgQuats ), READONLY, (char*)"grain average orientations"          },
	{(char*)"avgPhases", T_OBJECT, offsetof(PySegmentation, avgPhases), READONLY, (char*)"grain phases"                        },
	{(char*)"numPixels", T_OBJECT, offsetof(PySegmentation, numPixels), READONLY, (char*)"grain sizes"                         },
	{(char*)"numGrains", T_UINT  , offsetof(PySegmentation, numGrains), READONLY, (char*)"number of grains (including grain 0)"},
	{NULL}//sentinel
};

//type definition
static PyTypeObject SegmentationType = {
	PyVarObject_HEAD_INIT(NULL, 0)
	"Segmentation",                                         //tp_name
	sizeof(PySegmentation),                                 //tp_basicsize
	0,                                                      //tp_itemsize
	(destructor)Segmentation_dealloc,                       //tp_dealloc
	0,                                                      //tp_print
	0,                                                      //tp_getattr
	0,                                                      //tp_setattr
	0,                                                      //tp_reserved
	0,                                                      //tp_repr
	0,                                                      //tp_as_number
	0,                                                      //tp_as_sequence
	0,                                                      //tp_as_mapping
	0,                                                      //tp_hash 
	0,                                                      //tp_call
	0,                                                      //tp_str
	0,                                                      //tp_getattro
	0,                                                      //tp_setattro
	0,                                                      //tp_as_buffer
	Py_TPFLAGS_DEFAULT,                                     //tp_flags
	"Segmentation(orientation_map | tolerance, minQuality)",//tp_doc
	0,                                                      //tp_traverse
	0,                                                      //tp_clear
	0,                                                      //tp_richcompare
	0,                                                      //tp_weaklistoffset
	0,                                                      //tp_iter
	0,                                                      //tp_iternext
	Segmentation_methods,                                   //tp_methods
	Segmentation_members,                                   //tp_members
	0,                                                      //tp_getset
	&OrientationMapType,                                    //tp_base
	0,                                                      //tp_dict
	0,                                                      //tp_descr_get
	0,                                                      //tp_descr_set
	0,                                                      //tp_dictoffset
	(initproc)Segmentation_init,                            //tp_init
	0,                                                      //tp_alloc
	Segmentation_new,                                       //tp_new
};
static PySegmentation* Segmentation_create() {return (PySegmentation*)SegmentationType.tp_alloc(&SegmentationType, 0);}
static bool PySegmentation_Check(PyObject* obj) {return (PyTypeObject*)(obj->ob_type) == &SegmentationType;}

#endif//_orientation_map_wrapper_h_
