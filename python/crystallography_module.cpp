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
#include "quaternion_wrapper.hpp"
#include "rotations_wrapper.hpp"
#include "symmetry_wrapper.hpp"
#include "orientation_map_wrapper.hpp"

//method table
static PyMethodDef CrystallographyMethods[] = {//{function name in python module, c++ function, argument types, description}
	//other functions
	{"writeTif", (PyCFunction)Py_WriteTif, METH_VARARGS | METH_KEYWORDS, "writeTif(array, filename)\n\twrite array to tif (3d arrays will be interpreted as 2d images with multiple components per pixel"},

	//rotation functions
	{"setConvention",      (PyCFunction)setConvention      , METH_VARARGS                , "setConvention(convention)\n\tset the rotation convention (convention must be 'active' or 'passive')"},
	{"getConvention",      (PyCFunction)getConvention      , METH_NOARGS                 , "getConvention()\n\tget the rotation convention ('active' or 'passive')"                             },
	{"randomOrientations", (PyCFunction)random_orientations, METH_VARARGS | METH_KEYWORDS, "randomOrientations(count, representation = 'cu')\n\tgenerate random orientations"                   },

	//quaternion functions
	{"do2qu",  (PyCFunction)double2object, METH_VARARGS | METH_KEYWORDS, "do2qu(qu, axis = qu.ndim-1)\n\tconverts from array of floats array of Quaternion objects"  },
	{"qu2do",  (PyCFunction)object2double, METH_O                      , "qu2do(ob)\n\tconvert from array of Quaternion objects to array of floats"                  },

	//eu -->
	{"eu2om",  (PyCFunction)eu2om_wrapper, METH_VARARGS | METH_KEYWORDS, "eu2om(eu, axis = eu.ndim-1), see module level conversion routine documentation for details"},
	{"eu2ax",  (PyCFunction)eu2ax_wrapper, METH_VARARGS | METH_KEYWORDS, "eu2ax(eu, axis = eu.ndim-1), see module level conversion routine documentation for details"},
	{"eu2ro",  (PyCFunction)eu2ro_wrapper, METH_VARARGS | METH_KEYWORDS, "eu2ro(eu, axis = eu.ndim-1), see module level conversion routine documentation for details"},
	{"eu2qu",  (PyCFunction)eu2qu_wrapper, METH_VARARGS | METH_KEYWORDS, "eu2qu(eu, axis = eu.ndim-1), see module level conversion routine documentation for details"},
	{"eu2ho",  (PyCFunction)eu2ho_wrapper, METH_VARARGS | METH_KEYWORDS, "eu2ho(eu, axis = eu.ndim-1), see module level conversion routine documentation for details"},
	{"eu2cu",  (PyCFunction)eu2cu_wrapper, METH_VARARGS | METH_KEYWORDS, "eu2cu(eu, axis = eu.ndim-1), see module level conversion routine documentation for details"},

	//om -->
	{"om2eu",  (PyCFunction)om2eu_wrapper, METH_VARARGS | METH_KEYWORDS, "om2eu(om, axis = om.ndim-1), see module level conversion routine documentation for details"},
	{"om2ax",  (PyCFunction)om2ax_wrapper, METH_VARARGS | METH_KEYWORDS, "om2ax(om, axis = om.ndim-1), see module level conversion routine documentation for details"},
	{"om2ro",  (PyCFunction)om2ro_wrapper, METH_VARARGS | METH_KEYWORDS, "om2ro(om, axis = om.ndim-1), see module level conversion routine documentation for details"},
	{"om2qu",  (PyCFunction)om2qu_wrapper, METH_VARARGS | METH_KEYWORDS, "om2qu(om, axis = om.ndim-1), see module level conversion routine documentation for details"},
	{"om2ho",  (PyCFunction)om2ho_wrapper, METH_VARARGS | METH_KEYWORDS, "om2ho(om, axis = om.ndim-1), see module level conversion routine documentation for details"},
	{"om2cu",  (PyCFunction)om2cu_wrapper, METH_VARARGS | METH_KEYWORDS, "om2cu(om, axis = om.ndim-1), see module level conversion routine documentation for details"},

	//ax -->
	{"ax2eu",  (PyCFunction)ax2eu_wrapper, METH_VARARGS | METH_KEYWORDS, "ax2eu(ax, axis = ax.ndim-1), see module level conversion routine documentation for details"},
	{"ax2om",  (PyCFunction)ax2om_wrapper, METH_VARARGS | METH_KEYWORDS, "ax2om(ax, axis = ax.ndim-1), see module level conversion routine documentation for details"},
	{"ax2ro",  (PyCFunction)ax2ro_wrapper, METH_VARARGS | METH_KEYWORDS, "ax2ro(ax, axis = ax.ndim-1), see module level conversion routine documentation for details"},
	{"ax2qu",  (PyCFunction)ax2qu_wrapper, METH_VARARGS | METH_KEYWORDS, "ax2qu(ax, axis = ax.ndim-1), see module level conversion routine documentation for details"},
	{"ax2ho",  (PyCFunction)ax2ho_wrapper, METH_VARARGS | METH_KEYWORDS, "ax2ho(ax, axis = ax.ndim-1), see module level conversion routine documentation for details"},
	{"ax2cu",  (PyCFunction)ax2cu_wrapper, METH_VARARGS | METH_KEYWORDS, "ax2cu(ax, axis = ax.ndim-1), see module level conversion routine documentation for details"},

	//ro -->
	{"ro2eu",  (PyCFunction)ro2eu_wrapper, METH_VARARGS | METH_KEYWORDS, "ro2eu(ro, axis = ro.ndim-1), see module level conversion routine documentation for details"},
	{"ro2om",  (PyCFunction)ro2om_wrapper, METH_VARARGS | METH_KEYWORDS, "ro2om(ro, axis = ro.ndim-1), see module level conversion routine documentation for details"},
	{"ro2ax",  (PyCFunction)ro2ax_wrapper, METH_VARARGS | METH_KEYWORDS, "ro2ax(ro, axis = ro.ndim-1), see module level conversion routine documentation for details"},
	{"ro2qu",  (PyCFunction)ro2qu_wrapper, METH_VARARGS | METH_KEYWORDS, "ro2qu(ro, axis = ro.ndim-1), see module level conversion routine documentation for details"},
	{"ro2ho",  (PyCFunction)ro2ho_wrapper, METH_VARARGS | METH_KEYWORDS, "ro2ho(ro, axis = ro.ndim-1), see module level conversion routine documentation for details"},
	{"ro2cu",  (PyCFunction)ro2cu_wrapper, METH_VARARGS | METH_KEYWORDS, "ro2cu(ro, axis = ro.ndim-1), see module level conversion routine documentation for details"},

	//qu -->
	{"qu2eu",  (PyCFunction)qu2eu_wrapper, METH_O                      , "qu2eu(qu), see module level conversion routine documentation for details"                  },
	{"qu2om",  (PyCFunction)qu2om_wrapper, METH_O                      , "qu2om(qu), see module level conversion routine documentation for details"                  },
	{"qu2ax",  (PyCFunction)qu2ax_wrapper, METH_O                      , "qu2ax(qu), see module level conversion routine documentation for details"                  },
	{"qu2ro",  (PyCFunction)qu2ro_wrapper, METH_O                      , "qu2ro(qu), see module level conversion routine documentation for details"                  },
	{"qu2ho",  (PyCFunction)qu2ho_wrapper, METH_O                      , "qu2ho(qu), see module level conversion routine documentation for details"                  },
	{"qu2cu",  (PyCFunction)qu2cu_wrapper, METH_O                      , "qu2cu(qu), see module level conversion routine documentation for details"                  },

	//ho -->
	{"ho2eu",  (PyCFunction)ho2eu_wrapper, METH_VARARGS | METH_KEYWORDS, "ho2eu(ho, axis = ho.ndim-1), see module level conversion routine documentation for details"},
	{"ho2om",  (PyCFunction)ho2om_wrapper, METH_VARARGS | METH_KEYWORDS, "ho2om(ho, axis = ho.ndim-1), see module level conversion routine documentation for details"},
	{"ho2ax",  (PyCFunction)ho2ax_wrapper, METH_VARARGS | METH_KEYWORDS, "ho2ax(ho, axis = ho.ndim-1), see module level conversion routine documentation for details"},
	{"ho2ro",  (PyCFunction)ho2ro_wrapper, METH_VARARGS | METH_KEYWORDS, "ho2ro(ho, axis = ho.ndim-1), see module level conversion routine documentation for details"},
	{"ho2qu",  (PyCFunction)ho2qu_wrapper, METH_VARARGS | METH_KEYWORDS, "ho2qu(ho, axis = ho.ndim-1), see module level conversion routine documentation for details"},
	{"ho2cu",  (PyCFunction)ho2cu_wrapper, METH_VARARGS | METH_KEYWORDS, "ho2cu(ho, axis = ho.ndim-1), see module level conversion routine documentation for details"},

	//cu -->
	{"cu2eu",  (PyCFunction)cu2eu_wrapper, METH_VARARGS | METH_KEYWORDS, "cu2eu(cu, axis = cu.ndim-1), see module level conversion routine documentation for details"},
	{"cu2om",  (PyCFunction)cu2om_wrapper, METH_VARARGS | METH_KEYWORDS, "cu2om(cu, axis = cu.ndim-1), see module level conversion routine documentation for details"},
	{"cu2ax",  (PyCFunction)cu2ax_wrapper, METH_VARARGS | METH_KEYWORDS, "cu2ax(cu, axis = cu.ndim-1), see module level conversion routine documentation for details"},
	{"cu2ro",  (PyCFunction)cu2ro_wrapper, METH_VARARGS | METH_KEYWORDS, "cu2ro(cu, axis = cu.ndim-1), see module level conversion routine documentation for details"},
	{"cu2qu",  (PyCFunction)cu2qu_wrapper, METH_VARARGS | METH_KEYWORDS, "cu2qu(cu, axis = cu.ndim-1), see module level conversion routine documentation for details"},
	{"cu2ho",  (PyCFunction)cu2ho_wrapper, METH_VARARGS | METH_KEYWORDS, "cu2ho(cu, axis = cu.ndim-1), see module level conversion routine documentation for details"},
	{NULL, NULL, 0, NULL}//sentinel
};

//module definition
static struct PyModuleDef crystallographyModule = {
	PyModuleDef_HEAD_INIT,
	"crystallography", //module name
	//module documentation, may be NULL
	"\
orientation transform routines based on\n\
\n\
the following conventions are used:\n\
\t-rotation angle <= pi\n\
\t-rotation axis in positive z hemisphere for rotations of pi\n\
\t-rotation axis = [0, 0, 1] for rotations of 0\n\
\t-passive rotation convention (can be modified with 'setConvention'\n\
\n\
orientation conversion routine functions:\n\
\t-are based on:\n\
\t\t-Rowenhorst, David, et al. \"Consistent Representations of and Conversions Between 3D Rotations.\" Model. Simul. Mater. Sci. Eng. 23.8 (2015): 083501.\n\
\t\t-Roşca, D., et al. \"A New Method of Constructing a Grid in the Space of 3D rotations and its Applications to Texture Analysis.\" Model. Simul. Mater. Sci. Eng. 22.7 (2014): 075013.\n\
\t\t-fortran implementation of routines by Marc De Graef (https://github.com/marcdegraef/3Drotations)\n\
\t-are named __2__ (e.g. eu2om) with abbreviated representation names:\n\
\t\teu: euler angle triple (bunge convention, radians)\n\
\t\tom: orientation matrix (9 component vector or 3x3 matrix in row major order)\n\
\t\t\t-orientation matricies can be passed as 9 component vectors or 3x3 matricies\n\
\t\t\t\t-matricies are assumed by default with the second axis = axis+1\n\
\t\t\t\t-orientation matricies are assumed if there are at least 9 components in the axis or the the final axis is explicitly specified\
\t\tax: axis angle pair (nx, ny, nz, w)\n\
\t\tro: rodrigues vector (nx, ny, nz, tan(w/2))\n\
\t\tqu: unit quaternion (w, x, y, z) as a Quaternion object\n\
\t\tho: homochoric vector (x, y, z)\n\
\t\tcu: cubochoric vector (x, y ,z)\n\
\t-take a numpy array of doubles (or convertable argument) as the first argument and return a new numpy array\n\
\t-take the axis to convert as an optional second argument (defaults to the last axis or second to last for orientation matricies as 3x3 matricies)\n\
",
	-1,//size of per-interpreter state of the module or -1 if the module keeps state in global variables.
	CrystallographyMethods//method table
};

//module initialization function
PyMODINIT_FUNC PyInit_crystallography(void) {
	import_array();//import numpy
	import_umath();//import numpy umath functions

	//create module object
	PyObject* m = PyModule_Create(&crystallographyModule);
	if(m == NULL) return NULL;

	//register quaternion class
	QuaternionType.tp_name = "crystallography.Quaternion";//name must be preceded by modulename. to appear in module help
	// QuaternionType.tp_base = &PyGenericArrType_Type;//make quaternion type subclass of number array type
	if(PyType_Ready(&QuaternionType) < 0) return NULL;
	if(Quaternion_registerNumpy() < 0) return NULL;
	Py_INCREF(&QuaternionType);
	PyModule_AddObject(m, "quaternion", (PyObject *)&QuaternionType);

	//register symmetry class
	SymmetryType.tp_name = "crystallography.Symmetry";
	if(PyType_Ready(&SymmetryType) < 0) return NULL;
	Py_INCREF(&SymmetryType);
	PyModule_AddObject(m, "Symmetry", (PyObject*)&SymmetryType);

	//register orientation map class
	OrientationMapType.tp_name = "crystallography.OrientationMap";
	if(PyType_Ready(&OrientationMapType) < 0) return NULL;
	Py_INCREF(&OrientationMapType);
	PyModule_AddObject(m, "OrientationMap", (PyObject*)&OrientationMapType);

	//register segmentation class
	OrientationMapType.tp_name = "crystallography.Segmentation";
	if(PyType_Ready(&SegmentationType) < 0) return NULL;
	Py_INCREF(&SegmentationType);
	PyModule_AddObject(m, "Segmentation", (PyObject*)&SegmentationType);

	//add symmetry instances
	if(PyDict_SetItemString(SymmetryType.tp_dict, "Cubic"        , Symmetry_cubic        (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "CubicLow"     , Symmetry_cubiclow     (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "Hexagonal"    , Symmetry_hexagonal    (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "HexagonalLow" , Symmetry_hexagonallow (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "Tetragonal"   , Symmetry_tetragonal   (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "TetragonalLow", Symmetry_tetragonallow(NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "Trigonal"     , Symmetry_trigonal     (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "TrigonalLow"  , Symmetry_trigonallow  (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "Orthothombic" , Symmetry_orthorhombic (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "Monoclinic"   , Symmetry_monoclinic   (NULL)) < 0) return NULL;
	if(PyDict_SetItemString(SymmetryType.tp_dict, "Triclinic"    , Symmetry_triclinic    (NULL)) < 0) return NULL;
	
	return m;//return module object
}