# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                 #
# Copyright (c) 2017, William C. Lenthe                                           #
# All rights reserved.                                                            #
#                                                                                 #
# Redistribution and use in source and binary forms, with or without              #
# modification, are permitted provided that the following conditions are met:     #
#                                                                                 #
# 1. Redistributions of source code must retain the above copyright notice, this  #
#    list of conditions and the following disclaimer.                             #
#                                                                                 #
# 2. Redistributions in binary form must reproduce the above copyright notice,    #
#    this list of conditions and the following disclaimer in the documentation    #
#    and/or other materials provided with the distribution.                       #
#                                                                                 #
# 3. Neither the name of the copyright holder nor the names of its                #
#    contributors may be used to endorse or promote products derived from         #
#    this software without specific prior written permission.                     #
#                                                                                 #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"     #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE       #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE  #
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE    #
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL      #
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR      #
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER      #
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,   #
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE   #
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.            #
#                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
if(MINGW)
	add_definitions(-D_hypot=hypot) # if cmath is included on mingw there is a name conflict with python
endif()

# find python executable
find_package(PythonInterp 3 REQUIRED)

# find python include directory
execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "from distutils import sysconfig;print(sysconfig.get_python_inc(),end='')" OUTPUT_VARIABLE Python_Output ERROR_VARIABLE Python_Error)
STRING(REGEX REPLACE "\\\\" "/" Python_Output ${Python_Output}) # convert windows backslash to forward slash
if(NOT "${Python_Error}" STREQUAL "")
	message(FATAL_ERROR "failed to find python include directory: ${Python_Error}")
endif()
set(PYTHON_INCLUDE_DIR ${Python_Output} CACHE PATH "python include directories")
include_directories(${PYTHON_INCLUDE_DIR})

# find python library
if(WIN32)
	execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "from distutils import sysconfig;print(sysconfig.get_config_var(\"BINDIR\"),end='')" OUTPUT_VARIABLE Python_Output ERROR_VARIABLE Python_Error)
	STRING(REGEX REPLACE "\\\\" "/" Python_Output ${Python_Output})
	if(NOT "${Python_Error}" STREQUAL "")
		message(FATAL_ERROR "failed to find python include directory: ${Python_Error}")
	endif()
	set(PYTHON_LIBRARY_DIR ${Python_Output}/libs CACHE PATH "python library directories")

	if(MSVC)
		find_library(PYTHON_LIBRARY    python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} HINT ${PYTHON_LIBRARY_DIR} CACHE FILE)
	elseif(MINGW)
		find_library(PYTHON_LIBRARY libpython${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} HINT ${PYTHON_LIBRARY_DIR} CACHE FILE)
	endif()
else()
	execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "from distutils import sysconfig;print(sysconfig.get_config_var(\"LIBDIR\"),end='')" OUTPUT_VARIABLE Python_Output ERROR_VARIABLE Python_Error)
	if(NOT "${Python_Error}" STREQUAL "")
		message(FATAL_ERROR "failed to find python include directory: ${Python_Error}")
	endif()
	set(PYTHON_LIBRARY_DIR ${Python_Output} CACHE PATH "python library directories")
	find_library(PYTHON_LIBRARY libpython${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}m HINT ${PYTHON_LIBRARY_DIR}/x86_64-linux-gnu CACHE FILE)
endif()

# find numpy
execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "import numpy;print(numpy.get_include(),end='')" OUTPUT_VARIABLE Python_Output ERROR_VARIABLE Python_Error)
STRING(REGEX REPLACE "\\\\" "/" Python_Output ${Python_Output})
if(NOT "${Python_Error}" STREQUAL "")
	message(FATAL_ERROR "failed to find numpy include directory: ${Python_Error}")
endif()
find_path(NUMPY_INCLUDE_DIR numpy/arrayobject.h HINTS ${Python_Output})
get_filename_component(NUMPY_DIR ${NUMPY_INCLUDE_DIR} DIRECTORY)
include_directories(${NUMPY_INCLUDE_DIR})
find_library(NUMPY_LIBRARY NAMES npymath HINTS ${NUMPY_DIR}/lib)
get_filename_component(NUMPY_LIBRARY_DIR ${NUMPY_LIBRARY} DIRECTORY CACHE)

# add python extension
add_library(crystallography MODULE crystallography_module.cpp)
target_link_libraries(crystallography ${PYTHON_LIBRARY} ${NUMPY_LIBRARY})
set_property(TARGET crystallography PROPERTY CXX_STANDARD 11)
set_property(TARGET crystallography PROPERTY PREFIX "") # name crystallography instead of libcrystallography
if(WIN32)
	set_property(TARGET crystallography PROPERTY SUFFIX ".pyd") # name crystallography.pyd instead of crystallography.dll
else()
	set_property(TARGET crystallography PROPERTY SUFFIX ".so") # name crystallography.so
endif()

# add install rule to put python module in python path
if(WIN32)
	get_filename_component(PYTHON_DIR ${PYTHON_EXECUTABLE} DIRECTORY)
	install(TARGETS crystallography LIBRARY DESTINATION ${PYTHON_DIR}/Lib/site-packages)
else()
	install(TARGETS crystallography LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
endif()