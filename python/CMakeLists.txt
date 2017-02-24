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

# find python
find_package(PythonInterp REQUIRED) # should also define PYTHON_INCLUDE_DIRS and PYTHON_LIBRARIES but seems to fail
get_filename_component(Python_Dir "${PYTHON_EXECUTABLE}" DIRECTORY)

find_file(Python_Header "python.h" HINTS  ${Python_Dir})
get_filename_component(PYTHON_INCLUDE_DIR ${Python_Header} DIRECTORY CACHE)
include_directories(${PYTHON_INCLUDE_DIR})
unset(Python_Header CACHE)

if(MSVC)
	find_library(PYTHON_LIBRARY NAMES python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} HINTS ${Python_Dir}/libs)
else()
	find_library(PYTHON_LIBRARY NAMES libpython${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} HINTS ${Python_Dir}/libs)
endif()
get_filename_component(PYTHON_LIBRARY_DIR ${PYTHON_LIBRARY} DIRECTORY CACHE)

# find numpy
execute_process(COMMAND ${PYTHON_EXECUTABLE} "-c" "import numpy;print(numpy.get_include(),end='')" OUTPUT_VARIABLE Python_Output ERROR_VARIABLE Python_Error)
STRING(REGEX REPLACE "\\\\" "/" Python_Output ${Python_Output}) # handle windows \ -> /
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
set_property(TARGET crystallography PROPERTY SUFFIX ".pyd") # name crystallography.pyd instead of crystallography.(dll/dylib)