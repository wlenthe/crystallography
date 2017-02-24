# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                                 #
# Copyright (c) 2016 William Lenthe                                               #
# All rights reserved.                                                            #
#                                                                                 #
# Redistribution and use in source and binary forms, with or without              #
# modification, are permitted provided that the following conditions are met:     #
#     * Redistributions of source code must retain the above copyright            #
#       notice, this list of conditions and the following disclaimer.             #
#     * Redistributions in binary form must reproduce the above copyright         #
#       notice, this list of conditions and the following disclaimer in the       #
#       documentation and/or other materials provided with the distribution.      #
#     * Neither the name of the <organization> nor the                            #
#       names of its contributors may be used to endorse or promote products      #
#       derived from this software without specific prior written permission.     #
#                                                                                 #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND #
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          #
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY              #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      #
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    #
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     #
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      #
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   #
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    #
#                                                                                 #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import numpy as np
from tsl import OimScan
import h5py

sliceThickness = 0.75
firstSlice = 1
lastSlice = 10
sliceName = 'slice_%03i.hd5'
outputName = 'stack.h5ebsd'

# create h5ebsd file and write header
f = h5py.File(outputName, 'w')
f.create_dataset('EulerTransformationAngle', (1,), dtype='float32')[0] = 90
f.create_dataset('EulerTransformationAxis', (3,), dtype='float32')[:] = [0,0,1]
f.create_dataset('Index', (lastSlice-firstSlice+1,), dtype='int32')[:] = np.linspace(firstSlice, lastSlice, lastSlice-firstSlice+1)
f.create_dataset('Manufacturer', (1,), dtype='S4')[0] = b'TSL'
scan = OimScan(sliceName % lastSlice)
f.create_dataset('Max X Points', (1,), dtype='int64')[0] = scan.colsOdd
f.create_dataset('Max Y Points', (1,), dtype='int64')[0] = scan.rows
f.create_dataset('SampleTransformationAngle', (1,), dtype='float32')[0] = 180
f.create_dataset('SampleTransformationAxis', (3,), dtype='float32')[:] = [0,1,0]
stackOrder = f.create_dataset('Stacking Order', (1,), dtype='uint32')
stackOrder[0] = 1
stackOrder.attrs.create('Name', 'High To Low', dtype='S12')
f.create_dataset('X Resolution', (1,), dtype='float32')[0] = scan.xStep
f.create_dataset('Y Resolution', (1,), dtype='float32')[0] = scan.yStep
f.create_dataset('Z Resolution', (1,), dtype='float32')[0] = sliceThickness
f.create_dataset('ZEndIndex', (1,), dtype='int64')[0] = lastSlice
f.create_dataset('ZStartIndex', (1,), dtype='int64')[0] = firstSlice
f.attrs.create('FileVersion', np.array([5]), dtype = 'int32')

# read slices and add to h5ebsd file
for i in range(firstSlice, lastSlice+1):
	print('adding slice %i' % i)
	scan = OimScan(sliceName % i).writeToH5(f.create_group(str(i)))