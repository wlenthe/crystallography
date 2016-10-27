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

import re, enum, os
import numpy, h5py

class HKLFamily:
	def __init__(self):
		self.hkl = [0] * 3
		self.useInIndexing = 0
		self.diffractionIntensity = 0.0
		self.showBands = 0

class OimPhase:
	def __init__(self, number):
		self.number = number
		self.materialName = ''
		self.formula = ''
		self.info = ''
		self.symmetry = 0
		self.latticeConstants =[0] * 6 # a, b, c, alpha, beta, gamma (degrees)
		self.hklFamilies = []
		self.elasticConstants = []
		self.categories = []

class Grid(enum.Enum):
	Square = 'SqrGrid'
	Hexagonal = 'HexGrid'

class OimScan:
	def __init__(self, filename = None):
		self.clearHeader()
		self.allocate()
		if filename is not None:
			file, extension = os.path.splitext(filename)
			extension = extension.lower()
			if extension == '.ang':
				self.readAng(filename)
			if extension == '.h5' or extension == '.hdf' or extension == '.hdf5':
				self.readH5(filename)

	def clearHeader(self):
		self.xstar = 0
		self.ystar = 0
		self.zstar = 0
		self.workingDistance = 0
		self.phaseList = []
		self.gridType = Grid.Square
		self.colsOdd = 0
		self.colsEven = 0
		self.rows = 0
		self.xStep = 0.0
		self.yStep = 0.0
		self.operator = ''
		self.sampleId = ''
		self.scanId = ''

	def allocate(self, sem = False, fit = False):
		maxCols = max(self.colsOdd, self.colsEven)
		self.euler = numpy.zeros((self.rows, maxCols, 3))
		self.x = numpy.zeros((self.rows, maxCols))
		self.y = numpy.zeros((self.rows, maxCols))
		self.iq = numpy.zeros((self.rows, maxCols))
		self.ci = numpy.zeros((self.rows, maxCols))
		self.phase = numpy.zeros((self.rows, maxCols), dtype = 'int')
		self.sem = numpy.zeros((self.rows, maxCols)) if sem else None
		self.fit = numpy.zeros((self.rows, maxCols)) if fit else None

	def writeAng(self, filename):
		with open(filename, 'w') as file:
			# write header
			file.write('# TEM_PIXperUM 1.000000\n')
			file.write('# x-star ' + str(self.xstar) + '\n')
			file.write('# y-star ' + str(self.ystar) + '\n')
			file.write('# z-star ' + str(self.zstar) + '\n')
			file.write('# WorkingDistance ' + str(self.workingDistance) + '\n')
			file.write('#\n')
			for phase in self.phaseList:
				file.write('# Phase ' + str(phase.number) + '\n')
				file.write('# MaterialName ' + phase.materialName + '\n')
				file.write('# Formula ' + phase.formula + '\n')
				file.write('# Info ' + phase.info + '\n')
				file.write('# Symmetry ' + str(phase.symmetry) + '\n')
				constants = str(phase.latticeConstants[0]) + ' ' + str(phase.latticeConstants[1]) + ' ' + str(phase.latticeConstants[2]) + ' ' + str(phase.latticeConstants[3]) + ' ' + str(phase.latticeConstants[4]) + ' ' + str(phase.latticeConstants[5])
				file.write('# LatticeConstants ' + constants + '\n')
				file.write('# NumberFamilies ' + str(len(phase.hklFamilies)) + '\n')
				for family in phase.hklFamilies:
					file.write('# hklFamilies ' + str(family.hkl[0]) + ' ' + str(family.hkl[1]) + ' ' + str(family.hkl[2]) + ' ')
					file.write(str(family.useInIndexing) + ' ' + str(family.diffractionIntensity) + ' ' + str(family.showBands) + '\n')
				for line in phase.elasticConstants:
					constants = str(line[0]) + ' ' + str(line[1]) + ' ' + str(line[2]) + ' ' + str(line[3]) + ' ' + str(line[4]) + ' ' + str(line[5])
					file.write('# ElasticConstants ' + constants + '\n')
				file.write('# Categories')
				for category in phase.categories:
					file.write(' ' + str(category))
				file.write('\n#\n')
			file.write('# GRID: ' + self.gridType.value + '\n')
			file.write('# XSTEP: ' + str(self.xStep) + '\n')
			file.write('# YSTEP: ' + str(self.yStep) + '\n')
			file.write('# NCOLS_ODD: ' + str(self.colsOdd) + '\n')
			file.write('# NCOLS_EVEN: ' + str(self.colsEven) + '\n')
			file.write('# NROWS: ' + str(self.rows) + '\n')
			file.write('#\n')
			file.write('# OPERATOR: ' + self.operator + '\n')
			file.write('# SAMPLEID: ' + self.sampleId + '\n')
			file.write('# SCANID: ' + self.scanId + '\n')
			file.write('#\n')

			# write body
			oddRow = True
			for j in range(self.rows):
				for i in range(self.colsOdd if oddRow else self.colsEven):
					file.write(str(self.euler[j,i,0]))
					file.write(' ')
					file.write(str(self.euler[j,i,1]))
					file.write(' ')
					file.write(str(self.euler[j,i,2]))
					file.write(' ')
					file.write(str(self.x[j,i]))
					file.write(' ')
					file.write(str(self.y[j,i]))
					file.write(' ')
					file.write(str(self.iq[j,i]))
					file.write(' ')
					file.write(str(self.ci[j,i]))
					file.write(' ')
					file.write(str(self.phase[j,i]))
					if self.sem is not None:
						file.write(' ')
						file.write(str(self.sem[j,i]))
					if self.fit is not None:
						file.write(' ')
						file.write(str(self.fit[j,i]))
					file.write('\n')

	def readAng(self, filename):
		with open(filename, newline = '') as file:
			# parse header
			header = True
			while header:
				position = file.tell()
				line = file.readline()
				if not line:
					raise EOFError('ang file has no body')
				tokens = re.split('\s+', line.strip())
				header = tokens[0] == '#'
				if not header:
					file.seek(position)
				elif len(tokens) > 2:# '#', 'key', 'value'
						if tokens[1] == 'TEM_PIXperUM':
							self.pixPerUm = float(tokens[2])
						elif tokens[1] == 'x-star':
							self.xstar = float(tokens[2])
						elif tokens[1] == 'y-star':
							self.ystar = float(tokens[2])
						elif tokens[1] == 'z-star':
							self.zstar = float(tokens[2])
						elif tokens[1] == 'WorkingDistance':
							self.workingDistance = float(tokens[2])
						elif tokens[1] == 'Phase':
							self.phaseList.append(OimPhase(int(tokens[2])))
						elif tokens[1] == 'MaterialName':
							self.phaseList[-1].materialName = tokens[2]
						elif tokens[1] == 'Formula':
							self.phaseList[-1].formula = tokens[2]
						elif tokens[1] == 'Info':
							self.phaseList[-1].info = tokens[2]
						elif tokens[1] == 'Symmetry':
							self.phaseList[-1].symmetry = int(tokens[2])
						elif tokens[1] == 'NumberFamilies':
							pass
						elif tokens[1] == 'LatticeConstants':
							values = [0] * 6
							values[0] = float(tokens[2])
							values[1] = float(tokens[3])
							values[2] = float(tokens[4])
							values[3] = float(tokens[5])
							values[4] = float(tokens[6])
							values[5] = float(tokens[7])
							self.phaseList[-1].latticeConstants = values
						elif tokens[1] == 'hklFamilies':
							family = HKLFamily()
							family.hkl[0] = tokens[2]
							family.hkl[1] = tokens[3]
							family.hkl[2] = tokens[4]
							family.useInIndexing = int(tokens[5])
							family.diffractionIntensity = float(tokens[6])
							family.showBands = int(tokens[7])
							self.phaseList[-1].hklFamilies.append(family)
						elif tokens[1] == 'ElasticConstants':
							values = [0] * 6
							values[0] = float(tokens[2])
							values[1] = float(tokens[3])
							values[2] = float(tokens[4])
							values[3] = float(tokens[5])
							values[4] = float(tokens[6])
							values[5] = float(tokens[7])
							self.phaseList[-1].elasticConstants.append(values)
						elif tokens[1].startswith('Categories'): # tsl doesn't print space between categories and first number
							firstCategory = tokens[1][10:]
							values = [int(firstCategory)] if len(firstCategory) > 0 else []
							for i in range(len(tokens)-2):
								values.append(int(tokens[i+2]))
							self.phaseList[-1].categories = values
						elif tokens[1] == 'GRID:':
							self.gridType = Grid(tokens[2])
						elif tokens[1] == 'XSTEP:':
							self.xStep = float(tokens[2])
						elif tokens[1] == 'YSTEP:':
							self.yStep = float(tokens[2])
						elif tokens[1] == 'NCOLS_ODD:':
							self.colsOdd = int(tokens[2])
						elif tokens[1] == 'NCOLS_EVEN:':
							self.colsEven = int(tokens[2])
						elif tokens[1] == 'NROWS:':
							self.rows = int(tokens[2])
						elif tokens[1] == 'OPERATOR:':
							self.operator = tokens[2]
						elif tokens[1] == 'SAMPLEID:':
							self.sampleId = tokens[2]
						elif tokens[1] == 'SCANID:':
							self.scanId = tokens[2]

			# parse body lines: eu1, eu2, eu3, x, y, iq, ci, phase, sem, fit
			row = 0
			col = 0
			oddRow = True
			self.allocate(len(tokens) > 12, len(tokens) > 13)
			for line in file:
				tokens = re.split('\s+', line.strip())
				self.euler[row, col, 0] = float(tokens[0])
				self.euler[row, col, 1] = float(tokens[1])
				self.euler[row, col, 2] = float(tokens[2])
				self.x[row, col] = float(tokens[3])
				self.y[row, col] = float(tokens[4])
				self.iq[row, col] = float(tokens[5])
				self.ci[row, col] = float(tokens[6])
				self.phase[row, col] = int(tokens[7])
				if self.sem is not None:
					self.sem[row, col] = float(tokens[8])
				if self.fit is not None:
					self.fit[row, col] = float(tokens[9])
				col += 1
				if col == (self.colsOdd if oddRow else self.colsEven):
					col = 0
					row += 1
					oddRow = not oddRow

			#check for incomplete file
			if row != self.rows:
					raise EOFError('ang has incomplete body')

	def readH5(self, filename):
		with h5py.File(filename, 'r') as file:
			# get all scans in file
			scans = []
			for key in file.keys():
				if key == ' Manufacturer':
					pass
				elif key == ' Version':
					pass
				else:
					scans.append(key)

			# read first scan
			scan = file[scans[0]]['EBSD']
			header = scan['Header']
			body = scan['Data']

			# read header
			self.xstar = header['Pattern Center Calibration']['x-star'][0]
			self.ystar = header['Pattern Center Calibration']['y-star'][0]
			self.zstar = header['Pattern Center Calibration']['z-star'][0]
			self.workingDistance = header['Camera Elevation Angle'][0]
			self.gridType = Grid(header['Grid Type'][0].decode('utf-8'))
			for key in header['Phase'].keys():
				phase = OimPhase(int(key))
				phase.materialName = header['Phase'][key]['MaterialName'][0].decode('utf-8')
				phase.formula = header['Phase'][key]['Formula'][0].decode('utf-8')
				phase.info = header['Phase'][key]['Info'][0].decode('utf-8')
				phase.symmetry = header['Phase'][key]['Symmetry'][0]
				phase.latticeConstants[0] = header['Phase'][key]['Lattice Constant a'][0]
				phase.latticeConstants[1] = header['Phase'][key]['Lattice Constant b'][0]
				phase.latticeConstants[2] = header['Phase'][key]['Lattice Constant c'][0]
				phase.latticeConstants[3] = header['Phase'][key]['Lattice Constant alpha'][0]
				phase.latticeConstants[4] = header['Phase'][key]['Lattice Constant beta'][0]
				phase.latticeConstants[5] = header['Phase'][key]['Lattice Constant gamma'][0]
				for row in header['Phase'][key]['hkl Families']:
					family = HKLFamily()
					family.hkl = [row[0], row[1], row[2]]
					family.useInIndexing = row[4]
					family.diffractionIntensity = row[3]
					family.showBands = row[5]
					phase.hklFamilies.append(family)
				phase.elasticConstants = [[0] * 6 for i in range(6)]
				phase.categories = [0] * 5
				self.phaseList.append(phase)
			self.colsOdd = header['nColumns'][0]
			self.colsEven = self.colsOdd - 1
			self.rows = header['nRows'][0]
			self.xStep = header['Step X'][0]
			self.yStep = header['Step Y'][0]
			self.operator = header['Operator'][0].decode('utf-8')
			self.sampleId = header['Sample ID'][0].decode('utf-8')
			self.scanId = header['Scan ID'][0].decode('utf-8')

			# read data
			sem = 'SEM Signal' in body.keys()
			fit = 'Fit' in body.keys()
			self.allocate(sem, fit)
			maxCols = max(self.colsOdd, self.colsEven)
			self.euler[:,:,0] = numpy.reshape(body['Phi1'], (self.rows, maxCols))
			self.euler[:,:,1] = numpy.reshape(body['Phi'], (self.rows, maxCols))
			self.euler[:,:,2] = numpy.reshape(body['Phi2'], (self.rows, maxCols))
			self.x = numpy.reshape(body['X Position'], (self.rows, maxCols))
			self.y = numpy.reshape(body['Y Position'], (self.rows, maxCols))
			self.iq = numpy.reshape(body['IQ'], (self.rows, maxCols))
			self.ci = numpy.reshape(body['CI'], (self.rows, maxCols))
			self.phase = numpy.reshape(body['Phase'], (self.rows, maxCols))
			if sem:
				self.sem = numpy.reshape(body['SEM Signal'], (self.rows, maxCols))
			if fit:
				self.fit = numpy.reshape(body['Fit'], (self.rows, maxCols))