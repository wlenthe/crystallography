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

	@staticmethod
	def zeros_like(array, resolution = (1.0, 1.0), origin = (0.0, 0.0)):
		scan = OimScan()
		scan.rows = array.shape[0]
		scan.colsOdd = array.shape[1]
		scan.colsEven = array.shape[1]
		scan.gridType = Grid.Square
		scan.xStep = resolution[0]
		scan.yStep = resolution[1]
		scan.allocate()
		scan.ci.fill(-1)
		singleRow = numpy.linspace(0, scan.xStep * (scan.colsOdd-1), scan.colsOdd) + origin[0]
		for j in range(scan.rows):
		  scan.x[j] = singleRow
		singleCol = numpy.linspace(0, scan.yStep * (scan.rows-1), scan.rows) + origin[1]
		for i in range(scan.colsOdd):
		  scan.y[:,i] = singleCol
		return scan

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
			# compose header
			header = '# TEM_PIXperUM 1.000000\n'
			header += ('# TEM_PIXperUM 1.000000\n')
			header += ('# x-star ' + str(self.xstar) + '\n')
			header += ('# y-star ' + str(self.ystar) + '\n')
			header += ('# z-star ' + str(self.zstar) + '\n')
			header += ('# WorkingDistance ' + str(self.workingDistance) + '\n')
			header += ('#\n')
			for phase in self.phaseList:
				header += ('# Phase ' + str(phase.number) + '\n')
				header += ('# MaterialName ' + phase.materialName + '\n')
				header += ('# Formula ' + phase.formula + '\n')
				header += ('# Info ' + phase.info + '\n')
				header += ('# Symmetry ' + str(phase.symmetry) + '\n')
				constants = str(phase.latticeConstants[0]) + ' ' + str(phase.latticeConstants[1]) + ' ' + str(phase.latticeConstants[2]) + ' ' + str(phase.latticeConstants[3]) + ' ' + str(phase.latticeConstants[4]) + ' ' + str(phase.latticeConstants[5])
				header += ('# LatticeConstants ' + constants + '\n')
				header += ('# NumberFamilies ' + str(len(phase.hklFamilies)) + '\n')
				for family in phase.hklFamilies:
					header += ('# hklFamilies ' + str(family.hkl[0]) + ' ' + str(family.hkl[1]) + ' ' + str(family.hkl[2]) + ' ')
					header += (str(family.useInIndexing) + ' ' + str(family.diffractionIntensity) + ' ' + str(family.showBands) + '\n')
				for line in phase.elasticConstants:
					constants = str(line[0]) + ' ' + str(line[1]) + ' ' + str(line[2]) + ' ' + str(line[3]) + ' ' + str(line[4]) + ' ' + str(line[5])
					header += ('# ElasticConstants ' + constants + '\n')
				header += ('# Categories')
				for category in phase.categories:
					header += (' ' + str(category))
				header += ('\n#\n')
			header += ('# GRID: ' + self.gridType.value + '\n')
			header += ('# XSTEP: ' + str(self.xStep) + '\n')
			header += ('# YSTEP: ' + str(self.yStep) + '\n')
			header += ('# NCOLS_ODD: ' + str(self.colsOdd) + '\n')
			header += ('# NCOLS_EVEN: ' + str(self.colsEven) + '\n')
			header += ('# NROWS: ' + str(self.rows) + '\n')
			header += ('#\n')
			header += ('# OPERATOR: ' + self.operator + '\n')
			header += ('# SAMPLEID: ' + self.sampleId + '\n')
			header += ('# SCANID: ' + self.scanId + '\n')
			header += ('#\n')
			file.write(header)

			# write body
			if self.sem is not None and self.fit is not None:
				stackedData = numpy.stack((self.euler[:,:,0], self.euler[:,:,1], self.euler[:,:,2], self.x, self.y, self.iq, self.ci, self.phase, self.sem, self.fit))
			elif self.sem is not None:
				stackedData = numpy.stack((self.euler[:,:,0], self.euler[:,:,1], self.euler[:,:,2], self.x, self.y, self.iq, self.ci, self.phase, self.sem))
			elif self.fit is not None:
				stackedData = numpy.stack((self.euler[:,:,0], self.euler[:,:,1], self.euler[:,:,2], self.x, self.y, self.iq, self.ci, self.phase, self.fit))
			else:
				stackedData = numpy.stack((self.euler[:,:,0], self.euler[:,:,1], self.euler[:,:,2], self.x, self.y, self.iq, self.ci, self.phase))
			
			if self.gridType is Grid.Square:
				stackedData = numpy.reshape(stackedData, (stackedData.shape[0], stackedData.shape[1] * stackedData.shape[2])).T
				file.write('\n'.join(map(lambda line: ' '.join(map(str, line)), stackedData)))

			else:
				oddRow = True
				for j in range(self.rows):
					for i in range(self.colsOdd if oddRow else self.colsEven):
						tokens = stackedData[:,j,i]
						file.write(' '.join([str(item) for item in stackedData[:,j,i]]))
						file.write('\n')
					oddRow = not oddRow

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
			if self.gridType is Grid.Square:
				points = 0
				cols = max(self.colsOdd, self.colsEven)
				data = numpy.zeros((cols * self.rows, len(tokens)))
				for line in file:
					data[points] = numpy.fromstring(line, sep = ' ')
					points += 1
				self.euler = numpy.zeros((self.rows, cols, 3))
				self.euler[:,:,0] = numpy.reshape(data[:,0], (self.rows, cols))
				self.euler[:,:,1] = numpy.reshape(data[:,1], (self.rows, cols))
				self.euler[:,:,2] = numpy.reshape(data[:,2], (self.rows, cols))
				self.x = numpy.reshape(data[:,3], (self.rows, cols))
				self.y = numpy.reshape(data[:,4], (self.rows, cols))
				self.iq = numpy.reshape(data[:,5], (self.rows, cols))
				self.ci = numpy.reshape(data[:,6], (self.rows, cols))
				self.phase = numpy.reshape(data[:,7], (self.rows, cols)).astype('int')
				if len(tokens) > 12:
					self.sem = numpy.reshape(data[:,8], (self.rows, cols))
				else:
					self.sem = None
				if len(tokens) > 13:	
					self.fit = numpy.reshape(data[:,9], (self.rows, cols))
				else:
					self.fit = None
				complete = points == cols * self.rows

			else:
				row = 0
				col = 0
				oddRow = True
				end = self.colsOdd
				readSem = len(tokens) > 12
				readFit = len(tokens) > 13
				self.allocate(len(tokens) > 12, len(tokens) > 13)
				for line in file:
					tokens = numpy.fromstring(line, sep = ' ')
					self.euler[row, col, 0] = tokens[0]
					self.euler[row, col, 1] = tokens[1]
					self.euler[row, col, 2] = tokens[2]
					self.x[row, col] = tokens[3]
					self.y[row, col] = tokens[4]
					self.iq[row, col] = tokens[5]
					self.ci[row, col] = tokens[6]
					self.phase[row, col] = tokens[7]
					if readSem:
						self.sem[row, col] = tokens[8]
					if readFit:
						self.fit[row, col] = tokens[9]
					col += 1
					if col == end:
						col = 0
						row += 1
						oddRow = not oddRow
						end = self.colsOdd if oddRow else self.colsEven
				complete = row == self.rows

			#check for incomplete file
			if not complete:
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