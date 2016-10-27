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

import sys, math

class Quaternion:
	def __init__(self, wxyz = [1, 0, 0, 0]):
		if isinstance(wxyz, Quaternion):
			# construct from quaterion
			self.wxyz = wxyz.wxyz[:]
		else:
			# construct from list
			self.wxyz = wxyz[:4]

	# element accessors
	def w(self):
		return self.wxyz[0]
	def x(self):
		return self.wxyz[1]
	def y(self):
		return self.wxyz[2]
	def z(self):
		return self.wxyz[3]

	# iteration and operator[]
	def __getitem__(self, i):
		return self.wxyz[i]
	def __iter__(self):
		return self.wxyz.__iter__()

	# conversion to string representation
	def __str__(self):
		return str(self.wxyz)
	def __repr__(self):
		return str(self)

	# +-*/ 2 quaternions
	def __add__(self, other):
		return Quaternion([i + j for i, j in zip(self.wxyz, other.wxyz)])

	def __sub__(self, other):
		return Quaternion([i - j for i, j in zip(self.wxyz, other.wxyz)])
		
	def __mul__(self, other):
		w = self.w()*other.w() - self.x()*other.x() - self.y()*other.y() - self.z()*other.z()
		x = self.x()*other.w() + self.w()*other.x() - self.z()*other.y() + self.y()*other.z()
		y = self.y()*other.w() + self.z()*other.x() + self.w()*other.y() - self.x()*other.z()
		z = self.z()*other.w() - self.y()*other.x() + self.x()*other.y() + self.w()*other.z()
		return Quaternion([w, x, y, z])

	def __truediv__(self, other):
		n = other.norm2()
		w = (self.w()*other.w() + self.x()*other.x() + self.y()*other.y() + self.z()*other.z()) / n
		x = (self.x()*other.w() - self.w()*other.x() - self.z()*other.y() + self.y()*other.z()) / n
		y = (self.y()*other.w() + self.z()*other.x() - self.w()*other.y() - self.x()*other.z()) / n
		z = (self.z()*other.w() - self.y()*other.x() + self.x()*other.y() - self.w()*other.z()) / n
		return Quaternion([w, x, y, z])

	# +-*/ 2 quaternion with scalar
	def scalarAdd(self, s):
		return Quaternion([i + s for i in self.wxyz])

	def scalarSub(self, s):
		return Quaternion([i - s for i in self.wxyz])

	def scalarMul(self, s):
		return Quaternion([i * s for i in self.wxyz])

	def scalarDiv(self, s):
		return Quaternion([i / s for i in self.wxyz])

	# normalization
	def norm2(self):
		return sum([i * i for i in self.wxyz])

	def norm(self):
		return math.sqrt(self.norm2())

	def normalize(self):
		n = self.norm()
		return Quaternion([i / n for i in self.wxyz])

	# conjugate, inverse, and negative
	def conjugate(self):
		return Quaternion([self.w(), -self.x(), -self.y(), -self.z()])

	def negate(self):
		return Quaternion([-i for i in self.wxyz])

	def __neg__(self):
		return self.negate()

	def __abs__(self):
		return Quaternion([abs(i) for i in self.wxyz])

	def inverse(self):
		return self.conjugate().scalarDiv(self.norm2())

	# move rotation to 0->180
	def supplement(self):
		return self.negate() if self.w() < 0 else self

	# sort by rotation angle
	def __lt__(self, other):
		return self.w() > other.w()