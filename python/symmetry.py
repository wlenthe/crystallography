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

import enum, math, sys, itertools
import rotations
from quaternion import Quaternion

class Symmetry(enum.Enum):
	Cubic = 'O - 432 - m3m'
	CubicLow = 'T - 23 - m3'
	Hexagonal = 'D6 - 622 - 6/mmm'
	HexagonalLow = 'C6 - 6 - 6/m'
	Tetragonal = 'D4 - 422 - 4/mmm'
	TetragonalLow = 'C4 - 4 - 4/m'
	Orthorhombic = 'D2 - 222 - mmm'
	Trigonal = 'D3 - 32 - bar3m'
	TrigonalLow = 'C3 - 3 - bar3'
	Monoclinic = 'C2 - 2 2/m'
	Triclinic = 'C1 - 1 - bar1'

	# conversion from tsl symmetry numbers
	@staticmethod
	def FromTSL(number):
		return {
			43: Symmetry.Cubic,
			23: Symmetry.CubicLow,
			62: Symmetry.Hexagonal,
			6 : Symmetry.HexagonalLow,
			42: Symmetry.Tetragonal,
			4 : Symmetry.TetragonalLow,
			22: Symmetry.Orthorhombic,
			32: Symmetry.Trigonal,
			3 : Symmetry.TrigonalLow,
			2 : Symmetry.Monoclinic,
			20: Symmetry.Monoclinic,
			21: Symmetry.Monoclinic,
			1 : Symmetry.Triclinic
		}[number]

	# conversion to tsl symmetry numbers
	def toTSL(self):
		return {
			Symmetry.Cubic:         43,
			Symmetry.CubicLow:      23,
			Symmetry.Hexagonal:     62,
			Symmetry.HexagonalLow:   6,
			Symmetry.Tetragonal:    42,
			Symmetry.TetragonalLow:  4,
			Symmetry.Orthorhombic:  22,
			Symmetry.Trigonal:      32,
			Symmetry.TrigonalLow:    3,
			Symmetry.Monoclinic:     2,
			Symmetry.Monoclinic:    20,
			Symmetry.Monoclinic:    21,
			Symmetry.Triclinic:      1
		}[self]

	# symmetry operators as Quaternions
	def quOperators(self):
		a = 0.5
		b = 1.0 / math.sqrt(2.0)
		c = math.sqrt(3.0) / 2.0
		return {
			Symmetry.Cubic:         [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 1, 0, 0]),
			                         Quaternion([0, 0, 1, 0]),
			                         Quaternion([0, 0, 0, 1]),
			                         Quaternion([a, a, a, a]),
			                         Quaternion([a,-a, a, a]),
			                         Quaternion([a, a,-a, a]),
			                         Quaternion([a, a, a,-a]),
			                         Quaternion([a, a,-a,-a]),
			                         Quaternion([a,-a, a,-a]),
			                         Quaternion([a,-a,-a, a]),
			                         Quaternion([a,-a,-a,-a]),
			                         Quaternion([b, b, 0, 0]),
			                         Quaternion([b, 0, b, 0]),
			                         Quaternion([b, 0, 0, b]),
			                         Quaternion([0, b, b, 0]),
			                         Quaternion([0, b, 0, b]),
			                         Quaternion([0, 0, b, b]),
			                         Quaternion([b,-b, 0, 0]),
			                         Quaternion([b, 0,-b, 0]),
			                         Quaternion([b, 0, 0,-b]),
			                         Quaternion([0,-b, b, 0]),
			                         Quaternion([0,-b, 0, b]),
			                         Quaternion([0, 0,-b, b])],

			Symmetry.CubicLow:      [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 1, 0, 0]),
			                         Quaternion([0, 0, 1, 0]),
			                         Quaternion([0, 0, 0, 1]),
			                         Quaternion([a, a, a, a]),
			                         Quaternion([a,-a, a, a]),
			                         Quaternion([a, a,-a, a]),
			                         Quaternion([a, a, a,-a]),
			                         Quaternion([a, a,-a,-a]),
			                         Quaternion([a,-a, a,-a]),
			                         Quaternion([a,-a,-a, a]),
			                         Quaternion([a,-a,-a,-a])],

			Symmetry.Hexagonal:     [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 1, 0, 0]),
			                         Quaternion([0, 0, 1, 0]),
			                         Quaternion([0, 0, 0, 1]),
			                         Quaternion([a, 0, 0, c]),
			                         Quaternion([0, a, c, 0]),
			                         Quaternion([0, c, a, 0]),
			                         Quaternion([c, 0, 0, a]),
			                         Quaternion([a, 0, 0,-c]),
			                         Quaternion([0,-a, c, 0]),
			                         Quaternion([0,-c, a, 0]),
			                         Quaternion([c, 0, 0,-a])],

			Symmetry.HexagonalLow:  [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 0, 0, 1]),
			                         Quaternion([a, 0, 0, c]),
			                         Quaternion([c, 0, 0, a]),
			                         Quaternion([a, 0, 0,-c]),
			                         Quaternion([c, 0, 0,-a])],

			Symmetry.Tetragonal:    [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 1, 0, 0]),
			                         Quaternion([0, 0, 1, 0]),
			                         Quaternion([0, 0, 0, 1]),
			                         Quaternion([b, 0, 0, b]),
			                         Quaternion([0, b, b, 0]),
			                         Quaternion([b, 0, 0,-b]),
			                         Quaternion([0,-b, b, 0])],

			Symmetry.TetragonalLow: [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 0, 0, 1]),
			                         Quaternion([b, 0, 0, b]),
			                         Quaternion([b, 0, 0,-b])],

			Symmetry.Orthorhombic:  [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 1, 0, 0]),
			                         Quaternion([0, 0, 1, 0]),
			                         Quaternion([0, 0, 0, 1])],

			Symmetry.Trigonal:      [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 0, 0, 1]),
			                         Quaternion([a, 0, 0, c]),
			                         Quaternion([c, 0, 0, a]),
			                         Quaternion([a, 0, 0,-c]),
			                         Quaternion([c, 0, 0,-a])],

			Symmetry.TrigonalLow:   [Quaternion([1, 0, 0, 0]),
			                         Quaternion([a, 0, 0, c]),
			                         Quaternion([a, 0, 0,-c])],

			Symmetry.Monoclinic:    [Quaternion([1, 0, 0, 0]),
			                         Quaternion([0, 0, 0, 1])],

			Symmetry.Triclinic:     [Quaternion([1, 0, 0, 0])]
		}[self]

	@staticmethod
	def fzCyclic(ro, n):
		# top and bottom face at +/-tan(pi/2n)
		return abs(ro[0]) <= math.tan(math.pi / (2 * n))

	@staticmethod
	def fzDihedral(ro, n):
		# top and bottom face at +/-tan(pi/2n)
		t = math.tan(math.pi / (2 * n))
		if abs(ro[2]) > t:
			return False

		# 2n faces distance 1 from origin
		# y <= ((2+sqrt(2))*t - (1+sqrt(2))) * x + (1+sqrt(2))*(1-t)
		y, x = sorted([abs(ro[0]), abs(ro[1])])
		if x > 1:
			return False

		return {
			2: True,
			3: y / (1+math.sqrt(2)) + (1-math.sqrt(2/3)) * x < 1 - 1/math.sqrt(3),
			4: y + x < math.sqrt(2),
			6: y / (1+math.sqrt(2)) + (1-2*math.sqrt(2)+math.sqrt(6)) * x < math.sqrt(3)-1
		}[n]

	@staticmethod
	def fzT23(ro):
		return sum([abs(x) for x in ro]) <= 1.0

	@staticmethod
	def fzO432(ro):
		return Symmetry.fzT23(ro) and max([abs(x) for x in ro]) <= math.sqrt(2) - 1

	# check if a rodriguez vector is in the fundamental zone
	def roInFZ(self, ro):
		if len(ro) is 4:
			ro = [x * ro[3] for x in ro[:3]]
		return {
			Symmetry.Cubic:         Symmetry.fzO432(ro),
			Symmetry.CubicLow:      Symmetry.fzT23(ro),

			Symmetry.Hexagonal:     Symmetry.fzDihedral(ro, 6),
			Symmetry.Tetragonal:    Symmetry.fzDihedral(ro, 4),
			Symmetry.Trigonal:      Symmetry.fzDihedral(ro, 3),
			Symmetry.Orthorhombic:  Symmetry.fzDihedral(ro, 2),
			
			Symmetry.HexagonalLow:  Symmetry.fzCyclic(ro, 6),
			Symmetry.TetragonalLow: Symmetry.fzCyclic(ro, 4),
			Symmetry.TrigonalLow:   Symmetry.fzCyclic(ro, 3),
			Symmetry.Monoclinic:    Symmetry.fzCyclic(ro, 2),
			Symmetry.Triclinic:     True
		}[self]

	# check if a quaternion is in the fundamental zone
	def quInFZ(self, qu):
		return roInFZ(rotations.qu2ro(qu))

	# get the equivalent rotation in the fundamental zone
	def fzQu(self, qu, convention = rotations.convention):
		q = Quaternion(qu)
		assert(isinstance(convention, rotations.Convention))
		if convention is rotations.Convention.passive:
			qFZ = min([(sym * q).supplement() for sym in self.quOperators()])
		elif convention is rotations.Convention.active:
			qFZ = min([(q * sym).supplement() for sym in self.quOperators()])
		return qFZ if isinstance(qu, Quaternion) else qFZ.wxyz

	# computes disorientation between passive rotations qu1 and qu2 (active is from qu1 * sym * qu2.conjugate())
	def disoQuat(self, qu1, qu2, convention = rotations.convention):
		if self is not Symmetry.Cubic:
			if convention is rotations.Convention.passive:
				deltaQ = Quaternion(qu1) * Quaternion(qu2).conjugate()
				equivMisos = sorted([(q * j.conjugate()).supplement() for q, j in itertools.product([i * deltaQ for i in self.quOperators()], self.quOperators())]) # compute all equivilent misorientation representations and sort by misorientation angle
			elif convention is rotations.Convention.active:
				q1 = Quaternion(qu1)
				q2 = Quaternion(qu2)
				equivMisos = sorted([(q1 * i * q2).supplement() for i in self.quOperators()]) # compute all equivilent misorientation representations and sort by misorientation angle
			equivCount = 1
			while math.isclose(equivMisos[equivCount].w(), equivMisos[0].w(), abs_tol = 28 * sys.float_info.epsilon):
				equivCount += 1
			equivMisos = equivMisos[:equivCount] # keep all rotations with the same angle
			diso = sorted(equivMisos, key=lambda q: (q.z(), q.y(), q.x()), reverse = True)[0] # sort by z, then y, then x
		else:
			# cubic shortcut
			deltaQ = Quaternion(sorted(abs(Quaternion(qu1) * Quaternion(qu2).conjugate())))
			opertor4 = sum(deltaQ) / 2
			opertor17 = (deltaQ[2] + deltaQ[3]) / math.sqrt(2)
			if opertor4 > deltaQ[3] and opertor4 > opertor17:
				qu4 = Quaternion([0.5, 0.5, 0.5, 0.5])
				deltaQ = sorted(abs(deltaQ / qu4))
			elif opertor17 > deltaQ[3]:
				b = 1 / math.sqrt(2)
				qu17 = Quaternion([0, 0, b, b])
				deltaQ = sorted(abs(deltaQ / qu17))
			w = deltaQ[3]
			x = deltaQ[0]
			y = deltaQ[1]
			z = deltaQ[2]
			diso = Quaternion([w, x, y, z])
		return diso if isinstance(qu1, Quaternion) else diso.wxyz