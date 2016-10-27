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

# orientation transform routines based on
#  -Rowenhorst, David, et al. "Consistent Representations of and Conversions Between 3D Rotations." Model. Simul. Mater. Sci. Eng. 23.8 (2015): 083501.
#  -Roşca, D., et al. "A New Method of Constructing a Grid in the Space of 3D rotations and its Applications to Texture Analysis." Model. Simul. Mater. Sci. Eng. 22.7 (2014): 075013.
#  -fortran implementation of routines by Marc De Graef (https://github.com/marcdegraef/3Drotations)

# the following conventions are used:
#  -quaternions as [w, x, y, z]
#  -rotation angle <= pi
#  -rotation axis in positive z hemisphere for rotations of pi
#  -rotation axis = [0, 0, 1] for rotations of 0

import sys, enum, math

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Constants                                                                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# tolerance on zero 0
epsilon = sys.float_info.epsilon

# Passive / Active convention (passive by default)
class Convention(enum.Enum):
	passive = 1.0 # i*j*k = +1
	active = -1.0 # i*j*k = -1
convention = Convention.passive

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Helper Functions                                                                #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# returns upper hemisphere axis (for rotations of pi)
def orientAxis(n):
	if abs(n[2]) > epsilon:
		if n[2] < 0.0:
			return [-x for x in n] # [_, _, +z]
	elif abs(n[1]) > epsilon:
		if n[1] < 0.0:
			return [-x for x in n[:2]] + [0.0] # [_, +y, 0]
	elif n[0] < 0.0:
		return [-n[0]] + [0.0] * 2 # [+x, 0, 0]
	return [x for x in n]

# side length of cubochoric cube
def cuA():
	return math.pow(math.pi, 2.0 / 3.0)

# radius of homochoric sphere
def hoR():
	return math.pow(0.75 * math.pi, 1.0 / 3.0)

# inverse of ((3/4)*(x-sin(x)))^(2/3) via newton's method
def hoInv(y):
	if y < 0.0 or y > math.pow(0.75 * math.pi, 2.0 / 3.0) + 7.0 * epsilon:
		raise ValueError("domain error: '%f' is outside of [0, (3*pi/4)^(2/3)]" % y)
	x = 2.0 * math.acos(1.0 - y / 2.0) # initial guess from taylor expansion
	y = math.sqrt(y)
	prevErr = float('inf')
	for i in range(7): # converges within 6 calculation for all values tested within domain
		fx = math.pow(0.75 * (x - math.sin(x)), 1.0/3.0)
		delta = fx - y
		err = abs(delta)
		if 0.0 == delta or err == prevErr: # no error or flipping between +/- v
			return x
		x -= 4.0 * fx * fx * delta / (1.0 - math.cos(x))
		if err > prevErr: # flipping between +v / -2v (return )
			return x
		prevErr = err
	raise RuntimeError('failed to invert ((3/4)*(x-sin(x)))^(2/3) for %f' % y*y)

# return sextant type for cubochoric <-> homochoric transformation symmetry
def pyramidType(x):
	vMax = max([abs(i) for i in x])
	if abs(x[2]) == vMax:
		return 0
	if abs(x[0]) == vMax:
		return 1
	if abs(x[1]) == vMax:
		return 2
	raise RuntimeError('failed to find pyramid type for %f' % x)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Direct Conversions                                                              #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#A.1
def eu2om(eu):
	s = [math.sin(x) for x in eu]
	c = [math.cos(x) for x in eu]
	s = [0.0 if abs(x) <= epsilon else x for x in s]
	c = [0.0 if abs(x) <= epsilon else x for x in c]

	om = [[0.0] * 3 for i in range(3)]
	om[0][0] =  c[0] * c[2] - s[0] * c[1] * s[2]
	om[0][1] =  s[0] * c[2] + c[0] * c[1] * s[2]
	om[0][2] =  s[1] * s[2]
	om[1][0] = -c[0] * s[2] - s[0] * c[1] * c[2]
	om[1][1] = -s[0] * s[2] + c[0] * c[1] * c[2]
	om[1][2] =  s[1] * c[2]
	om[2][0] =  s[0] * s[1]
	om[2][1] = -c[0] * s[1]
	om[2][2] =  c[1]
	return om

#A.2
def eu2ax(eu):
	t = math.tan(eu[1] / 2.0)
	sigma = (eu[0] + eu[2]) / 2.0
	tau = math.sqrt(t * t + math.sin(sigma) * math.sin(sigma))
	if abs(tau) <= 2.0 * epsilon:
		return [0.0, 0.0, 1.0, 0.0] # handle 0 rotation
	delta = (eu[0] - eu[2]) / 2.0
	alpha = math.pi if abs(sigma - math.pi / 2.0) <= epsilon else 2.0 * math.atan(tau / math.cos(sigma))
	n = [-convention.value / math.copysign(tau, alpha)] * 3
	n[0] *= t * math.cos(delta)
	n[1] *= t * math.sin(delta)
	n[2] *= math.sin(sigma)

	# normalize
	mag = math.sqrt(math.fsum([x * x for x in n]))
	n = [x / mag for x in n]

	# handle ambiguous case (rotation angle of pi)
	alpha = abs(alpha)
	if alpha + epsilon >= math.pi:
		return orientAxis(n) + [math.pi]
	return n + [alpha]

#A.3
def eu2ro(eu):
	ax = eu2ax(eu)
	if math.pi == ax[3]:
		ax[3] = float('inf')
	elif abs(ax[3]) <= epsilon:
		ax = [0.0, 0.0, 1.0, 0.0]
	else:
		ax[3] = math.tan(ax[3] / 2.0)
	return ax

#A.4
def eu2qu(eu):
	eu = [x / 2.0 for x in eu]
	c = math.cos(eu[1])
	s = math.sin(eu[1])
	sigma = (eu[0] + eu[2])
	delta = (eu[0] - eu[2])
	qu = [-convention.value] * 4
	qu[0] = c * math.cos(sigma)
	qu[1] *= s * math.cos(delta)
	qu[2] *= s * math.sin(delta)
	qu[3] *= c * math.sin(sigma)
	if qu[0] < 0.0:
		qu = [-x for x in qu]

	# normalize
	mag = math.sqrt(math.fsum([x * x for x in qu]))
	qu = [x / mag for x in qu]

	# handle ambiguous case (rotation angle of pi)
	if abs(qu[0]) <= 3.0 * epsilon:
		return [0.0] + orientAxis(qu[1:])
	return qu

#A.5
def om2eu(om):
	eu = [0.0] * 3
	if abs(om[2][2]) >= 1.0 - 3.0 * epsilon:
		if om[2][2] > 0.0:
			eu[0] = math.atan2(om[0][1], om[0][0]) # eu = [_, 0, _]
		else:
			eu[0] = -math.atan2(-om[0][1], om[0][0]) # eu = [_, pi, _]
			eu[1] = math.pi
	else:
		eu[1] = math.acos(om[2][2])
		zeta = 1.0 / math.sqrt(1.0 - om[2][2] * om[2][2])
		eu[0] = math.atan2(om[2][0] * zeta, -om[2][1] * zeta)
		eu[2] = math.atan2(om[0][2] * zeta, om[1][2] * zeta)
	return [x + 2.0 * math.pi if x < 0.0 else x for x in eu]

#A.6
def om2ax(om):
	omega = (om[0][0] + om[1][1] + om[2][2] - 1.0) / 2.0
	if omega <= 2.0 * epsilon - 1.0:
		omega = math.copysign(1.0, omega)
	if 1.0 == omega:
		return [0.0, 0.0, 1.0, 0.0]

	# compute eigenvector for eigenvalue of 1 (cross product of 2 adjacent columns of A-y*I)
	om00 = om[0][0] - 1.0
	om11 = om[1][1] - 1.0
	om22 = om[2][2] - 1.0
	vecs = [[om[1][0]*om[2][1] - om[2][0]*  om11  , om[2][0]*om[0][1] -   om00  *om[2][1],   om00  *  om11   - om[1][0]*om[0][1]],
					[  om11  *  om22   - om[2][1]*om[1][2], om[2][1]*om[0][2] - om[0][1]*  om22  , om[0][1]*om[1][2] -   om11  *om[0][2]],
					[om[1][2]*om[2][0] -   om22  *om[1][0],   om22  *  om00   - om[0][2]*om[2][0], om[0][2]*om[1][0] - om[1][2]*  om00  ]]

	# select vector with largest magnitude
	mags = [math.fsum([x * x for x in v]) for v in vecs]
	mag = max(mags)
	i = mags.index(mag)
	mag = math.sqrt(mag)
	if mag <= epsilon:
		return [0.0, 0.0, 1.0, 0.0]
	n = [x / mag for x in vecs[i]]

	# check ambiguous case
	if omega <= epsilon - 1.0:
		return orientAxis(n) + [math.pi]
	
	# check axis sign
	n[0] = math.copysign(n[0], convention.value * (om[2][1] - om[1][2]))
	n[1] = math.copysign(n[1], convention.value * (om[0][2] - om[2][0]))
	n[2] = math.copysign(n[2], convention.value * (om[1][0] - om[0][1]))
	return n + [math.acos(omega)]

#A.7
def om2qu(om):
	qu = [1.0 + om[0][0] + om[1][1] + om[2][2],
	      1.0 + om[0][0] - om[1][1] - om[2][2],
	      1.0 - om[0][0] + om[1][1] - om[2][2],
	      1.0 - om[0][0] - om[1][1] + om[2][2]]

	# handle ambiguous case (rotation of pi)
	if abs(qu[0]) <= 3.0 * epsilon:
		return ax2qu(om2ax(om))

	# handle rotation of 0
	if qu[0] <= epsilon - 2.0:
		return [1.0, 0.0, 0.0, 0.0]

	qu = [0.0 if x < 0.0 else convention.value * math.sqrt(x) / 2.0 for x in qu]
	if convention.value * om[1][2] > convention.value * om[2][1]:
		qu[1] = -qu[1]
	if convention.value * om[2][0] > convention.value * om[0][2]:
		qu[2] = -qu[2]
	if convention.value * om[0][1] > convention.value * om[1][0]:
		qu[3] = -qu[3]

	# ensure rotation angle <= pi
	if qu[0] < 0.0:
		qu = [-x for x in qu]

	# normalize and return
	mag = math.sqrt(math.fsum([x * x for x in qu]))
	return [x / mag for x in qu]

#A.8
def ax2om(ax):
	c = math.cos(ax[3])
	s = math.sin(ax[3])
	omc = 1.0 - c
	om = [[0.0] * 3 for i in range(3)]
	om[0][0] = c + omc * ax[0] * ax[0]
	om[1][1] = c + omc * ax[1] * ax[1]
	om[2][2] = c + omc * ax[2] * ax[2]
	x = omc * ax[0] * ax[1]
	y = convention.value * s * ax[2]
	om[1][0] = x + y
	om[0][1] = x - y
	x = omc * ax[1] * ax[2]
	y = convention.value * s * ax[0]
	om[2][1] = x + y
	om[1][2] = x - y
	x = omc * ax[2] * ax[0]
	y = convention.value * s * ax[1]
	om[2][0] = x - y
	om[0][2] = x + y
	return om

#A.9
def ax2ro(ax):
	if abs(ax[3]) <= epsilon:
		return [0.0, 0.0, 1.0, 0.0]
	return ax[:3] + [float('inf') if abs(ax[3] - math.pi) <= epsilon else math.tan(ax[3] / 2.0)]

#A.10
def ax2qu(ax):
	if abs(ax[3]) <= epsilon:
		return [1.0, 0.0, 0.0, 0.0]
	s = math.sin(ax[3] / 2.0)
	qu = [math.cos(ax[3] / 2.0)] + [x * s for x in ax[:3]]
	mag = math.sqrt(math.fsum([x * x for x in qu]))
	return [x / mag for x in qu]

#A.11
def ax2ho(ax):
	k = math.pow(0.75 * ( ax[3] - math.sin(ax[3]) ), 1.0 / 3.0)
	return [x * k for x in ax[:3]]

#A.12
def ro2ax(ro):
	if abs(ro[3]) <= epsilon:
		return [0.0, 0.0, 1.0, 0.0]
	omega = 2.0 * math.atan(ro[3]) if ro[3] < float('inf') else math.pi
	return [x for x in ro[:3]] + [omega]

#A.13
def ro2ho(ro):
	t = 2.0 * math.atan(ro[3]) if ro[3] < float('inf') else math.pi
	f = math.pow(0.75 * (t - math.sin(t)), 1.0 / 3.0)
	return [x * f for x in ro[:3]]

#A.14
def qu2eu(qu):
	eu = [0.0] * 3
	qu0 = qu[0] * convention.value
	q03 = qu0 * qu0 + qu[3] * qu[3]
	q12 = qu[1] * qu[1] + qu[2] * qu[2]
	chi = math.sqrt(q03 * q12)
	if chi <= epsilon:
		if q12 <= epsilon:
			eu[0] = math.atan2(-2.0 * qu0 * qu[3], qu0 * qu0 - qu[3] * qu[3])
		else:
			eu[0] = math.atan2(2.0 * qu[1] * qu[2], qu[1] * qu[1] - qu[2] * qu[2])
			eu[1] = math.pi
	else:
		eu[0] = math.atan2((qu[1] * qu[3] - qu0 * qu[2]) / chi, (-qu[2] * qu[3] - qu0 * qu[1]) / chi)
		eu[1] = math.atan2(2.0 * chi, q03 - q12)
		eu[2] = math.atan2((qu[1] * qu[3] + qu0 * qu[2]) / chi, (qu[2] * qu[3] - qu0 * qu[1]) / chi)
	return [x + 2.0 * math.pi if x < 0.0 else x for x in eu]

#A.15
def qu2om(qu):
	om = [[0.0] * 3 for i in range(3)]
	qbar = qu[0] * qu[0] - math.fsum([x * x for x in qu[1:]])
	om[0][0] = qbar + 2.0 * qu[1] * qu[1]
	om[1][1] = qbar + 2.0 * qu[2] * qu[2]
	om[2][2] = qbar + 2.0 * qu[3] * qu[3]
	om[0][1] = 2.0 * (qu[1] * qu[2] - convention.value * qu[0] * qu[3])
	om[1][0] = 2.0 * (qu[2] * qu[1] + convention.value * qu[0] * qu[3])
	om[0][2] = 2.0 * (qu[1] * qu[3] + convention.value * qu[0] * qu[2])
	om[2][0] = 2.0 * (qu[3] * qu[1] - convention.value * qu[0] * qu[2])
	om[1][2] = 2.0 * (qu[2] * qu[3] - convention.value * qu[0] * qu[1])
	om[2][1] = 2.0 * (qu[3] * qu[2] + convention.value * qu[0] * qu[1])
	return om

#A.16
def qu2ax(qu):
	omega = 2.0 * math.acos(qu[0])
	if omega <= epsilon:
		return [0.0, 0.0, 1.0, 0.0]
	s = math.copysign(1.0 / math.sqrt(math.fsum([x * x for x in qu[1:]])), qu[0])
	return [s * n for n in qu[1:]] + [omega]

#A.17
def qu2ro(qu):
	if qu[0] <= epsilon:
		return qu[1:] + [float('inf')]
	s = math.sqrt(math.fsum([x * x for x in qu[1:]]))
	if s <= epsilon:
		return [0.0, 0.0, 1.0, 0.0]
	return [x / s for x in qu[1:]] + [math.tan(math.acos(qu[0]))]

#A.18
def qu2ho(qu):
	omega = 2.0 * math.acos(qu[0])
	if abs(omega) <= epsilon:
		return [0.0] * 3
	s = 1.0 / math.sqrt(math.fsum([x * x for x in qu[1:]]))
	f = math.pow(0.75 * (omega - math.sin(omega)), 1.0 / 3.0)
	return [s * x * f for x in qu[1:]]

#A.19
def ho2ax(ho):
	mag2 = math.fsum([x * x for x in ho])
	if mag2 <= epsilon:
		return [0.0, 0.0, 1.0, 0.0]
	theta = hoInv(mag2)
	ax = [x / math.sqrt(mag2) for x in ho]
	if theta >= math.pi - epsilon:
		ax = orientAxis(ax)
		theta = math.pi
	return ax + [theta]

def ho2cu(ho):
	# check bounds, get pyramid, and shuffle coordinates to +z pyramid
	rs = math.sqrt(math.fsum([x * x for x in ho])) # radius
	if rs > hoR() + 7.0 * epsilon:
		raise ValueError("'%f' lies outside the sphere of radius (3*pi/4)^(1/3)" % ho)
	p = pyramidType(ho)
	cu = ho[p:] + ho[:p]

	# handle origin
	if rs <= epsilon:
		return [0.0] * 3

	# invert operation M3
	cu[:2] = [x * math.sqrt(2.0 * rs / (rs + abs(cu[2]))) for x in cu[:2]]
	cu[2] = -rs * math.sqrt(math.pi / 6.0) if cu[2] < 0.0 else rs / math.sqrt(6.0 / math.pi)

	# invert operation M2
	sq = sorted([x * x for x in cu[:2]])
	mag = math.fsum(sq)
	if mag <= epsilon:
		cu = [0.0, 0.0, cu[2]]
	else:
		swapped = False
		if abs(cu[0]) > abs(cu[1]):
			swapped = True
			cu[0], cu[1] = cu[1], cu[0]
		k = math.sqrt((mag + sq[1]) * sq[1])
		sign = [-1 if x < 0.0 else 1 for x in cu[:2]]
		cu[:2] = [x * math.sqrt(math.pi / 3) * math.sqrt((mag + sq[1]) * mag / ((mag + sq[1]) - k)) / 2.0 for x in sign]
		k = (sq[0] + k) / (mag * math.sqrt(2.0))
		cu[0] *= 12.0 * math.acos(1.0 if 1.0 - k <= epsilon else k) / math.pi
		if swapped:
			cu[0], cu[1] = cu[1], cu[0]

	#invert operation M1, unshuffle coordinates, and return
	cu = [x / math.pow(math.pi / 6.0, 1.0 / 6.0) for x in cu]
	return cu[-p:] + cu[:-p]

def cu2ho(cu):
	# check bounds, get pyramid, and shuffle coordinates to +z pyramid
	if max([abs(i) for i in cu]) > cuA() / 2.0 + 7.0 * epsilon:
		raise ValueError("'%f' lies outside the cube of side length pi^(2/3)" % cu)
	p = pyramidType(cu)
	ho = cu[p:] + cu[:p]

	# handle origin
	if abs(ho[2]) <= epsilon:
		return [0.0] * 3

	# operation M1
	ho = [i * math.pow(math.pi / 6.0, 1.0 / 6.0) for i in ho]

	# operation M2
	if max([abs(i) for i in ho[:2]]) <= epsilon:
		ho = [0.0, 0.0, ho[2]] # handle points along z axis (to avoid divide by zero)
	else:
		swapped = False
		if abs(ho[0]) > abs(ho[1]):
			swapped = True
			ho[0], ho[1] = ho[1], ho[0]
		theta = (math.pi * ho[0]) / (12.0 * ho[1])
		k = math.sqrt(3.0 / math.pi) * math.pow(2.0, 0.75) * ho[1] / math.sqrt(math.sqrt(2.0) - math.cos(theta))
		ho[0] = math.sqrt(2.0) * math.sin(theta) * k
		ho[1] = (math.sqrt(2.0) * math.cos(theta) - 1.0) * k
		if swapped:
			ho[0], ho[1] = ho[1], ho[0]

	# operation M3
	k = ho[0]*ho[0] + ho[1]*ho[1]
	ho[:2] = [i * math.sqrt(1.0 - math.pi * k / (24.0 * ho[2]*ho[2])) for i in ho[:2]]
	ho[2] = math.sqrt(6.0 / math.pi) * ho[2] - k * math.sqrt(math.pi / 24) / ho[2]

	# unshuffle coordinates
	return ho[-p:] + ho[:-p]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Indirect Conversions                                                            #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def eu2ho(eu):
	return ax2ho(eu2ax(eu))
def eu2cu(eu):
	return ho2cu(eu2ho(eu))

def om2ro(om):
	return eu2ro(om2eu(om))
def om2ho(om):
	return ax2ho(om2ax(om))
def om2cu(om):
	return ho2cu(om2ho(om))

def ax2eu(ax):
	return om2eu(ax2om(ax))
def ax2cu(ax):
	return ho2cu(ax2ho(ax))

def ro2eu(ro):
	return om2eu(ro2om(ro))
def ro2om(ro):
	return ax2om(ro2ax(ro))
def ro2qu(ro):
	return ax2qu(ro2ax(ro))
def ro2cu(ro):
	return ho2cu(ro2ho(ro))

def qu2cu(qu):
	return ho2cu(qu2ho(qu))

def ho2eu(ho):
	return ax2eu(ho2ax(ho))
def ho2om(ho):
	return ax2om(ho2ax(ho))
def ho2ro(ho):
	return ax2ro(ho2ax(ho))
def ho2qu(ho):
	return ax2qu(ho2ax(ho))

def cu2eu(cu):
	return ho2eu(cu2ho(cu))
def cu2om(cu):
	return ho2om(cu2ho(cu))
def cu2ax(cu):
	return ho2ax(cu2ho(cu))
def cu2ro(cu):
	return ho2ro(cu2ho(cu))
def cu2qu(cu):
	return ho2qu(cu2ho(cu))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Testing                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
def defaultDist(a, b):
	return max([abs(i - j) for i, j in zip(a, b)])

def euDist(a, b):
	return defaultDist(eu2qu(a), eu2qu(b)) # infinite representations when 2nd angle is 0 or pi

def roDist(a, b):
	return defaultDist(ro2ax(a), ro2ax(b)) # tan(w/2) comparisons are bad for w near pi

def omDist(a, b):
	return max([defaultDist(i, j) for i, j in zip(a, b)]) # handle matrix layout

def test(n = 3, output = sys.stdout, verbose = False):
	# create test vector, total rotations = 4 * (n^3 - n^2) + n
	eulerList = []
	phi = [math.pi * i / (n - 1) for i in range((n - 1) * 2 + 1)]
	theta = [math.pi * i / (n - 1) for i in range(n)]
	for i in range(len(phi)):
		for j in range(len(theta)):
			for k in range(len(phi)):
				eulerList.append([phi[i], theta[j], phi[k]])

	# build matrix of functions
	names = ['eu', 'om', 'ax', 'ro', 'qu', 'ho', 'cu']
	representations = len(names)
	comparisons = [euDist, omDist, defaultDist, roDist, defaultDist, defaultDist, defaultDist]
	conversions = [[ None, eu2om, eu2ax, eu2ro, eu2qu, eu2ho, eu2cu],
	               [om2eu,  None, om2ax, om2ro, om2qu, om2ho, om2cu],
	               [ax2eu, ax2om,  None, ax2ro, ax2qu, ax2ho, ax2cu],
	               [ro2eu, ro2om, ro2ax,  None, ro2qu, ro2ho, ro2cu],
	               [qu2eu, qu2om, qu2ax, qu2ro,  None, qu2ho, qu2cu],
	               [ho2eu, ho2om, ho2ax, ho2ro, ho2qu,  None, ho2cu],
	               [cu2eu, cu2om, cu2ax, cu2ro, cu2qu, cu2ho,  None]]

	# check x = y2x(x2y(x)) 
	maxDiff = 0.0
	output.write('pairwise tests (' + convention.name + '):\n')
	for i in range(representations):
		if verbose:
			output.write(names[i] + ' test\n') 
		for j in range(representations):
			if i == j:
				continue
			for eu in eulerList:
				try:
					base = eu if 0 == i else conversions[0][i](eu)
					conv = conversions[j][i](conversions[i][j](base))
					diff = comparisons[i](conv, base)
					if diff > 1e-6 or verbose:
						output.write(names[i] + '2' + names[j] + ' max difference(' + str(base) + ') = ' + str(diff) + '\n')
					if diff > maxDiff:
						maxDiff = diff
						maxIndex = [i, j, base, eulerList.index(eu)]
				except ValueError as err:
					output.write(names[i] + '2' + names[j] + '[' + str(eulerList.index(eu)) + ']: ' + str(err) + '\n')
	output.write('max diff pairwise = ' + names[maxIndex[0]] + '2' + names[maxIndex[1]] + '(' + str(maxIndex[2]) + ') = ' + str(maxDiff) + '\n')

	# check x = z2x(y2z(x2y(x)))
	maxDiff = 0.0
	output.write('triplet tests (' + convention.name + '):\n')
	for i in range(representations):
		if verbose:
			output.write(names[i] + ' test\n')
		for j in range(representations):
			if i == j:
				continue
			for k in range(representations):
				if i == k or j == k:
					continue
				for eu in eulerList:
					try:
						base = eu if 0 == i else conversions[0][i](eu)
						conv = conversions[k][i](conversions[j][k](conversions[i][j](base)))
						diff = comparisons[i](conv, base)
						if diff > 1e-6 or verbose:
							output.write(names[k] + '2' + names[i] + '-' + names[j] + '2' + names[k] + '-' + names[i] + '2' + names[j] + ' max difference(' + str(base) + ') = ' + str(diff) + '\n')
						if diff > maxDiff:
							maxDiff = diff
							maxIndex = [i, j, k, base, eulerList.index(eu)]
					except ValueError as err:
						output.write(names[i] + '2' + names[j] + '[' + str(eulerList.index(eu)) + ']: ' + str(err) + '\n')
	output.write('max diff triplet = ' + names[maxIndex[2]] + '2' + names[maxIndex[0]] + '-' + names[maxIndex[1]] + '2' + names[maxIndex[2]] + '-' + names[maxIndex[0]] + '2' + names[maxIndex[1]] + '(' + str(maxIndex[3]) + ') = ' + str(maxDiff) + '\n')