import numpy, scipy.spatial
from PIL import Image
import matplotlib.pyplot as plt

import tsl, rotations
from symmetry import Symmetry

def cubicIPF(n):
	# move to stereographic triangle and manually handle 001 to avoid divide by zero
	n = sorted([abs(x) for x in n])
	if 0.0 == n[1]:
		return [255, 0, 0]

	# compute ipf color, rescale from 0-255, and return
	x = n[0] / n[1]
	y = min(numpy.math.acos(n[2]) / numpy.math.acos(numpy.math.sqrt(1 / (2 + x * x))), 1)
	x = min(4*numpy.math.atan(x)*y/numpy.pi, y)
	rgb = [1-y, y-x, x]
	rgb = [numpy.math.sqrt(x) for x in rgb]
	return [int(round(x * 255 / max(rgb))) for x in rgb]

# read an ang file, get size, and check if grid is square or hexagonal
scan = tsl.OimScan('scan.ang')
rows = scan.rows
cols = max(scan.colsOdd, scan.colsEven)
if scan.gridType is not tsl.Grid.Square:
	print('warning: assuming square grid for connectivity / image generation')

# show confidence index as image
plt.imshow(scan.ci)
plt.show()

# print information about phases
print(str(len(scan.phaseList)) + ' phases found:')
for phase in scan.phaseList:
	print('\tphase ' + str(phase.number) + ': ' + phase.materialName + ' - ' + str(Symmetry.FromTSL(phase.symmetry)))

# build mapping of phase number to symmetry
syms = [None]
if len(scan.phaseList) is 1:
	# for single phase scans phase is 0 everywhere
	syms[0] = Symmetry.FromTSL(scan.phaseList[0].symmetry)
else:
	# for multi phase scans phase is 1, 2, 3, ... n
	for phase in scan.phaseList:
			syms.append(Symmetry.FromTSL(phase.symmetry))

# compute ipf colors for cubic points and save as png
ipf = numpy.zeros(scan.euler.shape, dtype = 'uint8')
for j in range(rows):
	for i in range(cols):
		if scan.ci[j,i] >= 0 and syms[scan.phase[j,i]] is Symmetry.Cubic: # leave bad and non cubic pixels black
			om = numpy.array(rotations.eu2om(scan.euler[j,i]))
			ipf[j,i] = cubicIPF(om[:,2])
Image.fromarray(ipf).save('ipf.png')

# pad scan
top  = 10
bottom = 15
left  = 20
right = 25
scan.rows += top + bottom
scan.colsOdd += left + right
scan.colsEven += left + right

# tsl needs x/y values to be correct (can't just pad with zero)
xMin = scan.x[0,0] - left * scan.xStep
yMin = scan.y[0,0] - top * scan.yStep
xMax = scan.x[-1,-1] + right * scan.xStep
yMax = scan.y[-1,-1] + bottom * scan.yStep

# pad positions in x direction
scan.x = numpy.pad(scan.x, ((0, 0), (left, right)), 'linear_ramp', end_values = numpy.array([[0, 0], [xMin, xMax]]))
scan.y = numpy.pad(scan.y, ((0, 0), (left, right)), 'edge')

# pad positions in y direction
scan.x = numpy.pad(scan.x, ((top, bottom), (0, 0)), 'edge')
scan.y = numpy.pad(scan.y, ((top, bottom), (0, 0)), 'linear_ramp', end_values = numpy.array([[yMin, yMax], [0, 0]]))

# pad remaining arrays and write padded scan to new ang file
scan.euler = numpy.pad(scan.euler, ((top, bottom), (left, right), (0, 0)), 'constant', constant_values = 0)
scan.iq = numpy.pad(scan.iq, ((top, bottom), (left, right)), 'constant', constant_values = 0)
scan.ci = numpy.pad(scan.ci, ((top, bottom), (left, right)), 'constant', constant_values = -1)
scan.phase = numpy.pad(scan.phase, ((top, bottom), (left, right)), 'constant', constant_values = 0)
if scan.sem is not None:
	scan.sem = numpy.pad(scan.sem, ((top, bottom), (left, right)), 'constant', constant_values = 0)
if scan.fit is not None:
	scan.fit = numpy.pad(scan.fit, ((top, bottom), (left, right)), 'constant', constant_values = 0)
scan.writeAng('padded.ang')

# convert rotations from euler angles to quaternions
quats = numpy.zeros((rows, cols, 4))
for j in range(rows):
	for i in range(cols):
		quats[j,i,:] = rotations.eu2qu(scan.euler[j,i])

# segment grains with 5 degree misorientation tolerance (can be since slow disoQuat is an expensive operation for some symmetries)
grainIds = numpy.zeros((scan.euler.shape[0], scan.euler.shape[1]), dtype = 'int')
grainIds.fill(-1) # put every pixel in grain -1

numGrains = 0
pixelsBurned = 0
totalPixels = rows * cols
tolerance = numpy.math.cos((5.0 * numpy.pi / 180.0) / 2.0) # min w component of misorientation quaternion
for j in range(rows):
	for i in range(cols):
		pixel = (j,i)
		if grainIds[pixel] == -1: # skip pixels already assigned to a grain
			if scan.ci[pixel] > 0:
				numGrains += 1
				grainIds[pixel] = numGrains # seed new grain
				stack = [pixel] # add seed pixel to stack
				pixelsBurned += 1
				while len(stack) > 0:
					indJ,indI = stack.pop() # get pixel
					iMinus = indI - 1 if indI > 0 else None
					jMinus = indJ - 1 if indJ > 0 else None
					iPlus  = indI + 1 if indI + 1 < cols else None
					jPlus  = indJ + 1 if indJ + 1 < rows else None
					for neighbor in [(jMinus, indI), (jPlus, indI), (indJ, iMinus), (indJ, iPlus)]:
						if neighbor[0] is not None and neighbor[1] is not None: # check if neighboring pixel exists
							if grainIds[neighbor] == -1: # check if neighbor is already assigned to a grain
								if scan.phase[neighbor] == scan.phase[pixel]: # check if neighboring pixel is same phase
									misoW = syms[scan.phase[pixel]].disoQuat(quats[neighbor], quats[pixel])[0] # compute misorientation angle with neighbor
									if misoW >= tolerance:
										grainIds[neighbor] = numGrains # add to current grain
										stack.append(neighbor)
										pixelsBurned += 1
				progress = 100 * pixelsBurned / totalPixels
				print('segmenting: {0:.2f}'.format(progress) + '% complete', end = '\r')
			else:
				grainIds[j,i] = 0 # put bad pixels in grain 0
				pixelsBurned += 1
Image.fromarray(grainIds).save('grainIds.tif')

# compute the phase, average orientation, and number of pixels for each grain
pixels = numpy.zeros(numGrains)
avgCu = numpy.zeros((numGrains, 3))
phases = numpy.zeros(numGrains, dtype = 'int')
for j in range(rows):
	for i in range(cols):
		grain = grainIds[j,i]
		if grain > 0:
			pixels[grain-1] += 1
			phases[grain-1] = scan.phase[j,i]
			avgCu[grain-1] += rotations.qu2cu(syms[phases[grain-1]].fzQu(quats[j,i]))
avgCu[:,0] /= pixels
avgCu[:,1] /= pixels
avgCu[:,2] /= pixels

# generate a synthetic microstructure (voronoi tesselation) with 100 grains
print('generating voronoi tesselation')
seeds = 0
grains = 100
synRows = 256
synCols = 256
seedPoints = numpy.zeros((grains, 2))
syntheticGrains = numpy.zeros((synRows, synCols), dtype = 'int')
while seeds < grains:
	j = numpy.random.randint(0, synRows)
	i = numpy.random.randint(0, synCols)
	if 0 == syntheticGrains[j,i]:
		seedPoints[seeds] = [j,i]
		seeds += 1
		syntheticGrains[j,i] = seeds

tree = scipy.spatial.KDTree(seedPoints)
for j in range(synRows):
	for i in range(synCols):
		grain = tree.query(numpy.array([j,i]))[1]+1
		syntheticGrains[j,i] = grain

# generates a random number between -1 and 1 from the Epanechnikov distribution (for kernel density estimation)
def epanechnikov():
	return (2 * numpy.median(numpy.random.rand(3))) - 1

# generate orientations for each grain and convet to eulers
print('sampling orientations')
a = rotations.cuA();
a2 = a / 2;
avgEu = numpy.zeros((grains, 3))
for k in range(grains):
	if k < grains / 20:
		# first 5% of grains from random orientations
		cu = [numpy.random.uniform(-a2, a2), numpy.random.uniform(-a2, a2), numpy.random.uniform(-a2, a2)]
	else:
		# sample remainder from odf of input grains
		# select an orientation from the population
		cu = avgCu[numpy.random.randint(0, grains)]

		# apply spread (this effectively samples from the kernel density estimate)
		spread = numpy.array([epanechnikov(), epanechnikov(), epanechnikov()])
		spread *= 0.01 # ~1.5 degrees
		cu += spread

		#bring back within cube if needed
		for i in range(3):
			if cu[i] < -a2:
				cu[i] += a
			elif cu[i] > a2:
				cu[i] -= a
	avgEu[k] = numpy.array(rotations.cu2eu(list(cu)))

# create scan
scan = tsl.OimScan.zeros_like(syntheticGrains, resolution = (1.0, 1.0), origin = (0.0, 0.0))

# fill arrays
for j in range(synRows):
	for i in range(synCols):
		scan.euler[j,i] = avgEu[syntheticGrains[j,i]-1]
scan.iq.fill(1)
scan.ci.fill(1)

# fill phase list
scan.phaseList.append(tsl.OimPhase(1))
scan.phaseList[-1].materialName = 'Synthetic Phase'
scan.phaseList[-1].formula = 'Syn'
scan.phaseList[-1].symmetry = Symmetry.Cubic.toTSL()

# write to ang
scan.writeAng('synthetic.ang')
