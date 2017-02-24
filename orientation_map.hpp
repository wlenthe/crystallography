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
#ifndef _orientation_map_h_
#define _orientation_map_h_

#include <stack>
#include <vector>
#include <deque>
#include <set>
#include <memory>//shared_ptr
#include <algorithm>
#include <bitset>
#include <numeric>//iota
#include <limits>//infinity

#include "quaternion.hpp"
#include "rotations.hpp"
#include "symmetry.hpp"

#include "tsl.hpp"

template <typename T>
class OrientationMap {
	void readTSL(std::string fileName);
	public:
		//pixel array are stored in row major order with a normal cartesian coordinate system (not an image coordinate system)
		//keep shared pointers for pixel level arrays to avoid copying
		std::shared_ptr< std::vector< Quaternion<T> > > quats;//requires T be floating point type
		std::shared_ptr< std::vector<size_t> > phases;
		std::shared_ptr< std::vector<T> > quality;

		std::vector<Symmetry<T> const *> syms;
		size_t rows, cols;
		T xRes, yRes;

		OrientationMap(std::string fileName);//read an orientation map from a file
};

template <typename T>
class Segmentation : public OrientationMap<T> {
	static T DefaultThreshold() {return -std::numeric_limits<T>::infinity();}
	void segment(T tol, T qualityThreshold);
	public:
		std::vector<size_t> ids;
		std::vector< Quaternion<T> > avgQuats;
		std::vector<size_t> avgPhases;
		std::vector<size_t> numPixels;
		size_t numGrains;//includes grain 0, max grain id is numGrains - 1

		Segmentation(OrientationMap<T>& om, T tol = T(5), T thr = DefaultThreshold()) : OrientationMap<T>(om) {segment(tol, thr);}
		template<typename U> std::vector<U> mapGrainToPixel(const std::vector<U>& grainArray, const size_t numComp = 1) const;//copy a per grain array to a per pixel array
		std::vector< std::set<size_t> > findNeighbors() const;
		std::vector<std::uint8_t> mapColor();//6 color map coloring (nice alternative to unique id coloring) - linear 5 / quadratic 4 color algorithms exist but are more complex
};

////////////////////////////////////////////////////////////////////////////////////
//                        orientation map member functions                        //
////////////////////////////////////////////////////////////////////////////////////
template <typename T>
OrientationMap<T>::OrientationMap(std::string fileName) :
	quats(std::make_shared< std::vector< Quaternion<T> > >()),
	phases(std::make_shared< std::vector<size_t> >()),
	quality(std::make_shared< std::vector<T> >()),
	rows(0), cols(0), xRes(1), yRes(1) {
	if(TSL<T>::CanReadFile(fileName)) readTSL(fileName);
	else throw std::invalid_argument("no readers available for file type");
}

template <typename T>
void OrientationMap<T>::readTSL(std::string fileName) {
	//read tsl file and copy metadata
	TSL<T> file(fileName);
	if(tsl::GridType::Square != file.gridType) throw std::runtime_error("only square grids are supported");
	rows = file.nRows;
	cols = std::max(file.nColsOdd, file.nColsEven);
	xRes = file.xStep;
	yRes = file.yStep;
	syms = file.getSymmetries();

	//take scan data data
	phases->swap(file.phases);
	quality->swap(file.iq);

	//convert from eulers to quats
	quats->resize(phases->size());
	for(size_t i = 0; i < quats->size(); i++) rotations::eu2qu(file.eulers.data()+3*i, quats->operator[](i).data());
}

////////////////////////////////////////////////////////////////////////////////////
//                         segmentation member functions                          //
////////////////////////////////////////////////////////////////////////////////////
template <typename T>
void Segmentation<T>::segment(T tol, T qualityThreshold) {
	bool useQualityThreshold = DefaultThreshold() != qualityThreshold;
	tol = std::cos(tol * std::acos(T(0)) / T(180));//convert cutoff from degrees to cos(angle/2) once instead of acos on every misorientation
	size_t grains = 0;
	std::stack<size_t> s;
	ids.assign(OrientationMap<T>::rows * OrientationMap<T>::cols, 0);
	std::deque<size_t> seedPixels(1, 0);
	size_t* phases = OrientationMap<T>::phases->data();
	Quaternion<T>* quats = OrientationMap<T>::quats->data();
	Symmetry<T>const** syms = OrientationMap<T>::syms.data();
	T* quality = OrientationMap<T>::quality->data();
	for(size_t j = 0; j < OrientationMap<T>::rows; j++) {//loop over rows
		for(size_t i = 0; i < OrientationMap<T>::cols; i++) {//loop over columns
			bool passThreshold = useQualityThreshold ? quality[OrientationMap<T>::cols*j+i] > qualityThreshold : true;
			if(0 == ids[OrientationMap<T>::cols*j+i] && passThreshold) {//check if this pixel has already been assigned
				//create stack of pixels for flood filling
				s.push(OrientationMap<T>::cols*j+i);
				ids[s.top()] = ++grains;
				seedPixels.push_back(s.top());

				//flood through connected pixels
				while(!s.empty()) {
					//get top pixel and convert to index
					const size_t index = s.top();
					const size_t indI = index % OrientationMap<T>::cols;
					const size_t indJ = (index - indI) / OrientationMap<T>::cols;
					s.pop();

					//loop over neighbors (left, right, top, bottom) and check if connected (left, right, top, bottom)
					const bool exists[4] = {indI > 0, indI+1 < OrientationMap<T>::cols, indJ > 0, indJ+1 < OrientationMap<T>::rows};
					const size_t neighbor[4] = {index-1, index+1, index-OrientationMap<T>::cols, index+OrientationMap<T>::cols};
					for(size_t k = 0; k < 4; k++) {
						if(exists[k]) {//don't go off edge of scan
							passThreshold = useQualityThreshold ? quality[neighbor[k]] > qualityThreshold : true;
							if(phases[index] == phases[neighbor[k]] && 0 == ids[neighbor[k]] && passThreshold) {//check phases match and not already assigned
								if(syms[phases[index]]->disoQu(quats[index], quats[neighbor[k]]).w > tol) {//check disorientation angle
									ids[neighbor[k]] = grains;//assign to current grain
									s.push(neighbor[k]);//add pixel to stack
								}
							}
						}
					}
				}
			}
		}
	}
	numGrains = grains + 1;

	//determine grain phases
	avgPhases.reserve(numGrains);
	avgPhases.push_back(0);//grain 0
	for(size_t i = 1; i < numGrains; i++) avgPhases.push_back(phases[seedPixels[i]]);

	//replace the seed pixel with the highest quality pixel of that grain and move all seed pixels to fundamental zone 
	for(size_t i = 0; i < ids.size(); i++) {
		if(quality[i] > quality[seedPixels[ids[i]]]) seedPixels[ids[i]] = i;
	}
	for(size_t i = 1; i < numGrains; i++) quats[seedPixels[i]] = syms[avgPhases[i]]->fzQu(quats[seedPixels[i]]);

	//move all orientations to be the symmetrically equivalent representation closest to the seed pixel's orientation
	for(size_t i = 0; i < ids.size(); i++)
		quats[i] = syms[phases[i]]->nearbyQu(quats[seedPixels[ids[i]]], quats[i]);

	//convert all orientations to cubochoric space, average, and convert back to quaternions
	std::vector<T> avgCus(3 * numGrains, T(0));
	numPixels.assign(numGrains, 0);
	avgQuats.resize(numGrains);
	T cu[3];
	for(size_t i = 0; i < ids.size(); i++) {
		rotations::qu2cu(quats[i].data(), cu);
		std::transform(cu, cu+3, avgCus.data()+3*ids[i], avgCus.data()+3*ids[i], std::plus<T>());
		numPixels[ids[i]]++;
	}
	avgQuats[0] = Quaternion<T>::Zero();
	for(size_t i = 1; i < numGrains; i++) {
		avgCus[3*i+0] /= T(numPixels[i]);
		avgCus[3*i+1] /= T(numPixels[i]);
		avgCus[3*i+2] /= T(numPixels[i]);
		rotations::cu2qu(avgCus.data()+3*i, avgQuats[i].data());
		avgQuats[i] = syms[avgPhases[i]]->fzQu(avgQuats[i]);
	}
}

template <typename T>
template <typename U>
std::vector<U> Segmentation<T>::mapGrainToPixel(const std::vector<U>& grainArray, const size_t numComp) const {
	std::vector<U> pixelArray;
	size_t totalPixels = OrientationMap<T>::rows * OrientationMap<T>::cols;
	pixelArray.reserve(totalPixels * numComp);
	for(size_t i = 0; i < totalPixels; i++) {
		for(size_t j = 0; j < numComp; j++)
			pixelArray.push_back(grainArray[ids[i] * numComp + j]);
	}
	return pixelArray;
}

template <typename T>
std::vector< std::set<size_t> > Segmentation<T>::findNeighbors() const {
	std::vector< std::set<size_t> > neighbors(numGrains);
	for(size_t j = 1; j < OrientationMap<T>::rows; j++) {//loop over rows
		const size_t offset = j * OrientationMap<T>::cols;
		for(size_t i = 1; i < OrientationMap<T>::cols; i++) {//loop over columns
			const size_t id0 = ids[offset + i];
			const size_t idx = ids[offset + i - 1];
			const size_t idy = ids[offset + i - OrientationMap<T>::cols];
			if(id0 != idx && id0 > 0 && idx > 0) {
				neighbors[id0].insert(idx);
				neighbors[idx].insert(id0);
			}
			if(id0 != idy && id0 > 0 && idy > 0) {
				neighbors[id0].insert(idy);
				neighbors[idy].insert(id0);
			}
		}
	}
	return neighbors;
}

template <typename T>
std::vector<std::uint8_t> Segmentation<T>::mapColor() {
	std::vector<std::uint8_t> colors(numGrains, 0x00);
	if(numGrains <= 5) {
		std::iota(colors.begin(), colors.end(), 0x00);//handle trivial cases
	} else {
		//map coloring algorithm requires repeated removal of lowest degree node from graph
		//bucket based structure to store graph representation of neighbors (nodes of equal degree are fungible)
		//get lowest degree node: O(1), remove lowest degree node: ~O(ln(N)) for worst case (could use unordered_set to get O(1) typical, O(N) worst case)

		//compute degree of each node, determine min/max bucket, and fill buckets
		size_t maxDegree = 0;
		std::vector< size_t > bucketIds;//the current bucket for each node
		bucketIds.reserve(numGrains);
		bucketIds.push_back(0);//sk node 0
		std::vector< std::set<size_t> > neighbors = findNeighbors();//build adjacency list for each node
		for(size_t i = 1; i < numGrains; i++) {//skip node 0
			bucketIds.push_back(neighbors[i].size());//save degree of node
			if(bucketIds.back() > maxDegree) maxDegree = bucketIds.back();//track largest degree
		}
		std::vector< std::set<size_t> > buckets(maxDegree+1);//a bucket for each degree containing all nodes with that degree
		for(size_t i = 1; i < numGrains; i++) buckets[bucketIds[i]].insert(i);

		//remove node from lowest -> highest degree
		size_t minDegree = 0;
		std::vector<size_t> deletedNodes;
		deletedNodes.reserve(numGrains);
		for(size_t i = 1; i < numGrains; i++) {
			while(buckets[minDegree].empty()) ++minDegree;//if we've exhausted the lowest bucket find next non empty bucket
			size_t node = *buckets[minDegree].begin();//select any node from the lowest degree bucket
			buckets[minDegree].erase(buckets[minDegree].begin());//remove node from bucket
			for(std::set<size_t>::const_iterator iter = neighbors[node].cbegin(); iter != neighbors[node].cend(); ++iter) {//adjust node neighbors (maximum 5 for planar graph)
				neighbors[*iter].erase(node);//remove node from neighbor's neighbors
				buckets[bucketIds[*iter]].erase(*iter);//remove node from current bucket
				if(--bucketIds[*iter] < minDegree) --minDegree;//update neighbor bucket id, updating minimum degree if needed
				buckets[bucketIds[*iter]].insert(*iter);//add node to next bucket down
			}
			deletedNodes.push_back(node);//save the order that nodes are deleted
		}

		//color from highest -> lowest degree
		std::bitset<8> colorUsed;//can generalize to non-planar hueristic with vector<bool> and throw if more than 255 colors are used (no bound on chromatic color)
		for(size_t i = 1; i < numGrains; i++) {//color each node with lowest available color
			colorUsed.reset();
			for(std::set<size_t>::const_iterator iter = neighbors[deletedNodes.back()].cbegin(); iter != neighbors[deletedNodes.back()].cend(); ++iter)
				colorUsed[colors[*iter]] = true;
			     if(!colorUsed.test(1)) colors[deletedNodes.back()] = 0x01;
			else if(!colorUsed.test(2)) colors[deletedNodes.back()] = 0x02;
			else if(!colorUsed.test(3)) colors[deletedNodes.back()] = 0x03;
			else if(!colorUsed.test(4)) colors[deletedNodes.back()] = 0x04;
			else if(!colorUsed.test(5)) colors[deletedNodes.back()] = 0x05;
			else if(!colorUsed.test(6)) colors[deletedNodes.back()] = 0x06;
			else throw std::logic_error("insufficient colors for graph");//planar graphs are gaurenteed to be 6 colorable
			deletedNodes.pop_back();
		}
	}
	return mapGrainToPixel(colors);
}

#endif//_orientation_map_h_