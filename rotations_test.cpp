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
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>//transform
#include <numeric>//iota

#include "rotations.hpp"

template <typename T>
T defaultDist(T const * const a, T const * const b, const size_t i) {
	T maxDist = 0;
	for(size_t j = 0; j < i; j++) {
		const T dist = std::fabs(a[j] - b[j]);
		if(dist > maxDist) maxDist = dist;
	}
	return maxDist;
}

//distance for 2 euler angle triples (infinite representations when 2nd angle is 0 or pi)
template <typename T>
T euDist(T const * const a, T const * const b, const size_t i) {
	T qa[4], qb[4];
	rotations::eu2qu(a, qa);
	rotations::eu2qu(b, qb);
	return defaultDist(qa, qb, 4);
}

//distance for 2 rodriguez vectors (comparisons of w for rotations near pi are bad)
template <typename T>
T roDist(T const * const a, T const * const b, const size_t i) {
	T xa[4], xb[4];
	rotations::ro2ax(a, xa);
	rotations::ro2ax(b, xb);
	return defaultDist(xa, xb, 4);
}

//self consistency test
template <typename T>
void rotations_test(size_t n, bool threeWay, std::ostream& os = std::cout) {
	//build list of evenly spaced euler angles
	const T pi = rotations::Constants<T>::pi;
	std::vector<T> phi((n - 1) * 2 + 1);
	std::iota(phi.begin(), phi.end(), T(0));
	T n1 = (T) n-1; 
	std::transform(phi.begin(), phi.end(), phi.begin(), [&pi, &n1](const T i){return pi * i / n1;});
	std::vector<T> theta(phi.begin(), phi.begin() + n);
	//build tables of conversion / comparison functions
	std::string names[7] ={"eu", "om", "ax", "ro", "qu", "ho", "cu"};
	size_t lengths[7] = {3, 9, 4, 4, 4, 3, 3};
	typedef void (*conversionFunc)(T const * const, T * const);//signature for conversion function pointer
	typedef T (*comparisonFunc)(T const * const, T const * const, const size_t);//signature for comparison function pointer
	conversionFunc conversion[7][7] = {
		{                NULL, &rotations::eu2om<T>, &rotations::eu2ax<T>, &rotations::eu2ro<T>, &rotations::eu2qu<T>, &rotations::eu2ho<T>, &rotations::eu2cu<T>},
		{&rotations::om2eu<T>,                 NULL, &rotations::om2ax<T>, &rotations::om2ro<T>, &rotations::om2qu<T>, &rotations::om2ho<T>, &rotations::om2cu<T>},
		{&rotations::ax2eu<T>, &rotations::ax2om<T>,                 NULL, &rotations::ax2ro<T>, &rotations::ax2qu<T>, &rotations::ax2ho<T>, &rotations::ax2cu<T>},
		{&rotations::ro2eu<T>, &rotations::ro2om<T>, &rotations::ro2ax<T>,                 NULL, &rotations::ro2qu<T>, &rotations::ro2ho<T>, &rotations::ro2cu<T>},
		{&rotations::qu2eu<T>, &rotations::qu2om<T>, &rotations::qu2ax<T>, &rotations::qu2ro<T>,                 NULL, &rotations::qu2ho<T>, &rotations::qu2cu<T>},
		{&rotations::ho2eu<T>, &rotations::ho2om<T>, &rotations::ho2ax<T>, &rotations::ho2ro<T>, &rotations::ho2qu<T>,                 NULL, &rotations::ho2cu<T>},
		{&rotations::cu2eu<T>, &rotations::cu2om<T>, &rotations::cu2ax<T>, &rotations::cu2ro<T>, &rotations::cu2qu<T>, &rotations::cu2ho<T>,                 NULL}
	};
	comparisonFunc comparison[7] = {&euDist, &defaultDist, &defaultDist, &roDist, &defaultDist, &defaultDist, &defaultDist};

	//check x = y2x(x2y(x)) 
	T maxDiff = T(0);
	size_t maxI = 0, maxJ = 0, maxK = 0, maxM = 0, maxN = 0, maxP = 0;
	T eu[3], x[9], y[9], z[9], final[9];
	os << "pairwise tests\n";
	for(size_t i = 0; i < 7; i++) {
		for(size_t j = 0; j < 7; j++) {
			if(i == j) continue;
			for(size_t m = 0; m < phi.size(); m++) {
				eu[0] = phi[m];
				for(size_t n = 0; n < theta.size(); n++) {
					eu[1] = theta[n];
					for(size_t p = 0; p < phi.size(); p++) {
						eu[2] = phi[p];
						//get base orientation from eulers
						if(i == 0) std::copy(eu, eu+3, x);
						else conversion[0][i](eu, x);

						//do conversion
						conversion[i][j](x, y);//x2y
						conversion[j][i](y, final);//y2x

						//compute error
						T diff = comparison[i](x, final, lengths[i]);
						if(diff > maxDiff) {
							maxDiff = diff;
							maxI = i;
							maxJ = j;
							maxM = m;
							maxN = n;
							maxP = p;
						}
					}
				}
			}
		}
	}
	os << "max diff pairwise = " << names[maxI] << "2" << names[maxJ] << "-" << names[maxJ] << "2" << names[maxI] << ": " << maxDiff << "\n";
	os << "eu = (" << phi[maxM] << ", " << theta[maxN] << ", " << phi[maxP] << ")\n";

	//check x = z2x(y2z(x2y(x)))
	if(!threeWay) return;
	maxDiff = 0;
	maxI = maxJ = maxK = maxM = maxN = maxP = 0;
	os << "three way tests\n";
	for(size_t i = 0; i < 7; i++) {
		for(size_t j = 0; j < 7; j++) {
			if(i == j) continue;
				for(size_t k = 0; k < 7; k++) {
				if(j == k || k == i) continue;
				for(size_t m = 0; m < phi.size(); m++) {
					eu[0] = phi[m];
					for(size_t n = 0; n < theta.size(); n++) {
						eu[1] = theta[n];
						for(size_t p = 0; p < phi.size(); p++) {
							eu[2] = phi[p];
							//get base orientation from eulers
							if(i == 0) std::copy(eu, eu+3, x);
							else conversion[0][i](eu, x);

							//do conversion
							conversion[i][j](x, y);//x2y
							conversion[j][k](y, z);//y2z
							conversion[k][i](z, final);//z2x

							//compute error
							T diff = comparison[i](x, final, lengths[i]);
							if(diff > maxDiff) {
								maxDiff = diff;
								maxI = i;
								maxJ = j;
								maxK = k;
								maxM = m;
								maxN = n;
								maxP = p;
							}
						}
					}
				}
			}
		}
	}
	os << "max diff three way = " << names[maxI] << "2" << names[maxK] << "-" << names[maxK] << "2" << names[maxJ]  << "-" << names[maxJ] << "2" << names[maxI] << ": " << maxDiff << "\n";
	os << "eu = (" << phi[maxM] << ", " << theta[maxN] << ", " << phi[maxP] << ")\n";
}

int main() {
	try {
		//test parameters
		static const size_t n = 20;
		static const bool threeWay = true;
		std::ostream& os = std::cout;

		//self consistency test (passive)
		rotations_test<float >(n, threeWay, os);
		rotations_test<double>(n, threeWay, os);

		//switch to active convention
		rotations::Constants<float >::convention = rotations::Constants<float >::active;
		rotations::Constants<double>::convention = rotations::Constants<double>::active;

		//self consistency test (active)
		rotations_test<float >(n, threeWay, os);
		rotations_test<double>(n, threeWay, os);
	} catch (std::exception& e) {
		std::cout << e.what();
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}