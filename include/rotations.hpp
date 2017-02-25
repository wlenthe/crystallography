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

//orientation transform routines based on
// -Rowenhorst, David, et al. "Consistent Representations of and Conversions Between 3D Rotations." Model. Simul. Mater. Sci. Eng. 23.8 (2015): 083501.
// -Rosca, D., et al. "A New Method of Constructing a Grid in the Space of 3D rotations and its Applications to Texture Analysis." Model. Simul. Mater. Sci. Eng. 22.7 (2014): 075013.
// -fortran implementation of routines by Marc De Graef (https://github.com/marcdegraef/3Drotations)

//the following conventions are used:
// -quaternions as [w, x, y, z]
// -rotation angle <= pi
// -rotation axis in positive z hemisphere for rotations of pi
// -rotation axis = [0, 0, 1] for rotations of 0

#ifndef _rotations_h_
#define _rotations_h_

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <type_traits>//is_floating_point
#include <algorithm>
#include <numeric>//inner_product
#include <functional>//negate

namespace rotations{
	template <typename T>
	struct Constants {
		static_assert(std::is_floating_point<T>::value, "T must be a floating point type");//disallow integers
		static_assert(std::numeric_limits<T>::has_infinity, "T must have infinity");//must have ieee infinity
		static const T thr;
		static const T pi;
		static const T inf;
		static const T hoR;//radius of homochoric sphere
		static const T cuA;//side length of cubochoric cube
		static const T active, passive;
		static T convention;
	};

	template <typename T> const T Constants<T>::thr = T(10) * std::numeric_limits<T>::epsilon();
	template <typename T> const T Constants<T>::pi = T(3.14159265358979323846264338328);
	template <typename T> const T Constants<T>::inf = std::numeric_limits<T>::infinity();
	template <typename T> const T Constants<T>::hoR = std::pow(T(3) * pi / T(4), T(1) / T(3));
	template <typename T> const T Constants<T>::cuA = std::pow(pi, T(2) / T(3));
	template <typename T> const T Constants<T>::active = T(-1);
	template <typename T> const T Constants<T>::passive = T(1);
	template <typename T> T Constants<T>::convention = Constants<T>::passive;//default to passive convention

	//helper function to select equivalent axis for rotations of pi
	template <typename T> inline void orientAxis(T * const n) {
		if(std::fabs(n[2]) > Constants<T>::thr) {
			if(n[2] < T(0)) std::transform(n, n+3, n, std::negate<T>());//+z hemisphere
		} else {
			n[2] = T(0);
			if(std::fabs(n[1]) > Constants<T>::thr) {
				if(n[1] < T(0)) std::transform(n, n+2, n, std::negate<T>());//+y semicircle
			} else {
				if(n[0] < T(0)) n[0] = -n[0];//+x00
				n[1] = T(0);
			}
		}
	}

	//helpers to wrap 3x3 om functions as 9x1
	template<typename T> void xx2om(void(*conv)(T const * const, T * const * const), T const * const xx, T * const om) {
		T * const m[3] = {&om[0], &om[3], &om[6]};
		conv(xx, m);
	}
	template<typename T> void om2xx(void(*conv)(T const * const * const, T * const), T const * const om, T * const xx) {
		T const * const m[3] = {&om[0], &om[3], &om[6]};
		conv(m, xx);
	}

	////////////////////////////////////////////////////////////////////////////////////
	//                               direct conversions                               //
	////////////////////////////////////////////////////////////////////////////////////
	//A.1
	template<typename T> void eu2om(T const * const eu, T * const * const om) {
		const T c0 = std::cos(eu[0]);
		const T c1 = std::cos(eu[1]);
		const T c2 = std::cos(eu[2]);
		const T s0 = std::sin(eu[0]);
		const T s1 = std::sin(eu[1]);
		const T s2 = std::sin(eu[2]);

		om[0][0] =  c0 * c2 - s0 * c1 * s2;
		om[0][1] =  s0 * c2 + c0 * c1 * s2;
		om[0][2] =  s1 * s2;
		om[1][0] = -c0 * s2 - s0 * c1 * c2;
		om[1][1] = -s0 * s2 + c0 * c1 * c2;
		om[1][2] =  s1 * c2;
		om[2][0] =  s0 * s1;
		om[2][1] = -c0 * s1;
		om[2][2] =  c1;
		for(size_t i = 0; i < 3; i++) {
			for(size_t j = 0; j < 3; j++) {
				if(std::fabs(om[i][j]) <= Constants<T>::thr) om[i][j] = T(0);
			}
		}
	}
	template<typename T> void eu2om(T const * const eu, T * const om) {xx2om(&eu2om, eu, om);}

	//A.2
	template<typename T> void eu2ax(T const * const eu, T * const ax) {
		const T t = std::tan(eu[1] / T(2));
		const T sigma = (eu[0] + eu[2]) / T(2);
		const T s = std::sin(sigma);
		const T tau = std::sqrt(t * t + s * s);
		if(std::fabs(tau) <= T(2) * Constants<T>::thr) {
			std::fill(ax, ax+4, T(0));
			ax[2] = T(1);
			return;
		}
		const T delta = (eu[0] - eu[2]) / T(2);
		T alpha = std::fabs(sigma - Constants<T>::pi / T(2)) <= Constants<T>::thr ? Constants<T>::pi : T(2) * std::atan(tau / std::cos(sigma));
		const T k = -Constants<T>::convention / std::copysign(tau, alpha);
		ax[0] = k * std::cos(delta) * t;
		ax[1] = k * std::sin(delta) * t;
		ax[2] = k * std::sin(sigma);

		//normalize
		const T mag = std::sqrt(std::inner_product(ax, ax+3, ax, T(0)));
		std::transform(ax, ax+3, ax, [mag](const T i){return i/mag;});

		//handle ambiguous case (rotation angle of pi)
		alpha = std::fabs(alpha);
		if(alpha + Constants<T>::thr >= Constants<T>::pi) {
			orientAxis(ax);
			ax[3] = Constants<T>::pi;
		} else if(alpha <= Constants<T>::thr) {
			std::fill(ax, ax+4, T(0));
			ax[3] = T(1);
		} else {
			ax[3] = alpha;
		}
	}

	//A.3
	template<typename T> void eu2ro(T const * const eu, T * const ro) {
		eu2ax(eu, ro);
		ro[3] = ro[3] == Constants<T>::pi ? Constants<T>::inf : std::tan(ro[3] / T(2));
	}

	//A.4
	template<typename T> void eu2qu(T const * const eu, T * const qu) {
		T c = std::cos(eu[1] / T(2));
		T s = std::sin(eu[1] / T(2));
		T sigma = (eu[0] + eu[2]) / T(2);
		T delta = (eu[0] - eu[2]) / T(2);
		qu[0] = c * std::cos(sigma);
		qu[1] = s * std::cos(delta) * T(-Constants<T>::convention);
		qu[2] = s * std::sin(delta) * T(-Constants<T>::convention);
		qu[3] = c * std::sin(sigma) * T(-Constants<T>::convention);
		if(qu[0] < T(0)) std::transform(qu, qu+4, qu, std::negate<T>());

		//normalize
		T mag = std::sqrt(std::inner_product(qu, qu+4, qu, T(0)));
		std::transform(qu, qu+4, qu, [mag](const T i){return i/mag;});

		//handle ambiguous case (rotation angle of pi)
		if(std::abs(qu[0]) <= Constants<T>::thr) {
			orientAxis(&qu[1]);
			qu[0] = T(0);
		}
	}

	//A.5
	template<typename T> void om2eu(T const * const * const om, T * const eu) {
		if(std::fabs(om[2][2]) >= T(1) - Constants<T>::thr) {
			if(om[2][2] > T(0)) {
				eu[0] =  std::atan2( om[0][1], om[0][0]);//eu = [_, 0, _]
				eu[1] = T(0);
			} else {
				eu[0] = -std::atan2(-om[0][1], om[0][0]);//eu = [_, pi, _]
				eu[1] = Constants<T>::pi;
			}
			eu[2] = T(0);
		} else {
			eu[1] = std::acos(om[2][2]);
			const T zeta = T(1) / std::sqrt(T(1) - om[2][2] * om[2][2]);
			eu[0] = std::atan2(om[2][0] * zeta, -om[2][1] * zeta);
			eu[2] = std::atan2(om[0][2] * zeta,  om[1][2] * zeta);
		}
		std::for_each(eu, eu+3, [](T& i){if(i < T(0)) i += T(2) * Constants<T>::pi;});
	}
	template<typename T> void om2eu(T const * const om, T * const eu) {
		T const * const m[3] = {&om[0], &om[3], &om[6]};
		om2eu(m, eu);
	}

	//A.6
	template<typename T> void om2ax(T const * const * const om, T * const ax) {
		T omega = (om[0][0] + om[1][1] + om[2][2] - T(1)) / T(2);
		if(omega >= T(1) - Constants<T>::thr) {
			std::fill(ax, ax+4, T(0));
			ax[2] = T(1);
			return;
		}

		//compute eigenvector for eigenvalue of 1 (cross product of 2 adjacent columns of A-y*I)
		const T om00 = om[0][0] - T(1);
		const T om11 = om[1][1] - T(1);
		const T om22 = om[2][2] - T(1);
		const T vecs[3][3] = {
			{om[1][0]*om[2][1] - om[2][0]*  om11  , om[2][0]*om[0][1] -   om00  *om[2][1],   om00  *  om11   - om[1][0]*om[0][1]},
			{  om11  *  om22   - om[2][1]*om[1][2], om[2][1]*om[0][2] - om[0][1]*  om22  , om[0][1]*om[1][2] -   om11  *om[0][2]},
			{om[1][2]*om[2][0] -   om22  *om[1][0],   om22  *  om00   - om[0][2]*om[2][0], om[0][2]*om[1][0] - om[1][2]*  om00  }
		};

		//select vector with largest magnitude
		const T mags[3] = {
			std::sqrt(std::inner_product(vecs[0], vecs[0]+3, vecs[0], T(0))),
			std::sqrt(std::inner_product(vecs[1], vecs[1]+3, vecs[1], T(0))),
			std::sqrt(std::inner_product(vecs[2], vecs[2]+3, vecs[2], T(0)))
		};
		const size_t i = std::distance(mags, std::max_element(mags, mags+3));
		if(mags[i] <= Constants<T>::thr) {
			std::fill(ax, ax+4, T(0));
			ax[2] = T(1);
			return;
		}
		const T mag = mags[i];
		std::transform(vecs[i], vecs[i]+3, ax, [mag](const T i){return i/mag;});

		//handle ambiguous case (rotation of pi)
		if(omega <= Constants<T>::thr - T(1)) {
			orientAxis(ax);
			ax[3] = Constants<T>::pi;
			return;
		}

		//check axis sign
		ax[0] = std::copysign(ax[0], Constants<T>::convention * (om[2][1] - om[1][2]));
		ax[1] = std::copysign(ax[1], Constants<T>::convention * (om[0][2] - om[2][0]));
		ax[2] = std::copysign(ax[2], Constants<T>::convention * (om[1][0] - om[0][1]));
		ax[3] = std::acos(omega);
	}
	template<typename T> void om2ax(T const * const om, T * const ax) {om2xx(&om2ax, om, ax);}

	//A.7
	template<typename T> void ax2qu(T const * const ax, T * const qu);//required to handle ambiguous case
	template<typename T> void om2qu(T const * const * const om, T * const qu) {
		qu[0] = T(1) + om[0][0] + om[1][1] + om[2][2];
		qu[1] = T(1) + om[0][0] - om[1][1] - om[2][2];
		qu[2] = T(1) - om[0][0] + om[1][1] - om[2][2];
		qu[3] = T(1) - om[0][0] - om[1][1] + om[2][2];

		//handle ambiguous case (rotation of pi)
		if(std::fabs(qu[0]) <= Constants<T>::thr) {
			om2ax(om, qu);
			ax2qu(qu, qu);
			return;
		}

		//handle rotation of 0
		if(qu[0] <= Constants<T>::thr - T(2)) {
			std::fill(qu+1, qu+4, T(0));
			qu[0] = T(1);
			return;
		}

		std::transform(qu, qu+4, qu, [](const T i){return i <= Constants<T>::thr ? T(0) : Constants<T>::convention * std::sqrt(i) / T(2);});
		if(Constants<T>::convention * om[1][2] > Constants<T>::convention * om[2][1]) qu[1] = -qu[1];
		if(Constants<T>::convention * om[2][0] > Constants<T>::convention * om[0][2]) qu[2] = -qu[2];
		if(Constants<T>::convention * om[0][1] > Constants<T>::convention * om[1][0]) qu[3] = -qu[3];

		//ensure rotation angle <= pi
		if(qu[0] < T(0)) std::transform(qu, qu+4, qu, std::negate<T>());

		//normalize
		T mag = std::sqrt(std::inner_product(qu, qu+4, qu, T(0)));
		std::transform(qu, qu+4, qu, [mag](const T i){return i/mag;});
	}
	template<typename T> void om2qu(T const * const om, T * const qu) {om2xx(&om2qu, om, qu);}

	//A.8
	template<typename T> void ax2om(T const * const ax, T * const * const om) {
		if(std::fabs(ax[3]) <= Constants<T>::thr) {
			for(size_t i = 0; i < 3; i++) {
				std::fill(om[i], om[i]+3, T(0));
				om[i][i] = T(1);
			}
		} else if(std::fabs(ax[3]) >= Constants<T>::pi - Constants<T>::thr) {
			for(size_t i = 0; i < 3; i++)
				om[i][i] =T(2) * ax[i] * ax[i] - T(1);
			for(size_t i = 0; i < 3; i++) {
				const size_t j = (i+1)%3;
				const T x = T(2) * ax[i] * ax[j];
				om[i][j] = x;
				om[j][i] = x;
			}
		} else {
			const T c = std::cos(ax[3]);
			const T s = std::sin(ax[3]);
			const T omc = T(1) - c;
			for(size_t i = 0; i < 3; i++)
				om[i][i] = c + omc * ax[i] * ax[i];
			for(size_t i = 0; i < 3; i++) {
				const size_t j = (i+1)%3;
				const T x = omc * ax[i] * ax[j];
				const T y = Constants<T>::convention * s * ax[(i+2)%3];
				om[i][j] = x - y;
				om[j][i] = x + y;
			}
		}

	}
	template<typename T> void ax2om(T const * const ax, T * const om) {xx2om(&ax2om, ax, om);}

	//A.9
	template<typename T> void ax2ro(T const * const ax, T * const ro) {
		if(std::fabs(ax[3]) <= Constants<T>::thr) {
			std::fill(ro, ro+4, T(0));
			ro[2] = T(1);
		} else {
			std::copy(ax, ax+3, ro);
			ro[3] = std::fabs(ax[3] - Constants<T>::pi) <= Constants<T>::thr ? Constants<T>::inf : std::tan(ax[3] / T(2));
		}
	}

	//A.10
	template<typename T> void ax2qu(T const * const ax, T * const qu) {
		if(std::fabs(ax[3]) <= Constants<T>::thr) {
			std::fill(qu, qu+4, T(0));
			qu[0] = T(1);
		} else {
			if(ax == qu) std::rotate(qu, qu+3, qu+4);//rotate_copy doesn't work if source and destination are the same
			else std::rotate_copy(ax, ax+3, ax+4, qu);
			const T s = std::sin(qu[0] / T(2));
			qu[0] = std::cos(qu[0] / T(2));
			std::for_each(qu+1, qu+4, [s](T& i){i*=s;});

			//normalize
			T mag = std::sqrt(std::inner_product(qu, qu+4, qu, T(0)));
			std::for_each(qu, qu+4, [mag](T& i){i/=mag;});
		}
	}

	//A.11
	template<typename T> void ax2ho(T const * const ax, T * const ho) {
		if(std::fabs(ax[3]) <= Constants<T>::thr) {
			std::fill(ho, ho+3, T(0));
		} else {
			const T k = std::pow(T(3) / T(4) * ( ax[3] - std::sin(ax[3]) ), T(1) / T(3));
			std::transform(ax, ax+3, ho, [k](const T i){return i*k;});
		}
	}

	//A.12
	template<typename T> void ro2ax(T const * const ro, T * const ax) {
		if(std::fabs(ro[3]) <= Constants<T>::thr) {
			std::fill(ax, ax+4, T(0));
			ax[2] = T(1);
		} else {
			std::copy(ro, ro+3, ax);
			ax[3] = ro[3] < Constants<T>::inf ? T(2) * std::atan(ro[3]) : Constants<T>::pi;
		}
	}

	//A.13
	template<typename T> void ro2ho(T const * const ro, T * const ho) {
		const T t = ro[3] < Constants<T>::inf ? T(2) * std::atan(ro[3]) : Constants<T>::pi;
		if(std::fabs(t) <= Constants<T>::thr) {
			std::fill(ho, ho+3, T(0));
		} else {
			const T k = std::pow(T(3) * ( t - std::sin(t) ) / T(4), T(1) / T(3));
			std::transform(ro, ro+3, ho, [k](const T i){return i*k;});
		}
	}

	//A.13
	template<typename T> void qu2eu(T const * const qu, T * const eu) {
		const T qu0 = qu[0] * Constants<T>::convention;
		const T q03 = qu0 * qu0 + qu[3] * qu[3];
		const T q12 = qu[1] * qu[1] + qu[2] * qu[2];
		const T chi = std::sqrt(q03 * q12);
		if(chi <= Constants<T>::thr) {
			if(q12 <= Constants<T>::thr){
				eu[0] = std::atan2(-T(2) * qu0 * qu[3], qu0 * qu0 - qu[3] * qu[3]);
				eu[1] = T(0);
			} else {
				eu[0] = std::atan2(T(2) * qu[1] * qu[2], qu[1] * qu[1] - qu[2] * qu[2]);
				eu[1] = Constants<T>::pi;
			}
			eu[2] = T(0);
		} else {
			const T y1 = qu[1] * qu[3];//can divide by chi (but atan2 is magnitude independent)
			const T y2 = qu[2] * qu0;  //can divide by chi (but atan2 is magnitude independent)
			const T x1 =-qu[1] * qu0;  //can divide by chi (but atan2 is magnitude independent)
			const T x2 = qu[2] * qu[3];//can divide by chi (but atan2 is magnitude independent)
			eu[0] = std::atan2(y1 - y2, x1 - x2);
			eu[1] = std::atan2(T(2) * chi, q03 - q12);
			eu[2] = std::atan2(y1 + y2, x1 + x2);
		}
		std::for_each(eu, eu+3, [](T& i){if(i < T(0)) i += T(2) * Constants<T>::pi;});
	}

	//A.15
	template<typename T> void qu2om(T const * const qu, T * const * const om) {
		const T qbar = qu[0] * qu[0] - qu[1] * qu[1] - qu[2] * qu[2] - qu[3] * qu[3];
		for(size_t i = 0; i < 3; i++) {
			om[i][i] = qbar + T(2) * qu[i+1] * qu[i+1];

			const size_t j = (i+1)%3;
			const T x = T(2) * qu[i+1] * qu[j+1];
			const T y = T(2) * Constants<T>::convention * qu[0] * qu[(i+2)%3+1];
			om[i][j] = x - y;
			om[j][i] = x + y;
		}
	}
	template<typename T> void qu2om(T const * const qu, T * const om) {xx2om(&qu2om, qu, om);}

	//A.16
	template<typename T> void qu2ax(T const * const qu, T * const ax) {
		const T omega = T(2) * std::acos(qu[0]);
		if(omega <= Constants<T>::thr) {
			std::fill(ax, ax+4, T(0));
			ax[2] = T(1);
		} else {
			const T s = std::copysign(T(1) / std::sqrt(std::inner_product(qu+1, qu+4, qu+1, T(0))), qu[0]);
			std::transform(qu+1, qu+4, ax, [s](const T i){return i*s;});
			if(omega >= Constants<T>::pi - Constants<T>::thr) {
				orientAxis(ax);
				ax[3] = Constants<T>::pi;			
			} else {
				ax[3] = omega;
			}
		}
	}

	//A.17
	template<typename T> void qu2ro(T const * const qu, T * const ro) {
		if(qu[0] <= Constants<T>::thr) {
			std::copy(qu+1, qu+4, ro);
			ro[3] = Constants<T>::inf;
			return;
		}
		const T s = std::sqrt(std::inner_product(qu+1, qu+4, qu+1, T(0)));
		if(s <= Constants<T>::thr) {
			std::fill(ro, ro+4, T(0));
			ro[2] = T(1);
		} else {
			std::transform(qu+1, qu+4, ro, [s](const T i){return i/s;});
			ro[3] = std::tan(std::acos(qu[0]));
		}
	}

	//A.18
	template<typename T> void qu2ho(T const * const qu, T * const ho) {
		const T omega = T(2) * std::acos(qu[0]);
		if(std::fabs(omega) <= Constants<T>::thr) {
			std::fill(ho, ho+3, T(0));
		} else {
			const T s = T(1) / std::sqrt(std::inner_product(qu+1, qu+4, qu+1, T(0)));
			const T k = std::pow(T(3) * ( omega - std::sin(omega) ) / T(4), T(1) / T(3)) * s;
			std::transform(qu+1, qu+4, ho, [k](const T i){return i*k;});
		}
	}

	//inverse of ((3/4)*(x-sin(x)))^(2/3) via newton's method
	template <typename T> T hoInv(T y) {
	static const T hoR2 = std::pow(T(3) * Constants<T>::pi / T(4), T(2) / T(3));
		if(y < T(0) || y >= hoR2 + Constants<T>::thr) {
			std::stringstream ss;
			ss << "homochoric magnitude^2 " << y << " is outside of 0, (3*pi/4)^(2/3)]";
			throw std::domain_error(ss.str());
		}
		if(y >= hoR2 - Constants<T>::thr) return Constants<T>::pi;
		else if(y <= Constants<T>::thr) return T(0);
		T x = T(2) * std::acos(T(1) - y / T(2));//initial guess from small taylor expansion
		y = std::sqrt(y);
		T prevErr = Constants<T>::inf;
		for(size_t i = 0; i < 7; i++) {
			const T fx = std::pow(T(3) * (x - std::sin(x)) / T(4), T(1) / T(3));//compute forward value
			const T delta = fx - y;//compute error
			const T err = std::fabs(delta);
			if(0 == delta || err == prevErr) return x;//no error or flipping between +/- v
			x -= T(4) * fx * fx * delta / (T(1) - std::cos(x));//update
			if(err > prevErr) return x;//flipping between v / -2v
			prevErr = err;
		}
		std::stringstream ss;
		ss << "failed to invert ((3/4)*(x-sin(x)))^(2/3) for " << y;
		throw std::runtime_error(ss.str());
		return 0;
	}

	//A.19
	template<typename T> void ho2ax(T const * const ho, T * const ax) {
		const T mag2 = std::inner_product(ho, ho+3, ho, T(0));
		if(mag2 <= Constants<T>::thr) {
				std::fill(ax, ax+4, T(0));
				ax[2] = T(1);
		} else {
			ax[3] = hoInv(mag2);
			const T mag = std::sqrt(mag2);
			std::transform(ho, ho+3, ax, [mag](const T i){return i/mag;});
			if(ax[3] >= Constants<T>::pi - Constants<T>::thr) {
				orientAxis(ax);
				ax[3] = Constants<T>::pi;
			}
		}
	}

	//sextant type for cubochoric <-> homochoric transformation symmetry
	template<typename T> size_t pyramidType(T const * const v) {
		std::pair<T const *, T const *> minMax = std::minmax_element(v, v+3);
		const T minEl = std::fabs(*(minMax.first));
		const T maxEl =  *(minMax.second);
		return std::distance(v, minEl > maxEl ? minMax.first : minMax.second);
	}

	//forward and reverse shuffling for cubochoric <-> homochoric transformation symmetry
	template<typename T> void shufflePyramid(T * const v, const size_t p, bool unshuffle = false) {
		if(p != 2)	std::rotate(v, v + (unshuffle ? 2-p : p+1), v+3);
	}

	template<typename T> void cu2ho(T const * const cu, T * const ho) {
		//get pyramid type, shuffle coordinates to +z pyramid, and check bounds
		const size_t p = pyramidType(cu);
		std::copy(cu, cu+3, ho);
		shufflePyramid(ho, p);
		if(std::fabs(ho[2]) >= Constants<T>::cuA / T(2) + Constants<T>::thr) {
			std::stringstream ss;
			ss << "cubochoric component " << ho[2] << " is outside of +/-pi^(2/3)";
			throw std::domain_error(ss.str());
		}

		//handle origin
		if(std::fabs(ho[2]) <= Constants<T>::thr) {
			std::fill(ho, ho+3, T(0));
			return;
		}

		//operation M1
		static const T k1 = std::pow(Constants<T>::pi / T(6), T(1) / T(6));
		std::transform(ho, ho+3, ho, [](const T i){return i*k1;});

		//operation M2
		bool swapped = std::fabs(ho[0]) > std::fabs(ho[1]);
		if(swapped) std::swap(ho[0], ho[1]);
		if(std::fabs(ho[1]) >= Constants<T>::thr) {//skip points along z axis to avoid divide by zero
			static const T k2 = Constants<T>::pi / T(12);
			static const T k3 = std::sqrt(T(3) / Constants<T>::pi) * std::pow(T(2), T(3) / T(4));
			static const T k4 = std::sqrt(T(2));
			const T theta = k2 * ho[0] / ho[1];
			const T k = k3 * ho[1] / std::sqrt(k4 - std::cos(theta));
			ho[0] = k * k4 * std::sin(theta);
			ho[1] = k * k4 * std::cos(theta) - k;
			if(swapped) std::swap(ho[0], ho[1]);
		} else {
			std::fill(ho, ho+2, T(0));
		}

		//operation M3
		static const T k5 = std::sqrt(T(6) / Constants<T>::pi);
		static const T k6 = std::sqrt(Constants<T>::pi / T(24));
		const T ko = std::inner_product(ho, ho+2, ho, T(0));
		const T k = std::sqrt(T(1) - Constants<T>::pi * ko / (T(24) * ho[2]*ho[2]));
		std::transform(ho, ho+2, ho, [k](const T i){return i*k;});
		ho[2] = k5 * ho[2] - ko * k6 / ho[2];

		//unshuffle
		shufflePyramid(ho, p, true);
	}

	template<typename T> void ho2cu(T const * const ho, T * const cu) {
		//check bounds, get pyramid type, and shuffle coordinates to +z pyramid
		std::copy(ho, ho+3, cu);
		const T rs = std::sqrt(std::inner_product(cu, cu+3, cu, T(0)));
		if(rs >= Constants<T>::hoR + Constants<T>::thr) {
			std::stringstream ss;
			ss << "homochoric magnitude " << rs << " is outside of 0, (3*pi/4)^(2/3)]";
			throw std::domain_error(ss.str());
		}
		const size_t p = pyramidType(cu);
		shufflePyramid(cu, p);

		//handle origin
		if(rs <= Constants<T>::thr) {
			std::fill(cu, cu+3, T(0));
			return;
		}

		//invert operation M3
		static const T k5 = std::sqrt(T(6) / Constants<T>::pi);
		const T k = std::sqrt(T(2) * rs / (rs + std::fabs(cu[2])));
		std::transform(cu, cu+2, cu, [k](const T i){return i*k;});
		cu[2] = std::copysign(rs / k5, cu[2]);

		//invert operation M2
		T x2 = cu[0] * cu[0];
		T y2 = cu[1] * cu[1];
		const T mag2 = x2 + y2;
		if(mag2 >= Constants<T>::thr) {//skip points along z axis
			//only handle x <= y
			bool swapped = std::fabs(cu[0]) > std::fabs(cu[1]);
			if(swapped) {
				std::swap(cu[0], cu[1]);
				std::swap(x2, y2);
			}
			const T x2y = mag2 + y2;
			const T rx2y = std::sqrt(x2y);
			static const T k4 = std::sqrt(T(2));
			const T kxy = std::sqrt( (Constants<T>::pi / T(3)) * (x2y * mag2) / (x2y - std::fabs(cu[1]) * rx2y) ) / T(2);
			const T ckx = (x2 + std::fabs(cu[1]) * rx2y) / k4 / mag2;
			std::transform(cu, cu+2, cu, [kxy](const T i){return std::copysign(kxy, i);});
			if(ckx >= T(1))
				cu[0] = T(0);
			else if (ckx <= T(-1))
				cu[0] *= T(12);
			else
				cu[0] *= std::acos(ckx) / (Constants<T>::pi / T(12));
			if(swapped) std::swap(cu[0], cu[1]);
		} else {
			std::fill(cu, cu+2, T(0));
		}

		//invert operation M1
		static const T k1 = std::pow(Constants<T>::pi / T(6), T(1) / T(6));
		std::transform(cu, cu+3, cu, [](const T i){return i/k1;});

		//unshuffle
		shufflePyramid(cu, p, true);
	}

	////////////////////////////////////////////////////////////////////////////////////
	//                              indirect conversions                              //
	////////////////////////////////////////////////////////////////////////////////////
	//eu --> xx
	template<typename T> void eu2ho(T const * const eu, T * const ho) {
		T ax[4];
		eu2ax(eu, ax);
		ax2ho(ax, ho);
	}

	template<typename T> void eu2cu(T const * const eu, T * const cu) {
		eu2ho(eu, cu);
		ho2cu(cu, cu);
	}

	//om --> xx
	template<typename T> void om2ro(T const * const * const om, T * const ro) {
		om2eu(om, ro);
		eu2ro(ro, ro);
	}
	template<typename T> void om2ro(T const * const om, T * const ro) {om2xx(&om2ro, om, ro);}

	template<typename T> void om2ho(T const * const * const om, T * const ho) {
		T ax[4];
		om2ax(om, ax);
		ax2ho(ax, ho);
	}
	template<typename T> void om2ho(T const * const om, T * const ho) {om2xx(&om2ho, om, ho);}

	template<typename T> void om2cu(T const * const * const om, T * const cu) {
		om2ho(om, cu);
		ho2cu(cu, cu);
	}
	template<typename T> void om2cu(T const * const om, T * const cu) {om2xx(&om2cu, om, cu);}

	//ax --> xx
	template<typename T> void ax2eu(T const * const ax, T * const eu) {
		T om[9];
		ax2om(ax, om);
		om2eu(om, eu);
	}

	template<typename T> void ax2cu(T const * const ax, T * const cu) {
		ax2ho(ax, cu);
		ho2cu(cu, cu);
	}

	//ro --> xx
	template<typename T> void ro2om(T const * const ro, T * const om);
	template<typename T> void ro2eu(T const * const ro, T * const eu) {
		T om[9];
		ro2om(ro, om);
		om2eu(om, eu);
	}

	template<typename T> void ro2om(T const * const ro, T * const * const om) {
		T ax[4];
		ro2ax(ro, ax);
		ax2om(ax, om);
	}
	template<typename T> void ro2om(T const * const ro, T * const om) {xx2om(&ro2om, ro, om);}

	template<typename T> void ro2qu(T const * const ro, T * const qu) {
		ro2ax(ro, qu);
		ax2qu(qu, qu);
	}

	template<typename T> void ro2cu(T const * const ro, T * const cu) {
		ro2ho(ro, cu);
		ho2cu(cu, cu);
	}

	//qu --> xx
	template<typename T> void qu2cu(T const * const qu, T * const cu) {
		qu2ho(qu, cu);
		ho2cu(cu, cu);
	}

	//ho --> xx
	template<typename T> void ho2eu(T const * const ho, T * const eu) {
		T ax[4];
		ho2ax(ho, ax);
		ax2eu(ax, eu);
	}

	template<typename T> void ho2om(T const * const ho, T * const * const om) {
		T ax[4];
		ho2ax(ho, ax);
		ax2om(ax, om);
	}
	template<typename T> void ho2om(T const * const ho, T * const om) {xx2om(&ho2om, ho, om);}

	template<typename T> void ho2ro(T const * const ho, T * const ro) {
		ho2ax(ho, ro);
		ax2ro(ro, ro);
	}

	template<typename T> void ho2qu(T const * const ho, T * const qu) {
		ho2ax(ho, qu);
		ax2qu(qu, qu);
	}

	//cu --> xx
	template<typename T> void cu2eu(T const * const cu, T * const eu) {
		cu2ho(cu, eu);
		ho2eu(eu, eu);
	}

	template<typename T> void cu2om(T const * const cu, T * const * const om) {
		T ho[3];
		cu2ho(cu, ho);
		ho2om(ho, om);
	}
	template<typename T> void cu2om(T const * const cu, T * const om) {xx2om(&cu2om, cu, om);}

	template<typename T> void cu2ax(T const * const cu, T * const ax) {
		cu2ho(cu, ax);
		ho2ax(ax, ax);
	}

	template<typename T> void cu2ro(T const * const cu, T * const ro) {
		T ho[3];
		cu2ho(cu, ho);
		ho2ro(ho, ro);
	}

	template<typename T> void cu2qu(T const * const cu, T * const qu) {
		cu2ho(cu, qu);
		ho2qu(qu, qu);
	}
}

#endif//_rotations_h_