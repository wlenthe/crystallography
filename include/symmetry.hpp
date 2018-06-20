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
#ifndef _symmetry_h_
#define _symmetry_h_

#include <cmath>
#include <string>
#include <iterator>
#include <algorithm>

#include "rotations.hpp"
#include "quaternion.hpp"
#include "orientation_coloring.hpp"

//symmetry abstract base class
template <typename T>
class Symmetry {
	private:
		virtual T minFzW() const = 0;//cos(minimum angle from origin to FZ boundary/2)

	public:
		static Symmetry const * FromName(std::string name);
		virtual std::string name() const {return "unknown";}
		virtual bool roInFz(T const * const ro) const = 0;//check if a rodrigues vector is in the fundamental zone
		virtual std::vector<Quaternion<T> > const * const quOperators() const = 0;//get symmetry operators as quaternions

		//get the equivalent quaternion in the fundamental zone
		virtual void fzQu(T const * const qu, T * const fz) const;
		virtual Quaternion<T> fzQu(const Quaternion<T>& qu) const {
			Quaternion<T> fz;
			fzQu(qu.data(), fz.data());
			return fz;
		}

		//get the equivalent quaternion in the fundamental zone
		virtual void disoQu(T const * const qu1, T const * const qu2, T * const diso) const;
		Quaternion<T> disoQu(const Quaternion<T>& qu1, const Quaternion<T>& qu2) const {
			Quaternion<T> diso;
			disoQu(qu1.data(), qu2.data(), diso.data());
			return diso;
		}

		//compute the symmetric equivalent orientation of qu2 closest to qu1
		void nearbyQu(T const * const qu1, T const * const qu2, T * const equiv) const;
		virtual Quaternion<T> nearbyQu(const Quaternion<T>& qu1, const Quaternion<T>& qu2) const {
			Quaternion<T> equiv;
			nearbyQu(qu1.data(), qu2.data(), equiv.data());
			return equiv;
		}

		//compute ipf color (0-1 rgb) for orientation, assumes unit vector
		virtual void ipfColor(T const * const qu, T * const rgb, T const * const refDir) const = 0;
		virtual void ipfColor(const Quaternion<T>& qu, T * const rgb, T const * const refDir) const = 0;
};

//only centrosymmetric (Laue) groups are implemented
template <typename T> class CubicSymmetry;        //m-3m
template <typename T> class CubicLowSymmetry;     //m-3
template <typename T> class HexagonalSymmetry;    //6/mmm
template <typename T> class HexagonalLowSymmetry; //6/m
template <typename T> class TetragonalSymmetry;   //4/mmm
template <typename T> class TetragonalLowSymmetry;//4/m
template <typename T> class TrigonalSymmetry;     //-3m
template <typename T> class TrigonalLowSymmetry;  //-3
template <typename T> class OrthorhombicSymmetry; //mmm
template <typename T> class MonoclinicSymmetry;   //2/m
template <typename T> class TriclinicSymmetry;    //-1

namespace symmetry {
	template <typename T>
	struct Groups {
		static const CubicSymmetry<T>         Cubic;
		static const CubicLowSymmetry<T>      CubicLow;
		static const HexagonalSymmetry<T>     Hexagonal;
		static const HexagonalLowSymmetry<T>  HexagonalLow;
		static const TetragonalSymmetry<T>    Tetragonal;
		static const TetragonalLowSymmetry<T> TetragonalLow;
		static const TrigonalSymmetry<T>      Trigonal;
		static const TrigonalLowSymmetry<T>   TrigonalLow;
		static const OrthorhombicSymmetry<T>  Orthorhombic;
		static const MonoclinicSymmetry<T>    Monoclinic;
		static const TriclinicSymmetry<T>     Triclinic;
	};
	template <typename T> const CubicSymmetry<T>         Groups<T>::Cubic         = CubicSymmetry<T>        ();
	template <typename T> const CubicLowSymmetry<T>      Groups<T>::CubicLow      = CubicLowSymmetry<T>     ();
	template <typename T> const HexagonalSymmetry<T>     Groups<T>::Hexagonal     = HexagonalSymmetry<T>    ();
	template <typename T> const HexagonalLowSymmetry<T>  Groups<T>::HexagonalLow  = HexagonalLowSymmetry<T> ();
	template <typename T> const TetragonalSymmetry<T>    Groups<T>::Tetragonal    = TetragonalSymmetry<T>   ();
	template <typename T> const TetragonalLowSymmetry<T> Groups<T>::TetragonalLow = TetragonalLowSymmetry<T>();
	template <typename T> const TrigonalSymmetry<T>      Groups<T>::Trigonal      = TrigonalSymmetry<T>     ();
	template <typename T> const TrigonalLowSymmetry<T>   Groups<T>::TrigonalLow   = TrigonalLowSymmetry<T>  ();
	template <typename T> const OrthorhombicSymmetry<T>  Groups<T>::Orthorhombic  = OrthorhombicSymmetry<T> ();
	template <typename T> const MonoclinicSymmetry<T>    Groups<T>::Monoclinic    = MonoclinicSymmetry<T>   ();
	template <typename T> const TriclinicSymmetry<T>     Groups<T>::Triclinic     = TriclinicSymmetry<T>    ();

	template <typename T> struct operators {
		static const std::vector<Quaternion<T> > cubic;
		static const std::vector<Quaternion<T> > cubiclow;
		static const std::vector<Quaternion<T> > hexagonal;
		static const std::vector<Quaternion<T> > hexagonallow;
		static const std::vector<Quaternion<T> > tetragonal;
		static const std::vector<Quaternion<T> > tetragonallow;
		static const std::vector<Quaternion<T> > trigonal;
		static const std::vector<Quaternion<T> > trigonallow;
		static const std::vector<Quaternion<T> > orthorhombic;
		static const std::vector<Quaternion<T> > monoclinic;
		static const std::vector<Quaternion<T> > triclinic;
	};
}

////////////////////////////////////////////////////////////////////////////////////
//                      symmetry family intermediate classes                      //
////////////////////////////////////////////////////////////////////////////////////
template <typename T, size_t N>
class NFoldSymmetry : public Symmetry<T> {
	T minFzW() const {
		static const T kMw2 = T(1) / std::sqrt(T(2));                                 //cos(pi/4 )
		static const T kMw3 = std::sqrt(T(3)) / T(2);                                 //cos(pi/6 )
		static const T kMw4 = std::sqrt(T(2) + std::sqrt(T(2)) ) / T(2);              //cos(pi/8 )
		static const T kMw6 = ( T(1) + std::sqrt(T(3)) ) / ( T(2) * std::sqrt(T(2)) );//cos(pi/12)
		static const T kMw = (N == 2 ? kMw2 : (N == 3 ? kMw3 : (N == 4 ? kMw4 : (N == 6 ? kMw6 : 0))));//select constant based on N
		return kMw;
	}

	public:
		bool roInFz(T const * const ro) const {
			static const T kFz2 = T(1);                  //tan(pi/4 )
			static const T kFz3 = T(1) / std::sqrt(T(3));//tan(pi/6 )
			static const T kFz4 = std::sqrt(T(2)) - T(1);//tan(pi/8 )
			static const T kFz6 = T(2) - std::sqrt(T(3));//tan(pi/12)
			static const T kFz = (N == 2 ? kFz2 : (N == 3 ? kFz3 : (N == 4 ? kFz4 : (N == 6 ? kFz6 : 0))));//select constant based on N
			return std::fabs(ro[2] * ro[3]) <= kFz;//top and bottom faces z = +/-tan(pi/(2n))
		}
};

template <typename T, size_t N>
class CyclicSymmetry : public NFoldSymmetry<T, N> {
	public:
		virtual void ipfColor(T const * const qu, T * const rgb, T const * const refDir) const {
			T n[3];
			quaternion::rotateVector(qu, refDir, n);
			coloring::cyclicIpf<T,N>(n, rgb);
		}
		virtual void ipfColor(const Quaternion<T>& qu, T * const rgb, T const * const refDir) const {ipfColor(qu.data(), rgb, refDir);}
};

template <typename T, size_t N>
class DihedralSymmetry : public NFoldSymmetry<T, N> {
	public:
		bool roInFz(T const * const ro) const {
			if(!NFoldSymmetry<T, N>::roInFz(ro)) return false;

			//2n faces evenly spaced around z axis at distance 1
			T x = std::fabs(ro[0] * ro[3]);
			T y = std::fabs(ro[1] * ro[3]);
			if(y > x) std::swap(x, y);
			if(x > T(1)) return false;
			static const T kFz2  = T(0);
			static const T kFz3  = (T(1) + std::sqrt(T(2))) * (T(1) - T(1) / std::sqrt(T(3)));
			static const T kFz4  = std::sqrt(T(2));
			static const T kFz6  = (T(1) + std::sqrt(T(2))) * (std::sqrt(T(3)) - T(1));
			static const T kFzX2 = T(1);
			static const T kFzX3 = T(-1) - std::sqrt(T(2)) + ( T(2) + std::sqrt(T(2)) ) / std::sqrt(T(3));
			static const T kFzX4 = T(-1);
			static const T kFzX6 = std::sqrt(T(2)) + T(3) - std::sqrt(T(3)) * ( T(2) + std::sqrt(T(2)) );
			static const T kFz  = (N == 2 ? kFz2  : (N == 3 ? kFz3  : (N == 4 ? kFz4  : (N == 6 ? kFz6  : 0))));//select constant based on N
			static const T kFzX = (N == 2 ? kFzX2 : (N == 3 ? kFzX3 : (N == 4 ? kFzX4 : (N == 6 ? kFzX6 : 0))));//select constant based on N
			return y <= kFz + kFzX * x;//y <= (1+sqrt(2))*(1-t) + ((2+sqrt(2))*t - (1+sqrt(2))) * x
		}

		virtual void ipfColor(T const * const qu, T * const rgb, T const * const refDir) const {
			T n[3];
			quaternion::rotateVector(qu, refDir, n);
			coloring::dihedralIpf<T,N>(n, rgb);
		}
		virtual void ipfColor(const Quaternion<T>& qu, T * const rgb, T const * const refDir) const {ipfColor(qu.data(), rgb, refDir);}
};

template <typename T>
class PolyhedralSymmetry : public Symmetry<T> {
	T minFzW() const {static const T kMw = std::sqrt(T(2) + std::sqrt(T(2))) / T(2); return kMw;}
	public:
		bool roInFz(T const * const ro) const {return (ro[0] + ro[1] + ro[2]) * ro[3] <= T(1);}//regular octohedron with faces at tan(pi/6) away from origin
};

////////////////////////////////////////////////////////////////////////////////////
//                           specific symmetry classes                            //
////////////////////////////////////////////////////////////////////////////////////
template <typename T>
class CubicLowSymmetry : public PolyhedralSymmetry<T> {
	T minFzW() const {return T(1);}
	public:
		std::string name() const {return "Cubic Low (m-3m, Oh)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::cubiclow;}

		virtual void ipfColor(T const * const qu, T * const rgb, T const * const refDir) const {
			T n[3];
			quaternion::rotateVector(qu, refDir, n);
			coloring::cubicLowIpf(n, rgb);
		}
		virtual void ipfColor(const Quaternion<T>& qu, T * const rgb, T const * const refDir) const {ipfColor(qu.data(), rgb, refDir);}
};

template <typename T>
class CubicSymmetry : public PolyhedralSymmetry<T> {
	public:
		std::string name() const {return "Cubic (m-3, Th)";}
		bool roInFz(T const * const ro) const {
			//truncated cube (intersection of tetrahedral fz with cube of side length 2 * tan(pi/8))
			static const T kFz = std::sqrt(T(2)) - T(1);//tan(pi/8)
			if(!PolyhedralSymmetry<T>::roInFz(ro)) return false;
			return std::fabs(std::max(std::max(std::fabs(ro[0]), std::fabs(ro[1])), std::fabs(ro[2])) * ro[3]) <= kFz;
		}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::cubic;}

		//cubic shortcut for disorientation
		Quaternion<T> disoQu(const Quaternion<T>& qu1, const Quaternion<T>& qu2) const {return Symmetry<T>::disoQu(qu1, qu2);}//overloading -> base method hidden
		virtual void disoQu(T const * const qu1, T const * const qu2, T * const diso) const {
			T qu2conj[4];
			static const T kR2 = std::sqrt(T(2));
			quaternion::conjugate(qu2, qu2conj);
			quaternion::multiply(qu1, qu2conj, diso);
			quaternion::abs(diso, diso);
			std::sort(diso, diso+4, std::greater<T>());
			T operator11 = (diso[3] + diso[2] + diso[1] + diso[0]) / T(2);//add smallest -> largest for least round off
			T operator18 = (diso[0] + diso[1]) / kR2;
			if(operator11 >= operator18 && operator11 > diso[0]) {
				quaternion::multiply(diso, symmetry::operators<T>::cubic[11].data(), diso);
				quaternion::abs(diso, diso);
				std::sort(diso, diso+4, std::greater<T>());
			} else if(operator18 > operator11 && operator18 > diso[0]) {
				quaternion::multiply(diso, symmetry::operators<T>::cubic[18].data(), diso);
				quaternion::abs(diso, diso);
				std::sort(diso, diso+4, std::greater<T>());
			}
		}

		virtual void ipfColor(T const * const qu, T * const rgb, T const * const refDir) const {
			T n[3];
			quaternion::rotateVector(qu, refDir, n);
			coloring::cubicIpf(n, rgb);
		}
		virtual void ipfColor(const Quaternion<T>& qu, T * const rgb, T const * const refDir) const {ipfColor(qu.data(), rgb, refDir);}
};

template <typename T> class HexagonalSymmetry     : public DihedralSymmetry<T, 6> {
	public:
		std::string name() const {return "Hexagonal (6/mmm, D6h)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::hexagonal;}
};
template <typename T> class HexagonalLowSymmetry  : public CyclicSymmetry  <T, 6> {
	public:
		std::string name() const {return "Hexagonal Low (6/m, C6h)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::hexagonallow;}
};
template <typename T> class TetragonalSymmetry    : public DihedralSymmetry<T, 4> {
	public:
		std::string name() const {return "Tetragonal (4/mmm, D4h)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::tetragonal;}
};
template <typename T> class TetragonalLowSymmetry : public CyclicSymmetry  <T, 4> {
	public:
		std::string name() const {return "Tetragonal Low (4/m, C4h,)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::tetragonallow;}
};
template <typename T> class TrigonalSymmetry      : public DihedralSymmetry<T, 3> {
	public:
		std::string name() const {return "Trigonal (-3m, D3d)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::trigonal;}
};
template <typename T> class TrigonalLowSymmetry   : public CyclicSymmetry  <T, 3> {
	public:
		std::string name() const {return "Trigonal Low (-3, C3i)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::trigonallow;}
};
template <typename T> class OrthorhombicSymmetry  : public DihedralSymmetry<T, 2> {
	public:
		std::string name() const {return "Orthorhombic (mmm, D2h)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::orthorhombic;}
};
template <typename T> class MonoclinicSymmetry    : public CyclicSymmetry  <T, 2> {
	public:
		std::string name() const {return "Monoclinic (2/m, C2h)";}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::monoclinic;}
};

template <typename T>
class TriclinicSymmetry : public Symmetry<T> {
	T minFzW() const {return T(0);}
	public:
		std::string name() const {return "Triclinic (-1, Ci)";}
		bool roInFz(T const * const ro) const {return true;}
		std::vector<Quaternion<T> > const * const quOperators() const {return &symmetry::operators<T>::triclinic;}

		virtual void ipfColor(T const * const qu, T * const rgb, T const * const refDir) const {
			T n[3];
			quaternion::rotateVector(qu, refDir, n);
			coloring::hemiIpf(n, rgb);
		}
		virtual void ipfColor(const Quaternion<T>& qu, T * const rgb, T const * const refDir) const {ipfColor(qu.data(), rgb, refDir);}
};

////////////////////////////////////////////////////////////////////////////////////
//                      symmetry member function definitions                      //
////////////////////////////////////////////////////////////////////////////////////
template <typename T>
Symmetry<T> const * Symmetry<T>::FromName(std::string name){
	//remove spaces and convert to lowercase
	std::string::iterator end = std::remove_if(name.begin(), name.end(), [](char c){return std::isspace(c);});
	name.erase(end, name.end());
	std::transform(name.begin(), name.end(), name.begin(), [](char c){return std::tolower(c);});

	//get pointer
	if     (0 == name.compare(0,  5, "cubic"))         return &symmetry::Groups<T>::Cubic;
	else if(0 == name.compare(0,  8, "cubiclow")     ) return &symmetry::Groups<T>::CubicLow;
	else if(0 == name.compare(0,  9, "hexagonal")    ) return &symmetry::Groups<T>::Hexagonal;
	else if(0 == name.compare(0, 12, "hexagonallow") ) return &symmetry::Groups<T>::HexagonalLow;
	else if(0 == name.compare(0, 10, "tetragonal")   ) return &symmetry::Groups<T>::Tetragonal;
	else if(0 == name.compare(0, 13, "tetragonallow")) return &symmetry::Groups<T>::TetragonalLow;
	else if(0 == name.compare(0,  8, "trigonal")     ) return &symmetry::Groups<T>::Trigonal;
	else if(0 == name.compare(0, 11, "trigonallow")  ) return &symmetry::Groups<T>::TrigonalLow;
	else if(0 == name.compare(0, 12, "orthorhombic") ) return &symmetry::Groups<T>::Orthorhombic;
	else if(0 == name.compare(0, 10, "monoclinic")   ) return &symmetry::Groups<T>::Monoclinic;
	else if(0 == name.compare(0,  9, "triclinic")    ) return &symmetry::Groups<T>::Triclinic;
	else throw std::invalid_argument("uknown symmetry name " + name);
}

template <typename T>
void Symmetry<T>::fzQu(T const * const qu, T * const fz) const {
	T ro[4];
	std::vector<Quaternion<T> > const * const ops = quOperators();
	if(rotations::Constants<T>::convention == rotations::Constants<T>::passive) {
		//crystal symmetry * orientation * sample symmetry
		for(typename std::vector<Quaternion<T> >::const_iterator it = ops->begin(); it != ops->end(); ++it) {
			quaternion::multiply(it->data(), qu, fz);
			quaternion::explement(fz, fz);
			rotations::qu2ro(fz, ro);
			if(roInFz(ro)) break;
		}
	} else if(rotations::Constants<T>::convention == rotations::Constants<T>::active) {
		//sample symmetry * orientation * crystal symmetry
		for(typename std::vector<Quaternion<T> >::const_iterator it = ops->begin(); it != ops->end(); ++it) {
			quaternion::multiply(qu, it->data(), fz);
			quaternion::explement(fz, fz);
			rotations::qu2ro(fz, ro);
			if(roInFz(ro)) break;
		}
	} else {
		throw std::runtime_error("invalid convention state");
	}
}

template <typename T>
void Symmetry<T>::disoQu(T const * const qu1, T const * const qu2, T * const diso) const {
	T qu2conj[4], miso[4], minMiso[4];
	minMiso[0] = T(0);
	quaternion::conjugate(qu2, qu2conj);
	std::vector<Quaternion<T> > const * const ops = quOperators();
	if(rotations::Constants<T>::convention == rotations::Constants<T>::passive) {
		//crystal symmetry * orientation 1 * orientation 2 ^-1 * crystal symmetry ^-1
		T deltaQ[4], iDeltaQ[4];
		quaternion::multiply(qu1, qu2conj, deltaQ);
		for(typename std::vector<Quaternion<T> >::const_iterator qi = ops->begin(); qi != ops->end(); ++qi) {
			quaternion::multiply(qi->data(), deltaQ, iDeltaQ);
			for(typename std::vector<Quaternion<T> >::const_iterator qj = ops->begin(); qj != ops->end(); ++qj) {
				quaternion::multiply(iDeltaQ, qj->data(), miso);
				quaternion::explement(miso, miso);
				if(miso[0] > minMiso[0]) std::swap(miso, minMiso);
			}
		}
	} else if(rotations::Constants<T>::convention == rotations::Constants<T>::active) {
		//orientation 1 * crystal symmetry * orientation 2 ^-1
		T Iqu2Conj[4];
		quaternion::conjugate(qu2, qu2conj);
		for(typename std::vector<Quaternion<T> >::const_iterator qi = ops->begin(); qi != ops->end(); ++qi) {
			quaternion::multiply(qi->data(), qu2conj, Iqu2Conj);
			quaternion::multiply(qu1, Iqu2Conj, miso);;
			quaternion::explement(miso, miso);
			if(miso[0] > minMiso[0]) std::swap(miso, minMiso);
		}
	} else {
		throw std::runtime_error("invalid convention state");
	}
	std::copy(minMiso, minMiso+4, diso);
}

template <typename T>
void Symmetry<T>::nearbyQu(T const * const qu1, T const * const qu2, T * const equiv) const {
	T eq[4], wq[4];
	T maxDot = T(0);//dot product is w component of qu1 * qu2.conjugate()
	std::vector<Quaternion<T> > const * const ops = quOperators();
	if(rotations::Constants<T>::convention == rotations::Constants<T>::passive) {
		//crystal symmetry * orientation * sample symmetry
		for(typename std::vector<Quaternion<T> >::const_iterator qi = ops->begin(); qi != ops->end(); ++qi) {
			quaternion::multiply(qi->data(), qu2, wq);
			const T dot = std::fabs(quaternion::dot(qu1, wq));
			if(dot > maxDot) {
				maxDot = dot;
				std::swap(wq, eq);
				if(dot >= minFzW()) break;
			}
		}
	} else if(rotations::Constants<T>::convention == rotations::Constants<T>::active) {
		//sample symmetry * orientation * crystal symmetry
		for(typename std::vector<Quaternion<T> >::const_iterator qi = ops->begin(); qi != ops->end(); ++qi) {
			quaternion::multiply(qu2, qi->data(), wq);
			const T dot = std::fabs(quaternion::dot(qu1, wq));
			if(dot > maxDot) {
				maxDot = dot;
				std::swap(wq, eq);
				if(dot >= minFzW()) break;
			}
		}
	} else {
		throw std::runtime_error("invalid convention state");
	}
	quaternion::explement(eq, equiv);
}

////////////////////////////////////////////////////////////////////////////////////
//                                   constants                                    //
////////////////////////////////////////////////////////////////////////////////////
template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::cubic = {
	Quaternion<T>(                  T(1),                   T(0),                   T(0),                   T(0)),
	Quaternion<T>(                  T(0),                   T(0),                   T(0),                   T(1)),
	Quaternion<T>(                  T(0),                   T(1),                   T(0),                   T(0)),
	Quaternion<T>(                  T(0),                   T(0),                   T(1),                   T(0)),
	Quaternion<T>(           T(1) / T(2),            T(1) / T(2),            T(1) / T(2),            T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),           -T(1) / T(2),            T(1) / T(2),            T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),            T(1) / T(2),           -T(1) / T(2),            T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),            T(1) / T(2),            T(1) / T(2),           -T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),            T(1) / T(2),           -T(1) / T(2),           -T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),           -T(1) / T(2),            T(1) / T(2),           -T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),           -T(1) / T(2),           -T(1) / T(2),            T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),           -T(1) / T(2),           -T(1) / T(2),           -T(1) / T(2)),
	Quaternion<T>(T(1) / std::sqrt(T(2)), T(1) / std::sqrt(T(2)),                   T(0),                   T(0)),
	Quaternion<T>(T(1) / std::sqrt(T(2)),                   T(0), T(1) / std::sqrt(T(2)),                   T(0)),
	Quaternion<T>(T(1) / std::sqrt(T(2)),                   T(0),                   T(0), T(1) / std::sqrt(T(2))),
	Quaternion<T>(                  T(0), T(1) / std::sqrt(T(2)), T(1) / std::sqrt(T(2)),                   T(0)),
	Quaternion<T>(                  T(0), T(1) / std::sqrt(T(2)),                   T(0), T(1) / std::sqrt(T(2))),
	Quaternion<T>(                  T(0),                   T(0), T(1) / std::sqrt(T(2)), T(1) / std::sqrt(T(2))),
	Quaternion<T>(T(1) / std::sqrt(T(2)),-T(1) / std::sqrt(T(2)),                   T(0),                   T(0)),
	Quaternion<T>(T(1) / std::sqrt(T(2)),                   T(0),-T(1) / std::sqrt(T(2)),                   T(0)),
	Quaternion<T>(T(1) / std::sqrt(T(2)),                   T(0),                   T(0),-T(1) / std::sqrt(T(2))),
	Quaternion<T>(                  T(0),-T(1) / std::sqrt(T(2)), T(1) / std::sqrt(T(2)),                   T(0)),
	Quaternion<T>(                  T(0),-T(1) / std::sqrt(T(2)),                   T(0), T(1) / std::sqrt(T(2))),
	Quaternion<T>(                  T(0),                   T(0),-T(1) / std::sqrt(T(2)), T(1) / std::sqrt(T(2))),
};

template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::hexagonal = {
	Quaternion<T>(                  T(1),                   T(0),                   T(0),                   T(0)),
	Quaternion<T>(                  T(0),                   T(0),                   T(0),                   T(1)),
	Quaternion<T>(           T(1) / T(2),                   T(0),                   T(0), std::sqrt(T(3)) / T(2)),
	Quaternion<T>(std::sqrt(T(3)) / T(2),                   T(0),                   T(0),            T(1) / T(2)),
	Quaternion<T>(           T(1) / T(2),                   T(0),                   T(0),-std::sqrt(T(3)) / T(2)),
	Quaternion<T>(std::sqrt(T(3)) / T(2),                   T(0),                   T(0),           -T(1) / T(2)),
	Quaternion<T>(                  T(0),                   T(1),                   T(0),                   T(0)),
	Quaternion<T>(                  T(0),                   T(0),                   T(1),                   T(0)),
	Quaternion<T>(                  T(0),            T(1) / T(2), std::sqrt(T(3)) / T(2),                   T(0)),
	Quaternion<T>(                  T(0), std::sqrt(T(3)) / T(2),            T(1) / T(2),                   T(0)),
	Quaternion<T>(                  T(0),           -T(1) / T(2), std::sqrt(T(3)) / T(2),                   T(0)),
	Quaternion<T>(                  T(0),-std::sqrt(T(3)) / T(2),            T(1) / T(2),                   T(0)),
};

template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::tetragonal = {
	Quaternion<T>(                  T(1),                   T(0),                   T(0),                   T(0)),
	Quaternion<T>(                  T(0),                   T(0),                   T(0),                   T(1)),
	Quaternion<T>(T(1) / std::sqrt(T(2)),                   T(0),                   T(0), T(1) / std::sqrt(T(2))),
	Quaternion<T>(T(1) / std::sqrt(T(2)),                   T(0),                   T(0),-T(1) / std::sqrt(T(2))),
	Quaternion<T>(                  T(0),                   T(1),                   T(0),                   T(0)),
	Quaternion<T>(                  T(0),                   T(0),                   T(1),                   T(0)),
	Quaternion<T>(                  T(0), T(1) / std::sqrt(T(2)), T(1) / std::sqrt(T(2)),                   T(0)),
	Quaternion<T>(                  T(0),-T(1) / std::sqrt(T(2)), T(1) / std::sqrt(T(2)),                   T(0))
};

template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::trigonal = {
	Quaternion<T>(                  T(1),                   T(0),                   T(0),                   T(0)),
	Quaternion<T>(           T(1) / T(2),                   T(0),                   T(0), std::sqrt(T(3)) / T(2)),
	Quaternion<T>(           T(1) / T(2),                   T(0),                   T(0),-std::sqrt(T(3)) / T(2)),
	Quaternion<T>(                  T(0),                   T(0),                   T(0),                   T(1)),
	Quaternion<T>(std::sqrt(T(3)) / T(2),                   T(0),                   T(0),            T(1) / T(2)),
	Quaternion<T>(std::sqrt(T(3)) / T(2),                   T(0),                   T(0),           -T(1) / T(2))
};

template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::orthorhombic  = std::vector<Quaternion<T> >(       symmetry::operators<T>::cubic.begin(),        symmetry::operators<T>::cubic.begin() +  4);

template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::cubiclow      = std::vector<Quaternion<T> >(       symmetry::operators<T>::cubic.begin(),        symmetry::operators<T>::cubic.begin()  + 12);
template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::hexagonallow  = std::vector<Quaternion<T> >(   symmetry::operators<T>::hexagonal.begin(),    symmetry::operators<T>::hexagonal.begin()  + 6);
template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::tetragonallow = std::vector<Quaternion<T> >(  symmetry::operators<T>::tetragonal.begin(),   symmetry::operators<T>::tetragonal.begin()  + 4);
template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::trigonallow   = std::vector<Quaternion<T> >(    symmetry::operators<T>::trigonal.begin(),     symmetry::operators<T>::trigonal.begin()  + 4);
template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::monoclinic    = std::vector<Quaternion<T> >(symmetry::operators<T>::orthorhombic.begin(), symmetry::operators<T>::orthorhombic.begin()  + 2);

template <typename T> const std::vector<Quaternion<T> > symmetry::operators<T>::triclinic = {Quaternion<T>(T(1), T(0), T(0), T(0))};

#endif//_symmetry_h_