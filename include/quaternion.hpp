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
#ifndef _quaternion_h_
#define _quaternion_h_

#include <algorithm>//transform
#include <functional>//negate
#include <numeric>//inner_product
#include <cmath>//sqrt

#include "rotations.hpp"


////////////////////////////////////////////////////////////////////////////////////
//       underlying quaternion math functions (operate on pointer to wxyz)        //
////////////////////////////////////////////////////////////////////////////////////
namespace quaternion {
	//operations between a quaternion and a scalar
	template <typename T> void scalarAdd(T const * const qIn, const T s, T * const qOut) {std::transform(qIn, qIn+4, qOut, [&s](const T i){return i+s;});}
	template <typename T> void scalarSub(T const * const qIn, const T s, T * const qOut) {std::transform(qIn, qIn+4, qOut, [&s](const T i){return i-s;});}
	template <typename T> void scalarMul(T const * const qIn, const T s, T * const qOut) {std::transform(qIn, qIn+4, qOut, [&s](const T i){return i*s;});}
	template <typename T> void scalarDiv(T const * const qIn, const T s, T * const qOut) {std::transform(qIn, qIn+4, qOut, [&s](const T i){return i/s;});}

	//operations between 2 quaternions
	template <typename T> void elementAdd(T const * const qIn1, T const * const qIn2, T * const qOut) {std::transform(qIn1, qIn1+4, qIn2, qOut, [](const T i, const T j){return i+j;});}
	template <typename T> void elementSub(T const * const qIn1, T const * const qIn2, T * const qOut) {std::transform(qIn1, qIn1+4, qIn2, qOut, [](const T i, const T j){return i-j;});}

	//dot product / magnitude
	template <typename T> T dot(T const * const q1, T const * const q2) {return std::inner_product(q1, q1+4, q2, T(0));}
	template <typename T> T mag2(T const * const q) {return dot(q, q);}
	template <typename T> T mag(T const * const q) {return std::sqrt(mag2(q));}
	template <typename T> void normalize(T const * const qIn, T * const qOut) {scalarDiv(qIn, mag(qIn), qOut);}

	//conjugate, inverse, negative, abs, and explement
	template <typename T> void conjugate(T const * const qIn, T * const qOut) {
		qOut[0] = qIn[0];
		std::transform(qIn+1, qIn+4, qOut+1, std::negate<T>());
	}
	template <typename T> void inverse(T const * const qIn, T * const qOut) {
		conjugate(qIn, qOut);
		scalarDiv(qOut, mag2(qOut), qOut);
	}
	template <typename T> void negate(T const * const qIn, T * const qOut) {std::transform(qIn, qIn+4, qOut, std::negate<T>());}
	template <typename T> void abs(T const * const qIn, T * const qOut) {std::transform(qIn, qIn+4, qOut, static_cast<T(*)(T)>(&std::fabs));}
	template <typename T> void explement(T const * const qIn, T * const qOut) {if(qIn[0] < 0) negate(qIn, qOut); else std::copy(qIn, qIn+4, qOut);}//complementary angles sum to 90, supplementary to 180, and explementary to 360

	//quaternion multiplication
	template <typename T> void multiply(T const * const qIn1, T const * const qIn2, T * const qOut) {
		const T w = qIn1[0]*qIn2[0] - qIn1[1]*qIn2[1] - qIn1[2]*qIn2[2] - qIn1[3]*qIn2[3];
		const T x = qIn1[0]*qIn2[1] + qIn1[1]*qIn2[0] + rotations::Constants<T>::convention * (qIn1[2]*qIn2[3] - qIn1[3]*qIn2[2]);
		const T y = qIn1[0]*qIn2[2] + qIn1[2]*qIn2[0] + rotations::Constants<T>::convention * (qIn1[3]*qIn2[1] - qIn1[1]*qIn2[3]);
		const T z = qIn1[0]*qIn2[3] + qIn1[3]*qIn2[0] + rotations::Constants<T>::convention * (qIn1[1]*qIn2[2] - qIn1[2]*qIn2[1]);
		qOut[0] = w;
		qOut[1] = x;
		qOut[2] = y;
		qOut[3] = z;
	}

	//quaternion division: (qIn1 * qIn2.conjugate()).normalize()
	template <typename T> void divide(T const * const qIn1, T const * const qIn2, T * const qOut) {
		conjugate(qIn2, qOut);
		multiply(qIn1, qOut, qOut);
		normalize(qOut, qOut);
	}

	//rotate a vector by a quaternion (q * v * q.conjugate())
	template <typename T> void rotateVector(T const * const q, T const * const vIn, T * const vOut) {
		//q * v
		const T w = -(q[1]*vIn[0] + q[2]*vIn[1] + q[3]*vIn[2]);
		const T x = q[0]*vIn[0] + rotations::Constants<T>::convention * (q[2]*vIn[2] - q[3]*vIn[1]);
		const T y = q[0]*vIn[1] + rotations::Constants<T>::convention * (q[3]*vIn[0] - q[1]*vIn[2]);
		const T z = q[0]*vIn[2] + rotations::Constants<T>::convention * (q[1]*vIn[1] - q[2]*vIn[0]);

		//(q * v) * q.conjugate()
		vOut[0] = -w*q[1] + x*q[0] - rotations::Constants<T>::convention * (y*q[3] - z*q[2]);
		vOut[1] = -w*q[2] + y*q[0] - rotations::Constants<T>::convention * (z*q[1] - x*q[3]);
		vOut[2] = -w*q[3] + z*q[0] - rotations::Constants<T>::convention * (x*q[2] - y*q[1]);
	}
}

////////////////////////////////////////////////////////////////////////////////////
//    quaternion wrapper class (may be slightly less efficient due to copying)    //
////////////////////////////////////////////////////////////////////////////////////
//POD struct instead of class + packed check allows casting from array of doubles to array of quats ([w1,x1,y1,z1,w2,x2,y2,z2,w3,x3,y3,z3] -> [q1,q2,q3])
template<typename T>
struct Quaternion {
	T w, x, y, z;
	Quaternion() {
		static_assert(sizeof(Quaternion) == 4 * sizeof(T), "Quaternion struct must be packed");
		static_assert(std::is_floating_point<T>::value, "Quaternion must be templated on floating point type");
	}
	Quaternion(const T vw, const T vx, const T vy, const T vz) : w(vw), x(vx), y(vy), z(vz) {}
	static Quaternion Zero() {return Quaternion(T(0), T(0), T(0), T(0));}
	static Quaternion Identity() {return Quaternion(T(1), T(0), T(0), T(0));}

	T      *const data()       {return (T      *const)this;}
	T const*const data() const {return (T const*const)this;}

	//quaternion / scalar arithmetic
	Quaternion& operator+=(const T& s) {quaternion::scalarAdd(data(), s, data()); return *this;}
	Quaternion& operator-=(const T& s) {quaternion::scalarSub(data(), s, data()); return *this;}
	Quaternion& operator*=(const T& s) {quaternion::scalarMul(data(), s, data()); return *this;}
	Quaternion& operator/=(const T& s) {quaternion::scalarDiv(data(), s, data()); return *this;}
	Quaternion operator+(const T& s) const {Quaternion q; quaternion::scalarAdd(data(), s, q.data()); return q;}
	Quaternion operator-(const T& s) const {Quaternion q; quaternion::scalarSub(data(), s, q.data()); return q;}
	Quaternion operator*(const T& s) const {Quaternion q; quaternion::scalarMul(data(), s, q.data()); return q;}
	Quaternion operator/(const T& s) const {Quaternion q; quaternion::scalarDiv(data(), s, q.data()); return q;}

	//quaternion / quaternion arithmetic
	Quaternion& operator+=(const Quaternion& q) {quaternion::elementAdd(data(), q.data(), data()); return *this;}
	Quaternion& operator-=(const Quaternion& q) {quaternion::elementSub(data(), q.data(), data()); return *this;}
	Quaternion& operator*=(const Quaternion& q) {quaternion::  multiply(data(), q.data(), data()); return *this;}
	Quaternion& operator/=(const Quaternion& q) {quaternion::    divide(data(), q.data(), data()); return *this;}
	Quaternion operator+(const Quaternion& q) const {Quaternion r; quaternion::elementAdd(data(), q.data(), r.data()); return r;}
	Quaternion operator-(const Quaternion& q) const {Quaternion r; quaternion::elementSub(data(), q.data(), r.data()); return r;}
	Quaternion operator*(const Quaternion& q) const {Quaternion r; quaternion::  multiply(data(), q.data(), r.data()); return r;}
	Quaternion operator/(const Quaternion& q) const {Quaternion r; quaternion::    divide(data(), q.data(), r.data()); return r;}

	//urnary negate
	Quaternion operator-() const {return negate();}

	//comparison operators
	inline bool operator< (const Quaternion& rhs) const {
		if     (w > rhs.w) return true; //first sort by smallest angle
		else if(w < rhs.w) return false;
		//rotation angles are equal
		if     (z < rhs.z) return true; //next sort by z axis
		else if(z > rhs.z) return false;
		//z components are equal
		if     (y < rhs.y) return true; //next sort by y axis
		else if(y > rhs.y) return false;
		//y components are equal
		return  x < rhs.x;              //finally sort by x axis
	}
	inline bool operator> (const Quaternion& rhs) const {return rhs < *this;}
	inline bool operator<=(const Quaternion& rhs) const {return !operator>(rhs);}
	inline bool operator>=(const Quaternion& rhs) const {return !operator<(rhs);}
	inline bool operator==(const Quaternion& rhs) const {return w == rhs.w && x == rhs.x && y == rhs.y && z == rhs.z;}
	inline bool operator!=(const Quaternion& rhs) const {return !operator==(rhs);}

	//dot product / magnitude
	T dot(const Quaternion& q) const {return quaternion::dot(data(), q.data());}
	T mag () const {return quaternion::mag(data());}
	T mag2() const {return quaternion::mag2(data());}

	//normalize, conjugate, inverse, abs, and explement
	Quaternion normalize() const {Quaternion q; quaternion::normalize(data(), q.data()); return q;}
	Quaternion conjugate() const {Quaternion q; quaternion::conjugate(data(), q.data()); return q;}
	Quaternion   inverse() const {Quaternion q; quaternion::  inverse(data(), q.data()); return q;}
	Quaternion       abs() const {Quaternion q; quaternion::      abs(data(), q.data()); return q;}
	Quaternion    negate() const {Quaternion q; quaternion::   negate(data(), q.data()); return q;}
	Quaternion explement() const {Quaternion q; quaternion::explement(data(), q.data()); return q;}

	//vector rotation
	void rotateVector(T const * const vIn, T * const vOut) const {quaternion::rotateVector(data(), vIn, vOut);}
};

#endif//_quaternion_h_