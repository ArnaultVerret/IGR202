#ifndef _QUATERNION_HPP_
#define _QUATERNION_HPP_

#include "typedefs.hpp"
#include "Vector3.hpp"
#include "Matrix3x3.hpp"
#include <cmath>

class Quaternion {

public:

  // init by values
  Quaternion(Real a, Real b, Real c, Real d){
    v0 = a;
    v1 = b;
    v2 = c;
    v3 = d;
  }

  // init by rotation matrix
  Quaternion(Mat3f &m){
    Real r = sqrt(1 + m.v00 + m.v11 + m.v22) / 2;
    v0 = r;
    v1 = (m.v21 - m.v12) / 4 / r;
    v2 = (m.v02 - m.v20) / 4 / r;
    v3 = (m.v10 - m.v01) / 4 / r;
  }

  Mat3f ToMat(){
    return Mat3f(
      1 - 2*v2*v2 - 2*v3*v3, 2*v1*v2 - 2*v0*v3, 2*v1*v3 + 2*v0*v2,
      2*v1*v2 + 2*v0*v3, 1 - 2*v1*v1 - 2*v3*v3, 2*v2*v3 - 2*v0*v1,
      2*v1*v3 - 2*v0*v2, 2*v2*v3 + 2*v0*v1, 1 - 2*v1*v1 - 2*v2*v2
    );
  }

  Quaternion& inverse() { return *this = inversed();}
  Quaternion inversed() {
    return Quaternion(v0, -v1, -v2, -v3) / normSqr();
  }
  Quaternion& normalize() { return *this = normalized();}
  Quaternion normalized() {
    return Quaternion(*this) / norm();
  }

  Real normSqr() { // Return the norm squared
    return v0*v0 + v1*v1 + v2*v2 + v3*v3;
  }
  Real norm() { // return the norm
    return sqrt(normSqr());
  }
  Real getRealPart(){
    return v0;
  }
  Vec3f getUnitVector(){
    return Vec3f(v1, v2, v3).normalize();
  }

  // assignment operators
  Quaternion& operator+=(const Quaternion &q) {
    v0 += q.v0; v1 += q.v1; v2 += q.v2; v3 += q.v3;
    return *this;
  }
  Quaternion& operator-=(const Quaternion &q) {
    v0 -= q.v0; v1 -= q.v1; v2 -= q.v2; v3 -= q.v3;
    return *this;
  }
  Quaternion& operator*=(const Real &r) {
    v0 *= r; v1 *= r; v2 *= r; v3 *= r;
    return *this;
  }

  Quaternion& operator/=(const Real &r) {
    v0 /= r; v1 /= r; v2 /= r; v3 /= r;
    return *this;
  }
  Quaternion operator*=(const Quaternion &q) {
    return Quaternion(*this) * q;
  }
  Quaternion operator/=(Quaternion &q) {
    return Quaternion(*this) * q.inversed();
  }

  // binary operators
  Quaternion operator+(const Quaternion &q) const { return Quaternion(*this)+=q; }
  Quaternion operator-(const Quaternion &q) const { return Quaternion(*this)-=q; }
  Quaternion operator*(const Quaternion &q) const {
    return Quaternion(
      q.v0*v0 - q.v1*v1 - q.v2*v2 - q.v3*v3,
      q.v0*v1 + q.v1*v0 - q.v2*v3 + q.v3*v2,
      q.v0*v2 + q.v1*v3 + q.v2*v0 - q.v3*v1,
      q.v0*v3 - q.v1*v2 + q.v2*v1 - q.v3*v0
    );
  }
  Quaternion operator/(Quaternion &q) const {
    return Quaternion(*this)/=q;
  }
  Quaternion operator*(const Real &r) const { return Quaternion(*this)*=r; }
  Quaternion operator/(const Real &r) const { return Quaternion(*this)/=r; }



  Real v0, v1, v2, v3;
};


#endif  /* _MATRIX3X3_HPP_ */
