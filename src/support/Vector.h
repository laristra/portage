/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

///////////////////////////////////////////////////////////////////////////////
// Copyright 2016 Los Alamos National Laboratory                             //
//                                                                           //
// Original Author: Paul Henning, from the SCF library.                      //
//                                                                           //
// Modified for GK by: Brian Jean                                            //
//                     TSM, ASCI Problem Setup Team                          //
//                     Applied Physics Division                              //
//                     Los Alamos National Laboratory                        //
//                     505.665.6374                                          //
//                     baj@lanl.gov                                          //
//                                                                           //
// Modified for Portage by: Rao Garimella, rao@lanl.gov                      //
///////////////////////////////////////////////////////////////////////////////

#ifndef SRC_SUPPORT_VECTOR_H_
#define SRC_SUPPORT_VECTOR_H_

#include <math.h>
#include <assert.h>
#include <iostream>
#include <vector>

namespace Portage {

/*!
  @class Vector "Vector.h"
  @brief Represents a vector in N-dimensional space
  @tparam D Indicates the dimensionality of the Vector (this will generally be one
  of [1, 2, 3]).
*/
template <long D> class Vector {
 private:
  double m_comp[D];

 public:

  /// Default constructor - zero Vector in D-space.
  inline Vector() {
    for (int i = 0; i < D; i++)
      m_comp[i] = 0.0;
  }

  /*!
    @brief Specialized constructor for 2d Vectors.
    @param[in] xm_comp,ym_comp The (x,y) coordinate pair.
  */
  inline Vector(const double& xm_comp,
                const double& ym_comp) {
    assert(D == 2);
    m_comp[0] = xm_comp;
    m_comp[1] = ym_comp;
  }

  /*!
    @brief Specialized constructor for 3d Vectors.
    @param[in] xm_comp,ym_comp,zm_comp The (x,y,z) coordinate triple.
  */
  inline Vector(const double& xm_comp,
                const double& ym_comp,
                const double& zm_comp) {
    assert(D == 3);
    m_comp[0] = xm_comp;
    m_comp[1] = ym_comp;
    m_comp[2] = zm_comp;
  }

  /*!
    @brief Constructor from a std:vector
  */
  inline
  Vector(std::vector<double> const& invec) {
    assert(D == invec.size());
    for (int i = 0; i < D; i++) m_comp[i] = invec[i];
  }

  /// Return component @c i of the Vector.
  inline const double& operator[](const int& i) const {
    return m_comp[i];
  }

  /// Return component @c i of the Vector.
  inline double& operator[](const int& i) {
    return m_comp[i];
  }

  /// Add the Vector @c rhs to this Vector.
  Vector& operator+=(const Vector<D>& rhs) {
    for (int i = 0; i < D; i++) m_comp[i] += rhs.m_comp[i];
    return *this;
  }

  /// Subtract the Vector @c rhs from this vector.
  Vector& operator-=(const Vector<D>& rhs) {
    for (int i = 0; i < D; i++) m_comp[i] -= rhs.m_comp[i];
    return *this;
  }

  /// Scalar multiplication of this Vector by @c s.
  Vector& operator*=(const double& s) {
    for (int i = 0; i < D; i++) m_comp[i] *= s;
    return *this;
  }

  /// Scalar division of this Vector by @c s.
  Vector& operator/=(const double& s) {
    for (int i = 0; i < D; i++) m_comp[i] /= s;
    return *this;
  }

  /*!
    @brief Calculate the norm of a Vector.
    @param[in] doSqrt OPTIONAL: Return the square root of the norm, i.e. the
    magnitude of the Vector.
  */
  double norm(bool doSqrt = true) const {
    double result = 0.0;
    for (int i = 0; i < D; i++) result += (m_comp[i] * m_comp[i]);
    if (doSqrt) return sqrt(result);
    return result;
  }

  /// Convert this Vector into a unit Vector.
  void normalize() {
    double s = norm();
    *this /= s;
  }

  /// Convert this Vector into a zero Vector.
  void zero() {
    for (int i = 0; i < D; i++) m_comp[i] = 0;
  }

  /*!
    @brief Convenience method for constructing a unit Vector along a particular
    axis
    @param[in] nonZero The coordinate axis along which the Vector should point.
  */
  void axis(int nonZero) {
    zero();
    m_comp[nonZero] = 1;
  }

  /// Read in a Vector from an input stream.
  std::istream& readFromStream(std::istream& is) {
    for (int i = 0; i < D; i++)
      is >> m_comp[i];
    return is;
  }

  /// Pretty printing of a Vector to an output stream.
  std::ostream& writeToStream(std::ostream& os) {
    for (int i = 0; i < D; i++) {
      if (i > 0) os << ' ';
      os << m_comp[i];
    }
    return os;
  }
};

/// Alias for creating a 3D vector
typedef Vector<3> Vector3;

/// Alias for creating a 3D vector
typedef Vector<2> Vector2;

/// Dot product of two vectors, @f$\vec{a} \cdot \vec{b}@f$.
template<long D>
inline
double dot(const Vector<D>& a, const Vector<D>& b) {
  double r = 0.0;
  for (int i = 0; i < D; i++) r += a[i] * b[i];
  return r;
}

/// Add two vectors.
template<long D>
inline
const Vector<D> operator+(const Vector<D>& a, const Vector<D>& b) {
  return Vector<D>(a) += b;
}

/// Subtract two vectors.
template<long D>
inline
const Vector<D> operator-(const Vector<D>& a, const Vector<D>& b) {
  return Vector<D>(a) -= b;
}

/// Multiply a vector by a scalar, @f$ s \vec{a}@f$.
template<long D>
inline
const Vector<D> operator*(const Vector<D>& a, const double& s) {
  return Vector<D>(a) *= s;
}

/// Multiply a vector by a scalar, @f$ s \vec{a}@f$.
template<long D>
inline
const Vector<D> operator*(const double& s, const Vector<D>& a) {
  return Vector<D>(a) *= s;
}

/// Divide a vector by a scalar, @f$ \frac{1}{s} \vec{a}@f$.
template<long D>
inline
const Vector<D> operator/(const Vector<D>& a, const double& s) {
  return Vector<D>(a) /= s;
}

/// Pretty printing of a Vector to an output stream.
template<long D>
inline
std::ostream& operator<<(std::ostream& os, const Vector<D>& v) {
  return v.writeToStream(os);
}

/// Read in a Vector from an input stream.
template<long D> inline std::istream&
operator>>(std::istream& is, Vector<D>& v) {
  return v.readFromStream(is);
}

/// Cross product operator for two 2d vectors,  @f$\vec{a} \times \vec{b}@f$.

inline double cross(const Vector<2>& a, const Vector<2>& b) {
  return (a[0] * b[1] - a[1] * b[0]);
}

/// Cross product operator for two 3d vectors, @f$\vec{a} \times \vec{b}@f$.
inline const Vector<3> cross(const Vector<3>& a, const Vector<3>& b) {
  Vector<3> r;
  r[0] = a[1] * b[2] - a[2] * b[1];
  r[1] = a[2] * b[0] - a[0] * b[2];
  r[2] = a[0] * b[1] - a[1] * b[0];
  return r;
}

/*!
  @brief Obtain the value and index of the maximum component of a Vector.
  @param[in] v The input Vector.
  @param[in,out] icomp The index of the maximum component of @c v.
  @return The maximum component of @c v.
*/
template<long D>
inline double MaxComponent(const Vector<D>& v, long& icomp) {
  double max = v[0];
  icomp = 0;
  for (int i = 1; i < D; i++)
    if (max < v[i]) {
      max = v[i];
      icomp = i;
    }
  return max;
}

}  // namespace Portage


#endif  // SRC_SUPPORT_VECTOR_H_
