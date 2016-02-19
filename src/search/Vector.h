#ifndef __GK_VECTOR_H__ 
#define __GK_VECTOR_H__

/**
///////////////////////////////////////////////////////////////////////////////
// Original Author: Paul Henning, from the SCF library.                      //
//                                                                           //
// Modified for GK by: Brian Jean                                            //
//                     TSM, ASCI Problem Setup Team                          //
//                     Applied Physics Division                              //
//                     Los Alamos National Laboratory                        //
//                     505.665.6374                                          //
//                     baj@lanl.gov                                          //
///////////////////////////////////////////////////////////////////////////////
*/

#include <math.h>
#include <assert.h>
#include <iostream>

namespace gk
{

    /*!
      @class Vector "Vector.h"
      @brief Represents a vector in N-dimensional space
      @tparam D Indicates the dimensionality of the Vector (this will generally be one
      of [1, 2, 3]).
      */
    template <long D> class Vector
    {
        private: 
            double m_comp[D];

        public:

            /// Default constructor - zero Vector in D-space.
            inline Vector() {
                for (long i=0;i<D;i++)
                    m_comp[i] = 0.0;
            }

            /*!
              @brief Specialized constructor for 2d Vectors.
              @param[in] xm_comp,ym_comp The (x,y) coordinate pair.
              */
            inline Vector(const double& xm_comp, 
                          const double& ym_comp) {
                assert (D==2);
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
                assert (D==3);
                m_comp[0] = xm_comp;
                m_comp[1] = ym_comp;
                m_comp[2] = zm_comp;
            }

            /// Return component @c i of the Vector.
            inline const double& operator[](const long& i) const {
                return m_comp[i]; 
            }

            /// Return component @c i of the Vector.
            inline double& operator[](const long& i) {
                return m_comp[i]; 
            }

            /// Add the Vector @c rhs to this Vector.
            Vector& operator+=(const Vector<D>& rhs) {
                for(long i = 0; i < D; i++) m_comp[i] += rhs.m_comp[i];
                return *this;
            }

            /// Subtract the Vector @c rhs from this vector.
            Vector& operator-=(const Vector<D>& rhs) {
                for(long i = 0; i < D; i++) m_comp[i] -= rhs.m_comp[i];
                return *this;
            }

            /// Scalar multiplication of this Vector by @c s.
            Vector& operator*=(const double& s) {
                for(long i = 0; i < D; i++) m_comp[i] *= s;
                return *this;
            }

            /// Scalar division of this Vector by @c s.
            Vector& operator/=(const double& s) {
                for(long i = 0; i < D; i++) m_comp[i] /= s;
                return *this;
            }

            /*!
              @brief Calculate the norm of a Vector.
              @param[in] doSqrt OPTIONAL: Return the square root of the norm, i.e. the
              magnitude of the Vector.
              */
            double norm(bool doSqrt = true) const {
                double result = 0.0;
                for(long i = 0; i < D; i++) result += (m_comp[i] * m_comp[i]);
                if(doSqrt) return sqrt(result);
                return result;
            }

            /// Convert this Vector into a unit Vector.
            void normalize() {
                double s = norm();
                *this /= s;
            }

            /// Convert this Vector into a zero Vector.
            void zero() {
                for(long i = 0; i < D; i++) m_comp[i] = 0;
            }

            /*!
              @brief Convenience method for constructing a unit Vector along a particular
              axis
              @param[in] nonZero The coordinate axis along which the Vector should point.
              */
            void axis(long nonZero) {
                zero();
                m_comp[nonZero] = 1;
            }

            /// Read in a Vector from an input stream.
            std::istream& readFromStream(std::istream& is) {
                for (long i=0;i<D;i++)
                    is >> m_comp[i];
                return is;
            }

            /// Pretty printing of a Vector to an output stream.
            std::ostream& writeToStream(std::ostream& os) {
                for (long i=0;i<D;i++)
                {
                    if (i>0) os << ' ';
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
    template<long D> inline double
        dot(const Vector<D>& a, const Vector<D>& b)
        {
            double r = 0.0;
            for(long i = 0; i < D; i++) r += a[i] * b[i];
            return r;
        }

    /// Add two vectors.
    template<long D> inline const Vector<D>
        operator+(const Vector<D>& a, const Vector<D>& b) 
        {
            return Vector<D>(a) += b;
        }

    /// Subtract two vectors.
    template<long D> inline const Vector<D>
        operator-(const Vector<D>& a, const Vector<D>& b) 
        {
            return Vector<D>(a) -= b;
        }

    /// Multiply a vector by a scalar, @f$ s \vec{a}@f$.
    template<long D> inline const Vector<D>
        operator*(const Vector<D>& a, const double& s) 
        {
            return Vector<D>(a) *= s;
        }

    /// Multiply a vector by a scalar, @f$ s \vec{a}@f$. 
    template<long D> inline const Vector<D>
        operator*(const double& s, const Vector<D>& a) 
        {
            return Vector<D>(a) *= s;
        }

    /// Divide a vector by a scalar, @f$ \frac{1}{s} \vec{a}@f$.
    template<long D> inline const Vector<D>
        operator/(const Vector<D>& a, const double& s) 
        {
            return Vector<D>(a) /= s;
        }

    /// Pretty printing of a Vector to an output stream.
    template<long D> inline std::ostream&
        operator<<(std::ostream& os, const Vector<D>& v) 
        {
            return v.writeToStream(os);
        }

    /// Read in a Vector from an input stream.
    template<long D> inline std::istream&
        operator>>(std::istream& is, Vector<D>& v) 
        {
            return v.readFromStream(is);
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
            for (long i=1;i<D;i++)
                if (max < v[i]) {
                    max = v[i];
                    icomp = i;
                }
            return max;
        }

}

#endif
