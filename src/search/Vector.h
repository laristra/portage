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

/** Represent a vector in N-dimensional space.
*   
*   The template argument <D> indicates the dimensionality of the Vector
*   (this will generally be (1, 2, or 3).  
*/
template <long D> class Vector
{
  private: 
    double m_comp[D];

  public:

    inline Vector() {
      for (long i=0;i<D;i++)
        m_comp[i] = 0.0;
    }

    inline Vector(const double& xm_comp, 
                  const double& ym_comp)
      {
        assert (D==2);
        m_comp[0] = xm_comp;
        m_comp[1] = ym_comp;
      }

    inline Vector(const double& xm_comp, 
                  const double& ym_comp, 
                  const double& zm_comp)
      {
        assert (D==3);
        m_comp[0] = xm_comp;
        m_comp[1] = ym_comp;
        m_comp[2] = zm_comp;
      }

    inline const double& operator[](const long& i) const
      {
        return m_comp[i]; 
      }

    inline double& operator[](const long& i)
      {
        return m_comp[i]; 
      }

    ///
    Vector& operator+=(const Vector<D>& rhs) 
      {
        for(long i = 0; i < D; i++) m_comp[i] += rhs.m_comp[i];
        return *this;
      }

    ///
    Vector& operator-=(const Vector<D>& rhs) 
      {
        for(long i = 0; i < D; i++) m_comp[i] -= rhs.m_comp[i];
        return *this;
      }

    Vector& operator*=(const double& s) 
      {
        for(long i = 0; i < D; i++) m_comp[i] *= s;
        return *this;
      }

    Vector& operator/=(const double& s) 
      {
        for(long i = 0; i < D; i++) m_comp[i] /= s;
        return *this;
      }

    double norm(bool doSqrt = true) const 
      {
        double result = 0.0;
        for(long i = 0; i < D; i++) result += (m_comp[i] * m_comp[i]);
        if(doSqrt) return sqrt(result);
        return result;
      }

    ///
    void normalize() 
      {
        double s = norm();
        *this /= s;
      }

    ///
    void zero() 
      {
        for(long i = 0; i < D; i++) m_comp[i] = 0;
      }

    void axis(long nonZero) 
      {
        zero();
        m_comp[nonZero] = 1;
      }

    ///
    std::istream& readFromStream(std::istream& is)
      {
        for (long i=0;i<D;i++)
          is >> m_comp[i];
        return is;
      }

    ///
    std::ostream& writeToStream(std::ostream& os)
      {
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

template<long D> inline double
dot(const Vector<D>& a, const Vector<D>& b)
  {
    double r = 0.0;
    for(long i = 0; i < D; i++) r += a[i] * b[i];
    return r;
  }

template<long D> inline const Vector<D>
operator+(const Vector<D>& a, const Vector<D>& b) 
  {
    return Vector<D>(a) += b;
  }

template<long D> inline const Vector<D>
operator-(const Vector<D>& a, const Vector<D>& b) 
  {
    return Vector<D>(a) -= b;
  }

template<long D> inline const Vector<D>
operator*(const Vector<D>& a, const double& s) 
  {
    return Vector<D>(a) *= s;
  }

template<long D> inline const Vector<D>
operator*(const double& s, const Vector<D>& a) 
  {
    return Vector<D>(a) *= s;
  }

template<long D> inline const Vector<D>
operator/(const Vector<D>& a, const double& s) 
  {
    return Vector<D>(a) /= s;
  }

template<long D> inline std::ostream&
operator<<(std::ostream& os, const Vector<D>& v) 
  {
    return v.writeToStream(os);
  }

template<long D> inline std::istream&
operator>>(std::istream& is, Vector<D>& v) 
  {
    return v.readFromStream(is);
  }

/// Cross product operator for two vectors.
inline const Vector<3> cross(const Vector<3>& a, const Vector<3>& b)
  {
    Vector<3> r;
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
    return r;
  }

template<long D>
inline double MaxComponent(const Vector<D>& v, long& icomp)
{
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
