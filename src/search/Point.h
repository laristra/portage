///////////////////////////////////////////////////////////////////////////////
// Original Author: Paul Henning                                             //
//                  TSM, ASCI Problem Setup Team                             //
//                  Applied Physics Division                                 //
//                  Los Alamos National Laboratory                           //
//                  phenning@lanl.gov                                        //
//                                                                           //
// Modified for GK by: Brian Jean                                            //
//                     TSM, ASCI Problem Setup Team                          //
//                     Applied Physics Division                              //
//                     Los Alamos National Laboratory                        //
//                     505.665.6374                                          //
//                     baj@lanl.gov                                          //
///////////////////////////////////////////////////////////////////////////////

#ifndef __GK_POINT_H__
#define __GK_POINT_H__

#include <assert.h>
#include <iostream>
#include <vector>
#include "Vector.h"

namespace gk
{

const long X = 0;
const long Y = 1;
const long Z = 2;

/*!
  @class Point "Point.h"
  @brief Represents a point in an N-dimensional space.
  @tparam D Indicates the dimensionality of the Point (this will generally be one
  of [1,2,3]).
 */
template <long D> class Point
{
  private: 
    double m_loc[D];

  public:

    /// Default constructor - Point at origin in D-space.
    inline Point() {
      for (long i=0;i<D;i++)
        m_loc[i] = 0.0;
    }

    /*!
      @brief Specialized constructor from a std::vector of arbitary size.
      @param[in] v std::vector of coordinates.
     */
    inline Point(const std::vector<double> &v)
      {
          assert (v.size() == D);
          for (long i=0;i<D;i++)
            m_loc[i] = v[i];
      }

    /*!
      @brief Specialized constructor for Points in 3d.
      @param[in] x,y,z The (x,y,z) coordinates of this Point.
     */
    inline Point(const double& x, const double& y, const double& z)
      {
        assert (D==3);
        m_loc[0] = x;
        m_loc[1] = y;
        m_loc[2] = z;
      }

    /*!
      @brief Specialized constructor for Points in 2d.
      @param[in] x,y The (x,y) coordinates of this Point.
     */
    inline Point(const double& x, const double& y)
      {
        assert (D==2);
        m_loc[0] = x;
        m_loc[1] = y;
      }

    /// Convert a Vector to a Point.
    explicit Point(const Vector<D>& v) {
      for(long i = 0; i < D; i++) 
        m_loc[i] = v[i];
    }

    /// Copy constructor.
    inline Point(const Point<D>& rhs)
      {
          for (long i=0;i<D;i++)
            m_loc[i] = rhs[i];
      }

    /// Return component @c i of the Point.
    inline const double& operator[](const long& i) const
      {
        return m_loc[i]; 
      }

    /// Return component @c i of the Point.
    inline double& operator[](const long& i)
      {
        return m_loc[i]; 
      }

    /// Translate this Point along the Vector @c v.
    inline Point<D>& operator+=(const Vector<D>& v) 
      {
        for(long i = 0; i < D; i++) m_loc[i] += v[i];
        return *this;
      }

    /// Read in the coordinates a Point from an input stream.
    std::istream& readFromStream(std::istream& is)
      {
        for (long i=0;i<D;i++)
          is >> m_loc[i];
        return is;
      }

    /// Pretty printing of the coordinates of a Point to an output stream.
    std::ostream& writeToStream(std::ostream& os) const
      {
        for (long i=0;i<D;i++)
          {
            if (i>0) os << ' ';
            os << m_loc[i];
          }
        return os;
      }

    /// Convert a Point to a Vector from the origin to the coordinates of the Point.
    inline Vector<D> asV()
      {
        Vector<D> v;
        for (long i=0;i<D;i++)
          v[i] = m_loc[i];
        return v;
      }
    
    /// Convert a Point to a Vector from the origin to the coordinates of the Point.
    inline Vector<D> asV() const
      {
        Vector<D> v;
        for (long i=0;i<D;i++)
          v[i] = m_loc[i];
        return v;
      }
};

/// Alias for creating a Point in 3d.
typedef Point<3> Point3;
/// Alias for creating a Point in 2d.
typedef Point<2> Point2;

template <long D> inline std::ostream&
operator<<(std::ostream& os, const Point<D>& p)
{
  return p.writeToStream(os);
}


template <long D> inline const Point<D>
operator+(const Point<D>& p, const Vector<D>& v) 
  {
    return Point<D>(p)+=v;
  }

template <long D> inline const Vector<D>
operator-(const Point<D>& p1, const Point<D>& p2) 
  {
    Vector<D> v;
    for(long i = 0; i < D; i++) v[i] = p1[i] - p2[i];
    return v;
  }


template <long D> inline bool
approxEq(const Point<D>& p1, const Point<D>& p2, double tol = 1.0e-8) 
  {
    // This uses the L_infty norm to avoid nearby subtractions
    for(long i = 0; i < D; i++)
      if(p1[i] < (p2[i] - tol) || p1[i] > (p2[i] + tol))
        return false;
    return true;
  }

template <long D> bool
operator<(const Point<D>& p1, const Point<D>& p2) {
  if (approxEq (p1,p2))
    return false;
  else 
    for (long i=0;i<D;++i) {
      if (p1[i]<p2[i])
        return true;
      else if (p2[i]<p1[i])
        return false;
    }
  return false;
}


inline Point2 ToCylindrical(const Point3& p)
{
  Point2 result;

  result[0] = sqrt(p[0]*p[0]+p[1]*p[1]);
  result[1] = p[2]; 
  return result;
}

inline const Point3
createP3(double x, double y, double z) 
{
  Point3 p; 

  p[0] = x;
  p[1] = y; 
  p[2] = z; 

  return p;
}

inline const Point2
createP2(double x, double y) 
{
  Point2 p; 

  p[0] = x; 
  p[1] = y; 

  return p;
}

#if 0
template <typename D> 
class AbsPoint : public Point
{

};

template <typename D>
class RelPoint : public Point
{

};
#endif

}
#endif
