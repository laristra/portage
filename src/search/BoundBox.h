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

#ifndef __GK_BBOX__
#define __GK_BBOX__

#include <iostream>
#include <assert.h>
#include <limits>
#include "Point.h"

namespace gk {

  const double REAL_MAX = std::numeric_limits<double>::max();
  const double REAL_MIN = std::numeric_limits<double>::min();

  // An isothetic (axis-aligned) bounding box

  template<long D> class IsotheticBBox 
  {
    public:
      IsotheticBBox() { clearBox(); }
      // Size functions
      void clear() { clearBox(); }
      bool empty() const { return min[0] == REAL_MAX; }

      // Functions to add new elements to the box
      void add(const Point<D>& p) {
        if(empty()) {
          min = p; max = p;
        } 
        else {
          for(long i = 0; i < D; i++) {
            if(p[i] < min[i]) min[i] = p[i];
            else if(p[i] > max[i]) max[i] = p[i];
          }
        }
      }

      void add(const IsotheticBBox<D>& box) {
        add(box.getMin());
        add(box.getMax());
      }

      // Expand the box by an additive constant
      void bulge(double epsilon) {
        for(long i = 0; i < D; i++) {
          max[i] += epsilon;
          min[i] -= epsilon;
        }
      }

      // Calculate box feature points
      Point<D> center() const { return Point<D>((min.asV() + max.asV())/2); }
      double center(long axis) const { return (min[axis] + max[axis])/2; }

      double radius(bool doSqrt = true) const { 
        const Vector<D> r((max-min)/2);
        return r.norm(doSqrt); 
      }

      const Point<D>& getMin() const { return min; }
      double getMin(long i) const { assert(i < D); return min[i]; }

      const Point<D>& getMax() const { return max; }
      double getMax(long i) const { assert(i < D); return max[i]; }

      long longAxis() const {
        // Find a longest axis
        const Vector<D> diff = max-min;
        double maxVal = diff[0];
        long maxIdx = 0;
        for(long i = 1; i < D; i++) 
          if(diff[i] > maxVal) {
            maxVal = diff[i];
            maxIdx = i;
          }
        return maxIdx;
      }

      double volume() const { 
        double result = 1.0;
        for(long i = 0; i < D; i++)  result *= max[i] - min[i];
        return result;
      }

      // Intersection predicates
      bool intersect(const Point<D>& p) const {
        for(long i = 0; i < D; i++) 
          if(p[i] > max[i] || p[i] < min[i]) return false;
        return true;
      }

        bool intersect(const IsotheticBBox<D>& b) const {
          for(long i = 0; i < D; i++) 
            if(b.max[i] < min[i] || b.min[i] > max[i]) return false;
          return true;
        }

        private:
        Point<D> min,max;
        void clearBox() { 
          for (size_t i=0;i<D;++i) {
            min[i] = REAL_MAX; 
            max[i] = REAL_MIN; 
          }
        }

      };


      /// Determine if two IsotheticBBox elements are coincident in space.
      /**
       * approxEq() returns true if the corner points defining box1 and box2
       * are within tol of each other.  The comparison is made using the
       * approxEq() function for the Point class.
       */
      template <long D> bool
        approxEq(const IsotheticBBox<D>& box1,
            const IsotheticBBox<D>& box2,
            const double& tol) {
          const Point<D>& min1 = box1.getMin();
          const Point<D>& max1 = box1.getMax();
          const Point<D>& min2 = box2.getMin();
          const Point<D>& max2 = box2.getMax();

          return (approxEq(min1,min2,tol) && approxEq(max1,max2,tol));
        }

      template <long D> bool 
        interval(const IsotheticBBox<D>& box, const Point<D>& orig, 
            const Vector<D>& magdir, double& a, double& b) 
        {
          a = 0;
          b = 1;

          for(long i = 0; i < D; i++) {
            if(magdir[i] == 0) {
              if(orig[i] < box.getMin(i) || orig[i] > box.getMax(i)) return false;
            } else {
              register const double dMin = (box.getMin(i)-orig[i])/magdir[i];
              register const double dMax = (box.getMax(i)-orig[i])/magdir[i];
              if(dMax < dMin) {
                if(b < dMax || a > dMin) return false;
                if(a < dMax) a = dMax;
                if(b > dMin) b = dMin;
              } else {
                if(b < dMin || a > dMax) return false;
                if(a < dMin) a = dMin;
                if(b > dMax) b = dMax;
              }
            }
          }
          return true;
        }

      typedef IsotheticBBox<2> IsotheticBBox2;
      typedef IsotheticBBox<3> IsotheticBBox3;


      template<long D> inline std::ostream& 
        operator<<(std::ostream& os, const gk::IsotheticBBox<D>& box) {
          if(box.empty())
            os << "(empty)";
          else
            os << "min=" << box.getMin() << " max=" << box.getMax();
          return os;
        }

  }

#endif
