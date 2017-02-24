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
#include "portage/support/Point.h"
#include "portage/support/Vector.h"

namespace gk {

    using Portage::Point;
    using Portage::Vector;

    const double REAL_MAX = std::numeric_limits<double>::max();
    const double REAL_MIN = std::numeric_limits<double>::min();

    /*!
      @class IsotheticBBox "BoundBox.h"
      @brief An isothetic (axis-aligned) N-dimensional bounding box.
      @tparam D Dimension in which the bounding box lives.  Usually one of [1,2,3].
      */
    template<long D> class IsotheticBBox 
    {
        public:
            IsotheticBBox() { clearBox(); }

            /// Re-initialize box limits to default values
            void clear() { clearBox(); }

            /// Tests if the bounding box contains space.
            bool empty() const { return min[0] == REAL_MAX; }

            // Functions to add new elements to the box
            /// Update the bounding box by adding an additional Point.
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

            /// Update the bounding box by adding the extents of another IsotheticBBox
            void add(const IsotheticBBox<D>& box) {
                add(box.getMin());
                add(box.getMax());
            }

            /// Expand the box by an additive constant
            void bulge(double epsilon) {
                for(long i = 0; i < D; i++) {
                    max[i] += epsilon;
                    min[i] -= epsilon;
                }
            }

            // Calculate box feature points
            /// Calculate the Point of the center of the bounding box.
            Point<D> center() const { return Point<D>((min.asV() + max.asV())/2); }

            /// Calculate the center along the @c axis axis.
            double center(long axis) const { return (min[axis] + max[axis])/2; }

            /// Calculate distance from the origin to the center of the bounding box.
            double radius(bool doSqrt = true) const { 
                const Vector<D> r((max-min)/2);
                return r.norm(doSqrt); 
            }

            /// Return the minimum of the bounding box.
            const Point<D>& getMin() const { return min; }

            /// Return the minimum of the bonding box along axis @c i.
            double getMin(long i) const { assert(i < D); return min[i]; }

            /// Return the maximum of the bounding box.
            const Point<D>& getMax() const { return max; }

            /// Return the maximum of the bounding box along axis @c i.
            double getMax(long i) const { assert(i < D); return max[i]; }

            /// Find which axis is the longest of the bounding box.
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

            /// Calculate the volume of the bounding box.
            double volume() const { 
                double result = 1.0;
                for(long i = 0; i < D; i++)  result *= max[i] - min[i];
                return result;
            }

            // Intersection predicates
            /// Determine if the Point @c p is within the bounding box.
            bool intersect(const Point<D>& p) const {
                for(long i = 0; i < D; i++) 
                    if(p[i] > max[i] || p[i] < min[i]) return false;
                return true;
            }

            /// Determine if the IsotheticBBox @c b is within the bounding box.
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
    template <long D> bool approxEq(const IsotheticBBox<D>& box1,
                                    const IsotheticBBox<D>& box2,
                                    const double& tol) 
    {
        const Point<D>& min1 = box1.getMin();
        const Point<D>& max1 = box1.getMax();
        const Point<D>& min2 = box2.getMin();
        const Point<D>& max2 = box2.getMax();

        return (approxEq(min1,min2,tol) && approxEq(max1,max2,tol));
    }

    template <long D> bool interval(const IsotheticBBox<D>& box, 
                                    const Point<D>& orig, 
                                    const Vector<D>& magdir, 
                                    double& a, double& b) 
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


    template<long D> inline std::ostream& operator<<(std::ostream& os, 
                                                     const gk::IsotheticBBox<D>& box) 
    {
        if(box.empty())
            os << "(empty)";
        else
            os << "min=" << box.getMin() << " max=" << box.getMax();
        return os;
    }

}

#endif
