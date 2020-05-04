/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
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
//                                                                           //
// Modified for Portage by: Rao Garimella                                    //
//                          rao@lanl.gov                                     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PORTAGE_SEARCH_BOUNDBOX_H_
#define PORTAGE_SEARCH_BOUNDBOX_H_

#include <iostream>
#include <cassert>
#include <limits>

#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"
#include "portage/support/portage.h"

namespace Portage {
    using Wonton::Point; 
    using Wonton::Vector;

    const double real_max = std::numeric_limits<double>::max();
    const double real_min = std::numeric_limits<double>::min();

    /*!
      @class IsotheticBBox "BoundBox.h"
      @brief An isothetic (axis-aligned) N-dimensional bounding box.
      @tparam D Dimension in which the bounding box lives.  Usually one of [1,2,3].
      */
    template<int D> class IsotheticBBox
    {
        public:
            IsotheticBBox() { clearBox(); }

            /// Re-initialize box limits to default values
            void clear() { clearBox(); }

            /// Tests if the bounding box contains space.
            bool empty() const { return min[0] == real_max; }

            // Functions to add new elements to the box
            /// Update the bounding box by adding an additional Point.
            void add(const Point<D>& p) {
                if(empty()) {
                    min = p; max = p;
                }
                else {
                    for(int i = 0; i < D; i++) {
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
                for(int i = 0; i < D; i++) {
                    max[i] += epsilon;
                    min[i] -= epsilon;
                }
            }

            // Calculate box feature points
            /// Calculate the Point of the center of the bounding box.
            Point<D> center() const { return Point<D>((min.asV() + max.asV())/2); }

            /// Calculate the center aint the @c axis axis.
            double center(int axis) const { return (min[axis] + max[axis])/2; }

            /// Calculate distance from the origin to the center of the bounding box.
            double radius(bool doSqrt = true) const {
                const Vector<D> r((max-min)/2);
                return r.norm(doSqrt);
            }

            /// Return the minimum of the bounding box.
            const Point<D>& getMin() const { return min; }

            /// Return the minimum of the bonding box along axis @c i.
            double getMin(int i) const { assert(i < D); return min[i]; }

            /// Return the maximum of the bounding box.
            const Point<D>& getMax() const { return max; }

            /// Return the maximum of the bounding box along axis @c i.
            double getMax(int i) const { assert(i < D); return max[i]; }

            /// Find which axis is the longest of the bounding box.
            int longAxis() const {
                // Find a longest axis
                const Vector<D> diff = max-min;
                double maxVal = diff[0];
                int maxIdx = 0;
                for(int i = 1; i < D; i++)
                    if(diff[i] > maxVal) {
                        maxVal = diff[i];
                        maxIdx = i;
                    }
                return maxIdx;
            }

            /// Calculate the volume of the bounding box.
            double volume() const {
                double result = 1.0;
                for(int i = 0; i < D; i++)  result *= max[i] - min[i];
                return result;
            }

            // Intersection predicates
            /// Determine if the Point @c p is within the bounding box.
            bool intersect(const Point<D>& p) const {
                for(int i = 0; i < D; i++)
                    if(p[i] > max[i] || p[i] < min[i]) return false;
                return true;
            }

            /// Determine if the IsotheticBBox @c b is within the bounding box.
            bool intersect(const IsotheticBBox<D>& b) const {
                for(int i = 0; i < D; i++)
                    if(b.max[i] < min[i] || b.min[i] > max[i]) return false;
                return true;
            }

        private:
            Point<D> min,max;
            void clearBox() {
                for (size_t i=0;i<D;++i) {
                    min[i] = real_max;
                    max[i] = real_min;
                }
            }

    };


    /// Determine if two IsotheticBBox elements are coincident in space.
    /**
     * approxEq() returns true if the corner points defining box1 and box2
     * are within tol of each other.  The comparison is made using the
     * approxEq() function for the Point class.
     */
    template <int D> bool approxEq(const IsotheticBBox<D>& box1,
                                    const IsotheticBBox<D>& box2,
                                    const double& tol)
    {
        const Point<D>& min1 = box1.getMin();
        const Point<D>& max1 = box1.getMax();
        const Point<D>& min2 = box2.getMin();
        const Point<D>& max2 = box2.getMax();

        return (approxEq(min1,min2,tol) && approxEq(max1,max2,tol));
    }

    template <int D> bool interval(const IsotheticBBox<D>& box,
                                    const Point<D>& orig,
                                    const Vector<D>& magdir,
                                    double& a, double& b)
    {
        a = 0;
        b = 1;

        for(int i = 0; i < D; i++) {
            if(magdir[i] == 0) {
                if(orig[i] < box.getMin(i) || orig[i] > box.getMax(i)) return false;
            } else {
                const double dMin = (box.getMin(i)-orig[i])/magdir[i];
                const double dMax = (box.getMax(i)-orig[i])/magdir[i];
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


    template<int D> inline std::ostream& operator<<(std::ostream& os,
                                                    const Portage::IsotheticBBox<D>& box)
    {
        if(box.empty())
            os << "(empty)";
        else
            os << "min=" << box.getMin() << " max=" << box.getMax();
        return os;
    }

}  // namaespace Portage

#endif  // PORTAGE_SEARCH_BOUNDBOX_H_
