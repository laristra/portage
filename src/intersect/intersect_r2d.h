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



#ifndef INTERSECT_R2D_H
#define INTERSECT_R2D_H

#include <array>
#include <stdexcept>
#include <vector>

extern "C" {
#include "r2d.h"
}

#include "portage/support/Point.h"


namespace Portage {

/*!
 * \class IntersectR2D <typename C> 2-D intersection algorithm
 */

template <typename SourceMeshType, typename TargetMeshType=SourceMeshType>
class IntersectR2D {

public:

  /// Alias for a collection of Points.
  typedef std::vector<Portage::Point<2>> Poly; 
  /// Alias to provide volume and centroid
  typedef std::pair<double, Portage::Point<2>> Moment;

  /// Constructor taking a source mesh @c s and a target mesh @c t.
  IntersectR2D(const SourceMeshType &s, const TargetMeshType &t)
    : sourceMeshWrapper(s), targetMeshWrapper(t) {}

  /*! \brief Intersect two cells and return the first two moments.
   * \param[in] cellA first cell index to intersect
   * \param[in] cellB second cell index to intersect
   * \return list of moments; ret[0] == 0th moment; ret[1] == first moment
   */
  std::vector<std::vector<double>> operator() (const int cellA,
            const int cellB) const {
    Poly source_poly, target_poly;
    sourceMeshWrapper.cell_get_coordinates(cellA, &source_poly);
    targetMeshWrapper.cell_get_coordinates(cellB, &target_poly);

    std::vector<double> moments(3, 0);

    const int size1 = source_poly.size();
    std::vector<r2d_rvec2> verts1(size1);
    for (int i=0; i<size1; ++i) {
      verts1[i].xy[0] = source_poly[i][0];
      verts1[i].xy[1] = source_poly[i][1];
    }
    // Only do this check in Debug mode. This test is important,
    // otherwise R2D returns invalid results.
#ifdef DEBUG
    double volume1 = 0.;
    for (int i=2; i<size1; ++i)
      volume1 += r2d_orient(verts1[0], verts1[i-1], verts1[i]);
    if (volume1 < 0.)
      throw std::runtime_error("source_poly has negative volume");
#endif

    const int size2 = target_poly.size();
    std::vector<r2d_plane> faces(size2);
    std::vector<r2d_rvec2> verts2(size2);
    for (int i=0; i<size2; ++i) {
      verts2[i].xy[0] = target_poly[i][0];
      verts2[i].xy[1] = target_poly[i][1];
    }

    // check for convexity
    bool convex = true;
    const double eps=1e-14;
    for (int i=0; i<size2; ++i) {
      int iprev = (i == 0 ? size2-1 : i-1);
      int inext = (i+1 == size2 ? 0 : i+1);
      if (r2d_orient(verts2[iprev], verts2[i], verts2[inext]) < -eps) {
        convex = false;
        break;
      }
    }

    // case 1:  target_poly is convex
    // can simply use faces of target_poly as clip planes
    if (convex) {
      r2d_poly_faces_from_verts(&faces[0], &verts2[0], size2);

      // Only do this check in Debug mode. This test is important,
      // otherwise R2D returns invalid results.
#ifdef DEBUG
      double volume2 = 0.;
      for (int i=2; i<size2; ++i)
        volume2 += r2d_orient(verts2[0], verts2[i-1], verts2[i]);
      if (volume2 < 0.)
        throw std::runtime_error("target_poly has negative volume");
#endif

      r2d_poly poly;
      r2d_init_poly(&poly, &verts1[0], size1);
      // Only do this check in Debug mode:
#ifdef DEBUG
      if (r2d_is_good(&poly) == 0)
        throw std::runtime_error("source_poly: invalid poly");
#endif

      // clip the first poly against the faces of the second
      r2d_clip(&poly, &faces[0], size2);

      // find the moments (up to quadratic order) of the clipped poly
      const int POLY_ORDER=1;
      // Only do this check in Debug mode:
#ifdef DEBUG
      if (R2D_NUM_MOMENTS(POLY_ORDER) != 3)
        throw std::runtime_error("Invalid number of moments");
#endif
      r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
      r2d_reduce(&poly, om, POLY_ORDER);

      // Check that the returned volume is positive (if the volume is zero,
      // i.e. abs(om[0]) < eps, then it can sometimes be slightly negative,
      // like om[0] == -1.24811e-16. For this reason we use the condition
      // om[0] < -eps.
      if (om[0] < -eps) throw std::runtime_error("Negative volume");

      // Copy moments:
      for (int j=0; j<3; ++j) {
        moments[j] = om[j];
      }

    }  // if convex

    // case 2:  target_poly is non-convex
    // must divide target_poly into triangles for clipping
    if (!convex) {
      Point<2> centroid;
      targetMeshWrapper.cell_centroid(cellB, &centroid);
      faces.resize(3);
      verts2.resize(3);

      for (int i=0; i<size2; ++i) {
        verts2[0].xy[0] = centroid[0];
        verts2[0].xy[1] = centroid[1];
        verts2[1].xy[0] = target_poly[i][0];
        verts2[1].xy[1] = target_poly[i][1];
        int inext = (i+1 == size2 ? 0 : i+1);
        verts2[2].xy[0] = target_poly[inext][0];
        verts2[2].xy[1] = target_poly[inext][1];
        if (r2d_orient(verts2[0], verts2[1], verts2[2]) < 0.) {
#ifdef DEBUG
					std::cerr << "WARNING: non-convex polygon is incorrectly triangulated (see bug #430)." << std::endl;
#endif
          std::swap(verts2[1].xy[0], verts2[2].xy[0]);
          std::swap(verts2[1].xy[1], verts2[2].xy[1]);
        }

        // Only do this check in Debug mode. This test is important,
        // otherwise R2D returns invalid results.
#ifdef DEBUG
        if (r2d_orient(verts2[0], verts2[1], verts2[2]) < 0.)
          throw std::runtime_error("target_wedge has negative volume");
#endif
        r2d_poly_faces_from_verts(&faces[0], &verts2[0], 3);

        r2d_poly poly;
        r2d_init_poly(&poly, &verts1[0], size1);
        // Only do this check in Debug mode:
#ifdef DEBUG
        if (r2d_is_good(&poly) == 0)
          throw std::runtime_error("source_poly: invalid poly");
#endif

        // clip the first poly against the faces of the second
        r2d_clip(&poly, &faces[0], 3);

        // find the moments (up to quadratic order) of the clipped poly
        const int POLY_ORDER=1;
        // Only do this check in Debug mode:
#ifdef DEBUG
        if (R2D_NUM_MOMENTS(POLY_ORDER) != 3)
          throw std::runtime_error("Invalid number of moments");
#endif
        r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
        r2d_reduce(&poly, om, POLY_ORDER);

        // Check that the returned volume is positive (if the volume is zero,
        // i.e. abs(om[0]) < eps, then it can sometimes be slightly negative,
        // like om[0] == -1.24811e-16. For this reason we use the condition
        // om[0] < -eps.
        if (om[0] < -eps) throw std::runtime_error("Negative volume");

        // Accumulate moments:
        for (int j=0; j<3; ++j) {
          moments[j] += om[j];
        }

      } // for i
    } // if !convex

    std::vector<std::vector<double>> moments_all;
    moments_all.push_back(moments);
    return moments_all;
  }

  IntersectR2D() = delete;

  //! Copy constructor (disabled)
  IntersectR2D(const IntersectR2D &) = delete;

  //! Assignment operator (disabled)
  IntersectR2D & operator = (const IntersectR2D &) = delete;

private:
  const SourceMeshType &sourceMeshWrapper;
  const TargetMeshType &targetMeshWrapper;

}; // class IntersectR2D


} // namespace Portage

#endif // INTERSECT_R2D_H
