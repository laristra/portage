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



#ifndef INTERSECT_R3D_H
#define INTERSECT_R3D_H

#include <array>
#include <stdexcept>
#include <vector>

extern "C" {
#include "r3d.h"
}

#include "portage/support/Point.h"
#include "portage/support/portage.h"


namespace Portage {

/*!
 * \class IntersectR3D <typename C> 3-D intersection algorithm
 */

class IntersectR3D {
 public:
  IntersectR3D() = default;

  /*! \brief Intersect two cells and return the first two moments.
   * \param[in] cellA first cell index to intersect
   * \param[in] cellB second cell index to intersect
   * \return list of moments; ret[0] == 0th moment; ret[1] == first moment
   */

  Portage::vector<double>
  operator() (const Portage::vector<Point<3>> & source_coords,
              const Portage::vector<Point<3>> & target_coords) const {

    // Bounding box of the target cell - will be used to compute
    // epsilon for bounding box check. We could use the source cell
    // coordinates to find the bounds as well but its probably an
    // unnecessary computation (EVEN if the target cell is much
    // smaller than the source cell)

    double target_cell_bounds[6] = {1e99, -1e99, 1e99, -1e99, 1e99, -1e99};

    int nsrctets = source_coords.size()/4;
    int ntrgtets = target_coords.size()/4;

    Portage::vector<double> moments(4, 0);

    Portage::vector<r3d_plane> target_tetfaces_all;
    target_tetfaces_all.reserve(4*ntrgtets);
    Portage::vector<double> target_tetbounds_all;
    target_tetbounds_all.reserve(6*ntrgtets);

    for (int tidx = 0; tidx < ntrgtets; tidx++) {
      r3d_plane faces[4];
      double bounds[6] = {1e99, -1e99, 1e99, -1e99, 1e99, -1e99};
      r3d_rvec3 verts2[4];
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++) {
          verts2[i].xyz[j] = target_coords[4*tidx+i][j];
          if (bounds[2*j] > verts2[i].xyz[j])
            bounds[2*j] = verts2[i].xyz[j];
          if (bounds[2*j+1] < verts2[i].xyz[j])
            bounds[2*j+1] = verts2[i].xyz[j];
        }
      target_tetbounds_all.insert(target_tetbounds_all.end(),
                                  bounds, bounds + 6);

      for (int j = 0; j < 3; j++) {
        if (target_cell_bounds[2*j] > bounds[2*j])
          target_cell_bounds[2*j] = bounds[2*j];
        if (target_cell_bounds[2*j+1] < bounds[2*j+1])
          target_cell_bounds[2*j+1] = bounds[2*j+1];
      }

      // Only do this check in Debug mode. This test is important,
      // otherwise R3D returns invalid results.
#ifdef DEBUG
      if (r3d_orient(verts2) < 0)
        throw std::runtime_error("target_wedge has negative volume");
#endif

      r3d_tet_faces_from_verts(&faces[0], verts2);
      target_tetfaces_all.insert(target_tetfaces_all.end(), faces, faces+4);
    }  // for tidx

    double MAXLEN = -1e99;
    for (int j = 0; j < 3; j++) {
      double len = target_cell_bounds[2*j+1]-target_cell_bounds[2*j];
      if (MAXLEN < len) MAXLEN = len;
    }

    double bbeps = 1.0e-12*MAXLEN;  // used only for bounding box check
    //                            // not for intersections


    // Now intersect each tet from the source cell against each tet of
    // the target cell (IF their bounding boxes overlap)

    for (int sidx = 0; sidx < nsrctets; sidx++) {
      r3d_rvec3 verts1[4];
      double source_tetbounds[6] =  {1e99, -1e99, 1e99, -1e99, 1e99, -1e99};
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)  {
          verts1[i].xyz[j] = source_coords[4*sidx+i][j];
          if (source_tetbounds[2*j] > verts1[i].xyz[j])
            source_tetbounds[2*j] = verts1[i].xyz[j];
          if (source_tetbounds[2*j+1] < verts1[i].xyz[j])
            source_tetbounds[2*j+1] = verts1[i].xyz[j];
        }

      // Only do this check in Debug mode. This test is important,
      // otherwise R3D returns invalid results.
#ifdef DEBUG
      if (r3d_orient(verts1) < 0)
        throw std::runtime_error("source_wedge has negative volume");
#endif

      for (int tidx = 0; tidx < ntrgtets; tidx++) {
        r3d_plane *faces = &(target_tetfaces_all[4*tidx]);

        const double *target_tetbounds = &(target_tetbounds_all[6*tidx]);

        // Check if the target and source bounding boxes overlap - bbeps
        // is used to subject touching tets to the full intersection
        // (just in case)

        bool check_bb = true;
        if (check_bb) {
          bool disjoint = false;
          for (int j = 0; j < 3 && !disjoint; ++j)
            disjoint = (target_tetbounds[2*j] > source_tetbounds[2*j+1]+bbeps ||
                        target_tetbounds[2*j+1] < source_tetbounds[2*j]-bbeps);
          if (disjoint) continue;
        }

        r3d_poly poly;
        r3d_init_tet(&poly, verts1);
        // Only do this check in Debug mode:
#ifdef DEBUG
        if (r3d_is_good(&poly) == 0)
          throw std::runtime_error("source_wedge: invalid poly");
#endif

        // clip the first tet against the faces of the second
        r3d_clip(&poly, faces, 4);

        // find the moments (up to quadratic order) of the clipped poly
        const int POLY_ORDER = 1;
        // Only do this check in Debug mode:
#ifdef DEBUG
        if (R3D_NUM_MOMENTS(POLY_ORDER) != 4)
          throw std::runtime_error("Invalid number of moments");
#endif
        r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];
        r3d_reduce(&poly, om, POLY_ORDER);

        // Check that the returned volume is positive (if the volume is zero,
        // i.e. abs(om[0]) < eps, then it can sometimes be slightly negative,
        // like om[0] == -1.24811e-16. For this reason we use the condition
        // om[0] < -eps.
        const double eps = 1e-14;  // @todo Should multiply by domain or element size
        if (om[0] < -eps) throw std::runtime_error("Negative volume");

        // Accumulate moments:
        for (int i = 0; i < 4; i++)
          moments[i] += om[i];
      }  // for tidx
    }  // for sidx

    return moments;
  }

  //! Copy constructor (disabled)
  IntersectR3D(const IntersectR3D &) = delete;

  //! Assignment operator (disabled)
  IntersectR3D & operator = (const IntersectR3D &) = delete;
}; // class IntersectR3D


} // namespace Portage

#endif // INTERSECT_R3D_H
