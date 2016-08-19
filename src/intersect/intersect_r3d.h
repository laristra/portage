/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef INTERSECT_R3D_H
#define INTERSECT_R3D_H

#include <array>
#include <stdexcept>
#include <vector>

extern "C" {
#include "r3d.h"
}

#include "portage/support/Point.h"


namespace Portage {

/*!
 * \class IntersectR3D <typename C> 3-D intersection algorithm
 */

template <typename SourceMeshType, typename TargetMeshType=SourceMeshType>
class IntersectR3D {
public:
  IntersectR3D(const SourceMeshType &s, const TargetMeshType &t)
    : sourceMeshWrapper(s), targetMeshWrapper(t) {}

  /*! \brief Intersect two cells and return the first two moments.
   * \param[in] cellA first cell index to intersect
   * \param[in] cellB second cell index to intersect
   * \return list of moments; ret[0] == 0th moment; ret[1] == first moment
   */

  std::vector<std::vector<double>> operator() (const int cellA,
            const int cellB) const {
    std::vector<std::array<Portage::Point<3>, 4>> source_coords, target_coords;
    sourceMeshWrapper.decompose_cell_into_tets(cellA, &source_coords);
    targetMeshWrapper.decompose_cell_into_tets(cellB, &target_coords);

    // Bounding box of the target cell - will be used to compute
    // epsilon for bounding box check. We could use the source cell
    // coordinates to find the bounds as well but its probably an
    // unnecessary computation (EVEN if the target cell is much
    // smaller than the source cell)

    std::array<double, 6> target_cell_bounds =  {1e99, -1e99, 1e99, -1e99,
                                                 1e99, -1e99};

    std::vector<double> moments(4, 0);

    std::vector<std::array<r3d_plane, 4>> target_tetfaces_all;
    target_tetfaces_all.reserve(target_coords.size());
    std::vector<std::array<double, 6>> target_tetbounds_all;
    target_tetbounds_all.reserve(target_coords.size());

    for (const auto &target_cell_tet : target_coords) {
      std::array<r3d_plane, 4> faces;
      std::array<double, 6> bounds = {1e99, -1e99, 1e99, -1e99, 1e99, -1e99};
      r3d_rvec3 verts2[4];      
      for (int i=0; i<4; i++)
        for (int j=0; j<3; j++) {
          verts2[i].xyz[j] = target_cell_tet[i][j];
          if (bounds[2*j] > verts2[i].xyz[j])
            bounds[2*j] = verts2[i].xyz[j];
          if (bounds[2*j+1] < verts2[i].xyz[j])
            bounds[2*j+1] = verts2[i].xyz[j];
        }
      target_tetbounds_all.emplace_back(bounds);

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
      target_tetfaces_all.emplace_back(faces);      
    }

    double MAXLEN = -1e;
    for (int j = 0; j < 3; j++) {
      double len = target_cell_bounds[2*j+1]-target_cell_bounds[2*j];
      if (MAXLEN < len) MAXLEN = len;
    }

    double bbeps = 1.0e-12*MAXLEN;  // used only for bounding box check
    //                            // not for intersections


    // Now intersect each tet from the source cell against each tet of
    // the target cell (IF their bounding boxes overlap)

    for (const auto &source_cell_tet : source_coords) {
      r3d_rvec3 verts1[4];
      std::array<double, 6> source_tetbounds =  {1e99, -1e99, 1e99, -1e99, 1e99, -1e99};
      for (int i=0; i<4; i++)
        for (int j=0; j<3; j++)  {
          verts1[i].xyz[j] = source_cell_tet[i][j];
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

      int t = 0;
      for (auto &faces : target_tetfaces_all) {

        const std::array<double, 6>& target_tetbounds = target_tetbounds_all[t++];
        
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
        r3d_clip(&poly, &faces[0], 4);

        // find the moments (up to quadratic order) of the clipped poly
        const int POLY_ORDER=1;
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
        const double eps=1e-14;  // @todo Should multiply by domain or element size
        if (om[0] < -eps) throw std::runtime_error("Negative volume");

        // Accumulate moments:
        for (int i=0; i<4; i++) {
          moments[i] += om[i];
        }
      }
    }

    std::vector<std::vector<double>> moments_all;
    moments_all.push_back(moments);
    return moments_all;
  }

  IntersectR3D() = delete;

  //! Copy constructor (disabled)
  IntersectR3D(const IntersectR3D &) = delete;

  //! Assignment operator (disabled)
  IntersectR3D & operator = (const IntersectR3D &) = delete;
private:
  const SourceMeshType &sourceMeshWrapper;
  const TargetMeshType &targetMeshWrapper;
}; // class IntersectR3D


} // namespace Portage

#endif // INTERSECT_R3D_H
