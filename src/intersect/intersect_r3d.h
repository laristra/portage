/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef INTERSECT_R3D_H
#define INTERSECT_R3D_H

#include "search.h"
#include "clipper.hpp"
#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <array>

extern "C" {
#include "r3d.h"
}

namespace Portage {

/*!
 * \class IntersectR3D <typename C> 3-D intersection algorithm
 */

template <typename SourceMeshType, typename TargetMeshType=SourceMeshType>
class IntersectR3D {
public:
  typedef std::tuple<double, double, double> Point;
  typedef std::vector<Point> Poly;
  //Provide volume and centroid
  typedef std::pair<double, Point> Moment;

  IntersectR3D(const SourceMeshType &s, const TargetMeshType &t)
    : sourceMeshWrapper(s), targetMeshWrapper(t) {}

  /*! \brief Intersect two cells and return the first two moments.
   * \param[in] cellA first cell index to intersect
   * \param[in] cellB second cell index to intersect
   * \return list of moments; ret[0] == 0th moment; ret[1] == first moment
   */

  std::vector<std::vector<double>> operator() (const int cellA,
            const int cellB) const {
    std::vector<std::array<std::array<double, 3>, 4>> source_coords, target_coords;
    sourceMeshWrapper.wedges_get_coordinates(cellA, &source_coords);
    targetMeshWrapper.wedges_get_coordinates(cellB, &target_coords);
    /*
    std::cout << "Number of source wedges: " << source_coords.size() << std::endl;
    std::cout << "Number of target wedges: " << target_coords.size() << std::endl;
    */
    /*
    for (int i=0; i<source_coords.size(); i++) {
      std::cout << "wedge: " << i << std::endl
        << " " << source_coords[i][0][0] << " " << source_coords[i][0][1] << " "
          << source_coords[i][0][2] << std::endl
        << " " << source_coords[i][1][0] << " " << source_coords[i][1][1] << " "
          << source_coords[i][1][2] << std::endl
        << " " << source_coords[i][2][0] << " " << source_coords[i][2][1] << " "
          << source_coords[i][2][2] << std::endl
        << " " << source_coords[i][3][0] << " " << source_coords[i][3][1] << " "
          << source_coords[i][3][2] << std::endl
        << std::endl;
    }
    */

    std::vector<double> moments(4, 0);

    for (const auto &source_wedge : source_coords) {
      r3d_rvec3 verts1[4];
      for (int i=0; i<4; i++)
        for (int j=0; j<3; j++)
          verts1[i].xyz[j] = source_wedge[i][j];
      // TODO: Only do this check in Debug mode. This test is important,
      // otherwise R3D returns invalid results.
      if (r3d_orient(verts1) < 0)
        throw std::runtime_error("source_wedge has negative volume");

      for (const auto &target_wedge : target_coords) {
        r3d_poly poly;
        r3d_init_tet(&poly, verts1);
        // TODO: Only do this check in Debug mode:
        if (r3d_is_good(&poly) == 0)
          throw std::runtime_error("source_wedge: invalid poly");

        r3d_plane faces[4];
        r3d_rvec3 verts2[4];
        for (int i=0; i<4; i++)
          for (int j=0; j<3; j++)
            verts2[i].xyz[j] = target_wedge[i][j];
        // TODO: Only do this check in Debug mode. This test is important,
        // otherwise R3D returns invalid results.
        if (r3d_orient(verts2) < 0)
          throw std::runtime_error("target_wedge has negative volume");

        r3d_tet_faces_from_verts(faces, verts2);

        // clip the first tet against the faces of the second
        r3d_clip(&poly, faces, 4);

        // find the moments (up to quadratic order) of the clipped poly
        const int POLY_ORDER=1;
        // TODO: Only do this check in Debug mode:
        if (R3D_NUM_MOMENTS(POLY_ORDER) != 4)
          throw std::runtime_error("Invalid number of moments");
        r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];
        r3d_reduce(&poly, om, POLY_ORDER);

        // TODO: R3D sometimes returns a negative volume, even though both
        // r3d_orient() above are positive. We should figure out the ordering
        // so that R3D always returns a positive volume. For now we just take
        // an absolute value:
        om[0] = std::abs(om[0]);

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

  IntersectR3D() {}

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
