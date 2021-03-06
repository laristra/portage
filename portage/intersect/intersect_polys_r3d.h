/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef INTERSECT_POLYS_R3D_H
#define INTERSECT_POLYS_R3D_H

#include <array>
#include <stdexcept>
#include <vector>
#include <algorithm>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

extern "C" {
#include "r3d.h"
}

#include "portage/support/portage.h"
#ifdef PORTAGE_HAS_TANGRAM
#include "tangram/support/MatPoly.h"
#endif



namespace Portage {

using Wonton::Point;
using Wonton::Vector;

typedef
struct facetedpoly {
  std::vector<std::vector<int>> facetpoints;  // we can flatten this to a
  //                                          // 1D array if we assume only
  //                                          // triangular facets or we
  //                                          // include the number of points
  //                                          // in each facet
  std::vector<Point<3>> points;
} facetedpoly_t;

#ifdef PORTAGE_HAS_TANGRAM
/*!
 @brief Facetizes a 3D MatPoly and converts it to the facetedpoly_t structure
 @param matpoly  3D material polyhedron object to be converted
 @return  Corresponding facetedpoly_t structure
*/
inline
facetedpoly_t get_faceted_matpoly(const Tangram::MatPoly<3>& matpoly) {
  // facet the matpoly
  Tangram::MatPoly<3> faceted_matpoly;
  matpoly.faceted_matpoly(&faceted_matpoly);
  
  // initialize the result with the face vertices
  facetedpoly_t result{faceted_matpoly.face_vertices()};
  
  // convert tangram points to portage points for facetedpoly_t.points
  for (auto p : faceted_matpoly.points()) result.points.push_back(Point<3>(p));
  return result; 
}
#endif


// Intersect one source polyhedron (possibly non-convex but with
// triangular facets only) with a bunch of tets forming a target
// polyhedron

// NOTE: Given that this routine is R3D specific and is called by an
// R3D specific functor, we should just send in R3D-ized data
// structures for the source polyhedron and the target tet faces. It
// will cut down some calls to R3D initialization routines
// (particularly on the target mesh side)

std::vector<double>
inline
intersect_polys_r3d(const facetedpoly_t &srcpoly,
                    const std::vector<std::array<Point<3>, 4>> &target_tet_coords,
                    NumericTolerances_t num_tols) {

  // Bounding box of the target cell - will be used to compute
  // epsilon for bounding box check. We could use the source cell
  // coordinates to find the bounds as well but its probably an
  // unnecessary computation (EVEN if the target cell is much
  // smaller than the source cell)

  std::vector<double> source_cell_bounds =  {1e99, -1e99, 1e99, -1e99,
                                             1e99, -1e99};

  // Initialize the source polyhedron description in a form R3D wants
  // Simultaneously compute the bounding box
  r3d_poly src_r3dpoly;
  int num_verts = srcpoly.points.size();
  r3d_rvec3 *verts = new r3d_rvec3[num_verts];
  for (int i = 0; i < num_verts; i++) {
    for (int j = 0; j < 3; j++) {
      verts[i].xyz[j] = srcpoly.points[i][j];

      if (source_cell_bounds[2*j] > verts[i].xyz[j])
        source_cell_bounds[2*j] = verts[i].xyz[j];
      if (source_cell_bounds[2*j+1] < verts[i].xyz[j])
        source_cell_bounds[2*j+1] = verts[i].xyz[j];
    }
  }

  // used only for bounding box check not for intersections
  double bbeps = num_tols.min_absolute_distance;

  int num_faces = srcpoly.facetpoints.size();
  r3d_int *face_num_verts = new r3d_int[num_faces];
  for (int i = 0; i < num_faces; i++)
    face_num_verts[i] = srcpoly.facetpoints[i].size();

  r3d_int **face_vert_ids = new r3d_int *[num_faces];
  for (int i = 0; i < num_faces; i++)
    face_vert_ids[i] = new r3d_int[face_num_verts[i]];

  for (int i = 0; i < num_faces; i++)
    std::copy(srcpoly.facetpoints[i].begin(), srcpoly.facetpoints[i].end(),
              face_vert_ids[i]);


  int ok = r3d_init_poly(&src_r3dpoly, verts, num_verts, face_vert_ids,
                         face_num_verts, num_faces);
  if (!ok)
    throw std::runtime_error("intersect_polys_r3d.h: Failed to initialize R3D polyhedron");

  // Finished building source poly; now intersect with tets of target cell

  std::vector<double> moments(4, 0);
  for (auto const & target_cell_tet : target_tet_coords) {
    std::vector<r3d_plane> faces(4);

    std::vector<double> target_tet_bounds = {1e99, -1e99, 1e99,
                                             -1e99, 1e99, -1e99};
    r3d_rvec3 verts2[4];
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 3; j++) {
        verts2[i].xyz[j] = target_cell_tet[i][j];

        if (target_tet_bounds[2*j] > verts2[i].xyz[j])
          target_tet_bounds[2*j] = verts2[i].xyz[j];
        if (target_tet_bounds[2*j+1] < verts2[i].xyz[j])
          target_tet_bounds[2*j+1] = verts2[i].xyz[j];
      }

    // Check if the target and source bounding boxes overlap - bbeps
    // is used to subject touching tets to the full intersection
    // (just in case)

    bool check_bb = true;
    if (check_bb) {
      bool disjoint = false;
      for (int j = 0; j < 3 && !disjoint; ++j)
        disjoint = (target_tet_bounds[2*j] > source_cell_bounds[2*j+1]+bbeps ||
                    target_tet_bounds[2*j+1] < source_cell_bounds[2*j]-bbeps);
      if (disjoint) continue;
    }

    r3d_tet_faces_from_verts(&faces[0], verts2);

    // clip the source poly against the faces of the target tet - but
    // make a copy of src_r3dpoly first because it will get modified
    // in the process of clipping
    r3d_poly src_r3dpoly_copy = src_r3dpoly;
    ok = r3d_clip(&src_r3dpoly_copy, &faces[0], 4);
    if (!ok)
      throw std::runtime_error("intersect_polys_r3d.h: r3d clip failed");

    // find the moments (up to quadratic order) of the clipped poly
    const int POLY_ORDER = 1;
    r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];
    r3d_reduce(&src_r3dpoly_copy, om, POLY_ORDER);

    // Accumulate moments:
    for (int i = 0; i < 4; i++)
      moments[i] += om[i];
  }

  delete [] verts;
  delete [] face_num_verts;
  for (int i = 0; i < num_faces; i++)
    delete [] face_vert_ids[i];
  delete [] face_vert_ids;

  return moments;
}  // intersect_polys_3D

}  // namespace Portage

#endif  // INTERSECT_POLYS_R3D_H
