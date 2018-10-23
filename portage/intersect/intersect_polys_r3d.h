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

extern "C" {
#include "r3d.h"
}

#ifdef HAVE_TANGRAM
#include "tangram/support/MatPoly.h"
#endif


#include "portage/support/portage.h"


namespace Portage {

typedef
struct facetedpoly {
  std::vector<std::vector<int>> facetpoints;  // we can flatten this to a
  //                                          // 1D array if we assume only
  //                                          // triangular facets or we
  //                                          // include the number of points
  //                                          // in each facet
  std::vector<Point<3>> points;
} facetedpoly_t;

#ifdef HAVE_TANGRAM
/*!
 @brief Facetizes a 3D MatPoly and converts it to the facetedpoly_t structure
 @param matpoly  3D material polyhedron object to be converted
 @return  Corresponding facetedpoly_t structure
*/
facetedpoly_t get_faceted_matpoly(const Tangram::MatPoly<3>& matpoly) {
  // facet the matpoly
  Tangram::MatPoly<3> faceted_matpoly;
  matpoly.faceted_matpoly(&faceted_matpoly);
  
  // initialize the result with the face vertices
  facetedpoly_t result{faceted_matpoly.face_vertices()};
  
  // convert tangram points to portage points for facetedpoly_t.points
  for (auto p : faceted_matpoly.points()) result.points.push_back(Portage::Point<3>(p));
  return result; 
}
#endif


// Intersect one source polyhedron (possibly non-convex but with
// triangular facets only) with a bunch of tets forming a target
// polyhedron

std::vector<double>
intersect_polys_r3d(const facetedpoly_t &srcpoly,
                    const std::vector<std::array<Point<3>, 4>> &target_tet_coords) {

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

  double MAXLEN = -1e99;
  for (int j = 0; j < 3; j++) {
    double len = source_cell_bounds[2*j+1]-source_cell_bounds[2*j];
    if (MAXLEN < len) MAXLEN = len;
  }
  double bbeps = 1.0e-12*MAXLEN;  // used only for bounding box check
  //                            // not for intersections


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

#ifdef DEBUG
  // Lets check the volume of the source polygon - If its convex or
  // mildly non-convex, we should get a positive volume
  
  // First calculate a center point
  Point<3> cen(0.0, 0.0, 0.0);
  for (int i = 0; i < num_verts; i++)
    cen += srcpoly.points[i];
  cen /= num_verts;

  // Assume triangular facets only
  double polyvol = 0.0;
  for (int i = 0; i < num_faces; i++) {
    // p0, p1, p2 traversed in order form a triangle whose normal
    // points out of the source polyhedron
    const Point<3> &p0 = srcpoly.points[srcpoly.facetpoints[i][0]];
    const Point<3> &p1 = srcpoly.points[srcpoly.facetpoints[i][1]];
    const Point<3> &p2 = srcpoly.points[srcpoly.facetpoints[i][2]];

    Vector<3> v0 = p1-p0;
    Vector<3> v1 = p2-p0;
    Vector<3> outnormal = cross(v0, v1);
    Vector<3> v2 = cen-p0;
    double tetvol = -dot(v2, outnormal)/6.0;
    if (tetvol < 0.0)
      std::cerr << "Wrong orientation of facets or non-convex polyhedron" <<
          std::endl;
    polyvol += tetvol;
  }
  if (polyvol < 0.0)
    throw std::runtime_error("Source polyhedron has negative volume");
#endif

  r3d_init_poly(&src_r3dpoly, verts, num_verts, face_vert_ids, face_num_verts,
                num_faces);

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

#ifdef DEBUG
    if (r3d_orient(verts2) < 0)
      throw std::runtime_error("target_wedge has negative volume");
#endif

    r3d_tet_faces_from_verts(&faces[0], verts2);

    // clip the source poly against the faces of the target tet - but
    // make a copy of src_r3dpoly first because it will get modified
    // in the process of clipping
    r3d_poly src_r3dpoly_copy = src_r3dpoly;
    r3d_clip(&src_r3dpoly_copy, &faces[0], 4);

    // find the moments (up to quadratic order) of the clipped poly
    const int POLY_ORDER = 1;
    r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];
    r3d_reduce(&src_r3dpoly_copy, om, POLY_ORDER);

    // Check that the returned volume is positive (if the volume is
    // zero, i.e. abs(om[0]) < eps, then it can sometimes be
    // slightly negative, like om[0] == -1.24811e-16. For this
    // reason we use the condition om[0] < -eps.

    const double eps = 1e-14;  // @todo multiply by domain or element size
    if (om[0] < -eps) throw std::runtime_error("Negative volume");

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
