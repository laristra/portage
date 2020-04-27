/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef INTERSECT_POLYS_R2D_H
#define INTERSECT_POLYS_R2D_H

#include <array>
#include <stdexcept>
#include <vector>
#include <algorithm>

extern "C" {
#include "wonton/intersect/r3d/r2d.h"
}

#include "portage/support/portage.h"
#include "wonton/support/Point.h"

namespace Portage {

// intersect one source polygon (possibly non-convex) with a
// triangular decomposition of a target polygon

inline
std::vector<double>
intersect_polys_r2d(std::vector<Wonton::Point<2>> const & source_poly,
                    std::vector<Wonton::Point<2>> const & target_poly,
                    NumericTolerances_t num_tols) {

  std::vector<double> moments(3, 0);
  bool src_convex = true;
  bool trg_convex = true;

  const int POLY_ORDER = 1;  // max degree of moments to calculate

  // Initialize source polygon

  const int size1 = source_poly.size();
  const int size2 = target_poly.size();
  if (!size1 || !size2)
    return moments;  // could allow top level code to avoid an 'if' statement

  std::vector<r2d_rvec2> verts1(size1);
  for (int i = 0; i < size1; ++i) {
    verts1[i].xy[0] = source_poly[i][0];
    verts1[i].xy[1] = source_poly[i][1];
  }
#ifdef DEBUG
  double volume1 = 0.;
  for (int i = 1; i < size1-1; ++i)
    volume1 += r2d_orient(verts1[0], verts1[i], verts1[i+1]);
  if (volume1 < 0.)
    throw std::runtime_error("target_poly has negative volume");
#endif

  r2d_poly srcpoly_r2d;
  r2d_init_poly(&srcpoly_r2d, &verts1[0], size1);

#ifdef DEBUG
  if (r2d_is_good(&srcpoly_r2d) == 0)
    throw std::runtime_error("source_poly: invalid poly");
#endif


  // Initialize target polygon and check for convexity

  std::vector<r2d_plane> faces(size2);
  std::vector<r2d_rvec2> verts2(size2);
  for (int i = 0; i < size2; ++i) {
    verts2[i].xy[0] = target_poly[i][0];
    verts2[i].xy[1] = target_poly[i][1];
  }

  // check for convexity of target polygon
  for (int i = 0; i < size2; ++i) {
    //Compute distance from target_poly[(i+1)%size2] 
    //to segment (target_poly[i], target_poly[(i+2)%size2])
    int ifv = i, imv = (i+1)%size2, isv = (i+2)%size2;
    Wonton::Vector<2> normal( target_poly[isv][1] - target_poly[ifv][1], 
                     -target_poly[isv][0] + target_poly[ifv][0]);
    normal.normalize();
    Wonton::Vector<2> fv2mv = target_poly[imv] - target_poly[ifv];
    double dst = Wonton::dot(fv2mv, normal);

    if (dst <= -num_tols.min_absolute_distance) {
      trg_convex = false;
      break;
    }
  }
  // case 1:  target_poly is convex
  // can simply use faces of target_poly as clip planes
  if (trg_convex) {
    r2d_poly_faces_from_verts(&faces[0], &verts2[0], size2);

    // clip the first poly against the faces of the second
    r2d_clip(&srcpoly_r2d, &faces[0], size2);

    // find the moments (up to quadratic order) of the clipped poly
    r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
    r2d_reduce(&srcpoly_r2d, om, POLY_ORDER);

    // Check that the returned volume is positive (if the volume is zero,
    // i.e. abs(om[0]) < eps, then it can sometimes be slightly negative,
    // like om[0] == -1.24811e-16.
    if (om[0] < num_tols.minimal_intersection_volume)
       throw std::runtime_error("Negative volume");

    // Copy moments:
    for (int j = 0; j < 3; ++j)
      moments[j] = om[j];

  } else {  // case 2:  target_poly is non-convex

    // check for convexity of source polygon
    for (int i = 0; i < size1; ++i) {
      //Compute distance from source_poly[(i+1)%size1] 
      //to segment (source_poly[i], source_poly[(i+2)%size1])
      int ifv = i, imv = (i+1)%size1, isv = (i+2)%size1;
      Wonton::Vector<2> normal( source_poly[isv][1] - source_poly[ifv][1], 
                       -source_poly[isv][0] + source_poly[ifv][0]);
      normal.normalize();
      Wonton::Vector<2> fv2mv = source_poly[imv] - source_poly[ifv];
      double dst = Wonton::dot(fv2mv, normal);

      if (dst <= -num_tols.min_absolute_distance) {
        src_convex = false;
        break;
      }
    }
    // if source polygon is convex while target polygon is non-convex,
    // call the routine with the polygons reversed

    if (src_convex)
      return intersect_polys_r2d(target_poly, source_poly, num_tols);
    else {

      // Must divide target_poly into triangles for clipping.  Choice
      // of the central point is crucial. Try the centroid first -
      // computed by the area weighted sum of centroids of any
      // triangulation of the polygon
      bool center_point_ok = true;
      Wonton::Point<2> cen(0.0, 0.0);
      r2d_rvec2 cenr2d;
      cenr2d.xy[0] = 0.0; cenr2d.xy[1] = 0.0;
      double area_sum = 0.0;
      for (int i = 1; i < size2; ++i) {
        double area = r2d_orient(verts2[0], verts2[i], verts2[(i+1)%size2]);
        area_sum += area;
        Wonton::Point<2> tricen =
            (target_poly[0] + target_poly[i] + target_poly[(i+1)%size2])/3.0;
        cen += area*tricen;
      }
      cen /= area_sum;
      cenr2d.xy[0] = cen[0]; cenr2d.xy[1] = cen[1];

      for (int i = 0; i < size2; ++i)
        if (r2d_orient(cenr2d, verts2[i], verts2[(i+1)%size2]) < 0.)
          center_point_ok = false;

      if (!center_point_ok) {
        // If the centroid is not ok, we have to find the center of
        // the feasible set of the polygon. This means clipping the
        // target_poly with its own face planes/lines. For a
        // non-convex polygon, this will give a new polygon whose
        // interior (See Garimella/Shashkov/Pavel paper on untangling)

        r2d_poly fspoly;
        r2d_init_poly(&fspoly, &verts2[0], size2);

        r2d_poly_faces_from_verts(&faces[0], &verts2[0], size2);

        r2d_clip(&fspoly, &faces[0], size2);

        // If the resulting polygon is empty, we are out of luck
        if (fspoly.nverts == 0)
          std::runtime_error("Could not find a valid center point to triangulate non-convex polygon");

        // Have R2D compute first and second moments of polygon and
        // get its centroid from that

        r2d_real fspoly_moments[R2D_NUM_MOMENTS(POLY_ORDER)];
        r2d_reduce(&fspoly, fspoly_moments, POLY_ORDER);

        cen[0] = cenr2d.xy[0] = fspoly_moments[1]/fspoly_moments[0];
        cen[1] = cenr2d.xy[1] = fspoly_moments[2]/fspoly_moments[0];

        // Even if the resulting feasible set polygon has vertices,
        // maybe its degenerate. So we have to verify that its centroid
        // indeed is a point that will give valid triangles when
        // paired with the edges of the target polygon.

        for (int i = 0; i < size2; ++i)
          if (r2d_orient(cenr2d, verts2[i], verts2[(i+1)%size2]) < 0.)
            center_point_ok = false;
      }

      // If we still don't have a good center point, we are out of luck
      if (!center_point_ok)
        std::runtime_error("Could not find a valid center point to triangulate non-convex polygon");


      for (int i = 0; i < size2; ++i) {
        verts2[0].xy[0] = cen[0];
        verts2[0].xy[1] = cen[1];
        verts2[1].xy[0] = target_poly[i][0];
        verts2[1].xy[1] = target_poly[i][1];
        verts2[2].xy[0] = target_poly[(i+1)%size2][0];
        verts2[2].xy[1] = target_poly[(i+1)%size2][1];

        r2d_poly_faces_from_verts(&faces[0], &verts2[0], 3);

        r2d_poly srcpoly_r2d_copy = srcpoly_r2d;

        // clip the first poly against the faces of the second
        r2d_clip(&srcpoly_r2d_copy, &faces[0], 3);

        // find the moments (up to quadratic order) of the clipped poly
        // Only do this check in Debug mode:
#ifdef DEBUG
        if (R2D_NUM_MOMENTS(POLY_ORDER) != 3)
          throw std::runtime_error("Invalid number of moments");
#endif
        r2d_real om[R2D_NUM_MOMENTS(POLY_ORDER)];
        r2d_reduce(&srcpoly_r2d_copy, om, POLY_ORDER);

        // Check that the returned volume is positive (if the volume is zero,
        // i.e. abs(om[0]) < eps, then it can sometimes be slightly negative,
        // like om[0] == -1.24811e-16.
        if (om[0] < num_tols.minimal_intersection_volume)
          throw std::runtime_error("Negative volume for triangle of polygon");

        // Accumulate moments:
        for (int j = 0; j < 3; ++j)
          moments[j] += om[j];
      }  // for i
    }  // if (src_convex) ... else ...
  }  // if convex {} else {}

  return moments;
}

} // namespace Portage

#endif // INTERSECT_POLYS_R2D_H
