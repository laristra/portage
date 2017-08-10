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

#include "portage/support/portage.h"
#include "portage/support/Point.h"


namespace Portage {

/*!
 * \class IntersectR2D <typename C> 2-D intersection algorithm
 */

class IntersectR2D {

 public:

  /*! \brief Intersect two polygons and return the first two moments.
   * \param[in] cellA vector of points of first polygon (or points of triangles of the first polygon)
   * \param[in] cellB vector of points of second polygon (or points of triangles of the second polygon)
   * \return vector of moments; area, 1st moment x component, 1st moment y component
   
   The polygons may be described by a counterclockwise list of points
   (no repeated points) or by a list of points representing a
   triangular decomposition of the polygon. Thus if we send in a list
   of points p1, p2, p3, p4 we can assume it is quad but if we send in
   p1, p2, p3, p1, p3, p4 we assume that it is a triangular
   decomposition of the same quad. The straightforward description is
   useful for convex polygons and the triangular decomposition when one
   or both of the polygons are non-convex. NOTE: Both polygon
   descriptions must use the same convention.
  */
  
  Portage::vector<double>
  operator()(const Portage::vector<Point<2>> & source_poly,
             const Portage::vector<Point<2>> & target_poly) const {

    /* Check if any point is repeated. If it is, then we assume that we
       have been given a triangular decomposition of a polygon. We
       check if any of the first three points are repeated in the rest
       of the list */

    int npnts1 = source_poly.size();
    Point<2> p1, p2;
    bool found = false;
    for (int i = 0; i < 3 && !found; i++)
      for (int j = 4; j < npnts1 && !found; j++)
        if (approxEq(source_poly[i], source_poly[j], 1.0e-20))
          found = true;

    if (found)
      return intersect_polytris(source_poly, target_poly);
    else
      return intersect_polypnts(source_poly, target_poly);
  }



  /*! \brief Intersect two convex polygons described by counterclockwise list of points and return the first two moments.
   * \param[in] cellA vector of points of first polygon (or points of triangles of the first polygon)
   * \param[in] cellB vector of points of second polygon (or points of triangles of the second polygon)
   * \return vector of moments; area, 1st moment x component, 1st moment y component
   */
  Portage::vector<double>
  intersect_polypnts(const Portage::vector<Point<2>> & source_poly,
                     const Portage::vector<Point<2>> & target_poly) const {
    Portage::vector<double> moments(3, 0.0);
   
    const int size1 = source_poly.size();
    Portage::vector<r2d_rvec2> verts1(size1);
    for (int i=0; i<size1; ++i) {
      verts1[i].xy[0] = source_poly[i][0];
      verts1[i].xy[1] = source_poly[i][1];
    }

    // Only do this check in Debug mode. This test is important,
    // otherwise R2D returns invalid results.
   
    // THIS COMMENT IS MISLEADING - IF IT'S AFFECTING THE RESULT, IT'S
    // AN OPERATION NOT A CHECK. AND IF IT'S AN OPERATION THAT CAN GIVE
    // AN INCORRECT RESULT, IF OMITTED, SHOULDN'T WE ALWAYS DO IT?
   
#ifdef DEBUG
    double volume1 = 0.;
    for (int i=2; i<size1; ++i) {
      volume1 += r2d_orient(verts1[0], verts1[i-1], verts1[i]);
    }
    if (volume1 < 0.0)
      throw std::runtime_error("source polygon has negative volume");
#endif
   
    const int size2 = target_poly.size();
    Portage::vector<r2d_plane> faces(size2);
    Portage::vector<r2d_rvec2> verts2(size2);
    for (int i=0; i<size2; ++i) {
      verts2[i].xy[0] = target_poly[i][0];
      verts2[i].xy[1] = target_poly[i][1];
    }

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
   
    for (int j = 0; j < 3; ++j) {
      moments[j] = om[j];
    }
    return moments;
  }


  /*! \brief Intersect two possibly non-convex polygons described by a set of triangles and return the first two moments 
   * \param[in] cellA points of triangles of first polygon (arranged linearly) 
   * \param[in] cellB points of triangles of second polygon (arranged linearly)
   * \return vector of moments; area, 1st moment x component, 1st moment y component
   *
   * In principle, the intersection of a non-convex polygon with
   * another (convex or non-convex) polygon may return multiple
   * disjoint pieces but for now we aggregate their moments and return
   * it; otherwise, we would have to check with pieces are connected
   * and which are not or we would have to return a multitude of pieces
   * corresponding to each triangle-triangle intersection
   */
  Portage::vector<double>
  intersect_polytris(const Portage::vector<Point<2>> & source_poly,
                     const Portage::vector<Point<2>> & target_poly) const {
    Portage::vector<double> moments(3, 0.0);

    const int ntris1 = source_poly.size()/3;

    // Only do this check in Debug mode. This test is important,
    // otherwise R2D returns invalid results.
   
    // THIS COMMENT IS MISLEADING - IF IT'S AFFECTING THE RESULT, IT'S
    // AN OPERATION NOT A CHECK. AND IF IT'S AN OPERATION THAT CAN GIVE
    // AN INCORRECT RESULT, WHEN OMITTED, SHOULDN'T WE ALWAYS DO IT?
   
#ifdef DEBUG
    double trivolume1 = 0.;
    for (int i = 0; i < ntris1; ++i) {
      r2d_rvec2 verts1[3];
      verts1[0].xy[0] = source_poly[3*i][0];
      verts1[0].xy[1] = source_poly[3*i][1];
      verts1[1].xy[0] = source_poly[3*i+1][0];
      verts1[1].xy[1] = source_poly[3*i+1][1];
      verts1[2].xy[0] = source_poly[3*i+2][0];
      verts1[2].xy[1] = source_poly[3*i+2][1];
      trivolume1 = r2d_orient(verts1[0], verts1[1], verts1[2]);
      if (trivolume1 < 0.0)
        throw std::runtime_error("triangle of source polygon has negative volume");
    }
#endif

    const int ntris2 = target_poly.size()/3;
   
    for (int i = 0; i < ntris2; ++i) {
      r2d_rvec2 verts2[3];
      verts2[0].xy[0] = target_poly[3*i][0];
      verts2[0].xy[1] = target_poly[3*i][1];
      verts2[1].xy[0] = target_poly[3*i+1][0];
      verts2[1].xy[1] = target_poly[3*i+1][1];
      verts2[2].xy[0] = target_poly[3*i+2][0];
      verts2[2].xy[1] = target_poly[3*i+2][1];
     
      // Only do this check in Debug mode. This test is important,
      // otherwise R2D returns invalid results. SEE NOTE ABOVE
#ifdef DEBUG
      if (r2d_orient(verts2[0], verts2[1], verts2[2]) < 0.)
        throw std::runtime_error("target_poly triangle has negative volume");
#endif
      r2d_plane faces[3];
      r2d_poly_faces_from_verts(&faces[0], verts2, 3);

      for (int j = 0; j < ntris1; ++j) {
        r2d_rvec2 verts1[3];
        verts1[0].xy[0] = source_poly[3*j][0];
        verts1[0].xy[1] = source_poly[3*j][1];
        verts1[1].xy[0] = source_poly[3*j+1][0];
        verts1[1].xy[1] = source_poly[3*j+1][1];
        verts1[2].xy[0] = source_poly[3*j+2][0];
        verts1[2].xy[1] = source_poly[3*j+2][1];

        r2d_poly poly;
        r2d_init_poly(&poly, verts1, 3);
        // Only do this check in Debug mode:
#ifdef DEBUG
        if (r2d_is_good(&poly) == 0)
          throw std::runtime_error("source_poly: invalid triangle");
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
        for (int j = 0; j < 3; ++j) {
          moments[j] += om[j];
        }
      }  // for j
    }  // for i
   
    return moments;
  }

  //! Use default constructor
  IntersectR2D() = default;

  //! Copy constructor (disabled)
  IntersectR2D(const IntersectR2D &) = delete;
 
  //! Assignment operator (disabled)
  IntersectR2D & operator = (const IntersectR2D &) = delete;

 private:
  const double eps = 1e-14;

}; // class IntersectR2D


} // namespace Portage

#endif // INTERSECT_R2D_H
