/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

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
    // TODO:  check indents throughout!
    const Poly source_poly = sourceMeshWrapper.cellToXY(cellA);
    const Poly target_poly = targetMeshWrapper.cellToXY(cellB);

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
    for (int i=2; i<size1; ++i)
      if (r2d_orient(verts1[0], verts1[i-1], verts1[i]) < 0)
        throw std::runtime_error("source_poly has negative volume");
#endif

    const int size2 = target_poly.size();
    std::vector<r2d_plane> faces(size2);
    std::vector<r2d_rvec2> verts2(size2);
    for (int i=0; i<size2; ++i) {
      verts2[i].xy[0] = target_poly[i][0];
      verts2[i].xy[1] = target_poly[i][1];
    }
    // Only do this check in Debug mode. This test is important,
    // otherwise R2D returns invalid results.
#ifdef DEBUG
    for (int i=2; i<size2; ++i)
      if (r2d_orient(verts2[0], verts2[i-1], verts2[i]) < 0)
        throw std::runtime_error("target_poly has negative volume");
#endif
    r2d_poly_faces_from_verts(&faces[0], &verts2[0], size2);

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
    const double eps=1e-14;
    if (om[0] < -eps) throw std::runtime_error("Negative volume");

    // Copy moments:
    for (int i=0; i<3; ++i) {
      moments[i] = om[i];
    }

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
