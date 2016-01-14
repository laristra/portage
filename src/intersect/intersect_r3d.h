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
    Poly elA = sourceMeshWrapper.cellToXYZ(cellA);
    Poly elB = targetMeshWrapper.cellToXYZ(cellB);

    // TODO: Convert elA and elB to tets iterate and fill them into verts1 and
    // verts2 below:

    // variables: the polyhedra and their moments
#define POLY_ORDER 1
    r3d_poly poly;
    r3d_plane faces[4];
    r3d_real om[R3D_NUM_MOMENTS(POLY_ORDER)];
    //Intersect unit tet with itself (should be vol 1/6, centroid .25, .25,.25)
    r3d_rvec3 verts1[4] = // {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}
      {{2,2,2},{2,-2,-2},{-2,2,-2},{-2,-2,2}};
    r3d_rvec3 verts2[4] = // {{0,0,0}, {1,0,0}, {0,1,0}, {0,0,1}}
      {{-1,-1,-1},{-1,1,1},{1,-1,1},{1,1,-1}};
    for(int count=0; count< 10000; count++){
      r3d_init_tet(&poly, verts1);
      r3d_tet_faces_from_verts(faces, verts2);
      // clip the first tet against the faces of the second
      r3d_clip(&poly, faces, 4);
      // find the moments (up to quadratic order) of the clipped poly
      r3d_reduce(&poly, om, POLY_ORDER);
    }

    std::cout << "volume is " << om[0] << std::endl;
    for(int i=1;i<sizeof(om)/sizeof(om[0]);i++){
      std::cout << "centroid [i] is " << om[i]/om[0] << std::endl;
    }

    std::vector<std::vector<double>> moments;
    for(int i=0; i<R3D_NUM_MOMENTS(POLY_ORDER); i++) {
      om[i] = std::abs(om[i]);
    }
    moments = {{om[0], om[1]/om[0], om[2]/om[0], om[3]/om[0]}};
    return moments;
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
