/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERPOLATE_COORDINATE_SYSTEM_H_
#define PORTAGE_INTERPOLATE_COORDINATE_SYSTEM_H_

#include "wonton/support/CoordinateSystem.h"

template<int D, class Mesh>
Wonton::Point<D> cell_centroid_rz(const Mesh& mesh, int cell)
{
  assert(D == 2);

  double r(0.0), beta, a, b, c, nx, ny;
  Wonton::Point<D> centroid, xa, xb, xc;
  std::vector<int> faces, nodes, dirs;

  mesh.cell_get_faces_and_dirs(cell, &faces, &dirs);
  mesh.cell_centroid(cell, &xc);
  for (auto&& f : faces) {
    mesh.face_get_nodes(f, &nodes);
    mesh.node_get_coordinates(nodes[0], &xa);
    mesh.node_get_coordinates(nodes[1], &xb);

    a = xa[0];
    b = xb[0];
    c = (a + b) / 2;
    nx = xb[1] - xa[1];
    ny = xa[0] - xb[0];

    beta = a * nx + xa[1] * ny;
    if ((a - xc[0]) * nx + (xa[1] - xc[1]) * ny < 0) beta = -beta;
    r += beta * (a * a + b * b + 4 * c * c) / 24;
  }

  mesh.cell_centroid(cell, &centroid);
  double area = centroid[0] * mesh.cell_volume(cell);
  centroid[0] = r / area;
  return centroid;
}

#endif

