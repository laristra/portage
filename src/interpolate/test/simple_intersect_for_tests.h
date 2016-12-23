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



#include <vector>
#include "Point.hh"


namespace BOX_INTERSECT {

void bounding_box(std::vector<JaliGeometry::Point> coords,
                  JaliGeometry::Point *pmin, JaliGeometry::Point *pmax) {
  int dim = coords[0].dim();
  *pmin = coords[0];
  *pmax = coords[0];
  int np = coords.size();
  for (auto pcoord : coords) {
    for (int d = 0; d < dim; d++) {
      if (pcoord[d] < (*pmin)[d])
        (*pmin)[d] = pcoord[d];
      if (pcoord[d] > (*pmax)[d])
        (*pmax)[d] = pcoord[d];
    }
  }
}

bool intersect_boxes(JaliGeometry::Point min1, JaliGeometry::Point max1,
                     JaliGeometry::Point min2, JaliGeometry::Point max2,
                     std::vector<double> *xsect_moments) {

  assert(min1.dim() == max1.dim());
  assert(min1.dim() == min2.dim());
  assert(min1.dim() == max2.dim());

  int dim = min1.dim();
  JaliGeometry::Point intmin(dim), intmax(dim);

  for (int d = 0; d < dim; ++d) {
    // check for non-intersection in this dimension

    if (min1[d] > max2[d]) return false;
    if (min2[d] > max1[d]) return false;

    // sort the min max vals in this dimension
    double val[4];
    val[0] = min1[d]; val[1] = max1[d]; val[2] = min2[d]; val[3] = max2[d];
    for (int i = 0; i < 3; i++)
      for (int j = i+1; j < 4; j++)
        if (val[i] > val[j]) {
          double tmp = val[i];
          val[i] = val[j];
          val[j] = tmp;
        }

    // pick the middle two as the min max coordinates of intersection
    // box in this dimension

    intmin[d] = val[1]; intmax[d] = val[2];
  }

  // Calculate the volume

  double vol = 1.0;
  for (int d = 0; d < dim; d++)
    vol *= intmax[d]-intmin[d];

  // Sanity check

  assert(vol >= 0.0);

  if (vol == 0.0) return false;

  // Calculate the centroid

  JaliGeometry::Point centroid = (intmin+intmax)/2.0;

  // moments

  xsect_moments->clear();
  xsect_moments->push_back(vol);
  for (int d = 0; d < dim; d++)
    xsect_moments->push_back(centroid[d]*vol);

  return true;
}


void intersection_moments(std::vector<JaliGeometry::Point> cell_xyz,
                           std::vector<std::vector<JaliGeometry::Point>> candidate_cells_xyz,
                           std::vector<int> *xcells,
                           std::vector<std::vector<double>> *xwts) {

  int dim = cell_xyz[0].dim();
  int num_candidates = candidate_cells_xyz.size();

  xwts->clear();

  JaliGeometry::Point cmin(dim), cmax(dim);
  bounding_box(cell_xyz, &cmin, &cmax);

  for (int c = 0; c < num_candidates; ++c) {
    JaliGeometry::Point cmin2(dim), cmax2(dim);
    bounding_box(candidate_cells_xyz[c], &cmin2, &cmax2);

    std::vector<double> xsect_moments;
    if (intersect_boxes(cmin, cmax, cmin2, cmax2, &xsect_moments)) {
      xwts->push_back(xsect_moments);
      xcells->push_back(c);
    }
  }
}


}  // namespace BOX_INTERSECT

