/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <vector>
#include "Point.hh"


namespace {

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
  bounding_box(cell_xyz,&cmin,&cmax);
  
  for (int c = 0; c < num_candidates; ++c) {
    JaliGeometry::Point cmin2(dim), cmax2(dim);
    bounding_box(candidate_cells_xyz[c],&cmin2,&cmax2);
    
    std::vector<double> xsect_moments;
    if (intersect_boxes(cmin,cmax,cmin2,cmax2,&xsect_moments)) {
      xwts->push_back(xsect_moments);
      xcells->push_back(c);
    }
  }
  
}


} // namespace
  
