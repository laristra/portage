/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/
/*!
  \class Remap_1stOrder remap_1st_order.cc
  \brief Remap_1stOrder does a 1st order remap of scalars. 

  Viewed simply, the value at target cell is the weighted average of
  values on from source entities and therefore, this can work for
  cell-cell, particle-cell, cell-particle and particle-particle remap.

  In the context of remapping from one cell-based mesh to another (as
  opposed to particles to cells), it is assumed that scalars are
  provided at cell centers. A piecewise constant "reconstruction" of
  the quantity is assumed over cells which means that the integral
  value over the cell or any piece of the cell is just the average
  multiplied by the volume of the cell or its piece. Then the integral
  value over a target cell is the sum of the integral values of some
  source cells (or their pieces) and the weights to be specified in
  the call are the areas/volumes of the donor cells (pieces). If an
  exact intersection is performed between the target cell and the
  source mesh, these weights are the area/volumes of the intersection
  pieces. So, in this sense, this is the Cell-Intersection-Based
  Donor-Cell (CIB/DC) remap referred to in the Shashkov, Margolin
  paper [1]. This remap is 1st order accurate and positivity
  preserving (target cell values will be positive if the field is
  positive on the source mesh).

  Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

*/


#include "remap_1st_order.h"

#include "portage/state/state.h"
#include "portage/state/state_vector.h"

namespace Portage {

// Input is a pair containing the vector of contributing source
// entities and vector of contribution weights

double Remap_1stOrder::operator() (std::pair< std::vector<int> const & , std::vector<double> const & > cells_and_weights) const {

  std::vector<int> const & source_cells = cells_and_weights.first; 
  int nsrccells = source_cells.size();
  if (!nsrccells) {
    std::cerr << "ERROR: No source cells contribute to target cell?" << std::endl;
    return 0.0;
  }

  std::vector<double> const & weights = cells_and_weights.second;
  if (weights.size() != nsrccells) {
    std::cerr << "ERROR: Not enough weights provided for remapping " << std::endl;
    return 0.0;
  }

  // contribution of the source cell is its field value weighted by
  // its "weight" (in this case, its 0th moment/area/volume)

  double val = 0.0;
  double sumofweights = 0.0;
  for (int j = 0; j < nsrccells; ++j) {
    int srccell = source_cells[j];
    double srcval = (*source_var_ptr_)[srccell];

    val += srcval*weights[j];
    sumofweights += weights[j];
  }

  // Normalize the value by sum of all the weights

  val /= sumofweights;

  return val;
}

}
