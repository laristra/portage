/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef REMAP_1STORDER_H
#define REMAP_1STORDER_H

/*!
  \class Remap_1stOrder remap_1st_order.h
  \brief Remap_1stOrder does a 1st order remap of scalars

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

  [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

*/

#include <cassert>

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

namespace Portage {

// Template on variable type ??

template<typename MeshType, typename StateType, typename OnWhatType>
class Remap_1stOrder {
 public:
  
  Remap_1stOrder(MeshType const & source_mesh, StateType const & source_state,
                 OnWhatType const on_what, std::string const remap_var_name) :
      source_mesh_(source_mesh), 
      source_state_(source_state),
      on_what_(on_what),
      remap_var_name_(remap_var_name),
      source_vals_(NULL)
  {
    source_state.get_data(on_what,remap_var_name,&source_vals_);
  }


  //! Copy constructor (disabled)
  Remap_1stOrder(const Remap_1stOrder &) = delete;
  
  //! Assignment operator (disabled)
  Remap_1stOrder & operator = (const Remap_1stOrder &) = delete;

  //! Destructor
  ~Remap_1stOrder() {}

  
  // Remap functor - Need to make a pair from the vector of source
  // cells that contribute to the the target cell value and the
  // contribution weights associated with each source cell. Source
  // cells may be repeated in the list if the intersection of a target
  // cell and a source cell consists of two or more disjoint pieces

  double 
  operator() (std::pair<std::vector<int> const &, 
              std::vector< std::vector<double> > const &> cells_and_weights) const;

 private:

  MeshType const & source_mesh_;
  StateType const & source_state_;
  OnWhatType const on_what_;
  std::string const & remap_var_name_;
  double * source_vals_;

};


// Input is a pair containing the vector of contributing source
// entities and vector of contribution weights

template<typename MeshType, typename StateType, typename OnWhatType>
double Remap_1stOrder<MeshType,StateType,OnWhatType> :: operator() 
    (std::pair<std::vector<int> const &, 
     std::vector< std::vector<double> > const &> cells_and_weights) const {

  std::vector<int> const & source_cells = cells_and_weights.first; 
  int nsrccells = source_cells.size();
  if (!nsrccells) {
    std::cerr << "ERROR: No source cells contribute to target cell?" << 
        std::endl;
    return 0.0;
  }

  std::vector< std::vector<double> > const & weights = cells_and_weights.second;
  if (weights.size() < nsrccells) {
    std::cerr << "ERROR: Not enough weights provided for remapping " << 
        std::endl;
    return 0.0;
  }

  // contribution of the source cell is its field value weighted by
  // its "weight" (in this case, its 0th moment/area/volume)

  double val = 0.0;
  double sumofweights = 0.0;
  for (int j = 0; j < nsrccells; ++j) {
    int srccell = source_cells[j];
    std::vector<double> pair_weights = weights[j];

    val += source_vals_[srccell] * pair_weights[0];  // 1st order 
    sumofweights += pair_weights[0];
  }

  // Normalize the value by sum of all the weights

  val /= sumofweights;

  return val;
}


} // namespace portage
                     
                     
#endif
