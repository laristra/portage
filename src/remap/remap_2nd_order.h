/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef REMAP_2NDORDER_H
#define REMAP_2NDORDER_H

/*!
  \class Remap_2ndOrder remap_2nd_order.h
  \brief Remap_2ndOrder does a 2nd order remap of scalars


  [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

  [2] Dukowicz, J.K. and Kodis, J.W. "Accurate Conservative Remapping
  (Rezoning) for Arbitrary Lagrangian-Eulerian Computations," SIAM
  Journal on Scientific and Statistical Computing, Vol. 8, No. 3,
  pp. 305-321, 1987.

*/

#include <cassert>

#include "portage/support/portage.h"
#include "portage/remap/gradient.h"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

namespace Portage {

// Template on variable type ??

template<typename MeshType, typename StateType, typename OnWhatType>
class Remap_2ndOrder {
 public:
  
  Remap_2ndOrder(MeshType const & source_mesh, StateType const & source_state,
                 OnWhatType const on_what, std::string const remap_var_name,
                 LimiterType const limiter_type) :
      source_mesh_(source_mesh), 
      source_state_(source_state),
      on_what_(on_what),
      remap_var_name_(remap_var_name),
      source_vals_(NULL)
  {
    // Extract the field data from the statemanager

    source_state.get_data(on_what,remap_var_name,&source_vals_);

    // Compute the limited gradients for the field

    Limited_Gradient<MeshType,StateType,OnWhatType> 
        limgrad(source_mesh,source_state,on_what,remap_var_name,limiter_type);


    int nentities = source_mesh_.end(on_what)-source_mesh_.begin(on_what);
    gradients_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // ****************************
    // NOTE *** NOTE **** NOTE **** COMPILER IS NOT ABLE TO FIGURE OUT
    // WHETHER IT SHOULD USE THRUST::TRANSFORM OR STD::TRANSFORM AND
    // IS COMPLAINING. SO, FOR NOW I WILL JUST USE STD::TRANSFORM

    //    transform(source_mesh_.begin(on_what),source_mesh_.end(on_what),
    //              gradients_.begin(),limgrad);
    //*****************************

    std::transform(source_mesh_.begin(on_what),source_mesh_.end(on_what),
                   gradients_.begin(),limgrad);
  }


  //! Copy constructor (disabled)
  Remap_2ndOrder(const Remap_2ndOrder &) = delete;
  
  //! Assignment operator (disabled)
  Remap_2ndOrder & operator = (const Remap_2ndOrder &) = delete;

  //! Destructor
  ~Remap_2ndOrder() {}

  
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
  double * source_vals_;  // must remove assumption that field is scalar

  // Make each gradient vector of fixed length 3 even though we may
  // store only 1, 2 or 3 values in it depending on the spatial
  // dimension of the problem

  std::vector<std::vector<double>> gradients_; 

};


// Input is a pair containing the vector of contributing source
// entities and vector of contribution weights

template<typename MeshType, typename StateType, typename OnWhatType>
double Remap_2ndOrder<MeshType,StateType,OnWhatType> :: operator() 
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

  // dimension of mesh 

  int spdim = source_mesh_.space_dimension();

  double val = 0.0;
  double sumofweights = 0.0;
  if (on_what_ == CELL) {
    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)
    
    for (int j = 0; j < nsrccells; ++j) {
      int srccell = source_cells[j];
      std::vector<double> xsect_weights = weights[j];
      double xsect_volume = xsect_weights[0];
      
      std::vector<double> srccell_centroid(spdim);
      source_mesh_.cell_centroid(srccell,&srccell_centroid);
      
      std::vector<double> xsect_centroid(spdim);
      for (int i = 0; i < spdim; ++i) 
        xsect_centroid[i] = xsect_weights[1+i]/xsect_volume; // (1st moment)/(vol)
      
      val += source_vals_[srccell] * xsect_volume;
      for (int i = 0; i < spdim; ++i)
        val += gradients_[srccell][i] * (srccell_centroid[i] - xsect_centroid[i]);
      
      sumofweights += xsect_volume;
    }
  }
  else if (on_what_ == NODE) {
    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)
    
    for (int j = 0; j < nsrccells; ++j) {
      int srccell = source_cells[j];
      std::vector<double> xsect_weights = weights[j];
      double xsect_volume = xsect_weights[0];
      
      std::vector<double> srccell_centroid(spdim);
      source_mesh_.dual_cell_centroid(srccell,&srccell_centroid);
      
      std::vector<double> xsect_centroid(spdim);
      for (int i = 0; i < spdim; ++i) 
        xsect_centroid[i] = xsect_weights[1+i]/xsect_volume; // (1st moment)/(vol)
      
      val += source_vals_[srccell] * xsect_volume;
      for (int i = 0; i < spdim; ++i)
        val += gradients_[srccell][i] * (srccell_centroid[i] - xsect_centroid[i]);
      
      sumofweights += xsect_volume;
    }
  }

  // Normalize the value by sum of all the 0th weights (which is the
  // same as the total volume of the source cell)

  val /= sumofweights;

  return val;
}



} // namespace portage
                     
                     
#endif
