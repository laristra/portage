/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SRC_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_
#define SRC_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_

#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <iostream>
#include <utility>
#include <vector>

#include "portage/support/portage.h"
#include "portage/interpolate/gradient.h"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

namespace Portage {

/*!
  @class Interpolate_2ndOrder interpolate_2nd_order.h
  @brief Interpolate_2ndOrder does a 2nd order interpolation of scalars
  @tparam MeshType The type of the mesh wrapper used to access mesh info
  @tparam StateType The type of the state manager used to access data.
  @tparam OnWhatType The type of entity-based data we wish to interpolate;
  e.g. does it live on nodes, cells, edges, etc.

  [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

  [2] Dukowicz, J.K. and Kodis, J.W. "Accurate Conservative Remapping
  (Rezoning) for Arbitrary Lagrangian-Eulerian Computations," SIAM
  Journal on Scientific and Statistical Computing, Vol. 8, No. 3,
  pp. 305-321, 1987.

  @todo Template on variable type ??
*/


template<typename SourceMeshType, typename TargetMeshType, typename StateType,
    Entity_kind on_what>
class Interpolate_2ndOrder {
 public:
  /*!
    @brief Constructor
    @param[in] source_mesh The mesh wrapper used to query source mesh info
    @param[in] target_mesh The mesh wrapper used to query target mesh info
    @param[in] source_state The state-manager wrapper used to query field info
    @param[in] interp_var_name Name of the field to be remapped
    @param[in] limiter_type Gradient limiter type (see gradient.h)
  */

  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state,
                       std::string const interp_var_name,
                       LimiterType const limiter_type) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_(interp_var_name),
      source_vals_(NULL) {
    
    // Extract the field data from the statemanager

    source_state.get_data(on_what, interp_var_name, &source_vals_);

    // Compute the limited gradients for the field

    Limited_Gradient<SourceMeshType, StateType, on_what>
        limgrad(source_mesh, source_state, interp_var_name, limiter_type);


    int nentities = source_mesh_.end(on_what)-source_mesh_.begin(on_what);
    gradients_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform

    Portage::transform(source_mesh_.begin(on_what), source_mesh_.end(on_what),
                       gradients_.begin(), limgrad);
  }


  /// Copy constructor (disabled)
  Interpolate_2ndOrder(const Interpolate_2ndOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_2ndOrder & operator = (const Interpolate_2ndOrder &) = delete;

  /// Destructor
  ~Interpolate_2ndOrder() {}


  /*!
    @brief Functor to do the actual interpolate calculation
    @param[in] cells_and_weights A pair of two vectors
    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.
    @c sources_and_weights.second() is the vector of vector weights for each
    of the source mesh entities in @c sources_and_weights.first().  Each element
    of the weights vector is a moment of the source data over the target
    entity; for first order interpolate, only the first element (or zero'th moment)
    of the weights vector (i.e. the volume of intersection) is used. Source
    entities may be repeated in the list if the intersection of a target entity
    and a source entity consists of two or more disjoint pieces
    @param[in] targetCellId The index of the target cell.

    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.
    @todo must remove assumption that field is scalar
  */

  double
  operator() (std::pair<std::vector<int> const &,
              std::vector<std::vector<double>> const &> sources_and_weights,
              const int targetCellId)
      const {
    // not implemented for all types - see specialization for cells and nodes

    std::cerr << "Interpolation operator not implemented for this entity type"
              << std::endl;
  }

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string const & interp_var_name_;
  double * source_vals_;

  std::vector<std::vector<double>> gradients_;
};




//////////////////////////////////////////////////////////////////////////////
/*!
  @brief 2nd order interpolate class specialization for cells
  @param[in] cells_and_weights Pair containing vector of contributing source
  cells and vector of contribution weights
*/

template<typename SourceMeshType, typename TargetMeshType, typename StateType>
class Interpolate_2ndOrder<SourceMeshType, TargetMeshType, StateType, CELL> {
 public:
  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state,
                       std::string const interp_var_name,
                       LimiterType const limiter_type) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_(interp_var_name),
      limiter_type_(limiter_type),
      source_vals_(NULL) {
    // Extract the field data from the statemanager

    source_state_.get_data(CELL, interp_var_name_, &source_vals_);
  }

  void compute_gradient() {
    // Compute the limited gradients for the field

    int nentities = source_mesh_.end(CELL, Entity_type::PARALLEL_OWNED)-source_mesh_.begin(CELL);
    gradients_ = std::make_shared<std::vector<double>>();
    gradients_->resize(3*nentities);
   
    Limited_Gradient<SourceMeshType, StateType, CELL>
        limgrad(source_mesh_, source_state_, interp_var_name_, &((*gradients_)[0]), limiter_type_);
 
    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform
    std::for_each(source_mesh_.begin(CELL), source_mesh_.end(CELL, Entity_type::PARALLEL_OWNED),
                      limgrad);
  }


  /// Copy constructor (disabled)
  Interpolate_2ndOrder(const Interpolate_2ndOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_2ndOrder & operator = (const Interpolate_2ndOrder &) = delete;

  /// Destructor
  ~Interpolate_2ndOrder() {}


  /*!
    @brief   Functor to do the 2nd order interpolation of cell values
    @param[in] cells_and_weights           A pair of two vectors
    @c sources_and_weights.first() is the vector of source entity
    indices in the source mesh that will contribute to the current
    target mesh entity.  @c sources_and_weights.second() is the vector
    of vector weights for each of the source mesh entities in @c
    sources_and_weights.first().  Each element of the weights vector
    is a moment of the source data over the target entity; for first
    order interpolation, only the first element (or zero'th moment) of the
    weights vector (i.e. the volume of intersection) is used. Source
    entities may be repeated in the list if the intersection of a
    target entity and a source entity consists of two or more disjoint
    pieces
    @param[in] targetCellId The index of the target cell.

    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.

    @todo SHOULD WE USE boost::tuple FOR SENDING IN SOURCES_AND_WEIGHTS SO
    THAT WE CAN USE boost::zip_iterator TO ITERATOR OVER IT?

    @todo must remove assumption that field is scalar
  */

  double
  operator() (std::pair<std::vector<int> const &,
              std::vector<std::vector<double>> const &> sources_and_weights,
              const int targetCellId)
      const;

  void set_gradients(std::shared_ptr<std::vector<double>> gradients)
  {
    gradients_ = gradients;
  }

  std::shared_ptr<std::vector<double>> get_gradients()
  {
    return gradients_;
  }  

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string const & interp_var_name_;
  LimiterType const & limiter_type_;
  double * source_vals_;
  std::shared_ptr<std::vector<double>> gradients_; 
};

// Implementation of the () operator for 2nd order interpolation on cells


template<typename SourceMeshType, typename TargetMeshType, typename StateType>
double Interpolate_2ndOrder<SourceMeshType, TargetMeshType,
    StateType, CELL>::operator()
    (std::pair<std::vector<int> const &,
     std::vector< std::vector<double> > const &> cells_and_weights,
     const int targetCellId) const {
  std::vector<int> const & source_cells = cells_and_weights.first;
  int nsrccells = source_cells.size();
  if (!nsrccells) {
    std::cerr << "WARNING: No source cells contribute to target cell." <<
        std::endl;
    return 0.0;
  }

  std::vector< std::vector<double> > const & weights = cells_and_weights.second;
  if (weights.size() < nsrccells) {
    std::cerr << "ERROR: Not enough weights provided for interpolating " <<
        std::endl;
    return 0.0;
  }

  // dimension of mesh

  int spdim = source_mesh_.space_dimension();

  double totalval = 0.0;

  // contribution of the source cell is its field value weighted by
  // its "weight" (in this case, its 0th moment/area/volume)

  /// @todo Should use zip_iterator here but I am not sure I know how to

  for (int j = 0; j < nsrccells; ++j) {
    int srccell = source_cells[j];
    std::vector<double> xsect_weights = weights[j];
    double xsect_volume = xsect_weights[0];

    double eps = 1e-30;
    if (xsect_volume <= eps) continue;  // no intersection

    std::vector<double> srccell_centroid(spdim);
    source_mesh_.cell_centroid(srccell, &srccell_centroid);

    std::vector<double> xsect_centroid(spdim);
    for (int i = 0; i < spdim; ++i)
      xsect_centroid[i] = xsect_weights[1+i]/xsect_volume;  // (1st moment)/vol


    double val = source_vals_[srccell];
    for (int i = 0; i < spdim; ++i)
      val += (*gradients_)[srccell*3+i] * (xsect_centroid[i]-srccell_centroid[i]);
    val *= xsect_volume;
    totalval += val;
  }

  // Normalize the value by sum of all the 0th weights (which is the
  // same as the total volume of the source cell)

  totalval /= target_mesh_.cell_volume(targetCellId);

  return totalval;
}
//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
/*!
  @brief 2nd order interpolate class specialization for nodes
  @param[in] dualcells_and_weights Pair containing vector of contributing
  source nodes (dual cells) and vector of contribution weights
*/

template<typename SourceMeshType, typename TargetMeshType, typename StateType>
class Interpolate_2ndOrder<SourceMeshType, TargetMeshType, StateType, NODE> {
 public:
  Interpolate_2ndOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state,
                       std::string const interp_var_name,
                       LimiterType const limiter_type) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_(interp_var_name),
      source_vals_(NULL) {
  // Extract the field data from the statemanager

  source_state.get_data(NODE, interp_var_name, &source_vals_);

  // Compute the limited gradients for the field

  Limited_Gradient<SourceMeshType, StateType, NODE>
      limgrad(source_mesh, source_state, interp_var_name, limiter_type);

  int nentities = source_mesh_.end(NODE)-source_mesh_.begin(NODE);
  gradients_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform

    Portage::transform(source_mesh_.begin(NODE), source_mesh_.end(NODE),
                       gradients_.begin(), limgrad);
  }


  /// Copy constructor (disabled)
  Interpolate_2ndOrder(const Interpolate_2ndOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_2ndOrder & operator = (const Interpolate_2ndOrder &) = delete;

  /// Destructor
  ~Interpolate_2ndOrder() {}


  /*!
    @brief Functor to do the 2nd order interpolation of node values
    @param[in] sources_and_weights      A pair of two vectors

    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.

    @c sources_and_weights.second() is the vector of vector weights
    for each of the source mesh entities in @c
    sources_and_weights.first().  Each element of the weights vector
    is a moment of the source data over the target entity; for first
    order interpolation, only the first element (or zero'th moment) of
    the weights vector (i.e. the volume of intersection) is
    used. Source entities may be repeated in the list if the
    intersection of a target entity and a source entity consists of
    two or more disjoint pieces

    @param[in] targetCellId The index of the target cell.

    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.
    @todo must remove assumption that field is scalar
  */

  double
  operator() (std::pair<std::vector<int> const &,
              std::vector< std::vector<double> > const &> sources_and_weights,
              const int targetCellId)
      const;

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string const & interp_var_name_;
  double * source_vals_;

  std::vector<std::vector<double>> gradients_;
};

/// implementation of the () operator for 2nd order interpolate on nodes

template<typename SourceMeshType, typename TargetMeshType, typename StateType>
double Interpolate_2ndOrder<SourceMeshType, TargetMeshType,
    StateType, NODE> :: operator()
    (std::pair<std::vector<int> const &,
     std::vector< std::vector<double> > const &> dualcells_and_weights,
     const int targetCellId) const {
  std::vector<int> const & source_cells = dualcells_and_weights.first;
  int nsrccells = source_cells.size();
  if (!nsrccells) {
    std::cerr << "WARNING: No source cells contribute to target cell." <<
        std::endl;
    return 0.0;
  }

  std::vector<std::vector<double>> const & weights =
      dualcells_and_weights.second;
  if (weights.size() < nsrccells) {
    std::cerr << "ERROR: Not enough weights provided for interpolating " <<
        std::endl;
    return 0.0;
  }

  // dimension of mesh

  int spdim = source_mesh_.space_dimension();

  double totalval = 0.0;

  // contribution of the source cell is its field value weighted by
  // its "weight" (in this case, its 0th moment/area/volume)

  /// @todo Should use zip_iterator here but I am not sure I know how to

  for (int j = 0; j < nsrccells; ++j) {
    int srccell = source_cells[j];
    std::vector<double> xsect_weights = weights[j];
    double xsect_volume = xsect_weights[0];

    double eps = 1e-30;
    if (xsect_volume <= eps) continue;  // no intersection

    // note: here we are getting the node coord, not the centroid of
    // the dual cell
    std::vector<double> srccell_coord(spdim);
    if (spdim == 2) {
      Point<2> point;
      source_mesh_.node_get_coordinates(srccell, &point);
      srccell_coord[0] = point[0];
      srccell_coord[1] = point[1];
    } else if (spdim == 3) {
      Point<3> point;
      source_mesh_.node_get_coordinates(srccell, &point);
      srccell_coord[0] = point[0];
      srccell_coord[1] = point[1];
      srccell_coord[2] = point[2];
    }

    std::vector<double> xsect_centroid(spdim);
    for (int i = 0; i < spdim; ++i)
      // (1st moment)/(vol)
      xsect_centroid[i] = xsect_weights[1+i]/xsect_volume;

    double val = source_vals_[srccell];
    for (int i = 0; i < spdim; ++i)
      val += gradients_[srccell][i] * (xsect_centroid[i]-srccell_coord[i]);
    val *= xsect_volume;
    totalval += val;
  }

  // Normalize the value by sum of all the 0th weights (which is the
  // same as the total volume of the source cell)

  totalval /= target_mesh_.cell_volume(targetCellId);

  return totalval;
}


}  // namespace Portage

#endif  // SRC_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_
