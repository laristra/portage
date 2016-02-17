/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef REMAP_2NDORDER_H
#define REMAP_2NDORDER_H

#include <cassert>
#include <stdexcept>

#include "portage/support/portage.h"
#include "portage/remap/gradient.h"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

namespace Portage {

/*!
  @class Remap_2ndOrder remap_2nd_order.h
  @brief Remap_2ndOrder does a 2nd order remap of scalars
  @tparam MeshType The type of the mesh wrapper used to access mesh info
  @tparam StateType The type of the state manager used to access data.
  @tparam OnWhatType The type of entity-based data we wish to remap; e.g. does
  it live on nodes, cells, edges, etc.

  [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

  [2] Dukowicz, J.K. and Kodis, J.W. "Accurate Conservative Remapping
  (Rezoning) for Arbitrary Lagrangian-Eulerian Computations," SIAM
  Journal on Scientific and Statistical Computing, Vol. 8, No. 3,
  pp. 305-321, 1987.

  @todo Template on variable type ??
*/


template<typename MeshType, typename StateType, Entity_kind on_what>
class Remap_2ndOrder {
 public:
  
  /*! 
    @brief Constructor
    @param[in] source_mesh The mesh wrapper used to query source mesh info
    @param[in] source_state The state-manager wrapper used to query field info
    @param[in] remap_var_name Name of the field to be remapped
    @param[in] limiter_type Gradient limiter type (see gradient.h)
  */

  Remap_2ndOrder(MeshType const & source_mesh, StateType const & source_state,
                 std::string const remap_var_name,
                 LimiterType const limiter_type) :
      source_mesh_(source_mesh), 
      source_state_(source_state),
      remap_var_name_(remap_var_name),
      source_vals_(NULL)
  {
    // Extract the field data from the statemanager

    source_state.get_data(on_what,remap_var_name,&source_vals_);

    // Compute the limited gradients for the field

    Limited_Gradient<MeshType,StateType,on_what> 
        limgrad(source_mesh,source_state,remap_var_name,limiter_type);


    int nentities = source_mesh_.end(on_what)-source_mesh_.begin(on_what);
    gradients_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform

    Portage::transform(source_mesh_.begin(on_what),source_mesh_.end(on_what),
                       gradients_.begin(),limgrad);

  }


  /// Copy constructor (disabled)
  Remap_2ndOrder(const Remap_2ndOrder &) = delete;
  
  /// Assignment operator (disabled)
  Remap_2ndOrder & operator = (const Remap_2ndOrder &) = delete;

  /// Destructor
  ~Remap_2ndOrder() {}

  
  /*! 
    @brief Functor to do the actual remap calculation
    @param[in] cells_and_weights A pair of two vectors
    @c sources_and_weights.first() is the vector of source entity indices 
    in the source mesh that will contribute to the current target mesh entity.
    @c sources_and_weights.second() is the vector of vector weights for each
    of the source mesh entities in @c sources_and_weights.first().  Each element
    of the weights vector is a moment of the source data over the target
    entity; for first order remap, only the first element (or zero'th moment)
    of the weights vector (i.e. the volume of intersection) is used. Source 
    entities may be repeated in the list if the intersection of a target entity
    and a source entity consists of two or more disjoint pieces
    
    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.
    @todo must remove assumption that field is scalar
  */

  double 
  operator() (std::pair<std::vector<int> const &, 
              std::vector< std::vector<double> > const &> sources_and_weights) const {
    // not implemented for all types - see specialization for cells and nodes
    
    std::cerr << "Remap operator not implemented for this entity type\n";
  }

 private:

  MeshType const & source_mesh_;
  StateType const & source_state_;
  std::string const & remap_var_name_;
  double * source_vals_;  

  std::vector<std::vector<double>> gradients_; 

};




//////////////////////////////////////////////////////////////////////////////
/*! 
  @brief 2nd order remap class specialization for cells
  @param[in] cells_and_weights Pair containing vector of contributing source 
  cells and vector of contribution weights
*/

template<typename MeshType, typename StateType>
class Remap_2ndOrder<MeshType,StateType,CELL> {
 public:
  
  Remap_2ndOrder(MeshType const & source_mesh, StateType const & source_state,
                 std::string const remap_var_name,
                 LimiterType const limiter_type) :
      source_mesh_(source_mesh), 
      source_state_(source_state),
      remap_var_name_(remap_var_name),
      source_vals_(NULL)
  {
    // Extract the field data from the statemanager

    source_state.get_data(CELL,remap_var_name,&source_vals_);

    // Compute the limited gradients for the field

    Limited_Gradient<MeshType,StateType,CELL> 
        limgrad(source_mesh,source_state,remap_var_name,limiter_type);


    int nentities = source_mesh_.end(CELL)-source_mesh_.begin(CELL);
    gradients_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform

    Portage::transform(source_mesh_.begin(CELL),source_mesh_.end(CELL),
                       gradients_.begin(),limgrad);

  }


  /// Copy constructor (disabled)
  Remap_2ndOrder(const Remap_2ndOrder &) = delete;
  
  /// Assignment operator (disabled)
  Remap_2ndOrder & operator = (const Remap_2ndOrder &) = delete;

  /// Destructor
  ~Remap_2ndOrder() {}

  
  /*! 
    @brief Functor to do the actual remap calculation
    @param[in] cells_and_weights A pair of two vectors
    @c sources_and_weights.first() is the vector of source entity
    indices in the source mesh that will contribute to the current
    target mesh entity.  @c sources_and_weights.second() is the vector
    of vector weights for each of the source mesh entities in @c
    sources_and_weights.first().  Each element of the weights vector
    is a moment of the source data over the target entity; for first
    order remap, only the first element (or zero'th moment) of the
    weights vector (i.e. the volume of intersection) is used. Source
    entities may be repeated in the list if the intersection of a
    target entity and a source entity consists of two or more disjoint
    pieces
    
    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.

    @todo SHOULD WE USE boost::tuple FOR SENDING IN SOURCES_AND_WEIGHTS SO 
    THAT WE CAN USE boost::zip_iterator TO ITERATOR OVER IT? 

    @todo must remove assumption that field is scalar
  */

  double 
  operator() (std::pair<std::vector<int> const &, 
              std::vector< std::vector<double> > const &> sources_and_weights) const;

 private:

  MeshType const & source_mesh_;
  StateType const & source_state_;
  std::string const & remap_var_name_;
  double * source_vals_;  

  std::vector<std::vector<double>> gradients_; 

};

/// @brief implementation of the () operator for 2nd order remap on cells


template<typename MeshType, typename StateType>
double Remap_2ndOrder<MeshType,StateType,CELL> :: operator() 
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

  // contribution of the source cell is its field value weighted by
  // its "weight" (in this case, its 0th moment/area/volume)

  /// @todo Should use zip_iterator here but I am not sure I know how to correctly

  for (int j = 0; j < nsrccells; ++j) {
    int srccell = source_cells[j];
    std::vector<double> xsect_weights = weights[j];
    double xsect_volume = xsect_weights[0];

    double eps = 1e-30;
    if (xsect_volume <= eps) continue; // no intersection
    
    std::vector<double> srccell_centroid(spdim);
    source_mesh_.cell_centroid(srccell,&srccell_centroid);
    
    std::vector<double> xsect_centroid(spdim);
    for (int i = 0; i < spdim; ++i) 
      xsect_centroid[i] = xsect_weights[1+i]/xsect_volume; // (1st moment)/(vol)
    
    val += source_vals_[srccell] * xsect_volume;
    for (int i = 0; i < spdim; ++i)
      val += gradients_[srccell][i] * (xsect_centroid[i] - srccell_centroid[i]) * xsect_volume;
    
    sumofweights += xsect_volume;
  }

  // Normalize the value by sum of all the 0th weights (which is the
  // same as the total volume of the source cell)

  val /= sumofweights;

  return val;
}
//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
/*! 
  @brief 2nd order remap class specialization for nodes
  @param[in] dualcells_and_weights Pair containing vector of contributing 
  source nodes (dual cells) and vector of contribution weights
*/

template<typename MeshType, typename StateType>
class Remap_2ndOrder<MeshType,StateType,NODE> {
 public:
  
  Remap_2ndOrder(MeshType const & source_mesh, StateType const & source_state,
                 std::string const remap_var_name,
                 LimiterType const limiter_type) :
      source_mesh_(source_mesh), 
      source_state_(source_state),
      remap_var_name_(remap_var_name),
      source_vals_(NULL)
  {
    // Extract the field data from the statemanager

    source_state.get_data(NODE,remap_var_name,&source_vals_);

    // Compute the limited gradients for the field

    Limited_Gradient<MeshType,StateType,NODE> 
        limgrad(source_mesh,source_state,remap_var_name,limiter_type);


    int nentities = source_mesh_.end(NODE)-source_mesh_.begin(NODE);
    gradients_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" gradient of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Portage::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Portage::transform

    Portage::transform(source_mesh_.begin(NODE),source_mesh_.end(NODE),
                       gradients_.begin(),limgrad);

  }


  /// Copy constructor (disabled)
  Remap_2ndOrder(const Remap_2ndOrder &) = delete;
  
  /// Assignment operator (disabled)
  Remap_2ndOrder & operator = (const Remap_2ndOrder &) = delete;

  /// Destructor
  ~Remap_2ndOrder() {}

  
  /*! 
    @brief Functor to do the actual remap calculation
    @param[in] sources_and_weights A pair of two vectors
    @c sources_and_weights.first() is the vector of source entity indices 
    in the source mesh that will contribute to the current target mesh entity.
    @c sources_and_weights.second() is the vector of vector weights for each
    of the source mesh entities in @c sources_and_weights.first().  Each element
    of the weights vector is a moment of the source data over the target
    entity; for first order remap, only the first element (or zero'th moment)
    of the weights vector (i.e. the volume of intersection) is used. Source 
    entities may be repeated in the list if the intersection of a target entity
    and a source entity consists of two or more disjoint pieces
    
    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.
    @todo must remove assumption that field is scalar
  */

  double 
  operator() (std::pair<std::vector<int> const &, 
              std::vector< std::vector<double> > const &> sources_and_weights) const;

 private:

  MeshType const & source_mesh_;
  StateType const & source_state_;
  std::string const & remap_var_name_;
  double * source_vals_;  

  std::vector<std::vector<double>> gradients_; 

};

/// @brief implementation of the () operator for 2nd order remap on nodes

template<typename MeshType, typename StateType>
double Remap_2ndOrder<MeshType,StateType,NODE> :: operator() 
    (std::pair<std::vector<int> const &, 
     std::vector< std::vector<double> > const &> dualcells_and_weights) const {

  std::vector<int> const & source_cells = dualcells_and_weights.first; 
  int nsrccells = source_cells.size();
  if (!nsrccells) {
    std::cerr << "ERROR: No source cells contribute to target cell?" << 
        std::endl;
    return 0.0;
  }

  std::vector< std::vector<double> > const & weights = dualcells_and_weights.second;
  if (weights.size() < nsrccells) {
    std::cerr << "ERROR: Not enough weights provided for remapping " << 
        std::endl;
    return 0.0;
  }

  // dimension of mesh 

  int spdim = source_mesh_.space_dimension();

  double val = 0.0;
  double sumofweights = 0.0;

  // contribution of the source cell is its field value weighted by
  // its "weight" (in this case, its 0th moment/area/volume)

  /// @todo Should use zip_iterator here but I am not sure I know how to correctly
  
  for (int j = 0; j < nsrccells; ++j) {
    int srccell = source_cells[j];
    std::vector<double> xsect_weights = weights[j];
    double xsect_volume = xsect_weights[0];

    double eps = 1e-30;
    if (xsect_volume <= eps) continue; // no intersection
    
    // note: here we are getting the node coord, not the centroid of the dual cell
    std::vector<double> srccell_coord(spdim); 
    if (spdim == 2) {
      std::pair<double,double> coord_pair;
      source_mesh_.node_get_coordinates(srccell,&coord_pair);
      srccell_coord[0] = coord_pair.first;
      srccell_coord[1] = coord_pair.second;
    }
    else if (spdim == 3) {
      std::tuple<double,double,double> coord_tuple;
      source_mesh_.node_get_coordinates(srccell,&coord_tuple);
      srccell_coord[0] = std::get<0>(coord_tuple);
      srccell_coord[1] = std::get<1>(coord_tuple);
      srccell_coord[2] = std::get<2>(coord_tuple);
    }

    std::vector<double> xsect_centroid(spdim);
    for (int i = 0; i < spdim; ++i) 
      xsect_centroid[i] = xsect_weights[1+i]/xsect_volume; // (1st moment)/(vol)
    
    val += source_vals_[srccell] * xsect_volume;
    for (int i = 0; i < spdim; ++i)
      val += gradients_[srccell][i] * (xsect_centroid[i] - srccell_coord[i]) * xsect_volume;
    
    sumofweights += xsect_volume;
  }

  // Normalize the value by sum of all the 0th weights (which is the
  // same as the total volume of the source cell)

  val /= sumofweights;

  return val;
}


} // namespace portage
                     
                     
#endif
