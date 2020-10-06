/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_INTERPOLATE_INTERPOLATE_3RD_ORDER_H_
#define SRC_INTERPOLATE_INTERPOLATE_3RD_ORDER_H_

#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <iostream>
#include <utility>
#include <vector>

#include "wonton/support/wonton.h"

#include "portage/support/portage.h"
#include "portage/interpolate/quadfit.h"
#include "portage/driver/parts.h"

namespace Portage {

/*!
  @class Interpolate_3rdOrder interpolate_3rd_order.h
  @brief Interpolate_3rdOrder does a 3rd order interpolation of scalars
  @tparam MeshType The type of the mesh wrapper used to access mesh info
  @tparam SourceStateType The type of the state manager used to access data.
  @tparam OnWhatType The type of entity-based data we wish to interpolate;
  e.g. does it live on nodes, cells, edges, etc.

  [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

  [2] Dukowicz, J.K. and Kodis, J.W. "Accurate Conservative Remapping
  (Rezoning) for Arbitrary Lagrangian-Eulerian Computations," SIAM
  Journal on Scientific and Statistical Computing, Vol. 8, No. 3,
  pp. 305-321, 1987.

  @todo Template on variable type (YES)
*/


template<int D, Entity_kind on_what,
         typename SourceMeshType,
         typename TargetMeshType,
         typename SourceStateType,
         typename TargetStateType = SourceStateType>
class Interpolate_3rdOrder {

 public:
  /*!
    @brief Constructor
    @param[in] source_mesh The mesh wrapper used to query source mesh info
    @param[in] target_mesh The mesh wrapper used to query target mesh info
    @param[in] source_state The state-manager wrapper used to query field info
    @param[in] sources_and_weights Vector of source entities and their weights for each target entity
  */

  Interpolate_3rdOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols) {}

  /// Copy constructor (disabled)
  //  Interpolate_3rdOrder(const Interpolate_3rdOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_3rdOrder & operator = (const Interpolate_3rdOrder &) = delete;

  /// Destructor
  ~Interpolate_3rdOrder() = default;


  /// Set the name of the interpolation variable and the limiter type

  void set_interpolation_variable(std::string const & interp_var_name,
                                  Limiter_type limiter_type = NOLIMITER,
                                  Boundary_Limiter_type boundary_limiter_type = BND_NOLIMITER,
                                  Wonton::vector<Vector<D>>* gradients = nullptr) {
    interp_var_name_ = interp_var_name;

    // Extract the field data from the statemanager

    source_state_.mesh_get_data(on_what, interp_var_name, &source_vals_);

    // Compute the limited quadfits for the field

    Limited_Quadfit<D, on_what, SourceMeshType, SourceStateType>
        limqfit(source_mesh_, source_state_, interp_var_name, 
                limiter_type, boundary_limiter_type);


    int nentities = source_mesh_.end(on_what)-source_mesh_.begin(on_what);
    quadfits_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" quadfit of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Wonton::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Wonton::transform

    Wonton::transform(source_mesh_.begin(on_what), source_mesh_.end(on_what),
                       quadfits_.begin(), limqfit);
  }

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
    @param[in] targetCellID The index of the target cell.

    @todo Cleanup the datatype for sources_and_weights - it is somewhat confusing.
    @todo must remove assumption that field is scalar
  */

  double operator() (int const targetCellID,
                     std::vector<Weights_t> const & sources_and_weights) const {
    // not implemented for all types - see specialization for cells and nodes

    throw std::runtime_error("Interpolation operator not implemented for this entity type");
    return 0.;
  }

  constexpr static int order = 3;

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  SourceStateType const & source_state_;
  std::string interp_var_name_;
  double const * source_vals_;
  NumericTolerances_t num_tols_;

  // Wonton::vector is generalization of std::vector and
  // Wonton::Vector<D*(D+3)/2> is a geometric vector
  Wonton::vector<Vector<D*(D+3)/2>> quadfits_;
};




//////////////////////////////////////////////////////////////////////////////
/*!
  @brief 3rd order interpolate class specialization for cells
  @param[in] cells_and_weights Pair containing vector of contributing source
  cells and vector of contribution weights
*/

template<int D,
         typename SourceMeshType,
         typename TargetMeshType,
         typename SourceStateType,
         typename TargetStateType>
class Interpolate_3rdOrder<
  D, Entity_kind::CELL,
  SourceMeshType, TargetMeshType,
  SourceStateType, TargetStateType> {

  // useful aliases
  using Parts = PartPair<
    D, SourceMeshType, SourceStateType,
    TargetMeshType, TargetStateType
  >;

 public:
  Interpolate_3rdOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols,
                       const Parts* const parts = nullptr) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols),
      parts_(parts) {}


  /// Set the name of the interpolation variable and the limiter type

  void set_interpolation_variable(std::string const & interp_var_name,
                                  Limiter_type limiter_type = NOLIMITER,
                                  Boundary_Limiter_type boundary_limiter_type = BND_NOLIMITER,
                                  Wonton::vector<Vector<D>>* gradients = nullptr) {

    interp_var_name_ = interp_var_name;

    // Extract the field data from the statemanager

    source_state_.mesh_get_data(Entity_kind::CELL, interp_var_name, &source_vals_);

    // Compute the limited quadfits for the field

    Limited_Quadfit<D, Entity_kind::CELL, SourceMeshType, SourceStateType>
        limqfit(source_mesh_, source_state_, interp_var_name_, limiter_type, boundary_limiter_type);

    int nentities = source_mesh_.end(Entity_kind::CELL)-source_mesh_.begin(Entity_kind::CELL);
    quadfits_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" quadfit of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Wonton::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Wonton::transform

    Wonton::transform(source_mesh_.begin(Entity_kind::CELL), source_mesh_.end(Entity_kind::CELL),
                       quadfits_.begin(), limqfit);
  }

  /// Copy constructor (disabled)
  //  Interpolate_3rdOrder(const Interpolate_3rdOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_3rdOrder & operator = (const Interpolate_3rdOrder &) = delete;

  /// Destructor
  ~Interpolate_3rdOrder() = default;


  /*!
    @brief   Functor to do the 3rd order interpolation of cell values

    @param[in] targetCellID The index of the target cell.

    @param[in] sources_and_weights Vector of source mesh entities and
    corresponding weight vectors.  Each element of the weights vector
    is a moment of the source data over the target entity; for first
    order interpolation, only the first element (or zero'th moment) of
    the weights vector (i.e. the volume of intersection) is
    used. Source entities may be repeated in the list if the
    intersection of a target entity and a source entity consists of
    two or more disjoint pieces

    @todo must remove assumption that field is scalar
  */

  /*! Implementation of the () operator for 3rd order interpolation on cells
   *  Method: Uses an SVD decomposition to compute a quadratic
   *  multinomial fit to a given field, and approximates the
   *  field using the quadratic multinomial at points around
   *  the CELL center.
   */
  double operator() (int const targetCellID,
                     std::vector<Weights_t> const & sources_and_weights) const {

    int nsrccells = sources_and_weights.size();
    if (!nsrccells) {
#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
      std::cerr << "WARNING: No source cells contribute to target cell." <<
        std::endl;
#endif
      return 0.0;
    }

    double totalval = 0.0;

    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)

    /// @todo Should use zip_iterator here but I am not sure I know how to

    for (int j = 0; j < nsrccells; ++j) {
      int srccell = sources_and_weights[j].entityID;
      // int N = D*(D+3)/2;
      std::vector<double> xsect_weights = sources_and_weights[j].weights;
      double xsect_volume = xsect_weights[0];

      if (xsect_volume <= num_tols_.min_absolute_volume)
        continue;  // no intersection

      Point<D> srccell_centroid;
      source_mesh_.cell_centroid(srccell, &srccell_centroid);

      Point<D> xsect_centroid;
      for (int i = 0; i < D; ++i)
        xsect_centroid[i] = xsect_weights[1+i]/xsect_volume;  // (1st moment)/vol

      Vector<D*(D+3)/2> quadfit = quadfits_[srccell];
      Vector<D> vec = xsect_centroid - srccell_centroid;
      Vector<D*(D+3)/2> dvec;
      for (int j = 0; j < D; ++j) {
        int j1 = D;
        dvec[j] = vec[j];
        // Add the quadratic terms
        for (int k = 0; k <= j; ++k) {
          dvec[j1] = dvec[k]*dvec[j];
          j1 += 1;
          //dvec[D+j+k] = dvec[k]*dvec[j];
        }
      }
      double val = source_vals_[srccell] + dot(quadfit,dvec);
      val *= xsect_volume;
      totalval += val;
    }

    // Normalize the value by sum of all the 0th weights (which is the
    // same as the total volume of the source cell)

    totalval /= target_mesh_.cell_volume(targetCellID);

    return totalval;
  }

  constexpr static int order = 3;

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  SourceStateType const & source_state_;
  std::string interp_var_name_;
  double const * source_vals_;
  NumericTolerances_t num_tols_;
  Parts const* parts_;

  // Wonton::vector is generalization of std::vector and
  // Wonton::Vector<D> is a geometric vector
  Wonton::vector<Vector<D*(D+3)/2>> quadfits_;
};

//////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////
/*!
  @brief 3rd order interpolate class specialization for nodes
  @param[in] dualcells_and_weights Pair containing vector of contributing
  source nodes (dual cells) and vector of contribution weights
*/

template<int D,
         typename SourceMeshType,
         typename TargetMeshType,
         typename SourceStateType,
         typename TargetStateType>
class Interpolate_3rdOrder<
  D, Entity_kind::NODE,
  SourceMeshType, TargetMeshType,
  SourceStateType, TargetStateType> {

 public:
  Interpolate_3rdOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols) {}

  /// Copy constructor (disabled)
  //  Interpolate_3rdOrder(const Interpolate_3rdOrder &) = delete;

  /// Assignment operator (disabled)
  Interpolate_3rdOrder & operator = (const Interpolate_3rdOrder &) = delete;

  /// Destructor
  ~Interpolate_3rdOrder() = default;


  /// Set the name of the interpolation variable and the limiter type

  void set_interpolation_variable(std::string const & interp_var_name,
                                  Limiter_type limiter_type = NOLIMITER,
                                  Boundary_Limiter_type boundary_limiter_type = BND_NOLIMITER,
                                  Wonton::vector<Vector<D>>* gradients = nullptr) {

    interp_var_name_ = interp_var_name;

    // Extract the field data from the statemanager

    source_state_.mesh_get_data(Entity_kind::NODE, interp_var_name, &source_vals_);

    // Compute the limited quadfits for the field

    Limited_Quadfit<D, Entity_kind::NODE, SourceMeshType, SourceStateType>
        limqfit(source_mesh_, source_state_, interp_var_name, limiter_type, boundary_limiter_type);

    int nentities = source_mesh_.end(Entity_kind::NODE)-source_mesh_.begin(Entity_kind::NODE);
    quadfits_.resize(nentities);

    // call transform functor to take the values of the variable on
    // the cells and compute a "limited" quadfit of the field on the
    // cells (for transform definition, see portage.h)

    // Even though we defined Wonton::transform (to be
    // thrust::transform or boost::transform) in portage.h, the
    // compiler is not able to disambiguate this call and is getting
    // confused. So we will explicitly state that this is Wonton::transform

    Wonton::transform(source_mesh_.begin(Entity_kind::NODE), source_mesh_.end(Entity_kind::NODE),
                       quadfits_.begin(), limqfit);
  }


  /*!
    @brief Functor to do the 3rd order interpolation of node values
    @param[in] sources_and_weights      A pair of two vectors

    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.

    @c sources_and_weights.second() is the
    @param[in] targetCellID The index of the target cell.

    @param[in] sources_and_weights Vector of source mesh entities and
    corresponding weight vectors.  Each element of the weights vector
    is a moment of the source data over the target entity; for first
    order interpolation, only the first element (or zero'th moment) of
    the weights vector (i.e. the volume of intersection) is
    used. Source entities may be repeated in the list if the
    intersection of a target entity and a source entity consists of
    two or more disjoint pieces

    @todo must remove assumption that field is scalar
  */

  /*! implementation of the () operator for 3rd order interpolate on nodes
   *  Method: Uses an SVD decomposition to compute a quadratic
   *  multinomial fit to a given field, and approximates the
   *  field using the quadratic multinomial at points around
   *  the central NODE.
   */
  double operator() (const int targetNodeID,
                     std::vector<Weights_t> const & sources_and_weights) const {

    int nsrcnodes = sources_and_weights.size();
    if (!nsrcnodes) {
#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
      std::cerr << "WARNING: No source nodes contribute to target node." <<
        std::endl;
#endif
      return 0.0;
    }

    double totalval = 0.0;

    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)

    /// @todo Should use zip_iterator here but I am not sure I know how to

    for (int j = 0; j < nsrcnodes; ++j) {
      int srcnode = sources_and_weights[j].entityID;
      std::vector<double> xsect_weights = sources_and_weights[j].weights;
      double xsect_volume = xsect_weights[0];

      if (xsect_volume <= num_tols_.min_absolute_volume)
        continue;  // no intersection

      // note: here we are getting the node coord, not the centroid of
      // the dual cell
      Point<D> srcnode_coord;
      source_mesh_.node_get_coordinates(srcnode, &srcnode_coord);

      Point<D> xsect_centroid;
      for (int i = 0; i < D; ++i)
        // (1st moment)/(vol)
        xsect_centroid[i] = xsect_weights[1+i]/xsect_volume;

      Vector<D*(D+3)/2> quadfit = quadfits_[srcnode];
      Vector<D> vec = xsect_centroid - srcnode_coord;
      Vector<D*(D+3)/2> dvec;
      for (int j = 0; j < D; ++j) {
        int j1 = D;
        dvec[j] = vec[j];
        // Add the quadratic terms
        for (int k = 0; k <= j; ++k) {
          dvec[j1] = dvec[k]*dvec[j];
          j1 += 1;
          // dvec[D+j+k] = dvec[k]*dvec[j];
        }
      }
      double val = source_vals_[srcnode] + dot(quadfit,dvec);
      val *= xsect_volume;
      totalval += val;
    }

    // Normalize the value by volume of the target dual cell

    totalval /= target_mesh_.dual_cell_volume(targetNodeID);

    return totalval;
  }

  constexpr static int order = 3;

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  SourceStateType const & source_state_;
  std::string interp_var_name_;
  double const * source_vals_;
  NumericTolerances_t num_tols_;

  // Wonton::vector is generalization of std::vector and
  // Wonton::Vector<D> is a geometric vector
  Wonton::vector<Vector<D*(D+3)/2>> quadfits_;
};
}  // namespace Portage

#endif  // SRC_INTERPOLATE_INTERPOLATE_3RD_ORDER_H_
