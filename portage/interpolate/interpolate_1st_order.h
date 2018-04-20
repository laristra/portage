/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERPOLATE_INTERPOLATE_1ST_ORDER_H_
#define PORTAGE_INTERPOLATE_INTERPOLATE_1ST_ORDER_H_


#include <cassert>
#include <string>
#include <iostream>
#include <utility>
#include <vector>

#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/support/portage.h"


namespace Portage {

/*!
  @class Interpolate_1stOrder interpolate_1st_order.h
  @brief Interpolate_1stOrder does a 1st order interpolation of scalars
  @tparam MeshType The type of the mesh wrapper used to access mesh info
  @tparam StateType The type of the state manager used to access data.
  @tparam OnWhatType The type of entity-based data we wish to interpolate;
  e.g. does it live on nodes, cells, edges, etc.

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
  paper [1]. This interpolation is 1st order accurate and positivity
  preserving (target cell values will be positive if the field is
  positive on the source mesh).

  [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
  conservative interpolation (remapping) on general grids." Journal of
  Computational Physics, v 184, n 1, pp. 266-298, 2003.

  [2] Dukowicz, J.K. and Kodis, J.W. "Accurate Conservative Remapping
  (Rezoning) for Arbitrary Lagrangian-Eulerian Computations," SIAM
  Journal on Scientific and Statistical Computing, Vol. 8, No. 3,
  pp. 305-321, 1987.

  @todo Template on variable type??

*/

template<int D, 
         Entity_kind on_what,
         typename SourceMeshType, 
         typename TargetMeshType, 
         typename StateType,
         template<class, int> class InterfaceReconstructorType =
          DummyInterfaceReconstructor>
class Interpolate_1stOrder {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType>;
#endif

 public:

  //Constructor with interface reconstructor
#ifdef HAVE_TANGRAM
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state, 
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interface_reconstructor_(ir)
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr) {}
#endif 

  /*!
    @brief Constructor without interface reconstructor.
    @param[in] source_mesh The input mesh.
    @param[in] target_mesh The output mesh.
    @param[in] source_state The state manager for data on the input mesh.
    @param[in] on_what The location where the data lives; e.g. on cells, nodes,
    edges, etc.
    @param[in] interp_var_name The string name of the variable to interpolate.
    @param[in] sources_and_weights Vector of source entities and their weights for each target entity
  */
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr) {}


  /// Copy constructor (disabled)
  //  Interpolate_1stOrder(const Interpolate_1stOrder &) = delete;

  /// Assignment operator (disabled)
  //  Interpolate_1stOrder & operator = (const Interpolate_1stOrder &) = delete;

  /// Destructor
  ~Interpolate_1stOrder() {}

  /// Set the material we are operating on

  int set_material(int m) {
    matid_ = m;
  } //set_material

  /// Set the variable name to be interpolated.

  // Even though 1st order accurate interpolation does not require
  // limiters and only higher order ones do, we need to have the
  // limiter_type as an argument. This is because the templated driver
  // does not know whether its being called with 1st order or higher
  // order interpolators and so all interpolators need to have a
  // uniform interface

  void set_interpolation_variable(std::string const & interp_var_name,
                                  LimiterType limtype=NOLIMITER) {
    interp_var_name_ = interp_var_name;
    field_type_ = source_state_.field_type(CELL, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(CELL, interp_var_name,
                                  &source_vals_);
    else
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
  } //set_interpolation_variable

  /*!
    @brief Functor to do the actual interpolation.
    @param[in] sources_and_weights A pair of two vectors.
    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.
    @c sources_and_weights.second() is the vector of vector weights for each
    of the source mesh entities in @c sources_and_weights.first().  Each element
    of the weights vector is a moment of the source data over the target
    entity; for first order interpolation, only the first element (or zero'th
    moment) of the weights vector (i.e. the volume of intersection) is used.
    Source entities may be repeated in the list if the intersection of a target
    entity and a source entity consists of two or more disjoint pieces
    @param[in] targetCellId The index of the target cell.

  */

  double operator() (int const targetEntityId,
                     std::vector<Weights_t> const & sources_and_weights) const {
    std::cerr << "Interpolation operator not implemented for this entity type"
              << std::endl;
    return 0.0;
  }

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string interp_var_name_;
  double const * source_vals_;
  int matid_;
  Field_type field_type_;


#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
}; //interpolate_1st_order base



//////////////////////////////////////////////////////////////////////////////
/*!
  @brief Interpolate_1stOrder specialization for cells
*/

template<int D,
         typename SourceMeshType, 
         typename TargetMeshType, 
         typename StateType,
         template<class, int> class InterfaceReconstructorType>
class Interpolate_1stOrder<D, 
                           CELL, 
                           SourceMeshType, 
                           TargetMeshType, 
                           StateType,
                           InterfaceReconstructorType> {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType>;
#endif
                            
 public:

#ifdef HAVE_TANGRAM
   Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state,
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interface_reconstructor_(ir),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr) {}
#endif 
  /*!
    @brief Constructor.
    @param[in] source_mesh The input mesh.
    @param[in] target_mesh The output mesh.
    @param[in] source_state The state manager for data on the input mesh.
    @param[in] interp_var_name The string name of the variable to interpolate.
    @param[in] sources_and_weights Vector of source entities and their weights for each target entity
  */
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr) {}


  /// Copy constructor (disabled)
  //  Interpolate_1stOrder(const Interpolate_1stOrder &) = delete;

  /// Assignment operator (disabled)
  //  Interpolate_1stOrder & operator = (const Interpolate_1stOrder &) = delete;

  /// Destructor
  ~Interpolate_1stOrder() {}


  /// Set the material we are operating on

  int set_material(int m) {
    matid_ = m;
  } //set_material

  /// Set the variable name to be interpolated


  // Even though 1st order accurate interpolation does not require
  // limiters and only higher order ones do, we need to have the
  // limiter_type as an argument. This is because the templated driver
  // does not know whether its being called with 1st order or higher
  // order interpolators and so all interpolators need to have a
  // uniform interface

  void set_interpolation_variable(std::string const & interp_var_name, 
                                  LimiterType limtype = NOLIMITER) {
    interp_var_name_ = interp_var_name;
    field_type_ = source_state_.field_type(CELL, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(CELL, interp_var_name,
                                  &source_vals_);
    else
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
  } //set_interpolation_variable

  /*!
    @brief Functor to do the actual interpolation.
    @param[in] sources_and_weights A pair of two vectors.
    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.
    @c sources_and_weights.second() is the vector of vector weights for each
    of the source mesh entities in @c sources_and_weights.first().  Each element
    of the weights vector is a moment of the source data over the target
    entity; for first order interpolation, only the first element (or zero'th
    moment) of the weights vector (i.e. the volume of intersection) is used.
    Source entities may be repeated in the list if the intersection of a target
    entity and a source entity consists of two or more disjoint pieces
    @param[in] targetCellId The index of the target cell.

  */

  double operator() (int const targetCellID,
                     std::vector<Weights_t> const & sources_and_weights) const
  {
    int nsrccells = sources_and_weights.size();
    if (!nsrccells) {
#ifdef DEBUG
        std::cerr << "WARNING: No source cells contribute to target cell" <<
          targetCellID << std::endl;
#endif
      return 0.0;
    }

    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)
  
    double val = 0.0;
    double wtsum0 = 0.0;
    double vol = target_mesh_.cell_volume(targetCellID);

    int nsummed = 0;
    if (field_type_ == Field_type::MESH_FIELD) {
      for (auto const& wt : sources_and_weights) {
        int srccell = wt.entityID;
        std::vector<double> pair_weights = wt.weights;
        if (pair_weights[0]/vol < 1.0e-12) continue;  // skip small intersections
        val += source_vals_[srccell] * pair_weights[0];
        wtsum0 += pair_weights[0];
        nsummed++;
      }
    } else if (field_type_ == Field_type::MULTIMATERIAL_FIELD) {
      for (auto const& wt : sources_and_weights) {
        int srccell = wt.entityID;
        std::vector<double> pair_weights = wt.weights;
        if (pair_weights[0]/vol < 1.0e-12) continue;  // skip small intersections
        int matcell = source_state_.cell_index_in_material(srccell, matid_);
        val += source_vals_[matcell] * pair_weights[0];  // 1st order
        wtsum0 += pair_weights[0];
        nsummed++;
      }
    }
  
    // Normalize the value by the volume of the intersection of the target cells
    // with the source mesh. This will do the right thing for single-material
    // and multi-material remap (conservative and constant preserving) if there
    // is NO mismatch between source and target mesh boundaries. IF THERE IS A
    // MISMATCH, THIS WILL PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE. THEN
    // WE HAVE TO DO A SEMI-LOCAL OR GLOBAL REPAIR.
    if (nsummed)
      val /= wtsum0;
    else
      val = 0.0;

#ifdef DEBUG
    static bool first = true;
    if (first && fabs((vol-wtsum0)/vol) > 1.0e-10) {
      std::cerr <<
          "WARNING: Meshes may be mismatched in the neighborhood of cell " <<
          targetCellID << " in the target mesh (and maybe other places too)\n";
      first = false;
    }
#endif

    return val;
  } //operator()

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string interp_var_name_;
  double const * source_vals_;
  int matid_;
  Field_type field_type_;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
}; //interpolate_1st_order specialization for cell

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
/*!
  @brief Interpolate_1stOrder specialization for nodes
*/

template<int D,
         typename SourceMeshType, 
         typename TargetMeshType, 
         typename StateType,
         template<class, int> class InterfaceReconstructorType>
class Interpolate_1stOrder<D, 
                           NODE, 
                           SourceMeshType, 
                           TargetMeshType, 
                           StateType,
                           InterfaceReconstructorType> {

#ifdef HAVE_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType>;
#endif

 public:

  //Constructor with interface reconstructor
#ifdef HAVE_TANGRAM
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state,
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interface_reconstructor_(ir),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr) {}
#endif   

  /*!
    @brief Constructor without interface reconstructor
    @param[in] source_mesh The input mesh.
    @param[in] target_mesh The output mesh.
    @param[in] source_state The state manager for data on the input mesh.
    @param[in] interp_var_name The string name of the variable to interpolate.
    @param[in] sources_and_weights Vector of source entities and their weights for each target entity
  */
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       StateType const & source_state) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr) {}


  /// Copy constructor (disabled)
  //  Interpolate_1stOrder(const Interpolate_1stOrder &) = delete;

  /// Assignment operator (disabled)
  //  Interpolate_1stOrder & operator = (const Interpolate_1stOrder &) = delete;

  /// Destructor
  ~Interpolate_1stOrder() {}

  /// Set the material we are operating on

  int set_material(int m) {
    matid_ = m;
  } //set_material

  /// Set the variable name to be interpolated

  // Even though 1st order accurate interpolation does not require
  // limiters and only higher order ones do, we need to have the
  // limiter_type as an argument. This is because the templated driver
  // does not know whether its being called with 1st order or higher
  // order interpolators and so all interpolators need to have a
  // uniform interface

  void set_interpolation_variable(std::string const & interp_var_name,
                                  LimiterType limtype = NOLIMITER) {
    interp_var_name_ = interp_var_name;
    field_type_ = source_state_.field_type(NODE, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(NODE, interp_var_name,
                                  &source_vals_);
    else {
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
      std::cerr << "Cannot remap NODE-centered multi-material data" << "\n";
    }
  } //set_interpolation_variable

  /*!
    @brief Functor to do the actual interpolation.
    @param[in] sources_and_weights A pair of two vectors.
    @c sources_and_weights.first() is the vector of source entity indices
    in the source mesh that will contribute to the current target mesh entity.
    @c sources_and_weights.second() is the vector of vector weights for each
    of the source mesh entities in @c sources_and_weights.first().  Each element
    of the weights vector is a moment of the source data over the target
    entity; for first order interpolation, only the first element (or zero'th
    moment) of the weights vector (i.e. the volume of intersection) is used.
    Source entities may be repeated in the list if the intersection of a target
    entity and a source entity consists of two or more disjoint pieces
    @param[in] targetCellId The index of the target cell.

  */

  double operator() (int const targetNodeID,
                     std::vector<Weights_t> const & sources_and_weights) const
  {
    if (field_type_ != Field_type::MESH_FIELD) return 0.0;

    int nsrcdualcells = sources_and_weights.size();
    if (!nsrcdualcells) {
#ifdef DEBUG
      std::cerr << "WARNING: No source nodes contribute to target node." <<
          std::endl;
#endif
      return 0.0;
    }

    // contribution of the source node (dual cell) is its field value
    // weighted by its "weight" (in this case, the 0th
    // moment/area/volume of its intersection with the target dual cell)
    double vol = target_mesh_.dual_cell_volume(targetNodeID);

    double val = 0.0;
    double wtsum0 = 0.0;
    int nsummed = 0;
    for (auto const& wt : sources_and_weights) {
      int srcnode = wt.entityID;
      std::vector<double> pair_weights = wt.weights;
      if (pair_weights[0]/vol < 1.0e-16) continue;  // skip small intersections
      val += source_vals_[srcnode] * pair_weights[0];  // 1st order
      wtsum0 += pair_weights[0];
      nsummed++;
    }

    // Normalize the value by the volume of the intersection of the target cells
    // with the source mesh. This will do the right thing for single-material
    // and multi-material remap (conservative and constant preserving) if there
    // is NO mismatch between source and target mesh boundaries. IF THERE IS A
    // MISMATCH, THIS WILL PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE. THEN
    // WE HAVE TO DO A SEMI-LOCAL OR GLOBAL REPAIR.

    if (nsummed)
      val /= wtsum0;
    else
      val = 0.0;
  
#ifdef DEBUG
    static bool first = true;
    if (first && fabs((vol-wtsum0)/vol) > 1.0e-10) {
      std::cerr <<
          "WARNING: Meshes may be mismatched in the neighborhood of node " <<
          targetNodeID << " in the target mesh (and maybe other places too)\n";
      first = false;
    }
#endif

    return val;
  } //operator()

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  StateType const & source_state_;
  std::string interp_var_name_;
  double const * source_vals_;
  int matid_;
  Field_type field_type_;

#ifdef HAVE_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
}; //interpolate_1st_order specialization for nodes

//////////////////////////////////////////////////////////////////////////////


}  // namespace Portage

#endif  // PORTAGE_INTERPOLATE_INTERPOLATE_1ST_ORDER_H_
