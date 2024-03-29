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
#include <cmath>

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Vector.h"
#include "wonton/support/CoordinateSystem.h"

// portage includes
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/support/portage.h"
#include "portage/driver/parts.h"

#ifdef PORTAGE_HAS_TANGRAM
#include "tangram/driver/driver.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"
#endif

namespace Portage {

/*!
  @class Interpolate_1stOrder interpolate_1st_order.h
  @brief Interpolate_1stOrder does a 1st order interpolation of scalars

  @tparam D  spatial dimension of problem
  @tparam on_what The type of entity-based data we wish to interpolate;
  e.g. does it live on nodes, cells, edges, etc.
  @tparam SourceMeshType Mesh wrapper class used to access source mesh info
  @tparam TargetMeshType Mesh wrapper class used to access target mesh info
  @tparam SourceStateType State manager used to access source data.
  @tparam InterfaceReconstructorType Class for reconstructing material interfaces
  @tparam MatPoly_Splitter Class used for splitting material polygons
  @tparam MatPoly_Clipper Class used for clipping material polygons
  @tparam CoordSys  What coordinate system are we operating in?
  
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


  Type T must support +, += operators. It must also support *
  operators with a scalar operand. Finally, it must support
  initialization to null values using the syntax T(0.0)

  We could enforce these requirements using SFINAE

*/

template<int D,
         Entity_kind on_what,
         typename SourceMeshType,
         typename TargetMeshType,
         typename SourceStateType,
         typename TargetStateType,
         typename T,
         template<class, int, class, class>
           class InterfaceReconstructorType = DummyInterfaceReconstructor,
         class Matpoly_Splitter = void,
         class Matpoly_Clipper = void,
         class CoordSys = Wonton::DefaultCoordSys
         >
class Interpolate_1stOrder {

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

  // Constructor with interface reconstructor
#ifdef PORTAGE_HAS_TANGRAM
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols,
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols),
      interface_reconstructor_(ir) {}
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
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols) {}


  /// Copy constructor (disabled)
  //  Interpolate_1stOrder(const Interpolate_1stOrder &) = delete;

  /// Assignment operator (disabled)
  //  Interpolate_1stOrder & operator = (const Interpolate_1stOrder &) = delete;

  /// Destructor
  ~Interpolate_1stOrder() = default;

  /// Set the material we are operating on

  void set_material(int m) {
    matid_ = m;
  }  // set_material

  /// Set the variable name to be interpolated.

  // Even though 1st order accurate interpolation does not require
  // limiters and only higher order ones do, we need to have the
  // limiter_type as an argument. This is because the templated driver
  // does not know whether its being called with 1st order or higher
  // order interpolators and so all interpolators need to have a
  // uniform interface

  void set_interpolation_variable(std::string const & interp_var_name,
                                  Wonton::vector<Wonton::Vector<D>>* gradients = nullptr) {
    interp_var_name_ = interp_var_name;
    field_type_ = source_state_.field_type(Entity_kind::CELL, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(Entity_kind::CELL, interp_var_name,
                                  &source_vals_);
    else
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
  }  // set_interpolation_variable

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

  T operator() (int const targetEntityId,
                     std::vector<Weights_t> const & sources_and_weights) const {
    throw std::runtime_error("Interpolation operator not implemented for this entity type");
  }
  
  constexpr static int order = 1;

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  SourceStateType const & source_state_;
  std::string interp_var_name_;
  T const * source_vals_;
  int matid_ = 0;
  Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;
  NumericTolerances_t num_tols_;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
};  // interpolate_1st_order base



//////////////////////////////////////////////////////////////////////////////
/*!
  @brief Interpolate_1stOrder specialization for cells
*/

template<int D,
         typename SourceMeshType,
         typename TargetMeshType,
         typename SourceStateType,
         typename TargetStateType,
         typename T,
         template<class, int, class, class> class InterfaceReconstructorType,
         class Matpoly_Splitter,
         class Matpoly_Clipper,
         class CoordSys>
class Interpolate_1stOrder<
  D, Entity_kind::CELL,
  SourceMeshType, TargetMeshType,
  SourceStateType, TargetStateType,
  T,
  InterfaceReconstructorType,
  Matpoly_Splitter, Matpoly_Clipper, CoordSys> {

  // useful aliases
  using Parts = PartPair<D,
    SourceMeshType, SourceStateType,
    TargetMeshType, TargetStateType
  >;

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef PORTAGE_HAS_TANGRAM
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols,
                       std::shared_ptr<InterfaceReconstructor> ir,
                       const Parts* const parts = nullptr) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols),
      interface_reconstructor_(ir),
      parts_(parts) {}
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


  /// Copy constructor (disabled)
  //  Interpolate_1stOrder(const Interpolate_1stOrder &) = delete;

  /// Assignment operator (disabled)
  //  Interpolate_1stOrder & operator = (const Interpolate_1stOrder &) = delete;

  /// Destructor
  ~Interpolate_1stOrder() = default;


  /// Set the material we are operating on

  void set_material(int m) {
    matid_ = m;
  }  // set_material

  /// Set the variable name to be interpolated


  // Even though 1st order accurate interpolation does not require
  // limiters and only higher order ones do, we need to have the
  // limiter_type as an argument. This is because the templated driver
  // does not know whether its being called with 1st order or higher
  // order interpolators and so all interpolators need to have a
  // uniform interface

  void set_interpolation_variable(std::string const & interp_var_name,
                                  Wonton::vector<Wonton::Vector<D>>* gradients = nullptr) {
    interp_var_name_ = interp_var_name;
    field_type_ = source_state_.field_type(Entity_kind::CELL, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(Entity_kind::CELL, interp_var_name,
                                  &source_vals_);
    else
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
  }  // set_interpolation_variable


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

  T operator() (int const targetCellID,
                     std::vector<Weights_t> const & sources_and_weights) const
  {
    int nsrccells = sources_and_weights.size();
    if (!nsrccells) return T(0.0);

    // contribution of the source cell is its field value weighted by
    // its "weight" (in this case, its 0th moment/area/volume)

    T val(0.0);
    double wtsum0 = 0.0;

    int nsummed = 0;
    if (field_type_ == Field_type::MESH_FIELD) {
      for (auto const& wt : sources_and_weights) {
        int srccell = wt.entityID;
        std::vector<double> pair_weights = wt.weights;
        if (fabs(pair_weights[0]) < num_tols_.min_absolute_volume)
          continue;  // skip small intersections
        val += source_vals_[srccell] * pair_weights[0];
        wtsum0 += pair_weights[0];
        nsummed++;
      }
    } else if (field_type_ == Field_type::MULTIMATERIAL_FIELD) {
      for (auto const& wt : sources_and_weights) {
        int srccell = wt.entityID;
        std::vector<double> pair_weights = wt.weights;
        if (fabs(pair_weights[0]) < num_tols_.min_absolute_volume)
          continue;  // skip small intersections
        int matcell = source_state_.cell_index_in_material(srccell, matid_);
        val += source_vals_[matcell] * pair_weights[0];  // 1st order
        wtsum0 += pair_weights[0];
        nsummed++;
      }
    }

    // Normalize the value by the volume of the intersection of the
    // target cells with the source mesh. This will do the right thing
    // for single-material and multi-material remap (conservative and
    // constant preserving) if there is NO mismatch between source and
    // target mesh boundaries. IF THERE IS A MISMATCH, THIS WILL
    // PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE. THEN WE HAVE
    // TO DO A SEMI-LOCAL OR GLOBAL REPAIR.

    // We use the * operator instead of / so as to reduce the number
    // of requirements on generic variable type
    
    if (nsummed)
      val *= (1.0/wtsum0);

    return val;
  }  // operator()
  
  constexpr static int order = 1;

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  SourceStateType const & source_state_;
  std::string interp_var_name_;
  T const * source_vals_;
  int matid_ = 0;
  Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;
  NumericTolerances_t num_tols_;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
  Parts const* parts_;
};  // interpolate_1st_order specialization for cell

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////
/*!
  @brief Interpolate_1stOrder specialization for nodes
*/

template<int D,
         typename SourceMeshType,
         typename TargetMeshType,
         typename SourceStateType,
         typename TargetStateType,
         typename T,
         template<class, int, class, class> class InterfaceReconstructorType,
         class Matpoly_Splitter,
         class Matpoly_Clipper,
         class CoordSys>
class Interpolate_1stOrder<
  D, Entity_kind::NODE,
  SourceMeshType, TargetMeshType,
  SourceStateType, TargetStateType,
  T,
  InterfaceReconstructorType,
  Matpoly_Splitter, Matpoly_Clipper, CoordSys> {

  // useful aliases
#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor =
      Tangram::Driver<InterfaceReconstructorType, D, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

  // Constructor with interface reconstructor
#ifdef PORTAGE_HAS_TANGRAM
  Interpolate_1stOrder(SourceMeshType const & source_mesh,
                       TargetMeshType const & target_mesh,
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols,
                       std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols),
      interface_reconstructor_(ir) {}
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
                       SourceStateType const & source_state,
                       NumericTolerances_t num_tols) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interp_var_name_("VariableNameNotSet"),
      source_vals_(nullptr),
      num_tols_(num_tols) {}


  /// Copy constructor (disabled)
  //  Interpolate_1stOrder(const Interpolate_1stOrder &) = delete;

  /// Assignment operator (disabled)
  //  Interpolate_1stOrder & operator = (const Interpolate_1stOrder &) = delete;

  /// Destructor
  ~Interpolate_1stOrder() = default;

  /// Set the material we are operating on

  void set_material(int m) {
    matid_ = m;
  }  // set_material

  /// Set the variable name to be interpolated

  // Even though 1st order accurate interpolation does not require
  // limiters and only higher order ones do, we need to have the
  // limiter_type as an argument. This is because the templated driver
  // does not know whether its being called with 1st order or higher
  // order interpolators and so all interpolators need to have a
  // uniform interface

  void set_interpolation_variable(std::string const & interp_var_name,
                                  Wonton::vector<Wonton::Vector<D>>* gradients = nullptr) {
    interp_var_name_ = interp_var_name;
    field_type_ = source_state_.field_type(Entity_kind::NODE, interp_var_name);
    if (field_type_ == Field_type::MESH_FIELD)
      source_state_.mesh_get_data(Entity_kind::NODE, interp_var_name,
                                  &source_vals_);
    else {
      source_state_.mat_get_celldata(interp_var_name, matid_, &source_vals_);
      throw std::runtime_error("Cannot remap NODE-centered multi-material data");
    }
  }  // set_interpolation_variable


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

  T operator() (int const targetNodeID,
                std::vector<Weights_t> const & sources_and_weights) const
  {
    if (field_type_ != Field_type::MESH_FIELD) return T(0.0);

    int nsrcdualcells = sources_and_weights.size();
    if (!nsrcdualcells) return T(0.0);

    // contribution of the source node (dual cell) is its field value
    // weighted by its "weight" (in this case, the 0th
    // moment/area/volume of its intersection with the target dual cell)

    T val(0.0);
    double wtsum0 = 0.0;
    int nsummed = 0;
    for (auto const& wt : sources_and_weights) {
      int srcnode = wt.entityID;
      std::vector<double> pair_weights = wt.weights;
      if (fabs(pair_weights[0]) < num_tols_.min_absolute_volume)
        continue;  // skip small intersections
      val += source_vals_[srcnode] * pair_weights[0];  // 1st order
      wtsum0 += pair_weights[0];
      nsummed++;
    }

    // Normalize the value by the volume of the intersection of the
    // target cells with the source mesh. This will do the right thing
    // for single-material and multi-material remap (conservative and
    // constant preserving) if there is NO mismatch between source and
    // target mesh boundaries. IF THERE IS A MISMATCH, THIS WILL
    // PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE. THEN WE HAVE
    // TO DO A SEMI-LOCAL OR GLOBAL REPAIR.

    // We use the * operator instead of / so as to reduce the number
    // of requirements on generic variable type
    
    if (nsummed)
      val *= (1.0/wtsum0);

    return val;
  }  // operator()

  constexpr static int order = 1;

 private:
  SourceMeshType const & source_mesh_;
  TargetMeshType const & target_mesh_;
  SourceStateType const & source_state_;
  std::string interp_var_name_;
  T const * source_vals_;
  int matid_ = 0;
  Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;
  NumericTolerances_t num_tols_;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
};  // interpolate_1st_order specialization for nodes

//////////////////////////////////////////////////////////////////////////////


}  // namespace Portage

#endif  // PORTAGE_INTERPOLATE_INTERPOLATE_1ST_ORDER_H_
