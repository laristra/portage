/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_
#define PORTAGE_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_

#include <cassert>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <iostream>
#include <utility>
#include <vector>

#include "portage/support/portage.h"
#include "portage/interpolate/gradient.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/driver/fix_mismatch.h"
#include "portage/driver/parts.h"

#include "wonton/support/CoordinateSystem.h"

#ifdef HAVE_TANGRAM
  #include "tangram/driver/driver.h"
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/support/MatPoly.h"
#endif

namespace Portage {

  /**
   * @class Interpolate_2ndOrder interpolate_2nd_order.h
   * @brief Interpolate_2ndOrder does a 2nd order interpolation of scalars.
   *
   * @tparam D  spatial dimension of problem
   * @tparam on_what: The type of entity-based data we wish to interpolate;
   *                 e.g. does it live on nodes, cells, edges, etc.
   * @tparam SourceMeshType: Mesh wrapper class used to access source mesh info
   * @tparam TargetMeshType: Mesh wrapper class used to access target mesh info
   * @tparam SourceStateType: State manager used to access source data.
   * @tparam InterfaceReconstructorType: Class for reconstructing material interfaces
   * @tparam MatPoly_Splitter: Class used for splitting material polygons
   * @tparam MatPoly_Clipper: Class used for clipping material polygons
   * @tparam CoordSys:  What coordinate system are we operating in?
   *
   * [1] Margolin, L.G. and Shashkov, M.J. "Second-order sign-preserving
   * conservative interpolation (remapping) on general grids." Journal of
   * Computational Physics, v 184, n 1, pp. 266-298, 2003.
   *
   * [2] Dukowicz, J.K. and Kodis, J.W. "Accurate Conservative Remapping
   * (Rezoning) for Arbitrary Lagrangian-Eulerian Computations," SIAM
   * Journal on Scientific and Statistical Computing, Vol. 8, No. 3,
   * pp. 305-321, 1987.
   */
  template<
    int D, Entity_kind on_what,
    typename SourceMeshType,
    typename TargetMeshType,
    typename SourceStateType,
    typename TargetStateType,
    typename T,
    template<class, int, class, class>
      class InterfaceReconstructorType = DummyInterfaceReconstructor,
    class Matpoly_Splitter = void, class Matpoly_Clipper = void,
    class CoordSys = Wonton::DefaultCoordSys
  >
  class Interpolate_2ndOrder {

#ifdef HAVE_TANGRAM
    using InterfaceReconstructor =
      Tangram::Driver<
        InterfaceReconstructorType, D, SourceMeshType,
        Matpoly_Splitter, Matpoly_Clipper
      >;
#endif

  public:

    /**
     * @brief Constructor without interface reconstructor.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     */
    Interpolate_2ndOrder(SourceMeshType const& source_mesh,
                         TargetMeshType const& target_mesh,
                         SourceStateType const& source_state,
                         NumericTolerances_t num_tols)
    : source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      variable_name_("VariableNameNotSet"),
      source_values_(nullptr),
      num_tols_(num_tols) { CoordSys::template verify_coordinate_system<D>(); }

#ifdef HAVE_TANGRAM
    /**
     * @brief Constructor with interface reconstructor.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] num_tols: numerical tolerances.
     * @param[in] ir: interface reconstructor for querying matpolys on source mesh.
     */
    Interpolate_2ndOrder(SourceMeshType const& source_mesh,
                         TargetMeshType const& target_mesh,
                         SourceStateType const& source_state,
                         NumericTolerances_t num_tols,
                         std::shared_ptr<InterfaceReconstructor> ir)
      : source_mesh_(source_mesh),
        target_mesh_(target_mesh),
        source_state_(source_state),
        interface_reconstructor_(ir),
        variable_name_("VariableNameNotSet"),
        source_values_(nullptr),
        num_tols_(num_tols) { CoordSys::template verify_coordinate_system<D>(); }
#endif

    /**
     * @brief Assignment operator (disabled).
     *
     * @param[in] other the interpolator to copy
     * @return current interpolator reference
     */
    Interpolate_2ndOrder& operator = (const Interpolate_2ndOrder& other) = delete;


    /**
     * @brief Destructor.
     *
     */
    ~Interpolate_2ndOrder() = default;

    /**
     * @brief Set the material we are operating on.
     *
     * @param[in] m: the material ID.
     */
    void set_material(int m) { material_id_ = m; }


    /**
     * @brief Set the name of the interpolation variable and the gradient field.
     *
     * @param[in] variable_name: the variable name
     * @param[in] gradient_field: the gradient field
     */
    void set_interpolation_variable(std::string const& variable_name,
                                    const Portage::vector<Vector<D>>* gradient_field = nullptr) {
      std::cerr << "Interpolation is available for only cells and nodes";
      std::cerr << std::endl;
    }

    /**
     * @brief Functor to do the actual interpolate calculation.
     *
     * @param targetCellID: the target cell index.
     * @param sources_and_weights: a list of source cells and intersection weights.
     * Each element of the weights vector is a moment of the source data over
     * the target entity. For first-order interpolation, only the first element
     * (or zero'th moment) of the weights vector (the volume of intersection)
     * is used. Source entities may be repeated in the list if the intersection
     * of a target cell and a source cell consists of two or more disjoint pieces.
     *
     * @return the interpolated value.
     * @todo cleanup the datatype for sources_and_weights - it is somewhat confusing.
     * @todo must remove assumption that field is scalar.
     */
    double operator()(int const targetCellID,
                      std::vector<Weights_t> const& sources_and_weights) const {

      // not implemented for all types - see specialization for cells and nodes
      std::cerr << "Error: interpolation operator not implemented for this entity type";
      std::cerr << std::endl;
      return 0.;
    }

    constexpr static int order = 2;

  private:
    SourceMeshType const& source_mesh_;
    TargetMeshType const& target_mesh_;
    SourceStateType const& source_state_;
    std::string variable_name_ = "";
    T const* source_values_;
    NumericTolerances_t num_tols_;
    int material_id_ = 0;
    Portage::vector<Wonton::Vector<D>> const* gradients_;
    Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;
#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
  };

  /* ------------------------------------------------------------------------ */

  /**
   * @brief second-order interpolate class specialization for scalar
   * fields on cells
   *
   * @tparam D: spatial dimension of problem
   * @tparam SourceMeshType: mesh wrapper class used to access source mesh info
   * @tparam TargetMeshType: mesh wrapper class used to access target mesh info
   * @tparam SourceStateType: state manager used to access source data.
   * @tparam InterfaceReconstructorType: class for material interfaces reconstruction
   * @tparam MatPoly_Splitter: class used for splitting material polygons
   * @tparam MatPoly_Clipper: class used for clipping material polygons
   * @tparam CoordSys: what coordinate system are we operating in?
   */
  template<
    int D,
    typename SourceMeshType,
    typename TargetMeshType,
    typename SourceStateType,
    typename TargetStateType,
    template<class, int, class, class>
      class InterfaceReconstructorType,
    class Matpoly_Splitter, class Matpoly_Clipper, class CoordSys
  >
  class Interpolate_2ndOrder<
    D, Entity_kind::CELL,
    SourceMeshType, TargetMeshType,
    SourceStateType, TargetStateType,
    double,
    InterfaceReconstructorType,
    Matpoly_Splitter, Matpoly_Clipper, CoordSys> {

    // useful aliases
    using Parts = PartPair<
      D, SourceMeshType, SourceStateType,
      TargetMeshType, TargetStateType
    >;

    using Gradient = Limited_Gradient<
      D, Entity_kind::CELL,
      SourceMeshType, SourceStateType,
      InterfaceReconstructorType,
      Matpoly_Splitter, Matpoly_Clipper, CoordSys
    >;

#ifdef HAVE_TANGRAM
    using InterfaceReconstructor = Tangram::Driver<
      InterfaceReconstructorType, D, SourceMeshType,
      Matpoly_Splitter, Matpoly_Clipper
    >;
#endif

  public:
    /**
     * @brief Constructor without interface reconstructor.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] num_tols: numerical tolerances.
     */
    Interpolate_2ndOrder(SourceMeshType const& source_mesh,
                         TargetMeshType const& target_mesh,
                         SourceStateType const& source_state,
                         NumericTolerances_t num_tols,
                         const Parts* const parts = nullptr)
      : source_mesh_(source_mesh),
        target_mesh_(target_mesh),
        source_state_(source_state),
        variable_name_("VariableNameNotSet"),
        source_values_(nullptr),
        num_tols_(num_tols),
        parts_(parts) { CoordSys::template verify_coordinate_system<D>(); }

#ifdef HAVE_TANGRAM
    /**
     * @brief Constructor with interface reconstructor.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] num_tols: numerical tolerances.
     * @param[in] ir: interface reconstructor for querying matpolys on source mesh.
     */
    Interpolate_2ndOrder(SourceMeshType const& source_mesh,
                         TargetMeshType const& target_mesh,
                         SourceStateType const& source_state,
                         NumericTolerances_t num_tols,
                         std::shared_ptr<InterfaceReconstructor> ir,
                         const Parts* const parts = nullptr)
      : source_mesh_(source_mesh),
        target_mesh_(target_mesh),
        source_state_(source_state),
        interface_reconstructor_(ir),
        variable_name_("VariableNameNotSet"),
        source_values_(nullptr),
        num_tols_(num_tols),
        parts_(parts) { CoordSys::template verify_coordinate_system<D>(); }
#endif

    /**
     * @brief Assignment operator (disabled).
     *
     * @param[in] other the interpolator to copy
     * @return current interpolator reference
     */
    Interpolate_2ndOrder &operator=(const Interpolate_2ndOrder &) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~Interpolate_2ndOrder() = default;

    /**
     * @brief Set the material we are operating on.
     *
     * @param[in] m: the material ID.
     */
    void set_material(int m) { material_id_ = m; }


    /**
     * @brief Set the name of the interpolation variable and the gradient field.
     *
     * @param[in] variable_name: the variable name
     * @param[in] gradient_field: the gradient field
     */
    void set_interpolation_variable(std::string const& variable_name,
                                    const Portage::vector<Vector<D>>* gradient_field = nullptr) {

      variable_name_ = variable_name;
      gradients_ = gradient_field;
      field_type_ = source_state_.field_type(Entity_kind::CELL, variable_name);

      if (field_type_ == Field_type::MESH_FIELD) {
        source_state_.mesh_get_data(Entity_kind::CELL, variable_name, &source_values_);
      } else {
        source_state_.mat_get_celldata(variable_name, material_id_, &source_values_);
      }
    }


    /**
     * @brief Functor to compute the interpolation of cell values.
     *
     * @param[in] cell_id: target cell index
     * @param[in] sources_and_weights: list of source mesh entities and
     * corresponding weight vectors. Each element of the weights vector
     * is a moment of the source data over the target entity; for first
     * order interpolation, only the first element (or zero'th moment)
     * of the weights vector (volume of intersection) is used.
     * Source entities may be repeated in the list if the intersection of a
     * target cell and a source cell consists of two or more disjoint pieces.
     *
     * @return the interpolated value.
     * @todo must remove assumption that field is scalar.
     */
    double operator()(int cell_id,
                      std::vector<Weights_t> const& sources_and_weights) const {

      if (sources_and_weights.empty())
        return 0.;

      auto const& gradient_field = *gradients_;
      double total_value = 0.;
      double normalization = 0.;

      /*
       * contribution of the source cell is its field value weighted by
       * its "weight" (in this case, its 0th moment/area/volume).
       */
      int nb_summed = 0;

      // Loop over source cells
      for (auto&& current : sources_and_weights) {
        // Get source cell and the intersection weights
        int src_cell = current.entityID;
        auto intersect_weights = current.weights;
        double intersect_volume = intersect_weights[0];

        if (fabs(intersect_volume) <= num_tols_.min_absolute_volume)
          continue;  // no intersection

        // Obtain source cell centroid
        Point<D> source_centroid;
        if (field_type_ == Field_type::MESH_FIELD) {
          source_mesh_.cell_centroid(src_cell, &source_centroid);
        }
#ifdef HAVE_TANGRAM
        else if (field_type_ == Field_type::MULTIMATERIAL_FIELD) {
          int const nb_mats = source_state_.cell_get_num_mats(src_cell);
          std::vector<int> cellmats;
          source_state_.cell_get_mats(src_cell, &cellmats);

          bool is_pure_cell =
            (nb_mats == 0 or (nb_mats == 1 and cellmats[0] == material_id_));

          if (is_pure_cell) {
            source_mesh_.cell_centroid(src_cell, &source_centroid);
          } else /* multi-material cell */ {
            assert(interface_reconstructor_ != nullptr);  // must be defined

            auto pos = std::find(cellmats.begin(), cellmats.end(), material_id_);
            bool found_material = (pos != cellmats.end());

            if (found_material) /* mixed cell contains this material */ {

              // obtain matpoly's for this material
              auto const& cellmatpoly = interface_reconstructor_->cell_matpoly_data(src_cell);
              auto matpolys = cellmatpoly.get_matpolys(material_id_);

              int cnt = 0;
              for (int k = 0; k < D; k++)
                source_centroid[k] = 0;

              /*
               * compute centroid of all matpoly's by summing all the
               * first order moments first, and then dividing by the
               * total volume of all matpolys.
               */
              double mvol = 0.;
              for (auto&& poly : matpolys) {
                auto moments = poly.moments();
                mvol += moments[0];
                for (int k = 0; k < D; k++)
                  source_centroid[k] += moments[k + 1];
              }

              for (int k = 0; k < D; k++)
                source_centroid[k] /= mvol;
            }
          }
        }
#endif

        // compute intersection centroid
        Point<D> intersect_centroid;
        // first-moment / volume
        for (int k = 0; k < D; ++k)
          intersect_centroid[k] = intersect_weights[1 + k] / intersect_volume;

        // retrieve the correct source cell index
        int const source_index = (field_type_ == Field_type::MULTIMATERIAL_FIELD
          ? source_state_.cell_index_in_material(src_cell, material_id_)
          : src_cell);

        Vector<D> gradient = gradient_field[source_index];
        Vector<D> dr = intersect_centroid - source_centroid;
        CoordSys::modify_line_element(dr, source_centroid);

        double value = source_values_[source_index] + dot(gradient, dr);
        value *= intersect_volume;
        total_value += value;
        normalization += intersect_volume;
        nb_summed++;
      }

      /*
       * Normalize the value by the volume of the intersection of the target cells
       * with the source mesh. This will do the right thing for single-material
       * and multi-material remap (conservative and constant preserving) if there
       * is NO mismatch between source and target mesh boundaries. IF THERE IS A
       * MISMATCH, THIS WILL PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE.
       * THEN WE HAVE TO DO A SEMI-LOCAL OR GLOBAL REPAIR.
       */
      return nb_summed ? total_value / normalization : 0.;
    }
    
    constexpr static int order = 2;

  private:
    SourceMeshType const& source_mesh_;
    TargetMeshType const& target_mesh_;
    SourceStateType const& source_state_;
    std::string variable_name_;
    double const* source_values_;
    NumericTolerances_t num_tols_;
    int material_id_ = 0;
    Portage::vector<Wonton::Vector<D>> const* gradients_;
    Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;
#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
    Parts const* parts_;
  };

  /* ------------------------------------------------------------------------ */

  /**
   * @brief second-order interpolate class specialization for scalar fields on nodes.
   *
   * @tparam D: spatial dimension of problem
   * @tparam SourceMeshType: mesh wrapper class used to access source mesh info
   * @tparam TargetMeshType: mesh wrapper class used to access target mesh info
   * @tparam SourceStateType: state manager used to access source data.
   * @tparam InterfaceReconstructorType: class for material interfaces reconstruction
   * @tparam MatPoly_Splitter: class used for splitting material polygons
   * @tparam MatPoly_Clipper: class used for clipping material polygons
   * @tparam CoordSys: what coordinate system are we operating in?
   */
  template<
    int D,
    typename SourceMeshType,
    typename TargetMeshType,
    typename SourceStateType,
    typename TargetStateType,
    template<class, int, class, class>
      class InterfaceReconstructorType,
    class Matpoly_Splitter, class Matpoly_Clipper, class CoordSys
  >
  class Interpolate_2ndOrder<
    D, Entity_kind::NODE,
    SourceMeshType, TargetMeshType,
    SourceStateType, TargetStateType,
    double,
    InterfaceReconstructorType,
    Matpoly_Splitter, Matpoly_Clipper, CoordSys> {

    // useful aliases
    using Gradient = Limited_Gradient<
      D, Entity_kind::NODE,
      SourceMeshType, SourceStateType,
      InterfaceReconstructorType,
      Matpoly_Splitter, Matpoly_Clipper, CoordSys
    >;

#ifdef HAVE_TANGRAM
    using InterfaceReconstructor = Tangram::Driver<
      InterfaceReconstructorType, D, SourceMeshType,
      Matpoly_Splitter, Matpoly_Clipper
    >;
#endif

    // get rid of long namespaces
    static auto const Node = Entity_kind::NODE;

  public:

    /**
     * @brief Constructor without interface reconstructor.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] num_tols: numerical tolerances.
     */
    Interpolate_2ndOrder(SourceMeshType const& source_mesh,
                         TargetMeshType const& target_mesh,
                         SourceStateType const& source_state,
                         NumericTolerances_t num_tols) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      variable_name_("VariableNameNotSet"),
      source_values_(nullptr),
      num_tols_(num_tols) {}

#ifdef HAVE_TANGRAM
    /**
     * @brief Constructor with interface reconstructor.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] num_tols: numerical tolerances.
     * @param[in] ir: interface reconstructor for querying matpolys on source mesh.
     */
    Interpolate_2ndOrder(SourceMeshType const& source_mesh,
                         TargetMeshType const& target_mesh,
                         SourceStateType const& source_state,
                         NumericTolerances_t num_tols,
                         std::shared_ptr<InterfaceReconstructor> ir) :
      source_mesh_(source_mesh),
      target_mesh_(target_mesh),
      source_state_(source_state),
      interface_reconstructor_(ir),
      variable_name_("VariableNameNotSet"),
      source_values_(nullptr),
      num_tols_(num_tols) {}
#endif

    /**
     * @brief Assignment operator (disabled).
     *
     * @param[in] other the interpolator to copy
     * @return current interpolator reference
     */
    Interpolate_2ndOrder& operator = (const Interpolate_2ndOrder&) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~Interpolate_2ndOrder() = default;

    /**
     * @brief Set the material we are operating on.
     *
     * @param[in] m: the material ID.
     */
    void set_material(int m) { material_id_ = m; }

    /**
     * @brief Set the name of the interpolation variable and the gradient field.
     *
     * @param[in] variable_name: the variable name
     * @param[in] gradient_field: the gradient field to use
     */
    void set_interpolation_variable(std::string const variable_name,
                                    Portage::vector<Vector<D>>* gradient_field = nullptr) {

      variable_name_ = variable_name;
      gradients_ = gradient_field;

      // Extract the field data from the statemanager
      field_type_ = source_state_.field_type(Entity_kind::NODE, variable_name);

      if (field_type_ == Field_type::MESH_FIELD) {
        source_state_.mesh_get_data(Entity_kind::NODE, variable_name, &source_values_);
      } else {
        std::cerr << "Sorry: cannot remap node-centered multi-material data.";
        std::cerr << std::endl;
        return;
      }
    }


    /**
     * @brief Functor to compute the interpolation of node values.
     *
     * @param[in] node_id: target node index
     * @param[in] sources_and_weights: list of source mesh entities and
     * corresponding weight vectors. Each element of the weights vector
     * is a moment of the source data over the target entity; for first
     * order interpolation, only the first element (or zero'th moment)
     * of the weights vector (volume of intersection) is used.
     * Source entities may be repeated in the list if the intersection of a
     * target cell and a source cell consists of two or more disjoint pieces.
     *
     * @return the interpolated value.
     * @todo: must remove assumption that field is scalar.
     */
    double operator()(int node_id,
                      std::vector<Weights_t> const& sources_and_weights) const {

      if (sources_and_weights.empty())
        return 0.;

      auto const& gradient_field = *gradients_;
      int const nb_source_nodes = sources_and_weights.size();
      double total_value = 0.;
      double normalization = 0.;

      // contribution of the source cell is its field value weighted by
      // its "weight" (in this case, its 0th moment/area/volume)
      int nb_summed = 0;

      for (auto&& current : sources_and_weights) {
        int src_node = current.entityID;
        auto intersect_weights = current.weights;
        double intersect_volume = intersect_weights[0];

        if (fabs(intersect_volume) <= num_tols_.min_absolute_volume)
          continue;  // no intersection

        // note: here we are getting the node coord, not the centroid of
        // the dual cell
        Point<D> source_coord;
        source_mesh_.node_get_coordinates(src_node, &source_coord);

        Point<D> intersect_centroid;
        // first-moment / volume
        for (int k = 0; k < D; ++k)
          intersect_centroid[k] = intersect_weights[1 + k] / intersect_volume;

        Vector<D> gradient = gradient_field[src_node];
        Vector<D> dr = intersect_centroid - source_coord;
        CoordSys::modify_line_element(dr, source_coord);

        double value = source_values_[src_node] + dot(gradient, dr);
        value *= intersect_volume;
        total_value += value;
        normalization += intersect_volume;
        nb_summed++;
      }

      /*
       * Normalize the value by volume of the target dual cell.
       *
       * Normalize the value by the volume of the intersection of the target cells
       * with the source mesh. This will do the right thing for single-material
       * and multi-material remap (conservative and constant preserving) if there
       * is NO mismatch between source and target mesh boundaries. IF THERE IS A
       * MISMATCH, THIS WILL PRESERVE CONSTANT VALUES BUT NOT BE CONSERVATIVE.
       * THEN WE HAVE TO DO A SEMI-LOCAL OR GLOBAL REPAIR.
       */
      return nb_summed ? total_value / normalization : 0.;
    }

    constexpr static int order = 2;

  private:
    SourceMeshType const& source_mesh_;
    TargetMeshType const& target_mesh_;
    SourceStateType const& source_state_;
    std::string variable_name_;
    double const* source_values_;
    NumericTolerances_t num_tols_;
    int material_id_ = 0;
    Portage::vector<Vector<D>>* gradients_ = nullptr;
    Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;
#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
  };
  /* ------------------------------------------------------------------------ */
}  // namespace Portage

#endif  // PORTAGE_INTERPOLATE_INTERPOLATE_2ND_ORDER_H_
