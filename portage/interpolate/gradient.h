/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERPOLATE_GRADIENT_H_
#define PORTAGE_INTERPOLATE_GRADIENT_H_

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

#include "portage/support/portage.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/driver/fix_mismatch.h"
#include "portage/driver/parts.h"

#include "wonton/support/lsfits.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

#ifdef HAVE_TANGRAM
  #include "tangram/driver/driver.h"
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/support/MatPoly.h"
#endif

namespace Portage {

  using Wonton::Point;
  using Wonton::Vector;

/*! @class Limited_Gradient gradient.h
    @brief Compute limited gradient of a field or components of a field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
    @tparam on_what An enum type which indicates different entity types
*/

  template<
    int D, Entity_kind on_what, typename MeshType, typename StateType,
    template<class, int, class, class>
      class InterfaceReconstructorType = DummyInterfaceReconstructor,
    class Matpoly_Splitter = void,
    class Matpoly_Clipper = void,
    class CoordSys = Wonton::DefaultCoordSys
  >
  class Limited_Gradient {

    // useful aliases
    using Parts = PartPair<D, on_what, MeshType, StateType>;

#ifdef HAVE_TANGRAM
    using InterfaceReconstructor =
      Tangram::Driver<
        InterfaceReconstructorType, D, MeshType,
        Matpoly_Splitter, Matpoly_Clipper
      >;
#endif

  public:
    /*! @brief Constructor
        @param[in] mesh  Mesh class than one can query for mesh info
        @param[in] state A state manager class that one can query for field info
        @param[in] on_what An enum that indicates what type of entity the field is
        on
        @param[in] var_name Name of field for which the gradient is to be computed
        @param[in] limiter_type An enum indicating if the limiter type (none,
        Barth-Jespersen, Superbee etc)
        @param[in] Boundary_Limiter_type An enum indicating the limiter type on the boundary

        @todo must remove assumption that field is scalar
     */
    Limited_Gradient(MeshType const& mesh, StateType const& state,
                     std::string const var_name,
                     Limiter_type limiter_type,
                     Boundary_Limiter_type boundary_limiter_type,
                     const Parts* const parts = nullptr)
      : mesh_(mesh),
        state_(state),
        values_(nullptr),
        variable_name_(var_name),
        limiter_type_(limiter_type),
        boundary_limiter_type_(boundary_limiter_type),
        parts_(parts) {}

#ifdef HAVE_TANGRAM
    Limited_Gradient(MeshType const &mesh, StateType const &state,
                     std::string const var_name,
                     Limiter_type limiter_type,
                     Boundary_Limiter_type boundary_limiter_type,
                     std::shared_ptr<InterfaceReconstructor> ir,
                     const Parts* const parts = nullptr)
      : mesh_(mesh),
        state_(state),
        values_(nullptr),
        variable_name_(var_name),
        limiter_type_(limiter_type),
        boundary_limiter_type_(boundary_limiter_type),
        parts_(parts) {}
#endif

    // Assignment operator (disabled)
    Limited_Gradient& operator = (const Limited_Gradient&) = delete;

    // Destructor
    ~Limited_Gradient() = default;

    // Functor - not implemented for all types - see specialization for
    // cells, nodes
    Vector<D> operator()(int entity_id) {
      std::cerr << "Limited gradient not implemented for this entity kind";
      std::cerr << std::endl;
    }

  private:
    MeshType const& mesh_;
    StateType const& state_;
    double const* values_;
    std::string variable_name_ = "";
    Limiter_type limiter_type_ = DEFAULT_LIMITER;
    Boundary_Limiter_type boundary_limiter_type_ = DEFAULT_BND_LIMITER;
    Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;
    int material_id_ = 0;
    Parts const* parts_;
  };

  //////////////////////////////////////////////////////////////////////////////

  /*! @class Limited_Gradient<MeshType,StateType,CELL> gradient.h
    @brief Specialization of limited gradient class for @c cell-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
  */


  template<
    int D, typename MeshType, typename StateType,
    template<class, int, class, class>
      class InterfaceReconstructorType,
    class Matpoly_Splitter, class Matpoly_Clipper, class CoordSys
  >
  class Limited_Gradient<
    D, Entity_kind::CELL, MeshType, StateType,
    InterfaceReconstructorType,
    Matpoly_Splitter, Matpoly_Clipper, CoordSys
  > {

    // useful aliases
    using Parts = PartPair<D, Entity_kind::CELL, MeshType, StateType>;

#ifdef HAVE_TANGRAM
    using InterfaceReconstructor =
      Tangram::Driver<
        InterfaceReconstructorType, D, MeshType,
        Matpoly_Splitter, Matpoly_Clipper
      >;
#endif

  public:
    //Constructor for single material remap
    Limited_Gradient(MeshType const& mesh,
                     StateType const& state,
                     std::string const var_name,
                     Limiter_type limiter_type,
                     Boundary_Limiter_type boundary_limiter_type,
                     const Parts* const parts = nullptr)
      : mesh_(mesh),
        state_(state),
        values_(nullptr),
        variable_name_(var_name),
        limiter_type_(limiter_type),
        boundary_limiter_type_(boundary_limiter_type),
        parts_(parts) {

      // Collect and keep the list of neighbors for each CELL as it may
      // be expensive to go to the mesh layer and collect this data for
      // each cell during the actual gradient calculation
      int const nb_cells = mesh_.num_entities(Entity_kind::CELL);
      cell_neighbors_.resize(nb_cells);

      if (parts_ == nullptr) /* entire mesh */ {
        auto collect_neighbors = [this](int c) {
          mesh_.cell_get_node_adj_cells(c, Entity_type::ALL, &(cell_neighbors_[c]));
        };

        Portage::for_each(mesh_.begin(Entity_kind::CELL),
                          mesh_.end(Entity_kind::CELL), collect_neighbors);
      } else /* only on source part */ {
        auto filter_neighbors = [this](int c) {
          cell_neighbors_[c] = parts_->get_source_filtered_neighbors(c);
        };

        auto const& part_cells = parts_->get_source_entities();
        Portage::for_each(part_cells.begin(), part_cells.end(), filter_neighbors);
      }

      set_interpolation_variable(var_name, limiter_type, boundary_limiter_type);
    }

#ifdef HAVE_TANGRAM
    // Constructor with interface reconstructor for multimaterial remaps.
    Limited_Gradient(MeshType const& mesh,
                     StateType const& state,
                     std::string const& var_name,
                     Limiter_type limiter_type,
                     Boundary_Limiter_type boundary_limiter_type,
                     std::shared_ptr<InterfaceReconstructor> ir,
                     const Parts* const parts = nullptr)
      : mesh_(mesh),
        state_(state),
        values_(nullptr),
        variable_name_(var_name),
        limiter_type_(limiter_type),
        boundary_limiter_type_(boundary_limiter_type),
        interface_reconstructor_(ir),
        parts_(parts) {

      // Collect and keep the list of neighbors for each CELL as it may
      // be expensive to go to the mesh layer and collect this data for
      // each cell during the actual gradient calculation
      int const nb_cells = mesh_.num_entities(Entity_kind::CELL);
      cell_neighbors_.resize(nb_cells);

      if (parts_ == nullptr) /* entire mesh */ {
        auto collect_neighbors = [this](int c) {
          mesh_.cell_get_node_adj_cells(c, Entity_type::ALL, &(cell_neighbors_[c]));
        };

        Portage::for_each(mesh_.begin(Entity_kind::CELL),
                          mesh_.end(Entity_kind::CELL), collect_neighbors);

      } else /* only on source part */ {
        auto filter_neighbors = [this](int c) {
          cell_neighbors_[c] = parts_->get_source_filtered_neighbors(c);
        };

        auto const& part_cells = parts_->get_source_entities();
        Portage::for_each(part_cells.begin(), part_cells.end(), filter_neighbors);
      }

      //If the field type is a MESH_FIELD, then the corresponding data will
      //be stored. If the field_type is a MULTIMATERIAL_FIELD, then the constructor
      //only stores the variable name. The user code must make a call to the method
      //set_material to store the material-wise data.
      set_interpolation_variable(var_name, limiter_type, boundary_limiter_type);
    }
#endif

    //This method should be called by the user if the field type is MULTIMATERIAL_FIELD.
    //As the constructor with interface reconstructor only sets the variable name
    //for such fields, this method is needed to properly set the multimaterial data local
    //to the routine.
    void set_material(int matid) {
      material_id_ = matid;
      // Extract the field data from the state manager
      if (field_type_ != Field_type::MESH_FIELD) {
        state_.mat_get_celldata(variable_name_, material_id_, &values_);
      }
    }

    void set_interpolation_variable(std::string const variable_name,
                                    Limiter_type limiter_type,
                                    Boundary_Limiter_type boundary_limiter_type) {

      variable_name_ = variable_name;
      limiter_type_ = limiter_type;
      boundary_limiter_type_ = boundary_limiter_type;

      // Extract the field data from the state manager
      field_type_ = state_.field_type(Entity_kind::CELL, variable_name_);
      if (field_type_ == Field_type::MESH_FIELD) {
        state_.mesh_get_data(Entity_kind::CELL, variable_name_, &values_);
      }
    }

    // @brief Implementation of Limited_Gradient functor for CELLs
    Vector<D> operator()(int cellid) {

      assert(values_);

      double phi = 1.0;
      Vector<D> grad;

      // useful predicates
      bool is_boundary_cell = mesh_.on_exterior_boundary(Entity_kind::CELL, cellid);
      bool apply_limiter = limiter_type_ == BARTH_JESPERSEN and
                           (not is_boundary_cell
                            or boundary_limiter_type_ == BND_BARTH_JESPERSEN);

      // Limit the boundary gradient to enforce monotonicity preservation
      if (is_boundary_cell and boundary_limiter_type_ == BND_ZERO_GRADIENT) {
        grad.zero();
        return grad;
      }

      // Include cell where grad is needed as first element
      std::vector<int> neighbors{cellid};

      if (not cell_neighbors_.empty()) {
        neighbors.insert(std::end(neighbors),
                         std::begin(cell_neighbors_[cellid]),
                         std::end(cell_neighbors_[cellid]));
      }

      std::vector<Point<D>> list_coords;
      std::vector<double> list_values;

      // Loop over cell where grad is needed and its neighboring cells
      for (auto&& neigh_global : neighbors) {
#ifdef HAVE_TANGRAM
        // Field values for each cells in each material are stored according to
        // the material's cell list. So, get the local index of each neighbor cell
        // in the material cell list to access the correct field value
        int const neigh_local = (field_type_ == Field_type::MESH_FIELD)
          ? neigh_global
          : state_.cell_index_in_material(neigh_global, material_id_);

        // In case of material-data, we need to check that neighbors contain
        // material of interest (i.e. have valid local id);
        // in the case of mesh data, this is always true, since local cellid
        // (neigh_local) is equal to global cellid (neigh_global).
        // nota bene: cell_index_in_material can return -1.
        if (neigh_local >= 0) {
          std::vector<int> cell_mats;
          state_.cell_get_mats(neigh_global, &cell_mats);
          int const nb_mats = cell_mats.size();

          if (nb_mats > 1 and interface_reconstructor_) /* multi-material cell */ {
            // Get cell's cellmatpoly
            auto cell_matpoly = interface_reconstructor_->cell_matpoly_data(neigh_global);

            // Collect all the matpolys in this cell for the material of interest
            auto matpolys = cell_matpoly.get_matpolys(material_id_);

            // If there are multiple matpolys in this cell for the material of interest,
            // aggregate moments to compute new centroid
            double mvol = 0.0;
            Point<D> centroid;
            for (auto&& poly : matpolys) {
              auto moments = poly.moments();
              mvol += moments[0];
              for (int k = 0; k < D; k++)
                centroid[k] += moments[k + 1];
            }

            // There are cases where r3d returns a single zero volume poly which
            // ultimately produces a nan, so here we are explicitly filtering this
            // case.
            if (mvol==0.) continue;

            for (int k = 0; k < D; k++)
              centroid[k] /= mvol;

            list_coords.push_back(centroid);

            // Populate least squares vectors with centroid for material
            // of interest and field value in the current cell for that material
            list_values.push_back(values_[neigh_local]);
          } else if (nb_mats == 1) /* single-material cell */ {

            // Ensure that the single material is the material of interest
            if (cell_mats[0] == material_id_) {
              // Get the cell-centered value for this material
              Point<D> centroid;
              mesh_.cell_centroid(neigh_global, &centroid);
              list_coords.push_back(centroid);
              list_values.push_back(values_[neigh_local]);
            }
          }
        }
#endif
        // If we get here, we must have mesh data which is cell-centered
        // and not dependent on material, so just get the centroid and value
        if (field_type_ == Field_type::MESH_FIELD) {
          Point<D> centroid;
          mesh_.cell_centroid(neigh_global, &centroid);
          list_coords.push_back(centroid);
          list_values.push_back(values_[neigh_global]);
        }
      }

      grad = Wonton::ls_gradient<D, CoordSys>(list_coords, list_values);

      // Limit the gradient to enforce monotonicity preservation
      if (apply_limiter) {

        phi = 1.0;

        // Min and max vals of function (cell centered vals) among neighbors
        // and the cell itself
        /// @todo: must remove assumption the field is scalar
        double minval = list_values[0];
        double maxval = list_values[0];
        double cellcenval = list_values[0];

        // Find min and max values among all neighbors (exlude the first element
        // in nbrids because it corresponds to the cell itself, not a neighbor)
        int const nb_values = list_values.size();
        for (int i = 1; i < nb_values; ++i) {
          minval = std::min(list_values[i], minval);
          maxval = std::max(list_values[i], maxval);
        }

        /* Per page 278 of [Kucharik, M. and Shaskov, M, "Conservative
           Multi-material Remap for Staggered Multi-material Arbitrary
           Lagrangian-Eulerian Methods," Journal of Computational Physics,
           v 258, pp. 268-304, 2014], if a cell is a multimaterial cell and the
           field is a material field, then the min value (for density, at
           least) should be set to 0 and the max value to infinity (or a large
           number), so that we don't end up limiting the gradient to 0 and drop
           down to 1st order. But we don't know if a variable is density (don't
           want to do silly things like string comparison) or pressure or
           something else? What if we limit pressure to 0 and it should
           actually be allowed to go -ve? In the end, along with the limiter
           type, the application should be able to tell Portage the global
           bounds that a variable has to satisfy. Then we can impose the global
           limits at multi-material cells and boundary cells without limiting
           the gradient to 0. */

        // Find the min and max of the reconstructed function in the cell
        // Since the reconstruction is linear, this will occur at one of
        // the nodes of the cell. So find the values of the reconstructed
        // function at the nodes of the cell
        int dim = mesh_.space_dimension();
        std::vector<Point<D>> cellcoords;
        mesh_.cell_get_coordinates(cellid, &cellcoords);

        for (auto&& coord : cellcoords) {
          auto vec = coord - list_coords[0];
          double diff = dot(grad, vec);
          double extremeval = (diff > 0.) ? maxval : minval;
          double phi_new = (diff == 0. ? 1. : (extremeval - cellcenval) / diff);
          phi = std::min(phi_new, phi);
        }
      }

      // Limited gradient is phi*grad
      return phi * grad;
    }

  private:
    MeshType const& mesh_;
    StateType const& state_;
    double const* values_;
    std::string variable_name_ = "";
    Limiter_type limiter_type_ = DEFAULT_LIMITER;
    Boundary_Limiter_type boundary_limiter_type_ = DEFAULT_BND_LIMITER;
    Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;

    int material_id_ = 0;
    std::vector<int> cell_ids_;
    std::vector<std::vector<int>> cell_neighbors_;
#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor> interface_reconstructor_;
#endif
    Parts const* parts_;
  };

  ///////////////////////////////////////////////////////////////////////////////

  /*! @class Limited_Gradient<MeshType,StateType,NODE> gradient.h
    @brief Specialization of limited gradient class for @c node-centered field
    @tparam MeshType A mesh class that one can query for mesh info
    @tparam StateType A state manager class that one can query for field info
  */

  template<
    int D, typename MeshType, typename StateType,
    template<class, int, class, class>
      class InterfaceReconstructorType,
    class Matpoly_Splitter, class Matpoly_Clipper, class CoordSys
  >
  class Limited_Gradient<
    D, Entity_kind::NODE, MeshType, StateType,
    InterfaceReconstructorType,
    Matpoly_Splitter, Matpoly_Clipper, CoordSys
  > {

    using Parts = PartPair<D, Entity_kind::CELL, MeshType, StateType>;

  public:
    /*! @brief Constructor
      @param[in] mesh  Mesh class than one can query for mesh info
      @param[in] state A state manager class that one can query for field info
      @param[in] var_name Name of field for which the gradient is to be computed
      @param[in] limiter_type An enum indicating if the limiter type (none, Barth-Jespersen, Superbee etc)
      @param[in] boundary_limiter_type An enum indicating the limiter type on the boundary

      @todo must remove assumption that field is scalar
    */
    Limited_Gradient(MeshType const& mesh,
                     StateType const& state,
                     std::string const var_name,
                     Limiter_type limiter_type,
                     Boundary_Limiter_type boundary_limiter_type,
                     const Parts* const parts = nullptr)
      : mesh_(mesh),
        state_(state),
        values_(nullptr),
        variable_name_(var_name),
        limiter_type_(limiter_type),
        boundary_limiter_type_(boundary_limiter_type),
        parts_(parts) {

      auto collect_node_neighbors = [this](int n) {
        this->mesh_.dual_cell_get_node_adj_cells(
          n, Entity_type::ALL, &(node_neighbors_[n])
        );
      };

      int const nnodes = mesh_.num_entities(Entity_kind::NODE);
      node_neighbors_.resize(nnodes);
      Portage::for_each(mesh_.begin(Entity_kind::NODE),
                        mesh_.end(Entity_kind::NODE),
                        collect_node_neighbors);

      set_interpolation_variable(var_name, limiter_type, boundary_limiter_type);

      if (parts_ != nullptr) {
        std::cout << "Warning: part-by-part remap is only defined for cells. ";
        std::cout << "Source and target part pair will be ignored" << std::endl;
      }
    }

    void set_material(int material_id) { material_id_ = material_id; }

    void set_interpolation_variable(std::string const& variable_name,
                                    Limiter_type limiter_type,
                                    Boundary_Limiter_type boundary_limiter_type) {

      variable_name_ = variable_name;
      limiter_type_ = limiter_type;
      boundary_limiter_type_ = boundary_limiter_type;
      state_.mesh_get_data(Entity_kind::NODE, variable_name_, &values_);
    }

    // @brief Limited gradient functor implementation for NODE
    Vector<D> operator()(int nodeid) {

      assert(values_);

      double phi = 1.0;
      Vector<D> grad;

      bool is_boundary_node = mesh_.on_exterior_boundary(Entity_kind::NODE, nodeid);
      bool apply_limiter = limiter_type_ == BARTH_JESPERSEN and
                           (not is_boundary_node
                            or boundary_limiter_type_ == BND_BARTH_JESPERSEN);

      if (is_boundary_node and boundary_limiter_type_ == BND_ZERO_GRADIENT) {
        grad.zero();
        return grad;
      }

      auto const& neighbors = node_neighbors_[nodeid];
      std::vector<Point<D>> node_coords(neighbors.size() + 1);
      std::vector<double> node_values(neighbors.size() + 1);
      mesh_.node_get_coordinates(nodeid, &(node_coords[0]));
      node_values[0] = values_[nodeid];

      int i = 1;
      for (auto&& current : neighbors) {
        mesh_.node_get_coordinates(current, &node_coords[i]);
        node_values[i] = values_[current];
        i++;
      }

      grad = Wonton::ls_gradient<D, CoordSys>(node_coords, node_values);

      if (apply_limiter) {
        // Min and max vals of function (cell centered vals) among neighbors
        double minval = values_[nodeid];
        double maxval = values_[nodeid];

        for (auto const& val : node_values) {
          minval = std::min(val, minval);
          maxval = std::max(val, maxval);
        }

        // Find the min and max of the reconstructed function in the cell
        // Since the reconstruction is linear, this will occur at one of
        // the nodes of the cell. So find the values of the reconstructed
        // function at the nodes of the cell
        double nodeval = values_[nodeid];

        std::vector<Point<D>> dual_cell_coords;
        mesh_.dual_cell_get_coordinates(nodeid, &dual_cell_coords);

        for (auto&& coord : dual_cell_coords) {
          auto vec = coord - node_coords[0];
          double diff = dot(grad, vec);
          double extremeval = (diff > 0.0) ? maxval : minval;
          double phi_new = (diff == 0.0) ? 1 : (extremeval - nodeval) / diff;
          phi = std::min(phi_new, phi);
        }
      }

      // Limited gradient is phi*grad
      return phi * grad;
    }

  private:
    MeshType const& mesh_;
    StateType const& state_;
    double const* values_;
    std::string variable_name_ = "";
    Limiter_type limiter_type_ = DEFAULT_LIMITER;
    Boundary_Limiter_type boundary_limiter_type_ = DEFAULT_BND_LIMITER;
    Field_type field_type_ = Field_type::UNKNOWN_TYPE_FIELD;

    int material_id_ = 0;
    std::vector<int> cell_ids_;
    std::vector<std::vector<int>> node_neighbors_;
    Parts const* parts_;
  };
}  // namespace Portage

#endif  // PORTAGE_INTERPOLATE_GRADIENT_H_
