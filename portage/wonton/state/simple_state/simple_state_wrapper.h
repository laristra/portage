/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#ifndef SIMPLE_STATE_WRAPPER_H_
#define SIMPLE_STATE_WRAPPER_H_

#include <string>

#include "portage/simple_mesh/simple_state.h"

#include "portage/support/portage.h"

/*!
  @file simple_state_wrapper.h
  @brief Definitions for a wrapper to Simple_State
*/
namespace Wonton {

using namespace Portage;
  /*!
    @class Simple_State_Wrapper simple_state_wrapper.h
    @brief A thin wrapper that implements state methods for Simple_State needed
    by Portage.
   */
class Simple_State_Wrapper {
 public:
  /*!
    @brief Constructor for the state wrapper
   */
  explicit Simple_State_Wrapper(Simple_State & state) : state_(state) { }

  /// Assignment operator (disabled).
  Simple_State_Wrapper & operator=(Simple_State_Wrapper const &) = delete;

  /// Destructor.
  ~Simple_State_Wrapper() { }

  //  void init_from_mesh() { state_.init_from_mesh(); }
  //  void export_to_mesh() { state_.export_to_mesh(); }

  /*!
    @brief Number of materials in problem
  */

  int num_materials() const {
    return 0;
  }

  /*!
    @brief Name of material
  */

  std::string material_name(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return "UNKNOWN";
  }

  /*!
    @brief Get number of cells containing a particular material
    @param matid    Index of material (0, num_materials()-1)
    @return         Number of cells containing material 'matid'
  */

  int mat_get_num_cells(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return 0;
  }

  /*!
    @brief Get cell indices containing a particular material
    @param matid    Index of material (0, num_materials()-1)
    @param matcells Cells containing material 'matid'
  */

  void mat_get_cells(int matid, std::vector<int> *matcells) const {
    assert(matid >= 0 && matid < num_materials());
    matcells->clear();
  }

  /*!
    @brief Get number of materials contained in a cell
    @param cellid  Index of cell in mesh
    @return        Number of materials in cell
  */

  int cell_get_num_mats(int cellid) const {
    return 0;
  }

  /*!
    @brief Get the IDs of materials in a cell
    @param cellid    Index of cell in mesh
    @param cellmats  Indices of materials in cell
  */

  void cell_get_mats(int cellid, std::vector<int> *cellmats) const {
    cellmats->clear();
  }

  /*!
    @brief Get the local index of mesh cell in material cell list
    @param meshcell    Mesh cell ID
    @param matid       Material ID
    @return             Local cell index in material cell list
  */

  int cell_index_in_material(int meshcell, int matid) const {
    return -1;
  }
  
  /*!
    @brief Type of field (MESH_FIELD or MULTIMATERIAL_FIELD)
    @param onwhat    Entity_kind that field is defined on
    @param varname   Name of field
    @return          Field type
  */

  Field_type field_type(Entity_kind on_what, std::string const& var_name)
      const {
    return Field_type::MESH_FIELD;
  }

  /*!
    @brief Get a pointer to data from the state manager with a given
    variable @c name and on @c on_what mesh entities.
    @param[in] on_what The Entity_kind (e.g. CELL) on which the data
    lives.
    @param[in] name The name of the variable.
    @param data A pointer to the data array.  If the requested data is not
    found in the state manager, a @c nullptr is returned.
   */
  void mesh_get_data(Entity_kind on_what, std::string const& name,
                double **data) {
    auto it = state_.find(name, on_what);
    if (it != state_.end()) {
      (*data) = &(it->second[0]);
      return;
    }

    std::cerr << "get_data: Could not find state variable " << name
              << std::endl;
    (*data) = nullptr;
  }

  /*!
    @brief Get a pointer to data from the state manager with a given
    variable @c name and on @c on_what mesh entities.
    @param[in] on_what The Entity_kind (e.g. CELL) on which the data
    lives.
    @param[in] name The name of the variable.
    @param data A @c pointer to the const data array.  If the requested data is
    not found in the state manager, a @c nullptr is returned.
   */
  void mesh_get_data(Entity_kind on_what, std::string const& name,
                  double const **data) const {
    auto it = state_.find(name, on_what);
    if (it != state_.end()) {
      (*data) = &(it->second[0]);
      return;
    }

    std::cerr << "get_data: Could not find state variable " << name
              << std::endl;
    (*data) = nullptr;
  }

  /*!
    @brief Get pointer to read-only scalar cell data for a particular material
    @param[in] var_name The string name of the data field
    @param[in] matid   Index (not unique identifier) of the material
    @param[out] data   vector containing the values corresponding to cells in the material
   */

  void mat_get_celldata(std::string const& var_name, int matid,
                        double const **data) const {
  }

  // TEMPORARY: UNTIL THIS GETS TEMPLATED ON TYPE OF DATA
  void mat_get_celldata(std::string const& var_name, int matid,
                        Portage::Point<2> const **data) const {
  }
  void mat_get_celldata(std::string const& var_name, int matid,
                        Portage::Point<3> const **data) const {
  }


  /*!
    @brief Get pointer to read-write scalar data for a particular material
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in] matid   Index (not unique identifier) of the material
    @param[out] data   vector containing the values corresponding to cells in the material
   */

  void mat_get_celldata(std::string const& var_name, int matid, double **data) {
  }

  /*!
    @brief Get a pointer to data from the state manager with a given
    variable @c name and on @c on_what mesh entities.
    @param[in] on_what The Entity_kind (e.g. CELL) on which the data
    lives.
    @param[in] name The name of the variable.
    @param data A @c pointer to the const data array.  If the requested data is
    not found in the state manager, a @c nullptr is returned.
   */
  void mesh_add_data(Entity_kind on_what, std::string const& name,
                  double const **data) const {
    auto vec = state_.add(name, on_what, *data);
  }

  /*!
   @brief Add a scalar multi-valued data field on cells and initialize its
   material data to a single value
   @param[in] var_name The name of the data field
   @param[in] value Initialize with this value

   The 2D array will be read and values copied according to which materials
   are contained in which cells. If a material+cell combination is not active
   providing a value for the array will have no effect. 

  */

  void mat_add_celldata(std::string const& var_name, double value) {
  }


  /*!
   @brief Add a scalar multi-valued data field on cells and initialize its
   material data according to a 2D array
   @param[in] var_name The name of the data field
   @param[in] layout  Whether 2D array is laid out with first index being the cell (CELL_CENRIC) or material (MATERIAL CENTRIC)
   @param[in] value Initialize with this value

   The 2D array will be read and values copied according to which materials
   are contained in which cells. If a material+cell combination is not active
   providing a value for the array will have no effect. 
   */

  void mat_add_celldata(std::string const& var_name,
                        double const * const *values = nullptr,
                        Data_layout layout = Data_layout::MATERIAL_CENTRIC) {
  }


  /*!
   @brief Add a scalar multi-valued data field on cells and add data to one of
   its materials
   @param[in] var_name The name of the data field
   @param[in] matid  Index of material in the problem
   @param[in] layout Data layout - 
   @param[in] values Initialize with this array of values

   Subsequent calls to this function with the same name will find the added
   field and just add the data.
   */

  void mat_add_celldata(std::string const& var_name, int matid,
                        double const * values) {
  }


  /*!
   @brief Add a scalar multi-valued data field on cells and initialize one of
   its material data to a uniform value
   @param[in] var_name The name of the data field
   @param[in] matid Index of material in the problem
   @param[in] value Initialize with this value

   Subsequent calls to this function with the same name will find the added
   field and just add the data.
   */

  void mat_add_celldata(std::string const& var_name, int matid, double value) {
  }


  /*!
    @brief Add cells to material (or add material to cells)
    @param[in] matid  Material ID
    @param[in] newcells Vector of new cells in material
  */

  void mat_add_cells(int matid, std::vector<int> const& newcells) {
  }


  /*!
    @brief Remove cells from material (or remove material from cells)
    @param[in] matid  Material ID
    @param[in] matcells Vector of to be removed cells
  */
  
  void mat_rem_cells(int matid, std::vector<int> const& delcells) {
  }


  /*!
    @brief Add a material to state
    @param[in] matname  Name of material
    @param[in] matcells Cells containing the material
   */
  void add_material(std::string const& matname,
                    std::vector<int> const& matcells) {
  }

  /*!
    @brief Given a variable name, get where it lives on the mesh.
    @param[in] name The name of the variable to be found.
    @returns The Entity_kind (e.g. CELL) indicating where the variable @name
    lives.
    @throws std::runtime_exception Variable not found.

    @TODO This will only pick the first defined variable, if the user
    adds multiple fields with the same variable name, but on different
    Entity_kinds.  This should be fixed, or get_entity should be ammended.
   */
  Entity_kind get_entity(std::string const& name) const {
    // Find where the name is in the variable name part of the map's keys
    auto it = state_.begin();
    while (it != state_.end()) {
      if (it->first.first == name)
        break;
      else
        ++it;
    }
    // Now get the second part of the key, if it was found
    if (it != state_.end()) {
      return it->first.second;
    } else {
      // We don't know this variable, so bail.
      std::runtime_error("Requested variable not found.");
    }
  }

  /*!
    @brief Get the number of elements in a specific variable from the state
    manager.
    @param[in] on_what The Entity_kind (e.g. CELL) of the variable for which the
    size is requested.
    @param[in] name The name of the variable for which the size is requested.
    @returns The number of elements; for example, the number of cells, if the
    data were added with the CELLS Entity_kind.
   */
  int get_data_size(Entity_kind on_what, std::string const& name) const {
    auto it = state_.find(name, on_what);
    if (it != state_.end()) {
      return it->second.size();
    }

    std::cerr << "get_data_size: Could not find state variable " << name
              << std::endl;
    return 0;
  }


  /// @brief Get the data type of the given field
  /// @param[in] var_name The string name of the data field
  /// @return A reference to the type_info struct for the field's data type

  const std::type_info& get_data_type(std::string const& var_name) const {
    auto it = state_.find(var_name);
    if (it != state_.end())
      return typeid(double);  // Thats the only type we can handle for now
    else
      return typeid(void);
  }

  /// An iterator to the beginning of the vector of names in the state manager.
  Simple_State::name_vec_it names_begin() const {
    return state_.names_begin();
  }

  /// An iterator to the ending of the vector of names in the state manager.
  Simple_State::name_vec_it names_end() const {
    return state_.names_end();
  }

 private:
  /// The state to be wrapped.
  Simple_State & state_;
};  // class Simple_State_Wrapper
}  // namespace Wonton 

#endif  // SIMPLE_STATE_WRAPPER_H_
