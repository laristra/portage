/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
////////////////////////////////////////////////////////////////////////////////

#pragma once

// user includes

// library includes
#include <portage/support/portage.h>

// system includes
#include <utility>
#include <cstring>
#include <string>

//namespace flecsale {
namespace Wonton {

////////////////////////////////////////////////////////////////////////////////
/// \brief Provides access to data stored in Flecsi_State
////////////////////////////////////////////////////////////////////////////////
template< typename M >
class flecsi_state_t {

  //============================================================================
  // Typedefs
  //============================================================================

  //! \brief The mesh type
  using mesh_t = M;
  //! \brief the size type
  using size_t = typename mesh_t::size_t;
  //! \brief the real type
  using real_t = typename mesh_t::real_t;


  //! \brief The entity kind type
  using entity_kind_t = Portage::Entity_kind;
  //! \brief The entity type 
  using entity_type_t = Portage::Entity_type;

public:


  //============================================================================
  // Constructors
  //============================================================================

  //!  \brief Default constructor.
  //!  \param[in] mesh The minimum coordinates of the domain.
  explicit flecsi_state_t(mesh_t & mesh) : mesh_(&mesh)
  {}

  //! Default constructor deleted
  flecsi_state_t() = default;

  //! Default copy constructor
  flecsi_state_t(const flecsi_state_t &) = default;

  //! Default assignment operator
  flecsi_state_t & operator=(const flecsi_state_t &) = default;

  //============================================================================
  // Public Members
  //============================================================================


  //! \brief Number of materials in problem

  int num_materials() const {
    return 0;
  }

  //! \brief Name of material

  std::string material_name(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return "UNKNOWN";
  }

  //! \brief Get number of cells containing a particular material
  //! \param matid    Index of material (0, num_materials()-1)
  //! \return         Number of cells containing material 'matid'

  int mat_get_num_cells(int matid) const {
    assert(matid >= 0 && matid < num_materials());
    return 0;
  }

  //! \brief Get cell indices containing a particular material
  //! \param matid    Index of material (0, num_materials()-1)
  //! \param matcells Cells containing material 'matid'

  void mat_get_cells(int matid, std::vector<int> *matcells) const {
    assert(matid >= 0 && matid < num_materials());
    matcells->clear();
  }

  //! \brief Get number of materials contained in a cell
  //! \param cellid  Index of cell in mesh
  //! \return        Number of materials in cell

  int cell_get_num_mats(int cellid) const {
    return 0;
  }

  //! \brief Get the IDs of materials in a cell
  //! \param cellid    Index of cell in mesh
  //! \param cellmats  Indices of materials in cell

  void cell_get_mats(int cellid, std::vector<int> *cellmats) const {
    cellmats->clear();
  }

  //! \brief Get the local index of mesh cell in material cell list
  //! \param meshcell    Mesh cell ID
  //! \param matid       Material ID
  //! \return             Local cell index in material cell list

  int cell_index_in_material(int meshcell, int matid) const {
    return -1;
  }
  
  //! \brief Type of field (MESH_FIELD or MULTIMATERIAL_FIELD)
  //! \param[in] onwhat   Entity_kind that field is defined on
  //! \param[in] varname  Name of field
  //! \return             Field_type

  Field_type field_type(Entity_kind on_what, std::string const& var_name)
      const {
    return Field_type::MESH_FIELD;  // MULTI-MATERIAL FIELDS NOT ACCESSED YET
  }
  

  //! \brief Get the entity type on which the given field is defined
  //! \param[in] var_name The string name of the data field
  //! \return The Entity_kind enum for the entity type on which the field is defined
  //!
  //! \todo  THIS ASSUMES ONLY DOUBLE VECTORS - WE HAVE TO ACCOUNT FOR OTHER TYPES
  //!        OR WE HAVE TO GENERALIZE THE FIND FUNCTION!!!
  //! \todo  THIS ALSO DOES NOT CHECK FOR OTHER ENTITY TYPES LIKE EDGE, FACE,
  //!        SIDE, WEDGE AND CORNER
  entity_kind_t get_entity(std::string const& var_name) const 
  {
    auto cell_field_list = flecsi_get_accessors_all(
      *mesh_, real_t, dense, 0, flecsi_is_at(cells)
    );
    for (auto var : cell_field_list)
      if (var.label() == var_name) 
        return entity_kind_t::CELL;

    auto nodal_field_list = flecsi_get_accessors_all(
      *mesh_, real_t, dense, 0, flecsi_is_at(vertices)
    );
    for (auto var : nodal_field_list)
      if (var.label() == var_name) 
        return entity_kind_t::NODE;

    return entity_kind_t::UNKNOWN_KIND;
  }


  //! \brief Get pointer to scalar data
  //! \param[in] on_what The entity type on which to get the data
  //! \param[in] var_name The string name of the data field
  //! \param[in,out] data A pointer to an array of data
  template <class T>
  void mesh_get_data(entity_kind_t on_what, std::string const& var_name, 
                     T ** data) const {
    // Ignore on_what here - the state manager knows where it lives
    // based on its name

    // first check cells
    auto cell_field_list = flecsi_get_accessors_all(
      *mesh_, real_t, dense, 0, flecsi_is_at(cells)
    );
    for (auto var : cell_field_list)
      if (var.label() == var_name) {
        *data = &var[0];
        return;
      }

    // now check nodes
    auto nodal_field_list = flecsi_get_accessors_all(
      *mesh_, real_t, dense, 0, flecsi_is_at(vertices)
    );
    for (auto var : nodal_field_list)
      if (var.label() == var_name)  {
        *data = &var[0];
        return;
      }

    // if we got here, there is something wrong
    raise_runtime_error( "Could not find variable to ReMAP!" );
  }


  //! \brief Get pointer to read-only scalar cell data for a particular material
  //! \param[in] var_name The string name of the data field
  //! \param[in] matid   Index (not unique identifier) of the material
  //! \param[out] data   vector containing the values corresponding to cells in the material

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


  //! \brief Get pointer to read-write scalar data for a particular material
  //! \param[in] on_what The entity type on which to get the data
  //! \param[in] var_name The string name of the data field
  //! \param[in] matid   Index (not unique identifier) of the material
  //! \param[out] data   vector containing the values corresponding to cells in the material

  void mat_get_celldata(std::string const& var_name, int matid, double **data) {
  }

  //! \brief Get a pointer to data from the state manager with a given
  //! variable @c name and on @c on_what mesh entities.
  //! \param[in] on_what The Entity_kind (e.g. CELL) on which the data lives.
  //! \param[in] name The name of the variable.
  //! \param data A @c pointer to the const data array.  If the requested
  //! data is not found in the state manager, a @c nullptr is returned.

  void mesh_add_data(Entity_kind on_what, std::string const& name,
                  double const **data) const {
  }

  //! \brief Add a scalar multi-valued data field on cells and initialize its
  //! material data to a single value
  //! \param[in] var_name The name of the data field
  //! \param[in] value Initialize with this value

  //! The 2D array will be read and values copied according to which materials
  //! are contained in which cells. If a material+cell combination is not active
  //! providing a value for the array will have no effect. 

  void mat_add_celldata(std::string const& var_name, double value) {
  }


  //! \brief Add a scalar multi-valued data field on cells and initialize its
  //! material data according to a 2D array
  //! \param[in] var_name The name of the data field
  //! \param[in] layout  Whether 2D array is laid out with first index being
  //! the cell (CELL_CENRIC) or material (MATERIAL CENTRIC)
  //! \param[in] value Initialize with this value
  //!
  //! The 2D array will be read and values copied according to which
  //! materials are contained in which cells. If a material+cell
  //! combination is not active providing a value for the array will
  //! have no effect.

  void mat_add_celldata(std::string const& var_name,
                        double const * const *values = nullptr,
                        Data_layout layout = Data_layout::MATERIAL_CENTRIC) {
  }


  //! \brief Add a scalar multi-valued data field on cells and add
  //! data to one of its materials
  //! \param[in] var_name The name of the data field
  //! \param[in] matid  Index of material in the problem
  //! \param[in] layout Data layout - 
  //! \param[in] values Initialize with this array of values
  //!
  //! Subsequent calls to this function with the same name will find the added
  //! field and just add the data.

  void mat_add_celldata(std::string const& var_name, int matid,
                        double const * values) {
  }

  void mat_add_celldata(std::string const& var_name, int matid,
                        Portage::Point<2> const *values) {
  }

  void mat_add_celldata(std::string const& var_name, int matid,
                        Portage::Point<3> const *values) {
  }


  //! \brief Add a scalar multi-valued data field on cells and initialize one
  //! of its material data to a uniform value
  //! \param[in] var_name The name of the data field
  //! \param[in] matid Index of material in the problem
  //! \param[in] value Initialize with this value
  //! Subsequent calls to this function with the same name will find the added
  //! field and just add the data.

  void mat_add_celldata(std::string const& var_name, int matid, double value) {
  }


  //! \brief Add cells to material (or add material to cells)
  //! \param[in] matid  Material ID
  //! \param[in] newcells Vector of new cells in material

  void mat_add_cells(int matid, std::vector<int> const& newcells) {
  }


  //! \brief Remove cells from material (or remove material from cells)
  //! \param[in] matid  Material ID
  //! \param[in] matcells Vector of to be removed cells
  
  void mat_rem_cells(int matid, std::vector<int> const& delcells) {
  }


  //! \brief Add a material to state
  //! \param[in] matname  Name of material
  //! \param[in] matcells Cells containing the material

  void add_material(std::string const& matname,
                    std::vector<int> const& matcells) {
  }


  int get_data_size(entity_kind_t on_what, std::string const& var_name) const 
  {
    raise_runtime_error( "get_data_size not implemented yet!" );
  }



  //!
  //! @brief Get the data type of the given field
  //! @param[in] var_name The string name of the data field
  //! @return A reference to the type_info struct for the field's data type

  const std::type_info& get_data_type(std::string const& var_name) const {
    return typeid(double);  // thats the only type we can represent
  }
  
  /*!
    @brief  Vector of names
    @return vector of strings
   */
  std::vector<std::string> names() const { 
    return std::vector<std::string>{};
  }
  

#if 0

  //! \brief Get pointer to scalar data
  //! \param[in] on_what The entity type on which to get the data
  //! \param[in] var_name The string name of the data field
  //! \param[in,out] data A pointer to an array of data
  template <class T>
  void mesh_get_data(entity_kind_t on_what,
                     std::string const& var_name, T** data) {
    // Ignore on_what here - the state manager knows where it lives
    // based on its name
    // NOTE: access_state needs a const_string_t, which is a flecsi
    // datatype that apparently can't convert from a std::string...
    // I added to the case where it can create from a char*
    const const_string_t flecsiname(var_name.c_str());
    auto dat = access_state(flecsi_mesh_,
                            std::forward<const const_string_t>(flecsiname),
                            real_t);
    // note this might be fragile in the case of non-dense storage
    *data = &dat[0];
  }

  //! \brief Get the data size for the given field
  //! \param[in] on_what  The entity type on which the data field is defined
  //! \param[in] var_name The string name of the data field
  //! \return The data size for the field with the given name on the given entity type
  int get_data_size(entity_kind_t on_what, std::string const& var_name) const {

   const const_string_t flecsiname(var_name.c_str());
   auto dat = access_state(flecsi_mesh_,
                           std::forward<const const_string_t>(flecsiname),
                           real_t);
   return (dat.size());

  }

#endif 

private:

  //! \brief the flecsi mesh pointer
  mesh_t * mesh_ = nullptr;

};  // Flecsi_State_Wrapper

} // namespace wonton 
//} // namespace 
