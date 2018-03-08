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


  int get_data_size(entity_kind_t on_what, std::string const& var_name) const 
  {
    raise_runtime_error( "get_data_size not implemented yet!" );
  }
#if 0

  /*!
    @brief Get pointer to scalar data
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of data
  */
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

  /*!
    @brief Get the data size for the given field
    @param[in] on_what  The entity type on which the data field is defined
    @param[in] var_name The string name of the data field
    @return The data size for the field with the given name on the given entity type
   */
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
