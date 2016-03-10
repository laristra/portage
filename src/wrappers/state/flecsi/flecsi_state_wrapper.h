/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef FLECSI_STATE_WRAPPER_H_
#define FLECSI_STATE_WRAPPER_H_

#include <utility>
#include <cstring>

#include "flecsi/specializations/burton/burton.h"
#include "flecsi/utils/const_string.h"

#include "portage/support/portage.h"

using mesh_t = flecsi::burton_mesh_t;

using cell_t = flecsi::mesh_t::cell_t;
using vertex_t = flecsi::mesh_t::vertex_t;

using real_t = flecsi::mesh_t::real_t;

using const_string_t = flecsi::const_string_t;

/*!
  @file flecsi_state_wrapper.h
  @brief Wrapper for interfacing with the Flecsi state manager
 */

namespace Portage {

/*!
  @class Flecsi_State_Wrapper "flecsi_state_wrapper.h"
  @brief Provides access to data stored in Flecsi_State
*/
class Flecsi_State_Wrapper {
 public:

  /*!
    @brief Constructor of Flecsi_State_Wrapper
    @param[in] flecsi_mesh A reference to a flecsi::burton_mesh_t instance
   */
  Flecsi_State_Wrapper(mesh_t & flecsi_mesh) : flecsi_mesh_(flecsi_mesh) {}

  /*!
    @brief Copy constructor of Flecsi_State_Wrapper - not a deep copy
    @param[in] state A reference to another Flecsi_State_Wrapper instance
   */
  Flecsi_State_Wrapper(Flecsi_State_Wrapper & state)
      : flecsi_mesh_(state.flecsi_mesh_) {}

  /*!
    @brief Assignment operator (disabled) - don't know how to implement (RVG)
   */
  Flecsi_State_Wrapper & operator=(Flecsi_State_Wrapper const &) = delete;

  /*!
    @brief Empty destructor
   */
  ~Flecsi_State_Wrapper() {};

  /* /\*! */
  /*   @brief Initialize fields from mesh file */
  /*  *\/ */
  /* void init_from_mesh() { flecsi_mesh_.init_from_mesh(); } */

  /* /\*! */
  /*   @brief Export fields to mesh file */
  /*  *\/ */
  /* void export_to_mesh() {flecsi_mesh_.export_to_mesh(); } */

  /*!
    @brief Get pointer to scalar data
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of data
  */
  template <class T>
  void get_data(const Entity_kind on_what,
                const std::string var_name, T** const data) const {
    // Ignore on_what here - the state manager knows where it lives
    // based on its name
    // NOTE: access_state needs a const_string_t, which is a flecsi
    // datatype that apparently can't convert from a std::string...
    // I added to the case where it can create from a char*
    const const_string_t flecsiname(var_name.c_str());
    auto dat = access_state(flecsi_mesh_,
                            //                            flecsiname,
                            std::forward<const const_string_t>(flecsiname),
                            real_t);
    // NOTE: I added this method
    *data = dat.data();
  }

  /* /\*! */
  /*   @brief Get the entity type on which the given field is defined */
  /*   @param[in] var_name The string name of the data field */
  /*   @return The Entity_kind enum for the entity type on which the field is defined */
  /*  *\/ */
  /* Entity_kind get_entity(const std::string var_name) const { */

  /*   std::shared_ptr<Flecsi::BaseStateVector> vector =  */
  /*       *(flecsi_mesh_.find(var_name, Flecsi::ANY_KIND)); */
  /*   if (vector != 0) return (Portage::Entity_kind) vector->on_what(); */

  /*   return Portage::UNKNOWN_KIND; */

  /* } */

  /* /\*! */
  /*   @brief Get the data type of the given field */
  /*   @param[in] var_name The string name of the data field */
  /*   @return A reference to the type_info struct for the field's data type */
  /*  *\/  */
  /* const std::type_info& get_type(const std::string var_name) const { */
    
  /*   std::shared_ptr<Flecsi::BaseStateVector> vector =  */
  /*       *(flecsi_mesh_.find(var_name,Flecsi::ANY_KIND)); */
  /*   if (vector != 0) return vector->get_type(); */
  /*   return typeid(0); */

  /* } */

  /* /\*! */
  /*   @brief Begin iterator on vector names */
  /*   @return Begin iterator on vector of strings */
  /*  *\/ */
  /* std::vector<std::string>::iterator names_begin() const {  */
  /*   return flecsi_mesh_.names_begin();  */
  /* } */

  /* /\*! */
  /*   @brief End iterator on vector names */
  /*   @return End iterator on vector of strings */
  /*  *\/ */
  /* std::vector<std::string>::iterator names_end() const {  */
  /*   return flecsi_mesh_.names_end();  */
  /* } */

  /* /\*! */
  /*   @brief Typedef for permutation iterator on vector of strings */
  /*  *\/ */
  /* typedef Flecsi::State::string_permutation string_permutation; */

  /* /\*! */
  /*   @brief Begin iterator on vector names of specific entity type */
  /*   @param[in] on_what The desired entity type */
  /*   @return Permutation iterator to start of string vector */
  /*  *\/ */
  /* string_permutation names_entity_begin(Entity_kind const on_what) const {  */
  /*   return flecsi_mesh_.names_entity_begin((Flecsi::Entity_kind)on_what);  */
  /* } */

  /* /\*! */
  /*   @brief End iterator on vector of names of specific entity type */
  /*   @param[in] on_what The desired entity type */
  /*  *\/ */
  /* string_permutation names_entity_end(Entity_kind const on_what) const {  */
  /*   return flecsi_mesh_.names_entity_end((Flecsi::Entity_kind)on_what);  */
  /* } */

 private:

  mesh_t & flecsi_mesh_;

}; // Flecsi_State_Wrapper

} // namespace Portage

#endif // FLECSI_STATE_WRAPPER_H_
