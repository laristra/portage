/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef FLECSI_STATE_WRAPPER_H_
#define FLECSI_STATE_WRAPPER_H_

#include <utility>
#include <cstring>
#include <string>

#include "flecsi/specializations/burton/burton.h"
#include "flecsi/utils/const_string.h"

#include "portage/support/portage.h"

using mesh_t = flecsi::burton_mesh_t;

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
  explicit Flecsi_State_Wrapper(mesh_t & flecsi_mesh)
      : flecsi_mesh_(flecsi_mesh) {}

  /*!
    @brief Copy constructor of Flecsi_State_Wrapper - not a deep copy
    @param[in] state A reference to another Flecsi_State_Wrapper instance
   */
  Flecsi_State_Wrapper(Flecsi_State_Wrapper & state)
      : flecsi_mesh_(state.flecsi_mesh_) {}

  /*!
    @brief Assignment operator (disabled)
   */
  Flecsi_State_Wrapper & operator=(Flecsi_State_Wrapper const &) = delete;

  /*!
    @brief Empty destructor
   */
  ~Flecsi_State_Wrapper() {}

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
                            std::forward<const const_string_t>(flecsiname),
                            real_t);
    // note this might be fragile in the case of non-dense storage
    *data = &dat[0];
  }

 private:
  mesh_t & flecsi_mesh_;
};  // Flecsi_State_Wrapper

}  // namespace Portage

#endif // FLECSI_STATE_WRAPPER_H_
