/*---------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SIMPLE_STATE_WRAPPER_H_
#define SIMPLE_STATE_WRAPPER_H_

#include <string>

#include "portage/simple_mesh/simple_state.h"

#include "portage/support/portage.h"

/*!
  @file simple_state_wrapper.h
  @brief Definitions for a wrapper to Simple_State
*/

namespace Portage {

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
    @brief Get a pointer to data from the state manager with a given
    variable @c name and on @c on_what mesh entities.
    @param[in] on_what The Entity_kind (e.g. CELL) on which the data
    lives.
    @param[in] name The name of the variable.
    @param data A pointer to the data array.  If the requested data is not
    found in the state manager, a @c nullptr is returned.
   */
  void get_data(const Entity_kind on_what, const std::string name,
                  double **data) const {
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
    @brief Get the number of elements in a specific variable from the state
    manager.
    @param[in] on_what The Entity_kind (e.g. CELL) of the variable for which the
    size is requested.
    @param[in] name The name of the variable for which the size is requested.
    @returns The number of elements; for example, the number of cells, if the
    data were added with the CELLS Entity_kind.
   */
  int get_data_size(const Entity_kind on_what, const std::string name) const {
    auto it = state_.find(name, on_what);
    if (it != state_.end()) {
      return it->second.size();
    }

    std::cerr << "get_data_size: Could not find state variable " << name
              << std::endl;
    return 0;
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
}  // namespace Portage

#endif  // SIMPLE_STATE_WRAPPER_H_
