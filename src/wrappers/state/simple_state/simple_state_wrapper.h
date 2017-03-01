/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
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
  void get_data(const Entity_kind on_what, const std::string name,
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
    @brief Get a pointer to data from the state manager with a given
    variable @c name and on @c on_what mesh entities.
    @param[in] on_what The Entity_kind (e.g. CELL) on which the data
    lives.
    @param[in] name The name of the variable.
    @param data A @c pointer to the const data array.  If the requested data is
    not found in the state manager, a @c nullptr is returned.
   */
  void add_data(const Entity_kind on_what, const std::string name,
                  double const **data) const {
    auto vec = state_.add(name, on_what, *data);
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
  Entity_kind get_entity(std::string const name) const {
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
