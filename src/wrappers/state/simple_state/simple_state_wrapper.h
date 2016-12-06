/*---------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SIMPLE_STATE_WRAPPER_H_
#define SIMPLE_STATE_WRAPPER_H_

#include <string>

#include "portage/simple_mesh/simple_state.h"

#include "portage/support/portage.h"

namespace Portage {
class Simple_State_Wrapper {
 public:
  explicit Simple_State_Wrapper(Simple_State & state) : state_(state) { }
  Simple_State_Wrapper & operator=(Simple_State_Wrapper const &) = delete;
  ~Simple_State_Wrapper() { }

  //  void init_from_mesh() { state_.init_from_mesh(); }
  //  void export_to_mesh() { state_.export_to_mesh(); }

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

  int get_data_size(const Entity_kind on_what, const std::string name) const {
    auto it = state_.find(name, on_what);
    if (it != state_.end()) {
      return it->second.size();
    }

    std::cerr << "get_data_size: Could not find state variable " << name
              << std::endl;
    return 0;
  }

  Simple_State::name_vec_it names_begin() const {
    return state_.names_begin();
  }
  Simple_State::name_vec_it names_end() const {
    return state_.names_end();
  }

 private:
  Simple_State & state_;
};  // class Simple_State_Wrapper
}  // namespace Portage

#endif  // SIMPLE_STATE_WRAPPER_H_
