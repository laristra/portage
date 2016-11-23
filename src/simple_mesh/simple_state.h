/*----------------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------------*/

#ifndef SRC_SIMPLE_MESH_SIMPLE_STATE_H_
#define SRC_SIMPLE_MESH_SIMPLE_STATE_H_

#include <vector>
#include <memory>
#include <string>
#include <utility>

#include "portage/simple_mesh/simple_mesh.h"

#include "portage/support/portage.h"

namespace Portage {

class Simple_State {
 public:
  explicit Simple_State(const std::shared_ptr<Simple_Mesh> mesh) : mesh_(mesh)
  { }

  typedef std::pair<std::string, Entity_kind> field_ID;
  typedef std::vector<field_ID> field_ID_list;
  typedef field_ID_list::const_iterator field_ID_iterator;

  std::vector<double>& add(std::string const name,
                           Entity_kind const on_what,
                           double const * const data = nullptr) {
    auto it = find(name, on_what);
    if (it == known_fields_.end()) {
      // This wasn't found in our list; add it
      auto numEntities = mesh_->num_entities(on_what,
                                             Entity_type::PARALLEL_OWNED);
      if (data == nullptr)
        data_.emplace_back(numEntities);
      else
        data_.emplace_back(data, data + numEntities);
      known_fields_.emplace_back(name, on_what);
      return data_[data_.size()-1];
    } else {
      std::cerr << "Attempted to add duplicate state vector; ignoring..."
                << std::endl;
      return data_[it - known_fields_.begin()];
    }
  }

  int numVars() const { return known_fields_.size(); }

  double* get_data(std::string const name,
                   Entity_kind const on_what) {
    auto it = find(name, on_what);
    if (it == known_fields_.end()) {
      std::cerr << "Could not find state variable " << name << std::endl;
      return nullptr;
    } else {
      int idx = it - known_fields_.begin();
      return &(data_[idx].front());
    }
  }

 private:
  field_ID_list known_fields_;
  std::vector<std::vector<double>> data_;
  const std::shared_ptr<Simple_Mesh> mesh_;

  field_ID_iterator find(std::string const name,
                         Entity_kind const on_what) const {
    auto it = known_fields_.begin();
    while (it != known_fields_.end()) {
      if ((it->first == name) && (it->second == on_what))
        break;
      else
        ++it;
    }
    return it;
  }
};

}  // namespace Portage

#endif  // SRC_SIMPLE_MESH_SIMPLE_STATE_H_
