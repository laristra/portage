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
#include <map>
#include <stdexcept>

#include "portage/simple_mesh/simple_mesh.h"

#include "portage/support/portage.h"


namespace Portage {

class Simple_State {
 public:
  explicit Simple_State(std::shared_ptr<Simple_Mesh> mesh) : mesh_(mesh) { }
  ~Simple_State() { }

  int numVars() const { return names_.size(); }

  typedef std::vector<double> vec;
  typedef vec::iterator vec_it;
  typedef std::vector<std::string> name_vec;
  typedef name_vec::iterator name_vec_it;
  typedef std::pair<std::string, Entity_kind> key;
  typedef std::map<key, vec> mymap;
  typedef mymap::iterator map_it;

  map_it begin() { return state_vectors_.begin(); }
  map_it end() { return state_vectors_.end(); }

  name_vec_it names_begin() { return names_.begin(); }
  name_vec_it names_end() { names_.end(); }

  map_it find(std::string const name,
              Entity_kind const on_what = Entity_kind::ANY_KIND) {
    key testKey(name, on_what);
    auto it = begin();
    while (it != end()) {
      if ((it->first == testKey) ||
          ((it->first.first == name) &&
           (on_what == Entity_kind::ANY_KIND)))
        break;
      else
        ++it;
    }
    return it;
  }

  vec& add(std::string const name, Entity_kind const on_what,
           double const * const data = nullptr) {
    auto it = find(name, on_what);
    if (it == end()) {
      // This is a new variable, so store it.
      auto num_ent = mesh_->num_entities(on_what, Entity_type::PARALLEL_OWNED);
      key thisKey(name, on_what);
      if (data == nullptr) {
        state_vectors_.insert(std::make_pair(thisKey,
                                             vec(num_ent)));
      } else {
        state_vectors_.insert(std::make_pair(thisKey,
                                             vec(data, data+num_ent)));
      }
      names_.emplace_back(name);
      return state_vectors_[thisKey];
    } else {
      // This already exists!
      std::cerr << "Attempted to add duplicate vectors.  Ignoring."
                << std::endl;
      return it->second;
    }
  }

  vec& get(std::string const name, Entity_kind const on_what) {
    auto it = find(name, on_what);
    if (it != end()) {
      return it->second;
    } else {
      throw std::runtime_error("Couldn't find the requested variable.");
    }
  }

 private:
  std::shared_ptr<Simple_Mesh> mesh_;
  mymap state_vectors_;
  name_vec names_;
};  // Simple_State

}  // namespace Portage

#endif  // SRC_SIMPLE_MESH_SIMPLE_STATE_H_
