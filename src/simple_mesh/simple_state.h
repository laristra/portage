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

/*!
  @file simple_state.h
  @brief A very light-weight state manager for @c double data living atop a
  @c Simple_Mesh.
 */

namespace Portage {

  /*!
    @class Simple_State "simple_mesh.h"
    @brief A very light-weight state manager for a Simple_Mesh.

    This class stores @c double (and only @c double) data atop a Simple_Mesh.
    The data can be located at any mesh location (e.g. @c CELL or @c NODE) that
    Simple_Mesh understands.  The field variables are stored in a map and are
    identified by key composed of a string _name_ and an Entity_kind; for
    example, ```pair("pressure", CELL)```.  The actual data corresponding to
    any field is stored as a std::vector of @c double's.

    In addition to the map, this class stores (and exposes) a vector containing
    just the string variable names.
   */
class Simple_State {
 public:
  /*!
    @brief Constructor for creating a Simple_State
    @param[in] mesh A smart pointer to the Simple_Mesh where this data lives.
   */
  explicit Simple_State(std::shared_ptr<Simple_Mesh> mesh) : mesh_(mesh) { }

  /// Destructor
  ~Simple_State() { }

  /// The number of variables that live in the state manager.
  int numVars() const { return names_.size(); }

  /// Convenience types
  typedef std::vector<double> vec;
  typedef vec::iterator vec_it;
  typedef std::vector<std::string> name_vec;
  typedef name_vec::iterator name_vec_it;
  typedef std::pair<std::string, Entity_kind> key;
  typedef std::map<key, vec> mymap;
  typedef mymap::iterator map_it;

  /// An iterator to the beginning of the field map.
  map_it begin() { return state_vectors_.begin(); }
  /// An iterator to the ending of the field map.
  map_it end() { return state_vectors_.end(); }

  /// An iterator to the beginning of the variable names vector.
  name_vec_it names_begin() { return names_.begin(); }
  /// An iterator to the ending fo the variable names vector.
  name_vec_it names_end() { names_.end(); }

  /*!
    @brief Search for a specific key (name, Entity_kind) within the
    map of known fields and return an iterator to its location.
    @param[in] name The variable name we wish to find.
    @param[in] on_what The Entity_kind (e.g. CELL) where the data we wish to
    find lives.
    @returns An iterator to the location of the requested key in the map; if
    this is the same as Simple_State::end(), then the variable was not found in
    the map.
   */
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

  /*!
    @brief Add a field to the state manager.
    @param[in] name The variable name of the field to be added.
    @param[in] on_what The Entity_kind (e.g. CELL) on which the data should
    live.
    @param[in] data A pointer to some data to initialize.  This defaults to
    @c nullptr, which allows you to not initialize any of the data, but the
    underlying std::vector will be sized appropriately.
    @returns A reference to the underlying std::vector for further
    modification/use.

    If the requested variable (@c name, @c on_what) pair already exists,
    this will output a message to std::cerr, but will simply return the
    reference to the underlying std::vector.
   */
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

  /*!
    @brief Add a field to the state manager.
    @param[in] name The variable name of the field to be added.
    @param[in] on_what The Entity_kind (e.g. CELL) on which the data should
    live.
    @param[in] val A value with which to iniitalize the underlying std::vector
    data.
    @returns A reference to the underlying std::vector for further
    modification/use.

    If the requested variable (@c name, @c on_what) pair already exists,
    this will output a message to std::cerr, but will simply return the
    reference to the underlying std::vector.
  */
  vec& add(std::string const name, Entity_kind const on_what,
           double const val) {
    auto it = find(name, on_what);
    if (it == end()) {
      // This is a new variable, so store it.
      auto num_ent = mesh_->num_entities(on_what, Entity_type::PARALLEL_OWNED);
      key thisKey(name, on_what);
      state_vectors_.insert(std::make_pair(thisKey, vec(num_ent, val)));
      names_.emplace_back(name);
      return state_vectors_[thisKey];
    } else {
      // This already exists!
      std::cerr << "Attempted to add duplicate vectors.  Ignoring."
                << std::endl;
      return it->second;
    }
  }

  /*!
    @brief Get a specific variable from the state manager.
    @param[in] name The variable name of the field to be retrieved.
    @param[in] on_what The Entity_kind (e.g. CELL) on whih the requested data
    lives.
    @returns A reference to the underlying std::vector of data for further
    modification/use.
    @throws std::runtime_error The requested variable is not in the state
    manager.
   */
  vec& get(std::string const name, Entity_kind const on_what) {
    auto it = find(name, on_what);
    if (it != end()) {
      return it->second;
    } else {
      throw std::runtime_error("Couldn't find the requested variable.");
    }
  }

 private:
  /// The mesh on which this state manager operates.
  std::shared_ptr<Simple_Mesh> mesh_;
  /// The map of key (name, Entity_kind) and data values.
  mymap state_vectors_;
  /// Vector of variable names.
  name_vec names_;
};  // Simple_State

}  // namespace Portage

#endif  // SRC_SIMPLE_MESH_SIMPLE_STATE_H_
