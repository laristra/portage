/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef SWARM_STATE_H_INC_
#define SWARM_STATE_H_INC_

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <cassert>

#include "portage/swarm/swarm.h"
#include "portage/support/portage.h"

namespace Portage { namespace Meshfree {

/**
 * @class SwarmState
 *
 * @brief Particle field state class.
 *
 * @tparam dim: the dimension of the problem.
 */
template<int dim>
class SwarmState {
public:
  /**
   * @brief Create an empty state.
   *
   */
  SwarmState() = default;

  /**
   * @brief Initialize from a reference swarm.
   *
   * @param swarm: the swarm with which the field data are associated.
   */
  explicit SwarmState(Swarm<dim> const& swarm)
    : num_local_points_(swarm.num_owned_particles())
  {}

  /**
   * @brief Initialize with a field size.
   *
   * @param size: size of each field.
   */
  explicit SwarmState(int size) : num_local_points_(size) {}

  /**
   * @brief Create the swarm state from a given mesh state wrapper.
   *
   * Copies fields from mesh state wrapper to the current swarm state
   * of the same size. It works only for double precision floating-point fields.
   *
   * @tparam State: type of the mesh state wrapper.
   * @param state: the mesh state wrapper.
   * @param kind: the entity kind to consider (node, cell).
   */
  template<typename State>
  SwarmState(State const& state, Wonton::Entity_kind kind) {
    // create return value
    int num_entities = 0;
    auto const field_names = state.names();

    // retrieve number of entities
    for (auto&& name : field_names) {
      if (state.get_entity(name) == kind) {
        num_entities = state.get_data_size(kind, name);
        break;
      }
    }
    // set number of particles
    num_local_points_ = num_entities;

    // copy data
    for (auto&& name : field_names) {
      if (state.get_entity(name) == kind) {
        // retrieve field from mesh state
        const double* values;
        state.mesh_get_data(kind, name, &values);
        assert(values != nullptr);

        // perform a deep copy then
        auto& field = fields_dbl_[name];
        field.resize(num_entities);
        std::copy(values, values + num_entities, field.begin());
      }
    }
  }

  /**
   * @brief Create the swarm state from set of mesh states wrappers.
   *
   * Copies fields from mesh state wrappers to the current swarm state of the
   * same total size. It is useful for mapping from several meshes at once,
   * e.g. for analysis of multiple times, or multiple simulations at the same
   * set of spatial points, to obtain statistical properties like average,
   * confidence bounds, and uncertaintities, for example. It works only for
   * double precision floating-point fields.
   *
   * @tparam State: type of the mesh state wrapper.
   * @param state: a list of mesh states wrappers.
   * @param kind: the entity kind to consider (node, cell).
   */
  template<typename State>
  SwarmState(std::vector<State*> const& states, Wonton::Entity_kind kind) {

    auto const& names = states[0]->names();
    int const num_fields = names.size();
    int const num_states = states.size();

    // check all fields in all states
    for (auto&& state : states) {
      auto current = state->names();
      if (current.size() == unsigned(num_fields)) {
        for (int i = 0; i < num_fields; ++i)
          if (names[i] != current[i])
            throw std::runtime_error("field names do not match");
      } else
        throw std::runtime_error("field names do not match");
    }

    // get sizes of data fields that match entity on each wrapper
    std::vector<std::vector<int>> sizes(num_fields);

    for (int i = 0; i < num_fields; ++i) {
      auto const& name = names[i];
      for (int j = 0; j < num_states; ++j) {
        auto const& state = *(states[j]);
        if (state.get_entity(name) == kind) {
          sizes[i].emplace_back(state.get_data_size(kind, name));
        }
      }
    }

    // ensure sizes are same across names for each wrapper
    for (int i = 0; i < num_fields; ++i) {
      for (int j = 0; j < num_states; ++j) {
        assert(sizes[0][j] > 0);
        if (sizes[i][j] != sizes[0][j])
          throw std::runtime_error("field sizes do not match");
      }
    }

    // compute particle offsets per state
    int num_entities = 0;
    int offset[num_states];
    for (int i = 0; i < num_states; i++) {
      offset[i] = num_entities;
      num_entities += sizes[0][i];
    }

    num_local_points_ = num_entities;

    // copy data
    for (auto&& name : names) {
      // resize particle field array
      auto& field = fields_dbl_[name];
      field.resize(num_entities);

      for (int i = 0; i < num_states; ++i) {
        // retrieve field from mesh state
        double* values = nullptr;
        states[i]->mesh_get_data(kind, name, &values);
        assert(values != nullptr);
        // perform a deep copy then
        std::copy(values, values + sizes[0][i], field.begin() + offset[i]);
      }
    }
  }

  /**
   * @brief Destructor.
   *
   */
  ~SwarmState() = default;

  /**
   * @brief Set a field on the swarm.
   *
   * @tparam T: field values type (int, double)
   * @param name: field name.
   * @param value: field values list.
   */
  template<typename T = double>
  void add_field(std::string name, Portage::vector<T> const&  value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");
    // sizes should match
    assert(value.size() == unsigned(num_local_points_));

    if (std::is_integral<T>::value) {
      auto& field = fields_int_[name];
      field.resize(value.size());
      std::copy(value.begin(), value.end(), field.begin());
    } else {
      auto& field = fields_dbl_[name];
      field.resize(value.size());
      std::copy(value.begin(), value.end(), field.begin());
    }
  }

#ifdef PORTAGE_ENABLE_THRUST
  /**
   * @brief Set a field on the swarm.
   *
   * @tparam T: field values type (int, double)
   * @param name: field name.
   * @param value: field values list.
   */
  template<typename T = double>
  void add_field(std::string name, std::vector<T> const&  value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");
    // sizes should match
    assert(value.size() == unsigned(num_local_points_));

    if (std::is_integral<T>::value) {
      auto& field = fields_int_[name];
      field.resize(value.size());
      std::copy(value.begin(), value.end(), field.begin());
    } else {
      auto& field = fields_dbl_[name];
      field.resize(value.size());
      std::copy(value.begin(), value.end(), field.begin());
    }
  }
#endif

  /**
   * @brief Set a field on the swarm.
   *
   * @tparam T: field values type (int, double)
   * @param name: field name.
   * @param value: field values array.
   */
  template<typename T = double>
  void add_field(std::string name, const T* const value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");
    // assume value is of the right length - we can't check it
    assert(value != nullptr);
    // do a deep copy on related array
    if (std::is_integral<T>::value) {
      auto& field = fields_int_[name];
      field.resize(num_local_points_);
      std::copy(value, value + num_local_points_, field.begin());
    } else {
      auto& field = fields_dbl_[name];
      field.resize(num_local_points_);
      std::copy(value, value + num_local_points_, field.begin());
    }
  }

  /**
   * @brief Set an empty field on the swarm.
   *
   * @tparam T: field values type (int, double)
   * @param name: field name.
   * @param value: default field elements value.
   */
  template<typename T=double>
  void add_field(std::string name, T value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");

    if (std::is_integral<T>::value) {
      auto& field = fields_int_[name];
      field.resize(num_local_points_, value);
    } else {
      auto& field = fields_dbl_[name];
      field.resize(num_local_points_, value);
    }
  }

  /**
   * @brief Retrieve a specified integer field.
   *
   * One cannot specialize template method without
   * specializing the template class, so just revert
   * to plain untemplated methods.
   *
   * @param name the name of the integer field
   * @return value the values in the field
   */
  Portage::vector<int>& get_field_int(std::string name) const {
    assert(fields_int_.count(name));
    using T = Portage::vector<int>;
    return const_cast<T&>(fields_int_.at(name));
  }

  /**
   * @brief Retrieve a specified real field.
   *
   * One cannot specialize template method without
   * specializing the template class, so just revert
   * to plain untemplated methods.
   *
   * @param name the name of the integer field
   * @return value the values in the field
   */
  Portage::vector<double>& get_field_dbl(std::string name) const {
    assert(fields_dbl_.count(name));
    using T = Portage::vector<double>;
    return const_cast<T&>(fields_dbl_.at(name));
  }

  /**
   *
   * @param name
   * @return
   */
  Portage::vector<double>& get_field(std::string name) const {
    return get_field_dbl(name);
  }

  /**
   * @brief Retrieve the specified field.
   *
   * @param name field name.
   * @param value a pointer to field values.
   */
  template<typename T = double>
  void copy_field(std::string name, T* value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");
    assert(value != nullptr);

    if (std::is_integral<T>::value) {
      assert(fields_int_.count(name));
      auto& field = fields_int_[name];
      std::copy(field.begin(), field.end(), value);
    } else {
      assert(fields_dbl_.count(name));
      auto& field = fields_dbl_[name];
      std::copy(field.begin(), field.end(), value);
    }
  }

  /**
   * @brief Get number of particles.
   *
   * @return the number of particles of the swarm.
   */
  int get_size() { return num_local_points_; }

  /**
   * @brief Retrieve the list of field names.
   *
   * @tparam T: field values type (int or double)
   * @return the list of field names
   */
  template<typename T = double>
  std::vector<std::string> get_field_names() {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");

    std::vector<std::string> list;
    if (std::is_integral<T>::value) {
      for (auto&& field : fields_int_)
        list.emplace_back(field.first);
    } else {
      for (auto&& field : fields_dbl_)
        list.emplace_back(field.first);
    }

    return list;
  }

  /**
   * @brief Extend the field with the given subfield.
   *
   * @tparam T: field values type (int or double).
   * @param name: name of the field to update.
   * @param values: the subfield to add.
   */
  template<typename T>
  void extend_field(std::string name, Portage::vector<T> const& values) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");

    if (std::is_integral<T>::value) {
      assert(fields_int_.count(name));
      auto& field = fields_int_[name];
      field.insert(field.begin(), values.begin(), values.end());
    } else {
      assert(fields_dbl_.count(name));
      auto& field = fields_dbl_[name];
      field.insert(field.end(), values.begin(), values.end());
    }
  }

 private:
  /** owned particles count */
  int num_local_points_ = 0;

  /** data fields */
  std::map<std::string, Portage::vector<int>>    fields_int_ {};
  std::map<std::string, Portage::vector<double>> fields_dbl_ {};
};

}} //namespace Portage::MeshFree

#endif // SWARM_STATE_H_INC_