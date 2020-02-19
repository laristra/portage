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

/*!
 @class SwarmState "swarm_state.h"
 @brief Holds state data for particle Swarms.
 */
template<int dim>
class SwarmState {
public:
  /*! @brief Integer data type allowed on the swarm.  */
  using IntVec = Portage::vector<int>;

  /*! @brief Double data type allowed on the swarm.  */
  using DblVec = Portage::vector<double>;

  /*! @brief Pointer to integer data type allowed on the swarm.  */
  using IntVecPtr = std::shared_ptr<Portage::vector<int>>;

  /*! @brief Pointer to double data type allowed on the swarm.  */
  using DblVecPtr = std::shared_ptr<Portage::vector<double>>;

  /*! @brief Constructor provided a reference swarm.
   * @param swarm the swarm with which the field data are associated.
   */
  explicit SwarmState(Swarm<dim> const& swarm)
    : num_owned_points_(swarm.num_owned_particles())
  {}

  /*! @brief Constructor provided a size.
   * @param data_size the number of data elements
   */
  explicit SwarmState(int data_size)
    : num_owned_points_(data_size)
  {}

  /*! @brief Set an integer field on the swarm.
   * @param name the name of the integer field
   * @param value a shared pointer to the values in the field.
   * If field does not exist, create it.
   */
  template<typename T>
  void add_field(std::string name, Portage::vector<T> const& value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");
    // sizes should match
    assert (value.size() == num_owned_points_);

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

  /**
   * @brief Set an integer field on the swarm.
   *
   * @param name name of the integer field
   * @param value a shared pointer to the values in the field.
   * If field does not exist, create it.
   */
  template<typename T>
  void add_field(std::string name, std::vector<T> const&  value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");
    // sizes should match
    assert(value.size() == num_owned_points_);

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

  /**
  * @brief Set an integer field on the swarm.
  *
  * @param name name of the integer field
  * @param value a shared pointer to the values in the field.
  * If field does not exist, create it.
  */
  template<typename T>
  void add_field(std::string name, const T* const value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");
    // assume value is of the right length - we can't check it
    assert(value != nullptr);
    // do a deep copy on related array
    if (std::is_integral<T>::value) {
      auto& field = fields_int_[name];
      field.resize(num_owned_points_);
      std::copy(value, value + num_owned_points_, field.begin());
    } else {
      auto& field = fields_dbl_[name];
      field.resize(num_owned_points_);
      std::copy(value, value + num_owned_points_, field.begin());
    }
  }

  /**
   * @brief Retrieve an integer field.
   *
   * One cannot specialize template method without
   * specializing the template class, so just revert
   * to plain untemplated methods.
   *
   * @param name the name of the integer field
   * @return value the values in the field
   */
  Portage::vector<int> get_field_int(std::string name) const {
    assert(fields_int_.count("name"));
    return fields_int_.at(name);
  }

  /**
   * @brief Retrieve a real value field.
   *
   * One cannot specialize template method without
   * specializing the template class, so just revert
   * to plain untemplated methods.
   *
   * @param name the name of the integer field
   * @return value the values in the field
   */
  Portage::vector<double> get_field_double(std::string name) const {
    assert(fields_dbl_.count("name"));
    return fields_dbl_.at(name);
  }

  /**
   * @brief Get an integer field off the swarm.
   *
   * @param name the name of the integer field
   * @param value the values in the field
   */
  template<typename T>
  void get_field(std::string name, T const** value) const {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");

    if (std::is_integral<T>::value) {
      assert(fields_int_.count("name"));
      *value = fields_int_.at(name).data();
    }
    else {
      assert(fields_dbl_.count("name"));
      *value = fields_dbl_.at(name).data();
    }
  }

  /**
   * @brief Get an integer field off the swarm - add it if it does
   *  not exist
   *
   * @param name the name of the integer field
   * @param value the values in the field
   */
  template<typename T>
  void get_field(std::string name, T **value) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");

    if (std::is_integral<T>::value) {
      assert(fields_int_.count("name"));
      *value = fields_int_.at(name).data();
    } else {
      assert(fields_dbl_.count("name"));
      *value = fields_dbl_.at(name).data();
    }
  }

  /*! @brief Get number of particles
   * @return number of points
   */
  int get_size() { return num_owned_points_; }

  /*! @brief Get the names of all integer fields
   */
  template<typename T>
  std::vector<std::string> get_field_names() {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");

    std::vector<std::string> list;
    if (std::is_integral<T>::value) {
      assert(fields_int_.count("name"));
      for (auto&& field : fields_int_)
        list.emplace_back(field.first);
    } else {
      assert(fields_dbl_.count("name"));
      for (auto&& field : fields_dbl_)
        list.emplace_back(field.first);
    }

    return list;
  }

  /**
   *
   * @tparam T
   * @param name
   * @param values
   */
  template<typename T>
  void extend_field(std::string name, Portage::vector<T> const& values) {

    static_assert(std::is_arithmetic<T>::value, "only numeric fields");

    if (std::is_integral<T>::value) {
      assert(fields_int_.count(name));
      // resize and perform a deep copy
      int const old_size = fields_int_[name].size();
      int const new_size = old_size + values.size();
      fields_int_[name].resize(new_size);
      std::copy(values.begin(), values.end(), fields_int_[name].begin());
    } else {
      assert(fields_dbl_.count(name));
      // resize and perform a deep copy
      int const old_size = fields_dbl_[name].size();
      int const new_size = old_size + values.size();
      fields_dbl_[name].resize(new_size);
      std::copy(values.begin(), values.end(), fields_dbl_[name].begin());
    }
  }

 private:
  /** owned particles count */
  int num_owned_points_ = 0;

  /** data fields */
  std::map<std::string, Portage::vector<int>>    fields_int_;
  std::map<std::string, Portage::vector<double>> fields_dbl_;
};

/* -------------------------------------------------------------------------- */
/*! @brief SwarmState factory, given a mesh state wrapper. Only does double fields.
 * Copies fields from mesh state wrapper to a swarm state wrapper of the same size.
 * @param state the field data on the mesh
 * @param entity entity on which to get data (e.g. CELL, NODE, etc.)
 * @return shared pointer to the resultant swarm state
 */
template<size_t dim, class StateWrapper>
std::shared_ptr<SwarmState<dim>> SwarmStateFactory(
  const StateWrapper &state,
  const Portage::Entity_kind entity)
{
  // create return value
  size_t ndata=0;

  for (std::string name : state.names()) {
    // Simple_State does not store separte lists of names by entity, 
    // so we have to filter.
    if (state.get_entity(name) == entity) {
      ndata = state.get_data_size(entity, name);
      break;
    }
  }
  std::shared_ptr<SwarmState<dim>> result=std::make_shared<SwarmState<dim>>(ndata);

  // copy data
  for (std::string name : state.names()) {
    if (state.get_entity(name) != entity) continue;

    // make sure all fields have same size
    assert(state.get_data_size(entity, name) == ndata);

    const double *datap;
    state.mesh_get_data(entity, name, &datap);
    typename SwarmState<dim>::DblVecPtr data = std::make_shared<vector<double>>(ndata);
    
    for (size_t i=0; i<ndata; i++) (*data)[i] = datap[i];
    result->add_field(name, data);
  }

  return result;
}

/*! @brief SwarmState factory, given a set of mesh state wrappers. Only does double fields.
 * Copies fields from mesh state wrappers to a swarm state wrapper of the total size.
 * This factory is useful for mapping from several meshes at once, 
 * e.g. for analysis of multiple times, or multiple simulations at the same set of spatial points, 
 * to obtain statistical properties like average, confidence bounds, and uncertaintities, for example.
 * @param states the field data on the meshes
 * @param entity entity on which to get data (e.g. CELL, NODE, etc.)
 *
 * All fields must exist in all states. In each state, all fields must have same size.
 */
template<size_t dim, class StateWrapper>
std::shared_ptr<SwarmState<dim>> SwarmStateFactory
  (const std::vector<StateWrapper*> states,
   const Portage::Entity_kind entity)
{
  // check all fields in all states
  std::vector<std::string> names = states[0]->names();
  for (size_t wrap=1; wrap<states.size(); wrap++) {
    std::vector<std::string> these = states[wrap]->names();
    if (names != these) {
      throw std::runtime_error("field names don't match across state wrappers");
    }
  }

  // get sizes of data fields that match entity on each wrapper
  std::vector<std::vector<int>> sizes(names.size(), std::vector<int>(states.size(),0));
  for (size_t n=0; n<names.size(); n++) {
    for (size_t wrap=0; wrap<states.size(); wrap++) {
      StateWrapper &state=*states[wrap];
      std::string name=state.names()[n];
      if (state.get_entity(name) == entity) {
        sizes[n][wrap] = state.get_data_size(entity, name);
      }
    }
  }
  // ensure sizes are same across names for each wrapper
  for (size_t wrap=0; wrap<states.size(); wrap++) {
    assert(sizes[0][wrap]>0);
    for (size_t n=1; n<names.size(); n++) {
      if (sizes[n][wrap] != sizes[0][wrap]) {
	throw std::runtime_error("field sizes don't match across state names");
      }
    }
  }
  // get size of output swarm state
  size_t ndata=0;
  std::vector<int> offset(states.size(),0);
  for (size_t wrap=0; wrap<states.size(); wrap++) {
    offset[wrap] = ndata;
    ndata += sizes[0][wrap];
  }

  // create output swarm state
  std::shared_ptr<SwarmState<dim>> result = std::make_shared<SwarmState<dim>>(ndata);

  // copy data
  for (std::string name : states[0]->names()) {
    typename SwarmState<dim>::DblVecPtr data = std::make_shared<vector<double>>(ndata);
    for (size_t wrap=0; wrap<states.size(); wrap++) {
      StateWrapper &state=*states[wrap];

      const double *datap;
      state.mesh_get_data(entity, name, &datap);
    
      for (size_t i=0; i<sizes[0][wrap]; i++) (*data)[i+offset[wrap]] = datap[i];
    }
    result->add_field(name, data);
  }

  return result;
}

}} //namespace Portage::MeshFree

#endif // SWARM_STATE_H_INC_

