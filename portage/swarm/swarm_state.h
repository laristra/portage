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

namespace Portage {
namespace Meshfree {

using std::string;
using std::shared_ptr;
using std::make_shared;
using std::map;
using std::pair;

/*!
 @class SwarmState "swarm_state.h"
 @brief Holds state data for particle Swarms.
 */
template<size_t dim>
class SwarmState {
 public:
  /*! @brief Integer data type allowed on the swarm.  */
  using IntVec=vector<int>;

  /*! @brief Double data type allowed on the swarm.  */
  using DblVec=vector<double>;

  /*! @brief Pointer to integer data type allowed on the swarm.  */
  using IntVecPtr=shared_ptr<vector<int>>;

  /*! @brief Pointer to double data type allowed on the swarm.  */
  using DblVecPtr=shared_ptr<vector<double>>;

  /*! @brief Constructor provided a reference swarm.
   * @param swarm the swarm with which the field data are associated.
   */
 SwarmState(Swarm<dim>& swarmin): npoints_owned_(swarmin.num_owned_particles())
    {}

  /*! @brief Constructor provided a size.
   * @param data_size the number of data elements
   */
 SwarmState(size_t data_size): npoints_owned_(data_size)
    {}

  /*! @brief Set an integer field on the swarm.
   * @param name the name of the integer field
   * @param value the values in the field
   * If field does not exist, create it.
   */
  void add_field(const string name, IntVecPtr value);

  /*! @brief Set a double field on the swarm centers.
   * @param name the name of the double field
   * @param value the values in the field
   * If field does not exist, create it.
   */
  void add_field(const string name, DblVecPtr value);

  /*! @brief Get an integer field off the swarm - throw exception if
   *  it does not exist
   *
   * @param name the name of the integer field
   * @param value the values in the field
   */
  void get_field(const string name, IntVecPtr &value) const;

  /*! @brief Get a double field off the swarm centers - throw
   *  exception if it does not exist
   *
   * @param name the name of the double field
   * @param value the values in the field
   */
   void get_field(const string name, DblVecPtr &value) const;

  /*! @brief Get an integer field off the swarm - add it if it does
   *  not exist
   *
   * @param name the name of the integer field
   * @param value the values in the field
   */
  void get_field(const string name, IntVecPtr &value);

  /*! @brief Get a double field off the swarm centers - add it if it
   *  does not exist
   *
   * @param name the name of the double field
   * @param value the values in the field
   */
  void get_field(const string name, DblVecPtr &value);

  /*! @brief Get number of points in swarm
   * @return number of points
   */
  int get_size(){return npoints_owned_;}

  /*! @brief Get the names of all integer fields
   */
  std::vector<std::string> field_names_int() {
    std::vector<std::string> result;
    for (auto iter=int_field_map_.begin(); iter!=int_field_map_.end(); iter++) {
      result.push_back(iter->first);
    }
    return result;
  }

  /*! @brief Get the names of all double fields
   */
  std::vector<std::string> field_names_double() {
    std::vector<std::string> result;
    for (auto iter=dbl_field_map_.begin(); iter!=dbl_field_map_.end(); iter++) {
      result.push_back(iter->first);
    }
    return result;
  }

  void extend_field(const string name, IntVec new_value);
  void extend_field(const string name, DblVec new_value);

 private:
  /** number of owned particles */
  int npoints_owned_;

  /** integer data fields */
  map<string, IntVecPtr> int_field_map_;

  /** double data fields */
  map<string, DblVecPtr> dbl_field_map_;
};

//=======================================================================

template<size_t dim>
void SwarmState<dim>::add_field(const string name, IntVecPtr value) {
  // check size
  if (value->size() != npoints_owned_) {
    throw std::runtime_error(
      string("incorrect size when adding attempting to add int field ")+name);
  }

  // check duplicate
  auto checkdup = int_field_map_.find(name);
  if (checkdup != int_field_map_.end()) {
    assert(checkdup->first == name);
    throw std::runtime_error(string("tried to add int field ")+name+
                             "when it already existed");
  }

  // add it
  int_field_map_.insert(pair<string, IntVecPtr>(name, value));
}

template<size_t dim>
void SwarmState<dim>::add_field(const string name, DblVecPtr value) {
  // check size
  if (value->size() != npoints_owned_) {
    throw std::runtime_error(
      string("incorrect size when adding attempting to add double field ")+name);
  }

  // check duplicate
  auto checkdup = dbl_field_map_.find(name);
  if (checkdup != dbl_field_map_.end()) {
    assert(checkdup->first == name);
    throw std::runtime_error(string("tried to add double field ")+name+
                             " when it already existed");
  }

  // add it
  dbl_field_map_.insert(pair<string, DblVecPtr>(name, value));
}

// Const version of get_field for integer field - throws exception if
// field does not exist
template<size_t dim>
void SwarmState<dim>::get_field(const string name, IntVecPtr &value) const {
  value = int_field_map_.at(name);
}

// Const version of get_field for real field - throws exception if
// field does not exist
template<size_t dim>
void SwarmState<dim>::get_field(const string name, DblVecPtr &value) const {
  value = dbl_field_map_.at(name);
}

// Non-const version of get_field for integer field - inserts the field
// if it does not exist
template<size_t dim>
void SwarmState<dim>::get_field(const string name, IntVecPtr &value) {
  value = int_field_map_[name];
}

// Non-const version of get_field for real field - inserts the field
// if it does not exist
template<size_t dim>
void SwarmState<dim>::get_field(const string name, DblVecPtr &value) {
  value = dbl_field_map_[name];
}


template<size_t dim>
void SwarmState<dim>::extend_field(const string name, IntVec new_value)  
{
   // check if the field already exists
  auto check = int_field_map_.find(name);
  if (check == int_field_map_.end()) {
    throw std::runtime_error(string("tried to extend an int field that does not exist ")+name);
  }

  IntVecPtr val = int_field_map_.at(name);
  val->insert(val->end(), new_value.begin(), new_value.end());
}

template<size_t dim>
void SwarmState<dim>::extend_field(const string name, DblVec new_value)  
{
   // check if the field already exists
  auto check = dbl_field_map_.find(name);
  if (check == dbl_field_map_.end()) {
    throw std::runtime_error(string("tried to extend a double field that does not exist ")+name);
  }

  DblVecPtr val = dbl_field_map_.at(name);
  val->insert(val->end(), new_value.begin(), new_value.end());
}

/*! @brief SwarmState factory, given a mesh state wrapper.
 * Copies fields from mesh state wrapper to a swarm state wrapper of the 
 * same size.
 * @param mesh the mesh with which the field data are associated.
 * @param entity entity on which to get data (e.g. CELL, NODE, etc.)
 * @param state the field data on the mesh
 */
template<size_t dim, class StateWrapper>
shared_ptr<SwarmState<dim>> SwarmStateFactory(
  const StateWrapper &state,
  const Portage::Entity_kind entity)
{
  // create return value
  size_t ndata=0;

  for (const std::string name : state.names()) {
    // Simple_State does not store separte lists of names by entity, 
    // so we have to filter.
    if (state.get_entity(name) == entity) {
      ndata = state.get_data_size(entity, name);
      break;
    }
  }
  shared_ptr<SwarmState<dim>> result=make_shared<SwarmState<dim>>(ndata);

  // copy data
  for (const std::string name : state.names()) {
    if (state.get_entity(name) != entity) continue;

    // make sure all fields have same size
    assert(state.get_data_size(entity, name) == ndata);

    const double *datap;
    state.mesh_get_data(entity, name, &datap);
    typename SwarmState<dim>::DblVecPtr data = make_shared<vector<double>>(ndata);
    
    for (size_t i=0; i<ndata; i++) (*data)[i] = datap[i];
    result->add_field(name, data);
  }

  return result;
}

} //namespace MeshFree
} //namespace Portage

#endif // SWARM_STATE_H_INC_

