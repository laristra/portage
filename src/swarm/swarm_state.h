/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SWARM_STATE_H_INC_
#define SWARM_STATE_H_INC_

#include <vector>
#include <map>
#include <memory>
#include <string>
#include <cassert>

#include "portage/wrappers/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wrappers/state/flat/flat_state_wrapper.h"
#include "swarm.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;
using std::shared_ptr;
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
  SwarmState(Swarm<dim> const& swarmin): swarm_(swarmin){}

  /*! @brief Constructor provided a flat mesh and flat mesh state.
   * @param mesh the mesh with which the field data are associated.
   * @param entity entity on which to get data (e.g. CELL, NODE, etc.)
   * @param state the field data on the mesh
   */
  SwarmState(Portage::Flat_Mesh_Wrapper<double> &mesh,
	     Portage::Entity_kind entity, 
	     Portage::Flat_State_Wrapper<double> &state);

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
  int get_size(){return swarm_.num_owned_particles();}

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

 private:
  /** reference swarm state */
  Swarm<dim> const& swarm_;

  /** integer data fields */
  map<string, IntVecPtr> int_field_map_;

  /** double data fields */
  map<string, DblVecPtr> dbl_field_map_;
};

//=======================================================================

template<size_t dim>
SwarmState<dim>::SwarmState(Portage::Flat_Mesh_Wrapper<double> &mesh,
                                 Portage::Entity_kind entity,
				 Portage::Flat_State_Wrapper<double> &state)
  : swarm_(Swarm<dim>(mesh, entity))
{
  if (dim != mesh.space_dimension()) {
    throw std::runtime_error(string("dimension mismatch"));
  }

  std::vector<std::string> dnames;
  state.get_names(entity, dnames);

  for (auto iter=dnames.begin(); iter!=dnames.end(); iter++) {
    double *datap;
    state.get_data(entity, *iter, &datap);
    int npart = swarm_.num_owned_particles();
    DblVecPtr data(new vector<double>(npart));
    for (size_t i=0; i<npart; i++) (*data)[i] = datap[i];
    add_field(*iter, data);
  }
}

template<size_t dim>
void SwarmState<dim>::add_field(const string name, IntVecPtr value) {
  // check size
  if (value->size() != swarm_.num_owned_particles()) {
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
  if (value->size() != swarm_.num_owned_particles()) {
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

// Non-const version of get_filed for integer field - inserts the field
// if it does not exist
template<size_t dim>
void SwarmState<dim>::get_field(const string name, IntVecPtr &value) {
  value = int_field_map_[name];
}

// Non-const version of get_filed for real field - inserts the field
// if it does not exist
template<size_t dim>
void SwarmState<dim>::get_field(const string name, DblVecPtr &value) {
  value = dbl_field_map_[name];
}

}
}

#endif // SWARM_STATE_H_INC_

