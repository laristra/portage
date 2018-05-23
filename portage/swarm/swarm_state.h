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

#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wonton/state/flat/flat_state_wrapper.h"
#include "portage/swarm/swarm.h"
#include "portage/support/portage.h"

namespace Portage {
namespace Meshfree {

using std::string;
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
 SwarmState(Swarm<dim>& swarmin): npoints_owned_(swarmin.num_owned_particles())
    {}

  /*! @brief Constructor provided a flat mesh and flat mesh state.
   * @param mesh the mesh with which the field data are associated.
   * @param entity entity on which to get data (e.g. CELL, NODE, etc.)
   * @param state the field data on the mesh
   */
  SwarmState(Wonton::Flat_Mesh_Wrapper<double> &mesh,
	     Portage::Entity_kind entity, 
	     Wonton::Flat_State_Wrapper<double> &state);

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
SwarmState<dim>::SwarmState(Wonton::Flat_Mesh_Wrapper<double> &mesh,
                                 Portage::Entity_kind entity,
				 Wonton::Flat_State_Wrapper<double> &state)
  : npoints_owned_(0)
{
  if (dim != mesh.space_dimension()) {
    throw std::runtime_error(string("dimension mismatch"));
  }

  if (entity==CELL) {
    npoints_owned_ = mesh.num_owned_cells();
  } else if (entity==NODE) {
    npoints_owned_ = mesh.num_owned_nodes();
  }

  assert(state.get_entity_size(entity) == npoints_owned_);

  std::vector<std::string> dnames;
  state.get_names(entity, dnames);

  for (auto iter=dnames.begin(); iter!=dnames.end(); iter++) {
    double *datap;
    state.get_data(entity, *iter, &datap);
    DblVecPtr data = make_shared<vector<double>>(npoints_owned_);
    for (size_t i=0; i<npoints_owned_; i++) (*data)[i] = datap[i];
    add_field(*iter, data);
  }
}

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

} //namespace MeshFree
} //namespace Portage

#endif // SWARM_STATE_H_INC_

