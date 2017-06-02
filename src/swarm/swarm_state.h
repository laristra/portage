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

  /*! @brief Integer data type allowed on the swarm.  */
  using IntVecPtr=shared_ptr<vector<int>>;

  /*! @brief Double data type allowed on the swarm.  */
  using DblVecPtr=shared_ptr<vector<double>>;

  /*! @brief Constructor provides a reference swarm.
   * @param swarm the swarm with which the field data are associated.
   */
  SwarmState(Swarm<dim> const& swarmin): swarm_(swarmin){}

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

  /*! @brief Get an integer field off the swarm.
   * @param name the name of the integer field
   * @param value the values in the field
   */
  void get_field(const string name, IntVecPtr &value);

  /*! @brief Get a double field off the swarm centers.
   * @param name the name of the double field
   * @param value the values in the field
   */
  void get_field(const string name, DblVecPtr &value);

  /*! @brief Get number of points in swarm
   * @return number of points
   */
  size_t get_size(){return swarm_.num_owned_particles();}

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
void SwarmState<dim>::add_field(const string name, IntVecPtr value) {
  // check size
  assert(value->size() == swarm_.num_owned_particles());

  // check duplicate
  /* couldn't get this to work
  auto checkdup = int_field_map_.find(name);
  if (checkdup != int_field_map_.end()) {
    IntVecPtr fieldp = checkdup->second;
    IntVec field(*fieldp);
    throw std::runtime_error(string("tried to add int field ")+name+
                             "when it already existed");
  }
  */

  // add it
  int_field_map_.insert(pair<string, IntVecPtr>(name, value));
}

template<size_t dim>
void SwarmState<dim>::add_field(const string name, DblVecPtr value) {
  // check size
  assert(value->size() == swarm_.num_owned_particles());

  // check duplicate
  /* couldn't get this to work
  auto checkdup = dbl_field_map_.find(name);
  if (checkdup != dbl_field_map_.end()) {
    DblVecPtr fieldp = checkdup->second;
    DblVec field(*fieldp);
    throw std::runtime_error(string("tried to add double field ")+name+
                             " when it already existed");
  }
  */

  // add it
  dbl_field_map_.insert(pair<string, DblVecPtr>(name, value));
}

template<size_t dim>
void SwarmState<dim>::get_field(const string name, IntVecPtr &value) {
  value = int_field_map_[name];
}

template<size_t dim>
void SwarmState<dim>::get_field(const string name, DblVecPtr &value) {
  value = dbl_field_map_[name];
}

}
}

#endif // SWARM_STATE_H_INC_

