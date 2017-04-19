/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SWARM_H_INC_
#define SWARM_H_INC_

#include <vector>
#include <map>
#include <tr1/shared_ptr.h>
#include <string>
#include <cassert>

#include "swarm.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;
using std::shared_ptr;
using std::map;

using PointVecPtr = shared_ptr<vector<Point<dim>>>;
using PointVec = vector<Point<dim>>;

/*!
 @class SwarmState "swarm_state.h"
 @brief Holds state data for particle Swarms.
 */
class SwarmState {
 public:

  /*! @brief Integer data type allowed on the swarm.  */
  using IntVecPtr=shared_ptr<vector<int>>;

  /*! @brief Double data type allowed on the swarm.  */
  using DblVecPtr=shared_ptr<vector<double>>;

  /*! @brief Constructor provides a reference swarm.
   * @param swarm the swarm with which the field data are associated.
   */
  SwarmState(shared_ptr<Swarm>)

  /*! @brief Set an integer field on the swarm.
   * @param name the name of the integer field
   * @param value the values in the field
   * If field does not exist, create it.
   */
  void add_field(string name, IntVecPtr value);

  /*! @brief Set a double field on the swarm centers.
   * @param name the name of the double field
   * @param value the values in the field
   * If field does not exist, create it.
   */
  void add_field(string name, DblVecPtr value);

  /*! @brief Get an integer field off the swarm.
   * @param name the name of the integer field
   * @param value the values in the field
   */
  void get_field(string name, IntVecPtr value);

  /*! @brief Get a double field off the swarm centers.
   * @param name the name of the double field
   * @param value the values in the field
   */
  void get_field(string name, DblVecPtr value);

 private:

  /** integer data fields */
  map<string, IntVecPtr> int_field_map_;

  /** double data fields */
  map<string, DblVecPtr> dbl_field_map_;
};

}
}

#endif // SWARM_H_INC_

