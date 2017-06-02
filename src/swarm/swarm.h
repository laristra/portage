/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SWARM_H_INC_
#define SWARM_H_INC_

#include <vector>
#include <memory>
#include <cassert>
#include <cmath>

#include "portage/support/Point.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;
using std::shared_ptr;

/*!
 @class Swarm "swarm.h"
 @brief An effective "mesh" class for a collection disconnected points (particles).
 @tparam dim The dimension of the problem space.
 */
template<size_t dim> class Swarm {
 public:

  // When using PointVec and PointVecPtr in another file, they have to
  // be prefixed by the keyword typename
  
  //  SEE --- https://stackoverflow.com/questions/610245/where-and-why-do-i-have-to-put-the-template-and-typename-keywords

  using PointVecPtr =  shared_ptr<std::vector<Point<dim>>>;
  using PointVec = vector<Point<dim>>;
  
  
  /*!
   * @brief A particle has a center point and smoothing lengths in each dimension.
   * @param points center points of the particles
   * @param extents the widths of the particle bounding boxes in each dimension
   */
  Swarm(PointVecPtr points)
      : points_(points), npoints_(points_->size()) 
  {}

  /*! @brief return the nubmer of particles in the swarm.
   * @return the number of particles in the swarm
   */
  int num_owned_particles() const {
    return npoints_;
  }

  /*! @brief Get the coordinates of the particle,
   * @param index index of particle of interest
   * @return the particle coordinates
   */
  Point<dim> get_particle_coordinates(const size_t index) const {
    return (*points_)[index];
  }

  

 private:
  /** the centers of the particles */
  PointVecPtr points_;

  /** the number of particles in the swarm */
  const size_t npoints_;
};

}
}

#endif // SWARM_H_INC_

