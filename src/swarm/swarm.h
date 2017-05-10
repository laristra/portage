/*---------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SWARM_H_INC_
#define SWARM_H_INC_

#include <vector>
#include <memory>
#include <cassert>

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

  using PointVecPtr = shared_ptr<vector<Point<dim>>>;
  using PointVec = vector<Point<dim>>;

  /*!
   * @brief A particle has a center point and smoothing lengths in each dimension.
   * @param points center points of the particles
   * @param extents the widths of the particle bounding boxes in each dimension
   */
  Swarm(PointVecPtr points, PointVecPtr extents)
      : points_(points),
    extents_(points), npoints_(points_->size()) 
  {
    assert(extents_->size() == npoints_);
  }

  /*! @brief return the nubmer of particles in the swarm.
   * @return the number of particles in the swarm
   */
  int num_owned_cells() {
    return npoints_;
  }

  /*! @brief Get the coordinates of the particle bounding boxes.
   * @param c index of particle of interest
   * @param cell_coord the bounding box vertices
   */
  void cell_get_coordinates(int c, PointVec* cell_coord) {
    size_t ncorners=1; for (size_t i=0; i<dim; i++) ncorners*=2;
    for (size_t corner=0; corner<ncorners; corner++) {
      size_t offset[dim], mask=1;
      for (size_t j=0; j<dim; j++) {offset[j] = corner & mask; mask*=2;}
      for (size_t j=0; j<dim; j++) cell_coord[corner][j] =
          (*points_)[c][j]+offset[j]*extents_[c][j];
    }
  }

 private:
  /** the centers of the particles */
  PointVecPtr points_;

  /** the widths of the particle bounding box sides */
  PointVecPtr extents_;

  /** the number of particles in the swarm */
  const size_t npoints_;
};

}
}

#endif // SWARM_H_INC_

