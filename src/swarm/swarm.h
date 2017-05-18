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
      : points_(points), extents_(extents), npoints_(points_->size()) 
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

    // cell_coord->resize(ncorners);
    // for (size_t corner=0; corner<ncorners; corner++) {
    //   size_t offset[dim], mask=1;
    //   for (size_t j=0; j<dim; j++) {offset[j] = corner & mask; mask*=2;}
    //   for (size_t j=0; j<dim; j++) (*cell_coord)[corner][j] =
    //       (*points_)[c][j]+offset[j]*(*extents_)[c][j];
    // }

    // Alternative way - couldn't figure above method to make it work
    // OK to fix above code and blow this one away

    Point<dim>& size = (*extents_)[c];
    Point<dim> botcorner= Point<dim>((*points_)[c] - 2*size);

    Point<dim> corner = botcorner;
    for (size_t i = 0; i < 2; i++) {
      corner[0] += i*4*size[0];
      if (dim > 1) {
        for (size_t j = 0; j < 2; j++) {
          corner[1] += j*4*size[1];
          if (dim > 2) {
            for (size_t k = 0; k < 2; k++) {
              corner[2] += k*4*size[2];
              cell_coord->push_back(corner);
            }
            corner[2] = botcorner[2];  // reset z coordinate 
          }
          else
            cell_coord->push_back(corner);
        }
        corner[1] = botcorner[1];  // reset y coordinate
      }
      else
        cell_coord->push_back(corner);
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

