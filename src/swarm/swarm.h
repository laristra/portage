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
#include <ctime>
#include <algorithm>
#include <cstdlib>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

namespace Portage {
namespace Meshfree {

using std::string;
using std::vector;
using std::shared_ptr;
using std::make_shared;

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

  using PointVecPtr = shared_ptr<std::vector<Point<dim>>>;
  using PointVec = vector<Point<dim>>;
  
  
  /*!
   * @brief A particle has a center point and smoothing lengths in each dimension.
   * @param points center points of the particles
   * @param extents the widths of the particle bounding boxes in each dimension
   */
  Swarm(PointVecPtr points)
      : points_(points), npoints_(points_->size()) 
  {}

  /*! @brief Dimensionality of points */

  unsigned int space_dimension() const {
    return dim;
  }

  /*! @brief return the number of particles in the swarm owned by this processor
   * @return the number of owned particles in the swarm
   */
  int num_owned_particles() const {
    return npoints_;
  }

  /*! @brief return the number of ghost (not owned by this processor)
   * particles in the swarm.
   * @return the number of ghost particles in the swarm
   */
  int num_ghost_particles() const {
    return 0;
  }

  /*! @brief return the number of particles of given type
   * @return the number of owned + ghost particles in the swarm
   */
  int num_particles(Entity_type etype = ALL) const {
    switch (etype) {
      case PARALLEL_OWNED:
        return num_owned_particles();
      case PARALLEL_GHOST:
        return num_ghost_particles();
      case ALL:
        return num_owned_particles() + num_ghost_particles();
      default:
        return 0;
    }
    return 0;
  }


  /*! @brief Get the coordinates of the particle,
   * @param index index of particle of interest
   * @return the particle coordinates
   */
  Point<dim> get_particle_coordinates(const size_t index) const {
    return (*points_)[index];
  }

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity = Entity_kind::PARTICLE,
                          Entity_type const etype = Entity_type::ALL) const {
    assert(entity == PARTICLE);
    int start_index = 0;
    return make_counting_iterator(start_index);
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity = Entity_kind::PARTICLE,
                        Entity_type const etype = Entity_type::ALL) const {
    assert(entity == PARTICLE);
    int start_index = 0;
    return (make_counting_iterator(start_index) + num_particles(etype));
  }

  

 private:
  /** the centers of the particles */
  PointVecPtr points_;

  /** the number of particles in the swarm */
  const size_t npoints_;
};

// Factory for making swarms in 1 dimensions with random or uniform
// distribution of particles. Reason we are having to return a
// shared_ptr to the Swarm is because of how its used in the
// DriverTest

std::shared_ptr<Swarm<1>> SwarmFactory(double xmin, double xmax,
                                         unsigned int nparticles,
                                         unsigned int distribution) {

  auto pts = std::make_shared<typename Swarm<1>::PointVec>(nparticles);
  if (distribution == 0) {  // random distribution of particles
    srand(time(NULL));
    unsigned int rand_state;
    for (size_t i = 0; i < nparticles; i++)
      (*pts)[i][0] = xmin +
          (xmax-xmin)*(static_cast<double>(rand_r(&rand_state))/RAND_MAX);
  } else {
    double h = (xmax-xmin)/(nparticles-1);
    for (size_t i = 0; i < nparticles; i++)
      (*pts)[i][0] = xmin + i*h;
    
    if (distribution == 2) {
      srand(time(NULL));
      unsigned int rand_state;
      double h = (xmax-xmin)/(nparticles-1);
      for (size_t i = 0; i < nparticles; i++)
        (*pts)[i][0] +=
            0.25*h*((2*static_cast<double>(rand_r(&rand_state))/RAND_MAX)-1.0);
    }
  }
  
  std::shared_ptr<Swarm<1>> swarm = make_shared<Swarm<1>>(pts);
  return swarm;
}

// Factories for making swarms in 1/2/3 dimensions with random or
// uniform distribution of particles

std::shared_ptr<Swarm<2>> SwarmFactory(double xmin, double ymin,
                                       double xmax, double ymax,
                                       unsigned int nparticles,
                                       unsigned int distribution) {

  auto pts = std::make_shared<typename Swarm<2>::PointVec>(nparticles);
  if (distribution == 0) {  // random distribution of particles
    srand(time(NULL));
    unsigned int rand_state;
    for (size_t i = 0; i < nparticles; i++) {
      (*pts)[i][0] = xmin +
          (xmax-xmin)*(static_cast<double>(rand_r(&rand_state))/RAND_MAX);
      (*pts)[i][1] = ymin +
          (ymax-ymin)*(static_cast<double>(rand_r(&rand_state))/RAND_MAX);
    }
  } else {
    int npdim = sqrt(nparticles);
    if (npdim*npdim != nparticles) {
      std::cerr << "Requested number of particles not a perfect square\n";
      std::cerr << "Generating only " << npdim << "particles in each dimension\n";
    }
    double hx = (xmax-xmin)/(npdim-1);
    double hy = (ymax-ymin)/(npdim-1);
    int n = 0;
    for (size_t i = 0; i < npdim; i++) {
      for (size_t j = 0; j < npdim; j++) {
        (*pts)[n][0] = xmin + i*hx;
        (*pts)[n][1] = ymin + j*hy;
        n++;
      }
    }
    if (distribution == 2) {
      srand(time(NULL));
      unsigned int rand_state;
      for (size_t i = 0; i < nparticles; i++) {
        (*pts)[i][0] +=
            0.25*hx*((2*static_cast<double>(rand_r(&rand_state))/RAND_MAX)-1.0);
        (*pts)[i][1] +=
            0.25*hy*((2*static_cast<double>(rand_r(&rand_state))/RAND_MAX)-1.0);
      }
    }
  }
  
  std::shared_ptr<Swarm<2>> swarm = make_shared<Swarm<2>>(pts);
  return swarm;
}


std::shared_ptr<Swarm<3>> SwarmFactory(double xmin, double ymin, double zmin,
                                       double xmax, double ymax, double zmax,
                                       unsigned int nparticles,
                                       unsigned int distribution) {

  auto pts = std::make_shared<typename Swarm<3>::PointVec>(nparticles);
  if (distribution == 0) {  // Random distribution
    srand(time(NULL));
    unsigned int rand_state;
    for (size_t i = 0; i < nparticles; i++) {
      (*pts)[i][0] = xmin +
          (xmax-xmin)*(static_cast<double>(rand_r(&rand_state))/RAND_MAX);
      (*pts)[i][1] = ymin +
          (ymax-ymin)*(static_cast<double>(rand_r(&rand_state))/RAND_MAX);
      (*pts)[i][2] = zmin +
          (zmax-zmin)*(static_cast<double>(rand_r(&rand_state))/RAND_MAX);
    }
  } else {
    int npdim = pow(nparticles, 1.0/3.0);
    if (npdim*npdim*npdim != nparticles) {
      std::cerr << "Requested number of particles not a perfect cube\n";
      std::cerr << "Generating only " << npdim << "particles in each dimension\n";
    }
    double hx = (xmax-xmin)/npdim;
    double hy = (ymax-ymin)/npdim;
    double hz = (zmax-zmin)/npdim;
    int n = 0;
    for (size_t i = 0; i < npdim; i++) {
      for (size_t j = 0; j < npdim; j++) {
        for (size_t k = 0; k < npdim; k++) {
          (*pts)[n][0] = xmin + i*hx;
          (*pts)[n][1] = ymin + j*hy;
          (*pts)[n][2] = zmin + k*hz;
          n++;
        }
      }
    }
    if (distribution == 2) {
      srand(time(NULL));
      unsigned int rand_state;
      for (size_t i = 0; i < nparticles; i++) {
        (*pts)[i][0] +=
            0.25*hx*((2*static_cast<double>(rand_r(&rand_state))/RAND_MAX)-1.0);
        (*pts)[i][1] +=
            0.25*hy*((2*static_cast<double>(rand_r(&rand_state))/RAND_MAX)-1.0);
        (*pts)[i][2] +=
            0.25*hz*((2*static_cast<double>(rand_r(&rand_state))/RAND_MAX)-1.0);
      }
    }
  }
  
  std::shared_ptr<Swarm<3>> swarm = make_shared<Swarm<3>>(pts);
  return swarm;
}

}
}

#endif // SWARM_H_INC_

