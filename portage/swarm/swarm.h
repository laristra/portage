/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef PORTAGE_SWARM_SWARM_H_
#define PORTAGE_SWARM_SWARM_H_

#include <vector>
#include <memory>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <iostream>
#include <random>

// portage includes
#include "portage/support/portage.h"
#include "wonton/support/Point.h"

namespace Portage { namespace Meshfree {

/*!
 @class Swarm "swarm.h"
 @brief An effective "mesh" class for a collection disconnected points (particles).
 @tparam dim The dimension of the problem space.
 */
template<int dim>
class Swarm {
 public:

  // When using PointVec and PointVecPtr in another file, they have to
  // be prefixed by the keyword typename

  //  SEE --- https://stackoverflow.com/questions/610245/where-and-why-do-i-have-to-put-the-template-and-typename-keywords

  using PointVecPtr = std::shared_ptr<Portage::vector<Wonton::Point<dim>>>;
  using PointVec = vector<Wonton::Point<dim>>;

  /*!
   * @brief A particle has a center point and smoothing lengths in each dimension.
   * @param points center points of the particles
   * @param extents the widths of the particle bounding boxes in each dimension
   */

  /**
   * @brief Create a particle field from a given list.
   *
   * @param points
   */
  explicit Swarm(Portage::vector<Wonton::Point<dim>> const& points)
    : points_(points),
      num_owned_points_(points_.size())
  {}

  /**
   * @brief Generate a random or uniform particle field.
   *
   * @param num_particles: number of points.
   * @param p_min: lower bounds on points coordinates.
   * @param p_max: upper bounds on points coordinates.
   * @param distribution: the kind of random numbers distribution.
   * @param seed: a seed used by the random number generator.
   */
  Swarm(int num_particles, int distribution,
        unsigned user_seed = 0,
        double x_min = 0.0, double x_max = 0.0,
        double y_min = 0.0, double y_max = 0.0,
        double z_min = 0.0, double z_max = 0.0);

  /**
   * @brief Create a particle field from an arbitary mesh wrapper.
   *
   * @tparam Mesh: type of mesh wrapper.
   * @param mesh: a mesh.
   * @param entity: the entity kind.
   */
  template<typename Mesh>
  Swarm(Mesh const& mesh, Wonton::Entity_kind entity);

  /**
   * @brief Create a particle field from a set of arbitary mesh wrappers.
   *
   * @tparam Mesh: type of mesh wrappers.
   * @param meshes: a set of meshes.
   * @param entity: the entity kind.
   */
  template<typename Mesh>
  Swarm(std::vector<Mesh*> const& meshes, Wonton::Entity_kind entity);

  /*! @brief Dimensionality of points */
  unsigned int space_dimension() const {
    return dim;
  }

  /*! @brief return the number of particles in the swarm owned by this processor
   * @return the number of owned particles in the swarm
   */
  int num_owned_particles() const {
    return num_owned_points_;
  }

  /*! @brief return the number of ghost (not owned by this processor)
   * particles in the swarm.
   * @return the number of ghost particles in the swarm that are stored
   * at the end of the particle list i.e., between num_owned_particles
   * and num_owned_particles+num_ghost_particle.
   */
  int num_ghost_particles() const {
    return points_.size() - num_owned_points_;
  }

  /*! @brief return the number of particles of given type
   * @return the number of owned + ghost particles in the swarm
   */
  int num_particles(Entity_type etype = Entity_type::ALL) const {
    switch (etype) {
      case Entity_type::PARALLEL_OWNED:
        return num_owned_particles();
      case Entity_type::PARALLEL_GHOST:
        return num_ghost_particles();
      case Entity_type::ALL:
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
  Wonton::Point<dim> get_particle_coordinates(const size_t index) const {
    return points_[index];
  }

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity = Entity_kind::PARTICLE,
                          Entity_type const etype = Entity_type::ALL) const {
    assert(entity == Entity_kind::PARTICLE);
    int start_index = 0;
    return make_counting_iterator(start_index);
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity = Entity_kind::PARTICLE,
                        Entity_type const etype = Entity_type::ALL) const {
    assert(entity == Entity_kind::PARTICLE);
    int start_index = 0;
    return (make_counting_iterator(start_index) + num_particles(etype));
  }

  /*! @brief Add new particles to swarm
   * @return
   */
  void extend_particle_list(std::vector<Wonton::Point<dim>>& new_pts)
  {
    (*points_).insert((*points_).end(), new_pts.begin(), new_pts.end());
  }

 private:
  /** the centers of the particles */
  Portage::vector<Wonton::Point<dim>> points_;

  /** the number of owned particles in the swarm */
  size_t num_owned_points_;
};

// Get a reasonable random number seed if needed.
inline
unsigned long get_seed() {
  unsigned long seed;
  auto clock = std::chrono::high_resolution_clock::now();
  auto clockdur = clock.time_since_epoch();
  auto clockmicro = std::chrono::duration_cast<std::chrono::microseconds>(clockdur);
  seed = clockmicro.count();
  seed &= (1<<16)-1;
  return seed;
}

/* -------------------------------------------------------------------------- */
template<>
Swarm<1>::Swarm(int num_particles, int distribution, unsigned user_seed,
                double x_min, double x_max,
                double /* unused */, double /* unused */,
                double /* unused */, double /* unused */) {

  assert(num_particles > 0);
  assert(0 <= distribution and distribution <= 2);

  // resize field
  points_.resize(num_particles);
  num_owned_points_ = num_particles;

  // set the random engine and generator
  std::random_device device;
  std::mt19937 engine { user_seed ? user_seed : device() };
  std::uniform_real_distribution<> generator(0.0, 1.0);

  // update coordinates
  if (distribution == 0) {
    for (auto&& current : points_)
      current = Wonton::Point<1>(x_min + (x_max - x_min) * generator(engine));

  } else {
    double const h = (x_max - x_min) / (num_particles - 1);
    for (int i = 0; i < num_particles; i++)
      points_[i] = Wonton::Point<1>(x_min + i * h);

    if (distribution == 2) {
      for (int i = 0; i < num_particles; i++) {
        Wonton::Point<1> current = points_[i];
        current[0] += 0.25 * h * (2 * generator(engine) - 1);
        points_[i] = Wonton::Point<1>(std::max(x_min, std::min(x_max, current[0])));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
Swarm<2>::Swarm(int num_particles, int distribution, unsigned user_seed,
                double x_min, double x_max,
                double y_min, double y_max,
                double /* unused */, double /* unused */) {

  assert(num_particles > 0);
  assert(0 <= distribution and distribution <= 2);

  // set the random engine and generator
  std::random_device device;
  std::mt19937 engine { user_seed ? user_seed : device() };
  std::uniform_real_distribution<double> generator(0.0, 1.0);

  if (distribution == 0) {
    // resize field and update coordinates
    num_owned_points_ = num_particles;
    points_.resize(num_owned_points_);

    for (auto&& current : points_)
      current = Wonton::Point<2>(x_min + (x_max - x_min) * generator(engine),
                                 y_min + (y_max - y_min) * generator(engine));
  } else {
    // resize field
    int const num_per_dim = std::floor(std::sqrt(num_particles));
    num_owned_points_ = num_per_dim * num_per_dim;
    points_.resize(num_owned_points_);
    double const hx = (x_max - x_min) / (num_per_dim - 1);
    double const hy = (y_max - y_min) / (num_per_dim - 1);

    // update coordinates
    int k = 0;
    for (int i = 0; i < num_per_dim; i++)
      for (int j = 0; j < num_per_dim; ++j)
        points_[k++] = Wonton::Point<2>(x_min + i * hx, y_min + j * hy);

    if (distribution == 2) {
      for (auto&& current : points_) {
        Wonton::Point<2> copy = current;
        copy[0] += 0.25 * hx * (2 * generator(engine) - 1);
        copy[1] += 0.25 * hy * (2 * generator(engine) - 1);
        current = Wonton::Point<2>(std::max(x_min, std::min(x_max, copy[0])),
                                   std::max(y_min, std::min(y_max, copy[1])));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template<>
Swarm<3>::Swarm(int num_particles, int distribution, unsigned user_seed,
                double x_min, double x_max,
                double y_min, double y_max,
                double z_min, double z_max) {

  assert(num_particles > 0);
  assert(0 <= distribution and distribution <= 2);

  // set the random engine and generator
  std::random_device device;
  std::mt19937 engine { user_seed ? user_seed : device() };
  std::uniform_real_distribution<> generator(0.0, 1.0);

  if (distribution == 0) {
    // resize field and update coordinates
    num_owned_points_ = num_particles;
    points_.resize(num_owned_points_);

    for (int i = 0; i < num_particles; i++)
      points_[i] = Wonton::Point<3>(x_min + (x_max - x_min) * generator(engine),
                                    y_min + (y_max - y_min) * generator(engine),
                                    z_min + (z_max - z_min) * generator(engine));
  } else {
    // resize field
    auto const cubic_root = std::pow(num_particles, 1./3.);
    auto const num_per_dim = static_cast<int>(std::ceil(cubic_root));
    auto const hx = (x_max - x_min) / (cubic_root - 1);
    auto const hy = (y_max - y_min) / (cubic_root - 1);
    auto const hz = (z_max - z_min) / (cubic_root - 1);

//    std::printf("hx: %.3f, hy: %.3f, hz: %.3f\n", hx, hy, hz);
//    std::printf("cubic_root: %.3f, num_per_dim: %d", cubic_root, num_per_dim);
//    std::fflush(stdout);

    num_owned_points_ = num_per_dim * num_per_dim * num_per_dim;
    points_.resize(num_owned_points_);

    // update coordinates
    int n = 0;
    for (int i = 0; i < num_per_dim; i++)
      for (int j = 0; j < num_per_dim; ++j)
        for (int k = 0; k < num_per_dim; ++k)
          points_[n++] = Wonton::Point<3>(x_min + i * hx,
                                          y_min + j * hy,
                                          z_min + k * hz);

    if (distribution == 2) {
      for (auto&& current : points_) {
        Wonton::Point<3> copy = current;
        copy[0] += 0.25 * hx * (2 * generator(engine) - 1);
        copy[1] += 0.25 * hy * (2 * generator(engine) - 1);
        copy[2] += 0.25 * hz * (2 * generator(engine) - 1);
        current = Wonton::Point<3>(std::max(x_min, std::min(x_max, copy[0])),
                                   std::max(y_min, std::min(y_max, copy[1])),
                                   std::max(z_min, std::min(z_max, copy[2])));
      }
    }
  }
}

#if USE_FACTORIES
// Factories for making swarms in 1/2/3 dimensions with random or uniform
// distribution of particles. Reason we are having to return a
// shared_ptr to the Swarm is because of how its used in the
// DriverTest

inline
std::shared_ptr<Swarm<1>> SwarmFactory(double xmin, double xmax,
                                       unsigned int nparticles,
                                       unsigned int distribution,
                                       unsigned int rand_seed=0) {

  shared_ptr<typename Swarm<1>::PointVec> pts_sp = std::make_shared<typename Swarm<1>::PointVec>(nparticles);
  typename Swarm<1>::PointVec &pts(*pts_sp);

  if (distribution == 0 or distribution == 2) {  // random distribution of particles
    long int seed;
    if (rand_seed == 0) {
      seed = get_seed();
    } else {
      seed = rand_seed;
    }
    srand48(seed);
  }
  if (distribution == 0) {  // random distribution of particles
    for (size_t i=0; i<nparticles; i++) {
      // for some reason intel does not like these two lines combined when using thrust, i.e. pts[i][0]
      Point<1> pt=pts[i];
      double rnum = drand48();
      pt[0] = xmin + (xmax-xmin)*rnum;
      pts[i] = pt;
    }
  } else { // regular distribution
    double h = (xmax-xmin)/(nparticles-1);
    for (size_t i=0; i < nparticles; i++) {
      Point<1> pt=pts[i];
      pt[0] = xmin + i*h;
      pts[i] = pt;
    }

    if (distribution == 2) { // perturbed regular
      for (size_t i=0; i < nparticles; i++) {
        Point<1> pt=pts[i];
        double rnum = drand48();
        pt[0] += 0.25*h*(2*rnum-1.0);
	pt[0] = fmax(xmin,fmin(xmax,pt[0]));
	pts[i] = pt;
      }
    }
  }

  std::shared_ptr<Swarm<1>> swarm = make_shared<Swarm<1>>(pts_sp);
  return swarm;
}

inline
std::shared_ptr<Swarm<2>> SwarmFactory(double xmin, double ymin,
                                       double xmax, double ymax,
                                       unsigned int nparticles,
                                       unsigned int distribution,
                                       unsigned int rand_seed=0) {

  auto pts_sp = std::make_shared<typename Swarm<2>::PointVec>(nparticles);
  typename Swarm<2>::PointVec &pts(*pts_sp);

  if (distribution == 0 or distribution == 2) {  // random distribution of particles
    long int seed;
    if (rand_seed == 0) {
      seed = get_seed();

    } else {
      seed = rand_seed;
    }
    srand48(seed);
  }
  if (distribution == 0) {
    for (size_t i = 0; i < nparticles; i++) {
      Point<2> pt=pts[i];
      double rnum = drand48();
      pt[0] = xmin + (xmax-xmin)*rnum;
      rnum = drand48();
      pt[1] = ymin + (ymax-ymin)*rnum;
      pts[i] = pt;
    }
  } else {
    int npdim = sqrt(nparticles);
    double hx = (xmax-xmin)/(npdim-1);
    double hy = (ymax-ymin)/(npdim-1);
    int n = 0;
    for (size_t i = 0; i < npdim; i++) {
      for (size_t j = 0; j < npdim; j++) {
	Point<2> pt=pts[n];
        pt[0] = xmin + i*hx;
        pt[1] = ymin + j*hy;
	pts[n] = pt;
        n++;
      }
    }
    if (distribution == 2) {
      for (size_t i = 0; i < nparticles; i++) {
	Point<2> pt=pts[i];
        double rnum = drand48();
        pt[0] += 0.25*hx*(2*rnum-1.0);
        rnum = drand48();
        pt[1] += 0.25*hy*(2*rnum-1.0);
        pt[0] = fmax(xmin,fmin(xmax,pt[0]));
        pt[1] = fmax(ymin,fmin(ymax,pt[1]));
        pts[i] = pt;
      }
    }
  }

  std::shared_ptr<Swarm<2>> swarm = make_shared<Swarm<2>>(pts_sp);
  return swarm;
}

inline
std::shared_ptr<Swarm<3>> SwarmFactory(double xmin, double ymin, double zmin,
                                       double xmax, double ymax, double zmax,
                                       unsigned int nparticles,
                                       unsigned int distribution,
                                       unsigned int rand_seed=0) {

  auto pts_sp = std::make_shared<typename Swarm<3>::PointVec>(nparticles);
  typename Swarm<3>::PointVec &pts(*pts_sp);

  if (distribution == 0 or distribution == 2) {  // random distribution of particles
    long int seed;
    if (rand_seed == 0) {
      seed = get_seed();
    } else {
      seed = rand_seed;
    }
    srand48(seed);
  }
  if (distribution == 0) {
    for (size_t i = 0; i < nparticles; i++) {
      Point<3> pt=pts[i];
      double rnum = drand48();
      pt[0] = xmin + (xmax-xmin)*rnum;
      rnum = drand48();
      pt[1] = ymin + (ymax-ymin)*rnum;
      rnum = drand48();
      pt[2] = zmin + (zmax-zmin)*rnum;
      pts[i] = pt;
    }
  } else {
    int npdim = static_cast<int>(pow(nparticles, 1.0/3.0) + 0.5);
    if (npdim*npdim*npdim != nparticles) {
      std::cerr << "Requested number of particles not a perfect cube\n";
      std::cerr << "Generating only " << npdim << " particles in each dimension\n";
    }
    double hx = (xmax-xmin)/(npdim-1);
    double hy = (ymax-ymin)/(npdim-1);
    double hz = (zmax-zmin)/(npdim-1);
    int n = 0;
    for (size_t i = 0; i < npdim; i++) {
      for (size_t j = 0; j < npdim; j++) {
        for (size_t k = 0; k < npdim; k++) {
	        Point<3> pt=pts[n];
          pt[0] = xmin + i*hx;
          pt[1] = ymin + j*hy;
          pt[2] = zmin + k*hz;
	        pts[n] = pt;
          n++;
        }
      }
    }
    if (distribution == 2) {
      for (size_t i = 0; i < nparticles; i++) {
	      Point<3> pt=pts[i];
        double rnum = drand48();
        pt[0] += 0.25*hx*(2*rnum-1.0);
        rnum = drand48();
        pt[1] += 0.25*hy*(2*rnum-1.0);
        rnum = drand48();
        pt[2] += 0.25*hz*(2*rnum-1.0);
        pt[0] = fmax(xmin,fmin(xmax,pt[0]));
        pt[1] = fmax(ymin,fmin(ymax,pt[1]));
        pt[2] = fmax(zmin,fmin(zmax,pt[2]));
        pts[i] = pt;
      }
    }
  }

  std::shared_ptr<Swarm<3>> swarm = make_shared<Swarm<3>>(pts_sp);
  return swarm;
}


/*!
 * @brief Create a Swarm from an arbitary mesh wrapper.
 * @tparm dim Spatial dimension
 * @tparm MeshWrapper class representing a mesh wrapper
 * @param wrapper Input mesh wrapper
 * @param entity  Where the data is located in the cell
 */
template<size_t dim, class MeshWrapper>
std::shared_ptr<Swarm<dim>> SwarmFactory(MeshWrapper &wrapper, Portage::Entity_kind entity)
{
  size_t npoints_owned;
  typename Swarm<dim>::PointVecPtr points;
  if (entity != Entity_kind::NODE and entity != Entity_kind::CELL) {
    throw(std::runtime_error("only nodes and cells allowed"));
  }

  if (entity == Entity_kind::NODE) {
    npoints_owned = wrapper.num_owned_nodes();
    points = make_shared<typename Swarm<dim>::PointVec>(npoints_owned);
    Point<dim> node;
    for (size_t i=0; i<npoints_owned; i++) {
      wrapper.node_get_coordinates(i, &node);
      (*points)[i] = node;
    }
  } else if (entity == Entity_kind::CELL) {
    npoints_owned = wrapper.num_owned_cells();
    points = make_shared<typename Swarm<dim>::PointVec>(npoints_owned);
    Point<dim> centroid;
    for (size_t i=0; i<npoints_owned; i++) {
      wrapper.cell_centroid(i, &centroid);
      (*points)[i] = centroid;
    }
  }

  std::shared_ptr<Swarm<dim>> result = make_shared<Swarm<dim>>(points);
  return result;
}


/*!
 * @brief Create a Swarm from an set of arbitary mesh wrappers.
 * This factory is useful for mapping from several meshes at once,
 * e.g. for analysis of multiple times, or multiple simulations at the same set of spatial points,
 * to obtain statistical properties like average, confidence bounds, and uncertaintities, for example.
 * @tparm dim Spatial dimension
 * @tparm MeshWrapper class representing a mesh wrapper
 * @param wrappers Input mesh wrappers
 * @param entity  Where the data is located in the cell
 */
template<size_t dim, class MeshWrapper>
  std::shared_ptr<Swarm<dim>> SwarmFactory(std::vector<MeshWrapper*> wrappers, Portage::Entity_kind entity)
{
  size_t npoints_owned;
  typename Swarm<dim>::PointVecPtr points = make_shared<typename Swarm<dim>::PointVec>();
  if (entity != Entity_kind::NODE and entity != Entity_kind::CELL) {
    throw(std::runtime_error("only nodes and cells allowed"));
  }

  size_t total_npoints_owned=0;
  if (entity == Entity_kind::NODE) {
    for (auto wrap=wrappers.begin(); wrap!=wrappers.end(); wrap++) {
      MeshWrapper &wrapper = **wrap;
      size_t npoints_owned = wrapper.num_owned_nodes();
      total_npoints_owned += npoints_owned;
      Point<dim> node;
      for (size_t i=0; i<npoints_owned; i++) {
        wrapper.node_get_coordinates(i, &node);
        (*points).push_back(node);
      }
    }
  } else if (entity == Entity_kind::CELL) {
    for (auto wrap=wrappers.begin(); wrap!=wrappers.end(); wrap++) {
      MeshWrapper &wrapper = **wrap;
      size_t npoints_owned = wrapper.num_owned_cells();
      total_npoints_owned += npoints_owned;
      Point<dim> centroid;
      for (size_t i=0; i<npoints_owned; i++) {
        wrapper.cell_centroid(i, &centroid);
        (*points).push_back(centroid);
      }
    }
  }

  if (points->size() != total_npoints_owned) {
    throw(std::runtime_error("size discrepancy"));
  }

  std::shared_ptr<Swarm<dim>> result = make_shared<Swarm<dim>>(points);
  return result;
}
#endif


/* -------------------------------------------------------------------------- */
template<int dim>
template<typename Mesh>
Swarm<dim>::Swarm(Mesh const& mesh, Wonton::Entity_kind entity) {
  switch (entity) {
    case Wonton::NODE: {
      num_owned_points_ = mesh.num_owned_nodes();
      points_.resize(num_owned_points_);
      Wonton::Point<dim> coord;
      for (int i = 0; i < num_owned_points_; ++i) {
        mesh.node_get_coordinates(i, &coord);
        points_[i] = coord;
      }
    } break;
    case Wonton::CELL: {
      num_owned_points_ = mesh.num_owned_cells();
      points_.resize(num_owned_points_);
      Wonton::Point<dim> centroid;
      for (int i = 0; i < num_owned_points_; ++i) {
        mesh.cell_centroid(i, &centroid);
        points_[i] = centroid;
      }
    } break;
    default: throw std::runtime_error("unsupported entity");
  }
}

/* -------------------------------------------------------------------------- */
template<int dim>
template<typename Mesh>
Swarm<dim>::Swarm(std::vector<Mesh*> const& meshes, Wonton::Entity_kind entity) {

  // resize particle field
  int n = 0;
  switch (entity) {
    case Wonton::NODE: for (auto&& mesh : meshes) n += mesh->num_owned_nodes(); break;
    case Wonton::CELL: for (auto&& mesh : meshes) n += mesh->num_owned_cells(); break;
    default: throw std::runtime_error("unsupported entity");
  }

  num_owned_points_ = n;
  points_.resize(num_owned_points_);

  n = 0;
  switch (entity) {
    case Wonton::NODE: {
      for (auto&& mesh : meshes) {
        int const num_nodes = mesh->num_owned_nodes();
        Wonton::Point<dim> coord;
        for (int i = 0; i < num_nodes; i++) {
          mesh->node_get_coordinates(i, &coord);
          points_[n++] = coord;
        }
      }
    } break;
    case Wonton::CELL: {
      for (auto&& mesh : meshes) {
        int const num_nodes = mesh->num_owned_nodes();
        Wonton::Point<dim> centroid;
        for (int i = 0; i < num_nodes; i++) {
          mesh->cell_centroid(i, &centroid);
          points_[n++] = centroid;
        }
      }
    } break;
    default: break;
  }
}

/* -------------------------------------------------------------------------- */

}}  // namespace Portage::Meshfree

#endif  // PORTAGE_SWARM_SWARM_H_
