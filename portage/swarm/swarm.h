/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#pragma once

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

/**
 * @class Swarm
 *
 * @brief Particle field class.
 *
 * @tparam dim: the dimension of the problem.
 */
template<int dim>
class Swarm {
public:

  /**
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
   * @brief Create a particle field from a set of mesh wrappers.
   *
   * @tparam Mesh: type of mesh wrappers.
   * @param meshes: a set of meshes.
   * @param entity: the entity kind.
   */
  template<typename Mesh>
  Swarm(std::vector<Mesh*> const& meshes, Wonton::Entity_kind entity);

  /**
   * @brief Get the dimension.
   *
   * @return the dimension.
   */
  int space_dimension() const { return dim; }

  /**
   * @brief Get owned particles count.
   *
   * @return number of owned particles.
   */
  int num_owned_particles() const { return num_owned_points_; }

  /**
   * @brief Get ghost particles count.
   *
   * @return number of ghost particles.
   */
  int num_ghost_particles() const { return points_.size() - num_owned_points_; }

  /**
   * @brief Get particles count
   *
   * @param type: entity type (owned, ghost, all).
   * @return the particle count for the given entity type.
   */
  int num_particles(Wonton::Entity_type type = Wonton::ALL) const {
    switch (type) {
      case Wonton::PARALLEL_OWNED: return num_owned_particles();
      case Wonton::PARALLEL_GHOST: return num_ghost_particles();
      default: return num_owned_particles() + num_ghost_particles();;
    }
  }

  /**
   * @brief Get the coordinates of the particle.
   *
   * @param index: index of particle.
   * @return its coordinates.
   */
  Wonton::Point<dim> get_particle_coordinates(int index) const {
    return points_[index];
  }

  /**
   * @brief Get an iterator at first element of the particle field.
   *
   * @param kind: entity kind (kept for consistency with mesh iterators)
   * @param type: entity type (kept for consistency with mesh iterators)
   * @return an iterator pointing to the first element of the particle field.
   */
  counting_iterator begin(Entity_kind const kind = Wonton::PARTICLE,
                          Entity_type const type = Wonton::ALL) const {

    assert(kind == Wonton::PARTICLE);
    return make_counting_iterator(0);
  }

  /**
   * @brief Get an iterator at last element of the particle field.
   *
   * @param kind: entity kind (kept for consistency with mesh iterators)
   * @param type: entity type (kept for consistency with mesh iterators)
   * @return an iterator pointing to the last element of the particle field.
   */
  counting_iterator end(Entity_kind const kind = Wonton::PARTICLE,
                        Entity_type const type = Wonton::ALL) const {

    assert(kind == Wonton::PARTICLE);
    return (make_counting_iterator(0) + num_particles(type));
  }

  /**
   * @brief Extend the particle field.
   *
   */
  void extend_particle_list(std::vector<Wonton::Point<dim>> const& list) {
    int const old_size = points_.size();
    int const new_size = old_size + list.size();
    points_.resize(new_size);
    std::copy(list.begin(), list.end(), points_.begin() + old_size);
  }

private:
  /** the centers of the particles */
  Portage::vector<Wonton::Point<dim>> points_;

  /** the number of owned particles in the swarm */
  size_t num_owned_points_;
};

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
  std::uniform_real_distribution<double> generator(0.0, 1.0);

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
  std::uniform_real_distribution<double> generator(0.0, 1.0);

  if (distribution == 0) {
    // resize field and update coordinates
    num_owned_points_ = num_particles;
    points_.resize(num_owned_points_);

    for (int i = 0; i < num_particles; i++)
      points_[i] = Wonton::Point<3>(x_min + (x_max - x_min) * generator(engine),
                                    y_min + (y_max - y_min) * generator(engine),
                                    z_min + (z_max - z_min) * generator(engine));
  } else {
    auto const cubic_root = std::pow(num_particles, 1./3.);
    auto const num_per_dim = static_cast<int>(std::ceil(cubic_root));
    auto const hx = (x_max - x_min) / (cubic_root - 1);
    auto const hy = (y_max - y_min) / (cubic_root - 1);
    auto const hz = (z_max - z_min) / (cubic_root - 1);

    // resize field
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

