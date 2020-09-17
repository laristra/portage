/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */
#ifndef SEARCH_POINTS_BINS_H
#define SEARCH_POINTS_BINS_H

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "portage/support/portage.h"
#include "portage/accumulate/accumulate.h"

namespace Portage {

using Meshfree::WeightCenter;

template<int dim, class SourceSwarm, class TargetSwarm>
class SearchPointsBins {

public:
  /**
   * @brief
   */
  SearchPointsBins() = delete;

  /**
   * @brief
   *
   * @param source_swarm
   * @param target_swarm
   * @param source_extents
   * @param target_extents
   */
  SearchPointsBins(SourceSwarm const& source_swarm,
                   TargetSwarm const& target_swarm,
                   Wonton::vector<Wonton::Point<dim>> const& /* unused */,
                   Wonton::vector<Wonton::Point<dim>> const& target_radius,
                   Meshfree::WeightCenter const center = Meshfree::Gather)
    : source_swarm_(source_swarm),
      target_swarm_(target_swarm),
      target_radius_(target_radius) {

    static_assert(dim > 0 and dim < 4, "invalid dimension");
    
    if (center == Meshfree::Scatter)
      throw std::runtime_error("scatter weight form not supported");

    // build helper grid to bin source points
    // - build bounding box and discretize into cells
    // - compute cell indices by hashing


  }

  /**
   * @brief
   *
   * @param target_id
   * @return
   */
  std::vector<int> operator() (int id) const {

    auto const& p = target_swarm_.get_particle_coordinates(id);
    auto const& h = target_radius_[id];

    // step 1: build bounding box of the radius of 'j'
    Wonton::Point<dim> b_min, b_max;
    for (int d = 0; d < dim; ++d) {
      b_min[d] = std::max(p[d] - h[d], p_min_[d]);
      b_max[d] = std::min(p[d] + h[d], p_max_[d]);
    }

    // step 2: push cells overlapped by the bounding box
    std::vector<int> cells;
    auto i_min = retrieve_indices(b_min);
    auto i_max = retrieve_indices(b_max);

    if (dim == 1) {
      for (int i = i_min[0]; i < i_max[0]; ++i) {
        cells.emplace_back(i);
      }
    } else if (dim == 2) {
      for (int j = i_min[1]; j < i_max[1]; ++j) {
        for (int i = i_min[0]; i < i_max[0]; ++i) {
          cells.emplace_back(i + j * num_sides_);
        }
      }
    } else if (dim == 3) {
      for (int k = i_min[2]; k < i_max[2]; ++k) {
        for (int j = i_min[1]; j < i_max[1]; ++j) {
          for (int i = i_min[0]; i < i_max[0]; ++i) {
            cells.emplace_back(i + j * num_sides_ + k * num_sides_ * num_sides_);
          }
        }
      }
    }

    // step 3: scan cells and check distance of each included source point
    std::vector<int> neighbors;

    for (int const& cell : cells) {
      for (int const& s : bins_[cell]) {
        bool contained = true;
        auto const& q = source_swarm_.get_particle_coordinates(s);
        for (int d = 0; d < dim; ++d) {
          if (std::abs(q[d] - p[d]) > h[d]) {
            contained = false;
            break;
          }
        }
        if (contained) {
          neighbors.emplace_back(s);
        }
      }
    }
    return neighbors;
  }

private:

  /**
   * @brief
   *
   * @param h
   */
  void build_helper_grid(double h) {

    assert(h > 0.);
    int const num_source_points = source_swarm_.num_particles();

    // step 1: compute bounding box of source points
    for (int d = 0; d < dim; ++d) {
      p_min_[d] = std::numeric_limits<double>::max();
      p_max_[d] = std::numeric_limits<double>::lowest();
    }

    for (int i = 0; i < num_source_points; ++i) {
      auto const& p = target_swarm_.get_particle_coordinates(i);
      for (int d = 0; d < dim; ++d) {
        if (p[d] < p_min_[d]) { p_min_[d] = p[d]; }
        if (p[d] > p_max_[d]) { p_max_[d] = p[d]; }
      }
    }

    // step 2: discretize into cells
    num_sides_ = 0;
    for (int d = 0; d < dim; ++d) {
      double range = p_max_[d] - p_min_[d];
      num_sides_ = std::max(num_sides_, static_cast<int>(range / h));
    }

    // step 3: resize grid and push source points
    int const num_cells = static_cast<int>(std::pow(num_sides_, dim));
    bins_.resize(num_cells);

    for (int i = 0; i < num_source_points; ++i) {
      auto const& p = target_swarm_.get_particle_coordinates(i);
      int const j = retrieve_index(p);
      assert(j > -1);
      bins_[j].emplace_back(i);
    }
  }

  /**
   * @brief
   *
   * @param p
   * @return
   */
  std::array<int, dim> retrieve_indices(Wonton::Point<dim> const& p) const {
    assert(num_sides_);
    std::array<int, dim> indices;
    for (int d = 0; d < dim; ++d) {
      double const t = p[d] - p_min_[d];
      if (t >= 0.) {
        double range = p_max_[d] - p_min_[d];
        indices[d] = static_cast<int>(std::floor(t * num_sides_ / range));
      } else
        throw std::runtime_error("p[d] < p_min[d]");
    }
    return indices;
  }

  int retrieve_index(std::array<int, dim> const& indices) const {
    switch (dim) {
      case 1: return indices[0];
      case 2: return indices[0] + indices[1] * num_sides_;
      case 3: return indices[0] + indices[1] * num_sides_ + indices[2] * num_sides_ * num_sides_;
      default: return -1;
    }
  }

  /**
   * @brief
   *
   * @param p
   * @return
   */
  int retrieve_index(Wonton::Point<dim> const& p) const {
    // (x,y,z) to (i,j,k) to i'
    return retrieve_index(retrieve_indices(p));
  }


  SourceSwarm const& source_swarm_;
  TargetSwarm const& target_swarm_;
  std::vector<Wonton::Point<dim>> const& target_radius_;
  Wonton::Point<dim> p_min_, p_max_;
  std::vector<std::vector<int>> bins_; // flatten multi-dimensional array
  int num_sides_ = 0;
};

} // namespace Portage

#endif