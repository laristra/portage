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
                   Meshfree::WeightCenter center = Meshfree::Gather)
    : source_swarm_(source_swarm),
      target_swarm_(target_swarm),
      target_radius_(target_radius) {

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
  std::vector<int> operator() (int target_id) const {

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
      p_min[d] = std::numeric_limits<double>::max();
      p_max[d] = std::numeric_limits<double>::lowest();
    }

    for (int i = 0; i < num_source_points; ++i) {
      auto const& p = target_swarm_.get_particle_coordinates(i);
      for (int d = 0; d < dim; ++d) {
        if (p[d] < p_min[d]) { p_min[d] = p[d]; }
        if (p[d] > p_max[d]) { p_max[d] = p[d]; }
      }
    }

    // step 2: discretize into cells
    num_sides_ = 0;
    for (int d = 0; d < dim; ++d) {
      double range = p_max[d] - p_min[d];
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
  int retrieve_index(Wonton::Point<dim> const& p) const {

    assert(num_sides_);

    // step 1: (x,y,z) to (i,j,k)
    int index[dim];
    for (int d = 0; d < dim; ++d) {
      double t = p[d] - p_min[d];
      double range = p_max[d] - p_min[d];
      index[d] = static_cast<int>(std::floor(t * num_sides_ / range));
    }

    // step 2: (i,j,k) to flat array index
    for (int d = 1; d < dim; ++d) {
      index[0] += index[d] * static_cast<int>(std::pow(num_sides_, d));
    }
    return index[0];
  }


  SourceSwarm const& source_swarm_;
  TargetSwarm const& target_swarm_;
  Wonton::vector<Wonton::Point<dim>> const& target_radius_;
  Wonton::Point<dim> p_min, p_max;
  std::vector<std::vector<int>> bins_; // flatten multi-dimensional array
  int num_sides_ = 0;
};

} // namespace Portage

#endif