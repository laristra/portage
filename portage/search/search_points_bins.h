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

private:

  /**
   * @brief
   *
   * @param h
   */
  void build_helper_grid(double h) {

    assert(h > 0.);

    // step 1: compute bounding box of source points
    for (int d = 0; d < dim; ++d) {
      p_min[d] = std::numeric_limits<double>::max();
      p_max[d] = std::numeric_limits<double>::lowest();
    }

    int const num_source_points = source_swarm_.num_particles();

    for (int i = 0; i < num_source_points; ++i) {
      auto const& p = target_swarm_.get_particle_coordinates(i);
      for (int d = 0; d < dim; ++d) {
        if (p[d] < p_min[d]) { p_min[d] = p[d]; }
        if (p[d] > p_max[d]) { p_max[d] = p[d]; }
      }
    }

    // step 2: discretize into cells
    // how to choose the spatial step h
    num_sides_ = 0;
    for (int d = 0; d < dim; ++d) {
      double range = p_max[d] - p_min[d];
      num_sides_ = std::max(num_sides_, static_cast<int>(range / h));
    }

    cells_.resize(std::pow(num_sides_, dim));
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
    switch (dim) {
      case 1: return index[0];
      case 2: return index[0] + index[1] * num_sides_;
      case 3: return index[0] + index[1] * num_sides_ + index[2] * num_sides_ * num_sides_;
      default: return -1;
    }
  }

  SourceSwarm const& source_swarm_;
  TargetSwarm const& target_swarm_;
  Wonton::vector<Wonton::Point<dim>> const& target_radius_;
  Wonton::Point<dim> p_min, p_max;
  std::vector<int> cells_; // flatten multi-dimensional array
  int num_sides_ = 0;
};

} // namespace Portage

#endif