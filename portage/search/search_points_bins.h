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
   * @brief Disabled default constructor.
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

    // step 1: build bounding box of the radius of target point
    Wonton::Point<dim> box_min, box_max;
    for (int d = 0; d < dim; ++d) {
      box_min[d] = std::max(p[d] - h[d], p_min_[d]);
      box_max[d] = std::min(p[d] + h[d], p_max_[d]);
    }

    // step 2: filter cells overlapped by the bounding box
    std::vector<int> cells;
    int bin_min[dim];
    int bin_max[dim];
    retrieve_index(box_min, bin_min);
    retrieve_index(box_max, bin_max);

    if (dim == 1) {
      for (int i = bin_min[0]; i < bin_max[0]; ++i) {
        cells.emplace_back(i);
      }
    } else if (dim == 2) {
      for (int j = bin_min[1]; j < bin_max[1]; ++j) {
        for (int i = bin_min[0]; i < bin_max[0]; ++i) {
          cells.emplace_back(i + j * num_sides_);
        }
      }
    } else if (dim == 3) {
      for (int k = bin_min[2]; k < bin_max[2]; ++k) {
        for (int j = bin_min[1]; j < bin_max[1]; ++j) {
          for (int i = bin_min[0]; i < bin_max[0]; ++i) {
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
      int indices[dim];
      int const j = retrieve_index(p, indices);
      bins_[j].emplace_back(i);
    }
  }

  /**
   * @brief
   *
   * @param p
   * @param indices
   * @return
   */
  int retrieve_index(Wonton::Point<dim> const& p, int indices[dim]) const {

    assert(num_sides_ > 0);
    assert(indices != nullptr);

    // step 1: (x,y,z) to (i,j,k)
    for (int d = 0; d < dim; ++d) {
      double const t = p[d] - p_min_[d];
      if (t < 0) {
        throw std::runtime_error("outside bounding box");
      } else {
        double const range = p_max_[d] - p_min_[d];
        indices[d] = static_cast<int>(std::floor(t * num_sides_ / range));
      }
    }

    // step 2: (i,j,k) to i'
    int index = indices[0];
    for (int d = 1; d < dim; ++d) {
      index += indices[d] * static_cast<int>(std::pow(num_sides_, d));
    }
    return index;
  }

  /** source points */
  SourceSwarm const& source_swarm_;
  /** target points */
  TargetSwarm const& target_swarm_;
  /** search radius for each target point */
  std::vector<Wonton::Point<dim>> const& target_radius_;
  /** helper grid bounding box */
  Wonton::Point<dim> p_min_, p_max_;
  /** source point bins for each cell of helper grid */
  std::vector<std::vector<int>> bins_;
  /** number of sides of helper grid */
  int num_sides_ = 0;
};

} // namespace Portage

#endif