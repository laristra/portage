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

using namespace Meshfree;

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
                   WeightCenter const center = Gather,
                   double radius_scale = 1)
    : source_swarm_(source_swarm),
      target_swarm_(target_swarm),
      target_radius_(target_radius),
      radius_scale_(radius_scale) {

    static_assert(dim > 0 and dim < 4, "invalid dimension");

    if (center == Scatter)
      throw std::runtime_error("scatter weight form not supported");

    int const num_source_points = source_swarm_.num_particles();
    int const num_target_points = target_swarm_.num_particles();
    double radius[dim];
    int const num_max_sides = std::max(static_cast<int>(std::ceil(std::pow(0.1 * num_target_points, 1./dim))), 1);

    // step 1: compute bounding box of source points
    for (int d = 0; d < dim; ++d) {
      p_min_[d] = std::numeric_limits<double>::max();
      p_max_[d] = std::numeric_limits<double>::lowest();
      num_sides_[d] = 0;
      radius[d] = 0;
    }

    for (int t = 0; t < num_target_points; ++t) {
      auto const& p = target_swarm.get_particle_coordinates(t);
      Wonton::Point<dim> const& h = target_radius[t];
      for (int d = 0; d < dim; ++d) {
        p_min_[d] = std::min(p_min_[d], p[d] - radius_scale * h[d]);
        p_max_[d] = std::max(p_max_[d], p[d] + radius_scale * h[d]);
        radius[d] += std::abs(h[d]);
      }
    }

    std::cout << "bounding box computed" << p_min_ << " | "<< p_max_ << std::endl;

    // step 2: deduce discretization step from search radii
    int num_bins = 1;

    for (int d = 0; d < dim; ++d) {
      double const step = radius[d] * 2 * radius_scale / num_target_points;
      double const range = p_max_[d] - p_min_[d];
      num_sides_[d] = std::min(static_cast<int>(range / step), num_max_sides);
      num_bins *= num_sides_[d];
    }

    std::cout << "num bins: " << num_bins << std::endl;

    // step 3: push source points into bins
    bucket_.resize(num_bins);

    for (int s = 0; s < num_source_points; ++s) {
      auto const& p = source_swarm.get_particle_coordinates(s);
      if (not is_outside(p)) {
        int const i = deduce_bin_index(p);
        std::cout << "bin index: "<< i << std::endl;
        bucket_[i].emplace_back(s);
      }
    }
  }

  /**
   * @brief
   *
   * @param target_id
   * @return
   */
  std::vector<int> operator() (int id) const {

    Wonton::Point<dim> const& p = target_swarm_.get_particle_coordinates(id);
    Wonton::Point<dim> const& h = target_radius_[id];

    // step 1: build bounding box of the radius of target point
    Wonton::Point<dim> box_min, box_max;
    for (int d = 0; d < dim; ++d) {
      box_min[d] = p[d] - radius_scale_ * h[d];
      box_max[d] = p[d] + radius_scale_ * h[d];
    }

    // step 2: filter cells overlapped by the bounding box
    std::vector<int> cells;
    std::cout << " =========" << std::endl;
    std::cout << "compute box min, max" << std::endl;
    auto const first = deduce_cell_index(box_min);
    auto const last  = deduce_cell_index(box_max);

#ifndef NDEBUG
    for (int d = 0; d < dim; ++d) {
      assert(first[d] <= last[d]);
    }
#endif

    std::cout << " =========" << std::endl;
    std::cout << "box_min: " << box_min << ", box_max: " << box_max << std::endl;
    std::cout << "first: [" << first[0] << ", " << first[1] << "]" << std::endl;
    std::cout << "last:  [" <<  last[0] << ", " <<  last[1] << "]" << std::endl;

    if (dim == 1) {
      for (int i = first[0]; i <= last[0]; ++i) {
        cells.emplace_back(i);
      }
    } else if (dim == 2) {
      for (int j = first[1]; j <= last[1]; ++j) {
        for (int i = first[0]; i <= last[0]; ++i) {
          cells.emplace_back(i + j * num_sides_[0]);
        }
      }
    } else if (dim == 3) {
      for (int k = first[2]; k <= last[2]; ++k) {
        for (int j = first[1]; j <= last[1]; ++j) {
          for (int i = first[0]; i <= last[0]; ++i) {
            cells.emplace_back(i + j * num_sides_[0] + k * num_sides_[0] * num_sides_[1]);
          }
        }
      }
    }

    std::cout << "cells size: " << cells.size() << std::endl;

    // step 3: scan cells and check distance of each included source point
    std::vector<int> source_neighbors;
    for (int c : cells) {
      std::cout << "c: "<< c << ", bin size: "<< bucket_[c].size() << std::endl;
      for (int s : bucket_[c]) {
        bool contained = true;
        auto const& q = source_swarm_.get_particle_coordinates(s);
        for (int d = 0; d < dim; ++d) {
          if (std::abs(q[d] - p[d]) > radius_scale_ * h[d]) {
            contained = false;
            std::cout << "q[d]: " << q[d] <<", p[d]: "<< p[d] << ", extent: "<< radius_scale_ * h[d] << std::endl;
            break;
          }
        }
        if (contained) {
          source_neighbors.emplace_back(s);
        }
      }
    }
    std::sort(source_neighbors.begin(), source_neighbors.end());
    return source_neighbors;
  }

private:

  /**
   * @brief Check whether the given point is outside the grid.
   *
   * @param p: current point coordinates.
   * @return true if outside, false otherwise.
   */
  bool is_outside(Wonton::Point<dim> const& p) const {
    for (int d = 0; d < dim; ++d) {
      if (p[d] < p_min_[d] or p[d] > p_max_[d]) {
        return true;
      }
    }
    return false;
  }

  /**
   * @brief Deduce cell indices (i,j,k) from physical coordinates (x,y,z).
   *
   * @param p: current point coordinates.
   * @return index of the cell containing the point in helper grid.
   */
  std::array<int, dim> deduce_cell_index(Wonton::Point<dim> const& p) const {
//    assert(num_sides_ > 0);
    std::array<int, dim> indices;

    for (int d = 0; d < dim; ++d) {
      double const t = p[d] - p_min_[d];
      double const range = p_max_[d] - p_min_[d];
      assert(t >= 0 and t <= range);
//      std::cout << "p[d]: "<< p[d] <<", p_min[d]: "<< p_min_[d] << ", p_max[d]: " << p_max_[d] << std::endl;
      std::cout << "range: " << range << ", t: " << t << ", num sides: " << num_sides_[d] << std::endl;
      indices[d] = std::min(static_cast<int>(std::floor(t * num_sides_[d] / range)), num_sides_[d] - 1);
    }

    std::cout << "(i,j): ("<< indices[0] <<", "<< indices[1] << ")" << std::endl;
    return indices;
  }

  /**
   * @brief Deduce bin index i' from physical coordinates (x,y,z).
   *
   * @param p: current point coordinates.
   * @return index of the bin containing the point.
   */
  int deduce_bin_index(Wonton::Point<dim> const& p) const {
    // step 1: (x,y,z) to (i,j,k)
    auto const cell = deduce_cell_index(p);
    // step 2: (i,j,k) to i'
//    int index = cell[0];
//    for (int d = 1; d < dim; ++d) {
//      index += cell[d] * static_cast<int>(std::pow(num_sides_[d - 1], d));
//    }
    switch (dim) {
      case 1: return cell[0];
      case 2: return cell[0] + cell[1] * num_sides_[0];
      case 3: return cell[0] + cell[1] * num_sides_[0] + cell[2] * num_sides_[1] * num_sides_[1];
      default: return -1;
    }
  }

  /** source points */
  SourceSwarm const& source_swarm_;
  /** target points */
  TargetSwarm const& target_swarm_;
  /** search radius for each target point */
  Wonton::vector<Wonton::Point<dim>> const& target_radius_;
  /** search radius scaling factor */
  double radius_scale_ = 1;
  /** bounding box of helper grid */
  Wonton::Point<dim> p_min_, p_max_;
  /** bins of source points within helper grid */
  std::vector<std::vector<int>> bucket_;
  /** number of sides of helper grid */
  int num_sides_[dim];
};

} // namespace Portage

#endif