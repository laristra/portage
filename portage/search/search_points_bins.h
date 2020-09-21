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

/**
 * @class SearchPointsBins.
 * @brief A lightweight linear-time search kernel for particles based on
 *        source points binning using a helper grid. It is only valid
 *        for gather weight forms.
 *
 * @tparam dim: dimension of source/target points.
 * @tparam SourceSwarm: source point cloud type.
 * @tparam TargetSwarm: target point cloud type.
 */
template<int dim, class SourceSwarm, class TargetSwarm>
class SearchPointsBins {

public:
  /**
   * @brief Default constructor (disabled).
   */
  SearchPointsBins() = delete;

  /**
   * @brief Create an instance of the search kernel.
   *
   * @param source_swarm: the source points.
   * @param target_swarm: the target points.
   * @param source_extents: search radius of each source point (unused).
   * @param target_radius: search radius of each target point.
   * @param center: weight form (gather only)
   * @param radius_scale: scale factor for search radius.
   */
  SearchPointsBins(SourceSwarm const& source_swarm,
                   TargetSwarm const& target_swarm,
                   Wonton::vector<Wonton::Point<dim>> const& /* unused */,
                   Wonton::vector<Wonton::Point<dim>> const& target_radius,
                   WeightCenter const center = Gather,
                   double radius_scale = 1.)
    : source_swarm_(source_swarm),
      target_swarm_(target_swarm),
      search_radius_(target_radius),
      radius_scale_(radius_scale) {

    static_assert(dim > 0 and dim < 4, "invalid dimension");

    if (center == Scatter)
      throw std::runtime_error("scatter weight form not supported");

    int const num_source_points = source_swarm.num_particles();
    int const num_target_points = target_swarm.num_particles();
    double radius[dim];
    Wonton::Point<dim> p_min, p_max;

    // step 1: construct helper grid from target points
    for (int d = 0; d < dim; ++d) {
      p_min[d] = std::numeric_limits<double>::max();
      p_max[d] = std::numeric_limits<double>::lowest();
      radius[d] = 0.;
      sides_[d] = 0;
    }

    for (int t = 0; t < num_target_points; ++t) {
      auto const& p = target_swarm.get_particle_coordinates(t);
      auto const h = radius_scale * Wonton::Point<dim>(target_radius[t]);
      for (int d = 0; d < dim; ++d) {
        p_min[d] = std::min(p_min[d], p[d] - h[d]);
        p_max[d] = std::max(p_max[d], p[d] + h[d]);
        radius[d] += std::abs(h[d]);
      }
    }

    // step 2: deduce number of bins from search radii
    int const max_sides = get_max_sides(num_target_points);
    int num_bins = 1;

    for (int d = 0; d < dim; ++d) {
      orig_[d] = p_min[d];
      span_[d] = p_max[d] - p_min[d];
      double const step = 2 * radius[d] / num_target_points;
      sides_[d] = std::min(static_cast<int>(span_[d] / step), max_sides);
      num_bins *= sides_[d];
    }

    // step 3: push source points into bins
    bucket_.resize(num_bins);

    for (int s = 0; s < num_source_points; ++s) {
      auto const& p = source_swarm.get_particle_coordinates(s);
      bool inside = true;
      for (int d = 0; d < dim; ++d) {
        if (p[d] < p_min[d] or p[d] > p_max[d]) {
          inside = false;
          break;
        }
      }
      if (inside) {
        int const i = deduce_bin_index(p);
        bucket_[i].emplace_back(s);
      }
    }
  }

  /**
   * @brief Retrieve source neighbors within the radius of a given target point.
   *
   * @param target_id: the given target point.
   * @return the filtered list of source neighbors.
   */
  std::vector<int> operator() (int id) const {

    auto const p = target_swarm_.get_particle_coordinates(id);
    auto const h = radius_scale_ * Wonton::Point<dim>(search_radius_[id]);

    // step 1: build bounding box of the radius of target point
    Wonton::Point<dim> box_min, box_max;
    for (int d = 0; d < dim; ++d) {
      box_min[d] = p[d] - h[d];
      box_max[d] = p[d] + h[d];
    }

    // step 2: filter cells overlapped by the bounding box
    std::vector<int> cells;
    auto const first = deduce_cell_index(box_min);
    auto const last  = deduce_cell_index(box_max);

#ifndef NDEBUG
    for (int d = 0; d < dim; ++d) {
      assert(first[d] <= last[d]);
    }
#endif

    if (dim == 1) {
      for (int i = first[0]; i <= last[0]; ++i) {
        cells.emplace_back(i);
      }
    } else if (dim == 2) {
      for (int j = first[1]; j <= last[1]; ++j) {
        for (int i = first[0]; i <= last[0]; ++i) {
          cells.emplace_back(i + j * sides_[0]);
        }
      }
    } else if (dim == 3) {
      for (int k = first[2]; k <= last[2]; ++k) {
        for (int j = first[1]; j <= last[1]; ++j) {
          for (int i = first[0]; i <= last[0]; ++i) {
            cells.emplace_back(i + j * sides_[0] + k * sides_[0] * sides_[1]);
          }
        }
      }
    }

    // step 3: scan cells and check distance of each included source point
    std::vector<int> neighbors;
    for (int c : cells) {
      for (int s : bucket_[c]) {
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
   * @brief Retrieve the maximum number of sides for the helper grid.
   *
   * @param num_points: the number of points of the swarm.
   * @return the maximum allowed number of sides for the grid.
   */
  int get_max_sides(int num_points) const {
    double const factor = 0.1;
    double const num_per_axis = std::pow(factor * num_points, 1./dim);
    return std::max(static_cast<int>(std::ceil(num_per_axis)), 1);
  }

  /**
   * @brief Deduce cell indices (i,j,k) from physical coordinates (x,y,z).
   *
   * @param p: current point coordinates.
   * @return index of the cell containing the point in helper grid.
   */
  std::array<int, dim> deduce_cell_index(Wonton::Point<dim> const& p) const {
    std::array<int, dim> cell;

    for (int d = 0; d < dim; ++d) {
      assert(sides_[d] > 0);
      double const shift = p[d] - orig_[d];
      assert(shift >= 0 and shift <= span_[d]);
      cell[d] = std::min(static_cast<int>(std::floor(shift * sides_[d] / span_[d])), sides_[d] - 1);
    }

    return cell;
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
    switch (dim) {
      case 1: return cell[0];
      case 2: return cell[0] + cell[1] * sides_[0];
      case 3: return cell[0] + cell[1] * sides_[0] + cell[2] * sides_[0] * sides_[1];
      default: return -1;
    }
  }

  SourceSwarm const& source_swarm_;                         /** source points */
  TargetSwarm const& target_swarm_;                         /** target points */
  Wonton::vector<Wonton::Point<dim>> const& search_radius_; /** search radius per target point */
  double radius_scale_ = 1.;                                /** radii scale factor */
  Wonton::Point<dim> orig_, span_;                          /** helper grid */
  std::vector<std::vector<int>> bucket_;                    /** source points bins */
  int sides_[dim] {};                                       /** sides per axis */
};

} // namespace Portage

#endif