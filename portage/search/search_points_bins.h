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
 *
 * @brief A lightweight search algorithm in linearithmic time for particles.
 *
 * Principle:
 * It retrieves the source neighbors of each target point based on its
 * smoothing lengths. It is suited only for gather weight forms.
 *
 * Details:
 * During the initialization step, it creates a cartesian grid that encloses
 * the bounding boxes of all target points based on their smoothing lengths.
 * After that, it puts the source points into the grid cells that contain them
 * using coordinates hashing. To find the source neighbors of a target point,
 * it retrieves the cells that overlap the bounding box enclosing the target
 * point, and all points included in those cells are retrieved and filtered.
 *
 * Performance considerations:
 * - restrict to cartesian grid with fixed cell size per dimension (for GPU).
 * - disable spatial reordering: it is not a k-nearest neighbor search.
 * - minimal footprint and data operations: no internal data recopy.
 * - direct queries: no strides, no cell boundary check.
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
   * It constructs the helper grid based on the smoothing lengths
   * defined on target points. It determines the spatial extents of
   * the search by creating a bounding box enclosing the search
   * radii of all target points. It discretizes this bounding box
   * by deducing the number of edges per axis from the total radius
   * for that axis. Hence the number of bins is simply given by the
   * product of the number of edges. It finally puts each source
   * point to the corresponding bin by hashing its coordinates:
   *  (x,y,z) -> (i,j,k) -> i'.
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

    /*
     * verify that the correct weight form is supplied.
     * it only supports gather weight form since it filters
     * the source neighbors included in the support of each target point.
     */
    if (center == Scatter)
      throw std::runtime_error("scatter weight form not supported");

    int const num_source_points = source_swarm.num_particles();
    int const num_target_points = target_swarm.num_particles();
    double radius[dim];
    Wonton::Point<dim> p_min, p_max;

    /* ------------------------------------------------------------
     *  step 1: create search bounding box from target points.
     * ------------------------------------------------------------
     * It determines the spatial extents per axis by creating a
     * bounding box that encloses the support of each target point.
     * It deduces the total search radius per axis: sum_i=1^n |h_i|.
     */
    for (int d = 0; d < dim; ++d) {
      p_min[d] = std::numeric_limits<double>::max();
      p_max[d] = std::numeric_limits<double>::lowest();
      radius[d] = 0.;
      num_edges_[d] = 0;
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

    /* -------------------------------------------------------------------
     *  step 2: discretize it and deduce number of bins.
     * -------------------------------------------------------------------
     * It discretizes the bounding box by computing the spatial step 'ds'
     * which is a scaled mean radius per axis: ds = (r/n) sum_i=1^n |h_i|.
     * For a given axis, it determines the number of edges by dividing
     * the search radius by the spatial step: m = (p_max - p_min) / ds,
     * provided that it does not exceed the maximum number of edges
     * defaulted to max[1, ceil(pow(r n, 1/d))]. The number of bins is
     * simply deduced by the product of the number of edges per axis:
     * n_bins = m[0].m[1].m[2].
     * The scale factor 'r' is chosen for compatibility reasons with the
     * legacy search algorithm.
     */
    int const num_max_edges = get_max_edges(num_target_points);
    int num_bins = 1;

    for (int d = 0; d < dim; ++d) {
      orig_[d] = p_min[d];
      span_[d] = p_max[d] - p_min[d];
      double const step = 2 * radius[d] / num_target_points;
      num_edges_[d] = std::min(static_cast<int>(span_[d] / step), num_max_edges);
      num_bins *= num_edges_[d];
    }

    /* --------------------------------------------------------------
     *  step 3: filter source points and push them to the bins
     * --------------------------------------------------------------
     * After resizing the container, it puts each source point to the
     * correct bin by hashing its coordinates while discarding those
     * outside the helper grid extents.
     */
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
   * @brief Filter source neighbors within the radius of a given target point.
   *
   * It queries the source neighbors of a given target point using the
   * prebuilt helper grid. It first creates a bounding box that encloses the
   * local search area of the target point. It then filters the grid cells
   * that overlap this local bounding box. Afterwards, it scans those cells
   * while verifying that each included source point is actually within the
   * support of the given target point.
   *
   * @param target_id: the given target point.
   * @return the filtered list of source neighbors.
   */
  std::vector<int> operator() (int id) const {

    auto const p = target_swarm_.get_particle_coordinates(id);
    auto const h = radius_scale_ * Wonton::Point<dim>(search_radius_[id]);

    /* --------------------------------------------------------------
     *  step 1: build local bounding box centered at target point
     * --------------------------------------------------------------
     * It determines the extents of the local search area based on
     * the support of the target point.
     */
    Wonton::Point<dim> box_min, box_max;
    for (int d = 0; d < dim; ++d) {
      box_min[d] = p[d] - h[d];
      box_max[d] = p[d] + h[d];
    }

    /* --------------------------------------------------------------
     *  step 2: filter cells overlapped by the bounding box
     * --------------------------------------------------------------
     * It retrieves the cells of the helper grid that overlap with
     * the local bounding box. It first determines the indices of the
     * cells that overlap with the box extents. After that, it will
     * retrieve all cells included in that index range for each axis
     * in a local cell list.
     * Note that the helper grid cells are virtually stored in a flat
     * array: each cell can be accessed using the cached number of
     * edges as a stride.
     */
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
          cells.emplace_back(i + j * num_edges_[0]);
        }
      }
    } else if (dim == 3) {
      for (int k = first[2]; k <= last[2]; ++k) {
        for (int j = first[1]; j <= last[1]; ++j) {
          for (int i = first[0]; i <= last[0]; ++i) {
            cells.emplace_back(i + j * num_edges_[0] + k * num_edges_[0] * num_edges_[1]);
          }
        }
      }
    }

    /* --------------------------------------------------------------
     *  step 3: scan retrieved cells and verify each source point.
     * --------------------------------------------------------------
     * It finally scans the cell list and verify that each included
     * source point is actually inside the local search area of the
     * target point. If so, it will be added the filtered list of
     * neighbors.
     */
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
   * @brief Retrieve the maximum number of edges per axis for the helper grid.
   *
   * It defines the default maximum number of edges to limit the
   * number of cells of the helper grid. It corresponds to the
   * d^th root of the number of target points scaled by a factor.
   * The scaling factor is chosen for compatibility reasons with
   * the legacy search algorithm.
   *
   * @param num_points: the number of points of the swarm.
   * @return the maximum allowed number of edges for the grid.
   */
  int get_max_edges(int num_points) const {
    double const factor = 0.1;
    double const num_per_axis = std::pow(factor * num_points, 1./dim);
    return std::max(static_cast<int>(std::ceil(num_per_axis)), 1);
  }

  /**
   * @brief Deduce cell indices (i,j,k) from physical coordinates (x,y,z).
   *
   * It computes the index of the helper grid cell that contains the given
   * point. For each axis, it computes the relative position of the point
   * from the helper grid origin. The index of cell is then given by:
   * i =  ns * (p - p_min) / (p_max - p_min). After that, we just ensure
   * that the resulting index does not exceed the maximum for each axis.
   *
   * @param p: current point coordinates.
   * @return index of the cell containing the point in helper grid.
   */
  std::array<int, dim> deduce_cell_index(Wonton::Point<dim> const& p) const {
    std::array<int, dim> cell;

    for (int d = 0; d < dim; ++d) {
      assert(num_edges_[d] > 0);
      double const shift = p[d] - orig_[d];
      assert(shift >= 0 and shift <= span_[d]);
      cell[d] = std::min(static_cast<int>(std::floor(shift * num_edges_[d] / span_[d])), num_edges_[d] - 1);
    }

    return cell;
  }

  /**
   * @brief Deduce bin index i' from physical coordinates (x,y,z).
   *
   * It retrieve the index of bin that should contain the given source point
   * and is done in two steps. First, it computes the index (i,j,k) of the
   * helper grid cell that contains the given source point. It then deduces
   * the bin index i' in the flat array using the number of edges per axis
   * as strides: i' = i + j * ns[0] + k * ns[0] * ns[1].
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
      case 2: return cell[0] + cell[1] * num_edges_[0];
      case 3: return cell[0] + cell[1] * num_edges_[0] + cell[2] * num_edges_[0] * num_edges_[1];
      default: return -1;
    }
  }

  /** reference to source points */
  SourceSwarm const& source_swarm_;
  /** reference to target points */
  TargetSwarm const& target_swarm_;
  /** reference to search radii per target point */
  Wonton::vector<Wonton::Point<dim>> const& search_radius_;
  /** search radii scale factor */
  double radius_scale_ = 1.;
  /** helper grid extents */
  Wonton::Point<dim> orig_, span_;
  /** source points bins */
  std::vector<std::vector<int>> bucket_;
  /** number of edges per axis */
  int num_edges_[dim] {};
};

} // namespace Portage

#endif