/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef SEARCH_POINTS_BY_CELLS_H
#define SEARCH_POINTS_BY_CELLS_H

#include <memory>
#include <vector>
#include <cmath>
#include <list>

#include "portage/support/portage.h"

#include "portage/accumulate/accumulate.h"

#include "pile.hh"
#include "lretypes.hh"
#include "pairs.hh"

using Wonton::Point;

namespace Portage {

  /*!
  @class SearchPointsByCells "search_points_by_cells.h"
  @brief A simple, crude search algorithm that does a linear-time search
  over a swarm of points.
  @tparam SourceSwarm The swarm type of the input swarm.
  @tparam TargetSwarm The swarm type of the output swarm.
   */
template <int dim, class SourceSwarm, class TargetSwarm>
class SearchPointsByCells {
public:

  //! Default constructor (disabled)
  SearchPointsByCells() = delete;

  // Constructor with swarms and extents
  /*!
    @brief Builds the search structure for finding neighboring points.
    @param[in] source_swarm Pointer to source swarm info.
    @param[in] target_swarm Pointer to target swarm info.
    @param[in] source_extents Array of extents for source particles.
    @param[in] target_extents Array of extents for target particles.

    Constructor for search structure for finding points from a source
    swarm that are near points in the target swarm.
  */
  SearchPointsByCells(SourceSwarm const& source_swarm,
                      TargetSwarm const& target_swarm,
                      Portage::vector<Point<dim>> const& source_extents,
                      Portage::vector<Point<dim>> const& target_extents,
                      Meshfree::WeightCenter center = Meshfree::Scatter)
    : source_swarm_(source_swarm),
      target_swarm_(target_swarm),
      source_extents_(source_extents),
      target_extents_(target_extents),
      center_(center) {

    using namespace Meshfree;

    int const nb_source = source_swarm_.num_particles();
    int const nb_target = target_swarm_.num_particles();
    bool const do_scatter = (center == Scatter);

    // check sizes
    if (do_scatter) {
      assert(nb_source == source_extents_->size());
    } else {
      assert(nb_target == target_extents_->size());
    }

    // transpose geometry data to lre namespace structures
    Pairs::vpile source_vp(dim, nb_source);
    Pairs::vpile target_vp(dim, nb_target);
    Pairs::vpile extents_vp(dim, do_scatter ? nb_source : nb_target);

    for (int i = 0; i < nb_source; i++ ) {
      auto p = source_swarm_.get_particle_coordinates(i);
      for (int m = 0; m < dim; m++) {
        source_vp[m][i] = p[m];
      }
    }

    for (int i = 0; i < nb_target; i++ ) {
      auto p = target_swarm_.get_particle_coordinates(i);
      for (int m = 0; m < dim; m++) {
        target_vp[m][i] = p[m];
      }
    }

    if (do_scatter) {
      for (int i = 0; i < nb_source; i++ ) {
        for (int m = 0; m < dim; m++) {
          Point<dim> p = source_extents_[i];
          extents_vp[m][i] = p[m];
          source_extents_[i] = p;
        }
      }
    } else if (center_ == Gather) {
      for (int i = 0; i < nb_target; i++ ) {
        for (int m = 0; m < dim; m++) {
          Point<dim> p = target_extents_[i];
          extents_vp[m][i] = p[m];
          target_extents_[i] = p;
        }
      }
    }
    // h on source for scatter, on target for gather
    pair_finder_.init(source_vp, target_vp, extents_vp, do_scatter);
  }

  //! Copy constructor - use default - std::transform needs this
  SearchPointsByCells(SearchPointsByCells const&) = default;

  //! Assignment operator (disabled)
  SearchPointsByCells & operator = (SearchPointsByCells const&) = delete;

  //! Destructor
  ~SearchPointsByCells() = default;

  /*!
    @brief Find the source swarm points within an appropriate distance
    of a target point.
    @param[in] pointId The index of the point in the target swarm for
    which we wish to find the candidate neighbor points in the source
    swarm.
    @param[in,out] candidates Pointer to a vector of potential candidate
    points in the source swarm.
  */
  std::vector<int> operator() (int pointId) const {
    auto result = pair_finder_.find(pointId);
    return std::vector<int>(result.begin(), result.end());
  }

private:
  // Aggregate data members
  SourceSwarm const& source_swarm_;
  TargetSwarm const& target_swarm_;
  Portage::vector<Point<dim>> source_extents_;
  Portage::vector<Point<dim>> target_extents_;
  Meshfree::Pairs::CellPairFinder pair_finder_;
  Meshfree::WeightCenter center_ = Meshfree::Scatter;

}; // class SearchPointsByCells

} // namespace Portage

#endif // SEARCH_POINTS_BY_CELLS_H

