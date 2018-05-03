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
#include "portage/support/Point.h"
#include "portage/accumulate/accumulate.h"

#include "pile.hh"
#include "lretypes.hh"
#include "pairs.hh"

namespace Portage {

  /*!
  @class SearchPointsByCells "search_points_by_cells.h"
  @brief A simple, crude search algorithm that does a linear-time search
  over a swarm of points.
  @tparam SourceSwarmType The swarm type of the input swarm.
  @tparam TargetSwarmType The swarm type of the output swarm.
   */
template <int D, class SourceSwarmType, class TargetSwarmType>
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
  SearchPointsByCells(
      const SourceSwarmType & source_swarm,
      const TargetSwarmType & target_swarm,
      std::shared_ptr<vector<Point<D>>> source_extents,
      std::shared_ptr<vector<Point<D>>> target_extents,
      Meshfree::WeightCenter center=Meshfree::Scatter)
      : sourceSwarm_(source_swarm), targetSwarm_(target_swarm),
        sourceExtents_(source_extents), targetExtents_(target_extents),
        pair_finder_(NULL), center_(center)
  {
    // check sizes
    if (center == Meshfree::Scatter) {
      assert(sourceSwarm_.num_particles() == sourceExtents_->size());
    } else if (center == Meshfree::Gather) {
      assert(targetSwarm_.num_particles() == targetExtents_->size());
    }

    // transpose geometry data to lre namespace structures
    Meshfree::Pairs::vpile source_vp(D, sourceSwarm_.num_particles());
    Meshfree::Pairs::vpile target_vp(D, targetSwarm_.num_particles());
    Meshfree::Pairs::vpile source_extents_vp(D, sourceSwarm_.num_particles());
    Meshfree::Pairs::vpile target_extents_vp(D, targetSwarm_.num_particles());
    for (size_t i=0; i<sourceSwarm_.num_particles(); i++ ) {
      Point<D> pt = sourceSwarm_.get_particle_coordinates(i);
      for (size_t m=0; m<D; m++) {
	source_vp[m][i] = pt[m];
      }
    }
    for (size_t i=0; i<targetSwarm_.num_particles(); i++ ) {
      Point<D> pt = targetSwarm_.get_particle_coordinates(i);
      for (size_t m=0; m<D; m++) {
	target_vp[m][i] = pt[m];
      }
    }
    if (center == Meshfree::Scatter) {
      for (size_t i=0; i<sourceSwarm_.num_particles(); i++ ) {
	for (size_t m=0; m<D; m++) {
	  Point<D> pt = (*sourceExtents_)[i];
	  source_extents_vp[m][i] = pt[m];
	  (*sourceExtents_)[i] = pt;
	}
      }
    } else if (center_ == Meshfree::Gather) {
      for (size_t i=0; i<targetSwarm_.num_particles(); i++ ) {
	for (size_t m=0; m<D; m++) {
	  Point<D> pt = (*targetExtents_)[i];
	  target_extents_vp[m][i] = pt[m];
	  (*targetExtents_)[i] = pt;
	}
      }
    }

    // h on source
    if (center_ == Meshfree::Scatter) {
      pair_finder_ = std::make_shared<Meshfree::Pairs::CellPairFinder>(
          source_vp, target_vp, source_extents_vp, true);

    // h on target
    } else if (center_ == Meshfree::Gather) {
      pair_finder_ = std::make_shared<Meshfree::Pairs::CellPairFinder>(
          source_vp, target_vp, target_extents_vp, false);
    }
  } // SearchPointsByCells::SearchPointsByCells

  //! Copy constructor - use default - std::transform needs this
  //  SearchPointsByCells(const SearchPointsByCells &) = delete;

  //! Assignment operator (disabled)
  SearchPointsByCells & operator = (const SearchPointsByCells &) = delete;

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
  std::vector<unsigned int> operator() (const int pointId) const;

  private:

  // Aggregate data members
  const SourceSwarmType & sourceSwarm_;
  const TargetSwarmType & targetSwarm_;
  std::shared_ptr<vector<Point<D>>> sourceExtents_;
  std::shared_ptr<vector<Point<D>>> targetExtents_;
  std::shared_ptr<Meshfree::Pairs::CellPairFinder> pair_finder_;
  Meshfree::WeightCenter center_;

}; // class SearchPointsByCells


template<int D, class SourceSwarmType, class TargetSwarmType>
std::vector<unsigned int>
SearchPointsByCells<D, SourceSwarmType, TargetSwarmType>::
operator() (const int pointId) const {

  auto result = pair_finder_->find(pointId);
  return std::vector<unsigned int>(result.begin(), result.end());

} // SearchPointsByCells::operator()


} // namespace Portage

#endif // SEARCH_POINTS_BY_CELLS_H

