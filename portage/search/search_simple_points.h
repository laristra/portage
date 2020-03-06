/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef PORTAGE_SEARCH_SEARCH_SIMPLE_POINTS_H_
#define PORTAGE_SEARCH_SEARCH_SIMPLE_POINTS_H_

#include <memory>
#include <vector>
#include <cmath>


#include "portage/support/portage.h"
#include "portage/accumulate/accumulate.h"
#include "wonton/support/Point.h"

using Wonton::Point;

namespace Portage {

  /*!
  @class SearchSimplePoints "search_simple_points.h"
  @brief A simple, crude search algorithm that does a quadratic-time search
  over a swarm of points.
  @tparam SourceSwarm The swarm type of the input swarm.
  @tparam TargetSwarm The swarm type of the output swarm.
   */
template<int dim, class SourceSwarm, class TargetSwarm>
class SearchSimplePoints {
public:

  //! Default constructor (disabled)
  SearchSimplePoints() = delete;

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
  SearchSimplePoints(
    SourceSwarm const& source_swarm,
    TargetSwarm const& target_swarm,
    Portage::vector<Point<dim>> const& source_extents,
    Portage::vector<Point<dim>> const& target_extents,
    swarm::WeightCenter center = swarm::Scatter)
      : source_swarm_(source_swarm),
        target_swarm_(target_swarm),
        source_extents_(source_extents),
        target_extents_(target_extents),
        center_(center)
  {} // SearchSimplePoints::SearchSimplePoints

  //! Copy constructor - use default - std::transfor needs this
  //  SearchSimplePoints(const SearchSimplePoints &) = delete;

  //! Assignment operator (disabled)
  SearchSimplePoints& operator = (SearchSimplePoints const&) = delete;

  //! Destructor
  ~SearchSimplePoints() = default;

  /*!
    @brief Find the source swarm points within an appropriate distance
    of a target point.
    @param[in] point_id The index of the point in the target swarm for
    which we wish to find the candidate neighbor points in the source
    swarm.
    @param[in,out] candidates Pointer to a vector of potential candidate
    points in the source swarm.
  */
  std::vector<int> operator() (int point_id) const {

    std::vector<int> candidates;

    // find coordinates of target point
    Point<dim> target_coord = target_swarm_.get_particle_coordinates(point_id);

    // now see which source points are within an appropriate distance
    // of the target point
    // do a naive linear search
    int const nb_points = source_swarm_.num_particles(Wonton::ALL);
    for (int i = 0; i < nb_points; ++i) {
      Point<dim> source_coord = source_swarm_.get_particle_coordinates(i);
      bool contained = true;
      Point<dim> source_point_extent = source_extents_[i];
      Point<dim> target_point_extent = target_extents_[point_id];

      for (int d = 0; d < dim; ++d) {
        double maxdist = 0.;
        if (center_ == swarm::Scatter) {
          maxdist = 2. * source_point_extent[d];
        } else if (center_ == swarm::Gather) {
          maxdist = 2. * target_point_extent[d];
        }
        contained &= (std::abs(target_coord[d] - source_coord[d]) < maxdist);
        if (!contained)
          break;
      }
      if (contained)
        candidates.push_back(i);
    }

    return candidates;
  }

private:
  // Aggregate data members
  SourceSwarm const& source_swarm_;
  TargetSwarm const& target_swarm_;
  Portage::vector<Point<dim>> const& source_extents_;
  Portage::vector<Point<dim>> const& target_extents_;
  swarm::WeightCenter center_ = swarm::Scatter;

}; // class SearchSimplePoints

}  // namespace Portage

#endif  // PORTAGE_SEARCH_SEARCH_SIMPLE_POINTS_H_
