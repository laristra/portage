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


namespace Portage {

  /*!
  @class SearchSimplePoints "search_simple_points.h"
  @brief A simple, crude search algorithm that does a quadratic-time search
  over a swarm of points.
  @tparam SourceSwarmType The swarm type of the input swarm.
  @tparam TargetSwarmType The swarm type of the output swarm.
   */
template <int D, class SourceSwarmType, class TargetSwarmType>
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
      const SourceSwarmType & source_swarm,
      const TargetSwarmType & target_swarm,
      std::shared_ptr<vector<Point<D>>> source_extents,
      std::shared_ptr<vector<Point<D>>> target_extents,
      Meshfree::WeightCenter center=Meshfree::Scatter)
      : sourceSwarm_(source_swarm), targetSwarm_(target_swarm),
        sourceExtents_(source_extents), targetExtents_(target_extents),
        center_(center)
  {

    // currently no structure, just save the swarms and extents

  } // SearchSimplePoints::SearchSimplePoints

  //! Copy constructor - use default - std::transfor needs this
  //  SearchSimplePoints(const SearchSimplePoints &) = delete;

  //! Assignment operator (disabled)
  SearchSimplePoints & operator = (const SearchSimplePoints &) = delete;

  //! Destructor
  ~SearchSimplePoints() = default;

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
  Meshfree::WeightCenter center_;

}; // class SearchSimplePoints


template<int D, class SourceSwarmType, class TargetSwarmType>
std::vector<unsigned int>
SearchSimplePoints<D, SourceSwarmType, TargetSwarmType>::
operator() (const int pointId) const {

  using std::abs;
  std::vector<unsigned int> candidates;

  // find coordinates of target point
  Point<D> tpcoord = targetSwarm_.get_particle_coordinates(pointId);

  // now see which source points are within an appropriate distance
  // of the target point
  // do a naive linear search
  const int numPoints = sourceSwarm_.num_particles(Entity_type::ALL);
  for (int p = 0; p < numPoints; ++p) {
    Point<D> spcoord = sourceSwarm_.get_particle_coordinates(p);
    bool contained = true;
    Point<D> spt=(*sourceExtents_)[p];
    Point<D> tpt=(*targetExtents_)[pointId];
    for (int d = 0; d < D; ++d) {
      double maxdist;
      if (center_ == Meshfree::Scatter) {
        maxdist = 2.*spt[d];
      } else if (center_ == Meshfree::Gather) {
        maxdist = 2.*tpt[d];
      }
      contained = contained && (abs(tpcoord[d] - spcoord[d]) < maxdist);
      if (!contained) break;
    }
    if (contained) {
      candidates.push_back(p);
    }
  }

  return candidates;
}  // SearchSimplePoints::operator()


}  // namespace Portage

#endif  // PORTAGE_SEARCH_SEARCH_SIMPLE_POINTS_H_
