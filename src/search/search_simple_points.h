/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/


#ifndef SEARCH_SIMPLE_POINTS_H
#define SEARCH_SIMPLE_POINTS_H

#include <memory>
#include <vector>
#include <cmath>

#include "portage/support/Point.h"


namespace Portage {

  /*!
  @class SearchSimplePoints "search_simple_points.h"
  @brief A simple, crude search algorithm that does a linear search
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
      std::shared_ptr<std::vector<Point<D>>> source_extents,
      std::shared_ptr<std::vector<Point<D>>> target_extents)
      : sourceSwarm_(source_swarm), targetSwarm_(target_swarm),
        sourceExtents_(source_extents), targetExtents_(target_extents)  {

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
  std::shared_ptr<std::vector<Point<D>>> sourceExtents_;
  std::shared_ptr<std::vector<Point<D>>> targetExtents_;

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
  const int numPoints = sourceSwarm_.num_owned_particles();
  for (int p = 0; p < numPoints; ++p) {
    Point<D> spcoord = sourceSwarm_.get_particle_coordinates(p);
    bool overlap = true;
    for (int d = 0; d < D; ++d) {
      double maxdist = 2. *
          ((*targetExtents_)[pointId][d] + (*sourceExtents_)[p][d]);
      overlap = overlap && (abs(tpcoord[d] - spcoord[d]) < maxdist);
      if (!overlap) break;
    }
    if (overlap) {
      candidates.push_back(p);
    }
  }

  return candidates;
} // SearchSimplePoints::operator()


} // namespace Portage

#endif // SEARCH_SIMPLE_POINTS_H

