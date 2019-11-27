/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef SEARCH_KDTREE_NANOFLANN_H
#define SEARCH_KDTREE_NANOFLANN_H

#include <memory>
#include <vector>
#include <utility>
#include <cmath>

#ifdef HAVE_NANOFLANN
#include "nanoflann.hpp"  // KDTree search package NanoFlann

#include "portage/support/portage.h"

#include "portage/accumulate/accumulate.h"  // For Meshfree::WeightCenter
//                                          // Hope we can get rid of it


namespace Portage {
template <size_t D, class SwarmType>
class SwarmNanoflann {
 public:
  SwarmNanoflann(SwarmType const& swarm) : swarm_(swarm) {}

  // Nanoflann + C++ is requiring us to provide the assignment operator
  SwarmNanoflann & operator=(SwarmNanoflann const& inswarm) {
    swarm_ = inswarm.swarm_;
    return *this;
  }

  //
  // METHODS REQUIRED BY NANOFLANN
  //

  // How many points in the swarm?
  inline size_t kdtree_get_point_count() const {
    return swarm_.num_particles(Entity_type::ALL);
  }

  // Square of Euclidean distance between query point and a point in the swarm
  inline double kdtree_distance(double const *p1, size_t const idx, size_t /* size */) const {    
    double dist2 = 0.0;
    Point<D> const& p2 = swarm_.get_particle_coordinates(idx);
    for (int i = 0; i < D; i++)
      dist2 += (p1[i]-p2[i])*(p1[i]-p2[i]);
    return dist2;
  }

  // get the d'th coordinate of a point 
  inline double kdtree_get_pt(size_t const idx, int d) const {
    assert(d >= 0 && d < D);
    Point<D> const& p = swarm_.get_particle_coordinates(idx);
    return p[d];
  }

  // Optional bounding-box computation; return false to default to a
  // standard bbox computation loop
  template<class BBOX>
  bool kdtree_get_bbox(BBOX& /*bb*/) const {return false;}

  
  //
  // MORE EFFICIENT METHOD FOR GETTING POINT COORDINATES FOR OUR USE
  //

  Point<D> get_point(size_t const idx) const {
    return swarm_.get_particle_coordinates(idx);
  }

 private:
  SwarmType const& swarm_;

};  // class SwarmNanoflann



/*!
  @class Search_KDTree_Nanoflann
  @brief A kd-tree based search algorithm for nearest neighbors based on the 'nanoflann' package
  @tparam SourceSwarmType The swarm type of the input swarm.
  @tparam TargetSwarmType The swarm type of the output swarm.
*/
template <int D, class SourceSwarmType, class TargetSwarmType>
class Search_KDTree_Nanoflann {
 public:

  // Particular kd-tree type based on nanoflann's offerings
  typedef
  nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<double, SwarmNanoflann<D, SourceSwarmType>>,
   SwarmNanoflann<D, SourceSwarmType>,
   D
   > kdtree_t;
  
  //! Default constructor (disabled)
  Search_KDTree_Nanoflann() = delete;

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
  Search_KDTree_Nanoflann(SourceSwarmType const& source_swarm,
                          TargetSwarmType const& target_swarm,
                          std::shared_ptr<std::vector<Point<D>>> source_extents,
                          std::shared_ptr<std::vector<Point<D>>> target_extents,
                          Meshfree::WeightCenter center = Meshfree::Gather) :
      sourceSwarm_(source_swarm), targetSwarm_(target_swarm) {
    
    assert(center == Meshfree::Gather);  // Our KDTree is built on source points

    int ntgt = target_swarm.num_owned_particles();
    search_radii_.resize(ntgt);
    for (int i = 0; i < ntgt; i++) {
      double len = 0.0;
      for (int d = 0; d < D; d++)
        len += (*target_extents)[i][d]*(*target_extents)[i][d];
      search_radii_[i] = sqrt(len);
    }

    kdtree_ = std::make_shared<kdtree_t>(D, sourceSwarm_,
                nanoflann::KDTreeSingleIndexAdaptorParams(10 /* maxleaf */));
    kdtree_->buildIndex();
  }

  /*!
    @brief Find the source swarm points within an appropriate distance
    of a target point.
    @param[in] pointId The index of the point in the target swarm for
    which we wish to find the candidate neighbor points in the source
    swarm.
    @param[in,out] candidates Pointer to a vector of potential candidate
    points in the source swarm.
  */
  std::vector<unsigned int> operator() (const size_t pointId) const;

 private:
  SwarmNanoflann<D, SourceSwarmType> const sourceSwarm_;
  SwarmNanoflann<D, TargetSwarmType> const targetSwarm_;
  std::vector<double> search_radii_;
  std::shared_ptr<kdtree_t> kdtree_ = nullptr;
};  // class Search_KDTree_Nanoflann


template<int D, class SourceSwarmType, class TargetSwarmType>
std::vector<unsigned int>
Search_KDTree_Nanoflann<D, SourceSwarmType, TargetSwarmType>::
operator() (const size_t pointId) const {
  std::vector<unsigned int> candidates;

  // find coordinates of target point
  Point<D> tpcoord = targetSwarm_.get_point(pointId);

  double p[D];
  for (int i = 0; i < D; i++)
    p[i] = tpcoord[i];

  std::vector<std::pair<size_t, double>> matches;

  nanoflann::SearchParams params;
  // params.sorted = false;        // one could make this true

  double radius = search_radii_[pointId];
  const size_t nMatches = kdtree_->radiusSearch(p, radius, matches, params);

  for (auto const &idx_dist_pair : matches)
    candidates.push_back(static_cast<int>(idx_dist_pair.first));

  return candidates;
}  // Search_KDTree_Nanoflann::operator()

}  // namespace Portage

#endif  // HAVE_NANOFLANN


#endif  // SEARCH_KDTREE_NANOFLANN_H
