/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/swarm/swarm.h"

#include "portage/search/search_points_bins.h"

TEST(search_points_bins, gather_2d) {

  using Wonton::Swarm;

  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points
  Wonton::vector<Wonton::Point<2>> source_points(16);
  Wonton::vector<Wonton::Point<2>> source_extent(16);
  Wonton::vector<Wonton::Point<2>> target_points(9);
  Wonton::vector<Wonton::Point<2>> target_extent(9);

  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      double const x = (i + 0.5);
      double const y = (j + 0.5);
      double const ext = 0.375;
      source_points[i + j * 4] = Wonton::Point<2>(x, y);
      source_extent[i + j * 4] = Wonton::Point<2>(ext, ext);
    }
  }

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      double const x = (i + 1.0);
      double const y = (j + 1.0);
      double const ext = 0.375;
      target_points[i + j * 3] = Wonton::Point<2>(x, y);
      target_extent[i + j * 3] = Wonton::Point<2>(ext, ext);
    }
  }

  Swarm<2> source_swarm(source_points);
  Swarm<2> target_swarm(target_points);

  Portage::SearchPointsBins<2, Swarm<2>, Swarm<2>>
    search(source_swarm, target_swarm,
           source_extent, target_extent, Portage::Meshfree::Gather, 2);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int tp = i + j * 3;
      auto candidates = search(tp);

      // there should be four candidate source points, in a square
      ASSERT_EQ(unsigned(4), candidates.size());
      // compute source_point_base = index of lower left source point
      const int source_point_base = i + j * 4;
      ASSERT_EQ(source_point_base,     candidates[0]);
      ASSERT_EQ(source_point_base + 1, candidates[1]);
      ASSERT_EQ(source_point_base + 4, candidates[2]);
      ASSERT_EQ(source_point_base + 5, candidates[3]);
    }
  }

}