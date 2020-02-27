/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "portage/support/portage.h"
#include "portage/search/search_simple_points.h"
#include "portage/swarm/swarm.h"
#include "portage/accumulate/accumulate.h"
#include "wonton/support/Point.h"

TEST(search_simple_points, scatter_2d) {

  using Portage::Meshfree::Swarm;

  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points
  Portage::vector<Wonton::Point<2>> source_points(16);
  Portage::vector<Wonton::Point<2>> source_extent(16);
  Portage::vector<Wonton::Point<2>> target_points(9);
  Portage::vector<Wonton::Point<2>> target_extent(9);

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

  Portage::SearchSimplePoints<2, Swarm<2>, Swarm<2>>
    search(source_swarm, target_swarm,
           source_extent, target_extent);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int target_point = i + j * 3;
      auto candidates = search(target_point);

      // there should be four candidate source points, in a square
      ASSERT_EQ(4, candidates.size());
      // compute source_point_base = index of lower left source point
      const int source_point_base = i + j * 4;
      ASSERT_EQ(source_point_base,     candidates[0]);
      ASSERT_EQ(source_point_base + 1, candidates[1]);
      ASSERT_EQ(source_point_base + 4, candidates[2]);
      ASSERT_EQ(source_point_base + 5, candidates[3]);
    }
  }

} // TEST(search_simple_points, scatter_2d)


TEST(search_simple_points, scatter_3d) {

  using Portage::Meshfree::Swarm;

  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
  // each target point should have eight candidate source points
  Portage::vector<Wonton::Point<3>> source_points(27);
  Portage::vector<Wonton::Point<3>> source_extent(27);
  Portage::vector<Wonton::Point<3>> target_points(8);
  Portage::vector<Wonton::Point<3>> target_extent(8);

  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double const x = (i + 0.5);
        double const y = (j + 0.5);
        double const z = (k + 0.5);
        double const ext = 0.375;
        int const index = i + 3 * j + 9 * k;
        source_points[index] = Wonton::Point<3>(x, y, z);
        source_extent[index] = Wonton::Point<3>(ext, ext, ext);
      }
    }
  }

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        double const x = (i + 1.0);
        double const y = (j + 1.0);
        double const z = (k + 1.0);
        double const ext = 0.375;
        int const index = i + 2 * j + 4 * k;
        target_points[index] = Wonton::Point<3>(x, y, z);
        target_extent[index] = Wonton::Point<3>(ext, ext, ext);
      }
    }
  }

  Swarm<3> source_swarm(source_points);
  Swarm<3> target_swarm(target_points);

  Portage::SearchSimplePoints<3, Swarm<3>, Swarm<3>>
    search(source_swarm, target_swarm,
           source_extent, target_extent);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        int const target_point_index = i + j * 2 + k * 4;
        auto candidates = search(target_point_index);
        // there should be eight candidate source points, in a cube
        ASSERT_EQ(8, candidates.size());

        // compute source_point_base = index of lower left source point
        int const source_point_base = i + j * 3 + k * 9;
        ASSERT_EQ(source_point_base + 0, candidates[0]);
        ASSERT_EQ(source_point_base + 1, candidates[1]);
        ASSERT_EQ(source_point_base + 3, candidates[2]);
        ASSERT_EQ(source_point_base + 4, candidates[3]);
        ASSERT_EQ(source_point_base + 9, candidates[4]);
        ASSERT_EQ(source_point_base +10, candidates[5]);
        ASSERT_EQ(source_point_base +12, candidates[6]);
        ASSERT_EQ(source_point_base +13, candidates[7]);
      }
    }
  }
} // TEST(search_simple_points, scatter_3d)


TEST(search_simple_points, gather_2d) {

  using Portage::Meshfree::Swarm;

  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points
  Portage::vector<Wonton::Point<2>> source_points(16);
  Portage::vector<Wonton::Point<2>> source_extent(16);
  Portage::vector<Wonton::Point<2>> target_points(9);
  Portage::vector<Wonton::Point<2>> target_extent(9);

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

  Portage::SearchSimplePoints<2, Swarm<2>, Swarm<2>>
    search(source_swarm, target_swarm,
           source_extent, target_extent, Portage::Meshfree::Gather);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int tp = i + j * 3;
      auto candidates = search(tp);

      // there should be four candidate source points, in a square
      ASSERT_EQ(4, candidates.size());
      // compute source_point_base = index of lower left source point
      const int source_point_base = i + j * 4;
      ASSERT_EQ(source_point_base,     candidates[0]);
      ASSERT_EQ(source_point_base + 1, candidates[1]);
      ASSERT_EQ(source_point_base + 4, candidates[2]);
      ASSERT_EQ(source_point_base + 5, candidates[3]);
    }
  }
} // TEST(search_simple_points, gather_2d)


TEST(search_simple_points, gather_3d) {

  using Portage::Meshfree::Swarm;

  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
  // each target point should have eight candidate source points
  Portage::vector<Wonton::Point<3>> source_points(27);
  Portage::vector<Wonton::Point<3>> source_extent(27);
  Portage::vector<Wonton::Point<3>> target_points(8);
  Portage::vector<Wonton::Point<3>> target_extent(8);

  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double const x = (i + 0.5);
        double const y = (j + 0.5);
        double const z = (k + 0.5);
        double const ext = 0.375;
        int const index = i + 3 * j + 9 * k;
        source_points[index] = Wonton::Point<3>(x, y, z);
        source_extent[index] = Wonton::Point<3>(ext, ext, ext);
      }
    }
  }

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        double const x = (i + 1.0);
        double const y = (j + 1.0);
        double const z = (k + 1.0);
        double const ext = 0.375;
        int const index = i + 2 * j + 4 * k;
        target_points[index] = Wonton::Point<3>(x, y, z);
        target_extent[index] = Wonton::Point<3>(ext, ext, ext);
      }
    }
  }

  Swarm<3> source_swarm(source_points);
  Swarm<3> target_swarm(target_points);

  Portage::SearchSimplePoints<3, Swarm<3>, Swarm<3>>
    search(source_swarm, target_swarm,
           source_extent, target_extent, Portage::Meshfree::Gather);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        int const target_point_index = i + j * 2 + k * 4;
        auto candidates = search(target_point_index);
        // there should be eight candidate source points, in a cube
        ASSERT_EQ(8, candidates.size());

        // compute source_point_base = index of lower left source point
        int const source_point_base = i + j * 3 + k * 9;
        ASSERT_EQ(source_point_base + 0, candidates[0]);
        ASSERT_EQ(source_point_base + 1, candidates[1]);
        ASSERT_EQ(source_point_base + 3, candidates[2]);
        ASSERT_EQ(source_point_base + 4, candidates[3]);
        ASSERT_EQ(source_point_base + 9, candidates[4]);
        ASSERT_EQ(source_point_base +10, candidates[5]);
        ASSERT_EQ(source_point_base +12, candidates[6]);
        ASSERT_EQ(source_point_base +13, candidates[7]);
      }
    }
  }
}  // TEST(search_simple_points, gather_3d)


