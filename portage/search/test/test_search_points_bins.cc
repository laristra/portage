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
#include "portage/search/search_points_by_cells.h"

using Portage::Meshfree::Gather;

class SearchBins : public ::testing::Test {
protected:
  /**
   * @brief Verify results on random 2D point clouds and compare to legacy kernel.
   *
   * @param num_source_points: number of source points.
   * @param num_target_points: number of target points.
   */
  static void test_random_2d(int num_source_points, int num_target_points) {

    using SearchCell = Portage::SearchPointsByCells<2, Wonton::Swarm<2>, Wonton::Swarm<2>>;
    using SearchBins = Portage::SearchPointsBins<2, Wonton::Swarm<2>, Wonton::Swarm<2>>;

    Wonton::vector<Wonton::Point<2>> source_points(num_source_points);
    Wonton::vector<Wonton::Point<2>> source_extent(num_source_points);
    Wonton::vector<Wonton::Point<2>> target_points(num_target_points);
    Wonton::vector<Wonton::Point<2>> target_extent(num_target_points);

    double const radius[] = { 1./std::sqrt(num_source_points),
                              1./std::sqrt(num_target_points) };

    std::random_device device;
    std::mt19937 engine { device() };
    std::uniform_real_distribution<double> generator(0.0, 1.0);

    for (int j = 0; j < num_source_points; ++j) {
      double const x = 1.2 * generator(engine) - 0.1;
      double const y = 1.2 * generator(engine) - 0.1;
      source_points[j] = { x, y };
      source_extent[j] = { radius[0], radius[0] };
    }

    for (int j = 0; j < num_target_points; ++j) {
      target_points[j] = { generator(engine), generator(engine) };
      target_extent[j] = { radius[1], radius[1] };
    }

    Wonton::Swarm<2> source_swarm(source_points);
    Wonton::Swarm<2> target_swarm(target_points);

    SearchCell oldsearch(source_swarm, target_swarm, source_extent, target_extent, Gather);
    SearchBins newsearch(source_swarm, target_swarm, source_extent, target_extent, Gather, 2);

    for (int i = 0; i < num_target_points; i++) {
      auto neighbors = std::make_pair(oldsearch(i), newsearch(i));
      std::sort(neighbors.first.begin(), neighbors.first.end());
      std::sort(neighbors.second.begin(), neighbors.second.end());

      ASSERT_EQ(neighbors.first.size(), neighbors.second.size());
      int const num_neigh = neighbors.first.size();
      for (int j =  0; j < num_neigh; j++) {
        ASSERT_EQ(neighbors.first[j], neighbors.second[j]);
      }
    }
  }

  /**
   * @brief Verify results on random 3D point clouds and compare to legacy kernel.
   *
   * @param num_source_points: number of source points.
   * @param num_target_points: number of target points.
   */
  static void test_random_3d(int num_source_points, int num_target_points) {

    using SearchCell = Portage::SearchPointsByCells<3, Wonton::Swarm<3>, Wonton::Swarm<3>>;
    using SearchBins = Portage::SearchPointsBins<3, Wonton::Swarm<3>, Wonton::Swarm<3>>;

    Wonton::vector<Wonton::Point<3>> source_points(num_source_points);
    Wonton::vector<Wonton::Point<3>> source_extent(num_source_points);
    Wonton::vector<Wonton::Point<3>> target_points(num_target_points);
    Wonton::vector<Wonton::Point<3>> target_extent(num_target_points);

    double const one_third = 1./3.;
    double const radius[] = { 1./std::pow(num_source_points, one_third),
                              1./std::pow(num_target_points, one_third) };

    std::random_device device;
    std::mt19937 engine { device() };
    std::uniform_real_distribution<double> generator(0.0, 1.0);

    for (int j = 0; j < num_source_points; ++j) {
      double const x = 1.2 * generator(engine) - 0.1;
      double const y = 1.2 * generator(engine) - 0.1;
      double const z = 1.2 * generator(engine) - 0.1;
      source_points[j] = { x, y, z };
      source_extent[j] = { radius[0], radius[0], radius[0] };
    }

    for (int j = 0; j < num_target_points; ++j) {
      target_points[j] = { generator(engine), generator(engine), generator(engine) };
      target_extent[j] = { radius[1], radius[1], radius[1] };
    }

    Wonton::Swarm<3> source_swarm(source_points);
    Wonton::Swarm<3> target_swarm(target_points);

    SearchCell oldsearch(source_swarm, target_swarm, source_extent, target_extent, Gather);
    SearchBins newsearch(source_swarm, target_swarm, source_extent, target_extent, Gather, 2);

    for (int i = 0; i < num_target_points; i++) {
      auto neighbors = std::make_pair(oldsearch(i), newsearch(i));
      std::sort(neighbors.first.begin(), neighbors.first.end());
      std::sort(neighbors.second.begin(), neighbors.second.end());

      ASSERT_EQ(neighbors.first.size(), neighbors.second.size());
      int const num_neigh = neighbors.first.size();
      for (int j =  0; j < num_neigh; j++) {
        ASSERT_EQ(neighbors.first[j], neighbors.second[j]);
      }
    }
  }
};

TEST_F(SearchBins, random_2d_case1) { test_random_2d(128, 128); }

TEST_F(SearchBins, random_2d_case2) { test_random_2d(256, 128); }

TEST_F(SearchBins, random_2d_case3) { test_random_2d(128, 256); }

TEST_F(SearchBins, random_3d_case1) { test_random_3d(1000, 1000); }

TEST_F(SearchBins, random_3d_case2) { test_random_3d(1000, 729); }

TEST_F(SearchBins, random_3d_case3) { test_random_3d(729, 1000); }

TEST_F(SearchBins, regular_2d) {

  using Search = Portage::SearchPointsBins<2, Wonton::Swarm<2>, Wonton::Swarm<2>>;

  // overlay a 3x3 target swarm on a 4x4 source swarm,
  // each target point should have four source neighbors
  Wonton::vector<Wonton::Point<2>> source_points(16);
  Wonton::vector<Wonton::Point<2>> source_extent(16);
  Wonton::vector<Wonton::Point<2>> target_points(9);
  Wonton::vector<Wonton::Point<2>> target_extent(9);

  double const radius = 0.375;

  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      source_points[i + j * 4] = { i + 0.5, j + 0.5 };
      source_extent[i + j * 4] = { radius, radius };
    }
  }

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      target_points[i + j * 3] = { i + 1., j + 1. };
      target_extent[i + j * 3] = { radius, radius };
    }
  }

  Wonton::Swarm<2> source_swarm(source_points);
  Wonton::Swarm<2> target_swarm(target_points);

  Search search(source_swarm, target_swarm, source_extent, target_extent, Gather, 2);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      int const target_point = i + j * 3;
      auto neighbors = search(target_point);
      std::sort(neighbors.begin(), neighbors.end());

      // there should be four candidate source points, in a square
      ASSERT_EQ(unsigned(4), neighbors.size());
      // compute source_point = index of lower left source point
      int const source_point = i + j * 4;
      ASSERT_EQ(source_point, neighbors[0]);
      ASSERT_EQ(source_point + 1, neighbors[1]);
      ASSERT_EQ(source_point + 4, neighbors[2]);
      ASSERT_EQ(source_point + 5, neighbors[3]);
    }
  }
}

TEST_F(SearchBins, regular_3d) {

  using Search = Portage::SearchPointsBins<3, Wonton::Swarm<3>, Wonton::Swarm<3>>;

  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm,
  // each target point should have eight source neighbors.
  Wonton::vector<Wonton::Point<3>> source_points(27);
  Wonton::vector<Wonton::Point<3>> source_extent(27);
  Wonton::vector<Wonton::Point<3>> target_points(8);
  Wonton::vector<Wonton::Point<3>> target_extent(8);

  double const radius = 0.375;

  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        int const index = i + j * 3 + k * 9;
        source_points[index] = { i + 0.5, j + 0.5, k + 0.5 };
        source_extent[index] = { radius, radius, radius };
      }
    }
  }

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        int const index = i + j * 2 + k * 4;
        target_points[index] = { i + 1., j + 1., k + 1. };
        target_extent[index] = { radius, radius, radius };
      }
    }
  }

  Wonton::Swarm<3> source_swarm(source_points);
  Wonton::Swarm<3> target_swarm(target_points);

  Search search(source_swarm, target_swarm, source_extent, target_extent, Gather, 2);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        int const target_point = i + j * 2 + k * 4;
        auto neighbors = search(target_point);
        std::sort(neighbors.begin(), neighbors.end());
        // there should be eight candidate source points, in a cube
        ASSERT_EQ(unsigned(8), neighbors.size());

        // compute source_point_base = index of lower left source point
        int const source_point = i + j * 3 + k * 9;
        ASSERT_EQ(source_point + 0, neighbors[0]);
        ASSERT_EQ(source_point + 1, neighbors[1]);
        ASSERT_EQ(source_point + 3, neighbors[2]);
        ASSERT_EQ(source_point + 4, neighbors[3]);
        ASSERT_EQ(source_point + 9, neighbors[4]);
        ASSERT_EQ(source_point + 10, neighbors[5]);
        ASSERT_EQ(source_point + 12, neighbors[6]);
        ASSERT_EQ(source_point + 13, neighbors[7]);
      }
    }
  }
}