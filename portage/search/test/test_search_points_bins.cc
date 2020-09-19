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
#include "portage/search/search_simple_points.h"

using Wonton::Swarm;
using Portage::Meshfree::Gather;

class SearchBins : public ::testing::Test {
protected:

  static void test_regular_2d() {
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
             source_extent, target_extent, Gather, 2);

    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        const int tp = i + j * 3;
        auto candidates = search(tp);
        std::sort(candidates.begin(), candidates.end());

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

  static void test_regular_3d() {
    // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
    // each target point should have eight candidate source points
    Wonton::vector<Wonton::Point<3>> source_points(27);
    Wonton::vector<Wonton::Point<3>> source_extent(27);
    Wonton::vector<Wonton::Point<3>> target_points(8);
    Wonton::vector<Wonton::Point<3>> target_extent(8);

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

    Portage::SearchPointsBins<3, Swarm<3>, Swarm<3>>
      search(source_swarm, target_swarm,
             source_extent, target_extent, Gather, 2);

    for (int k = 0; k < 2; ++k) {
      for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
          int const target_point_index = i + j * 2 + k * 4;
          auto candidates = search(target_point_index);
          std::sort(candidates.begin(), candidates.end());
          // there should be eight candidate source points, in a cube
          ASSERT_EQ(unsigned(8), candidates.size());

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
  }

  static void test_random_2d(int num_source_points, int num_target_points) {

    // random point sets and test against SearchSimplePoints
    Wonton::vector<Wonton::Point<2>> source_points(num_source_points);
    Wonton::vector<Wonton::Point<2>> source_extent(num_source_points);
    Wonton::vector<Wonton::Point<2>> target_points(num_target_points);
    Wonton::vector<Wonton::Point<2>> target_extent(num_target_points);

    for (int j = 0; j < num_source_points; ++j) {
      double x = 1.2 * rand()/RAND_MAX-.1;
      double y = 1.2 * rand()/RAND_MAX-.1;
      double ext = 1./sqrt(num_source_points * 1.0);
      source_points[j] = Wonton::Point<2>(x, y);
      source_extent[j] = Wonton::Point<2>(ext, ext);
    }

    for (int j = 0; j < num_target_points; ++j) {
      double x = 1.0 * rand()/RAND_MAX;
      double y = 1.0 * rand()/RAND_MAX;
      double ext = 1./sqrt(num_target_points * 1.0);
      target_points[j] = Wonton::Point<2>(x, y);
      target_extent[j] = Wonton::Point<2>(ext, ext);
    }

    Swarm<2> source_swarm(source_points);
    Swarm<2> target_swarm(target_points);

    Portage::SearchPointsBins<2, Swarm<2>, Swarm<2>>
      cellsearch(source_swarm, target_swarm, source_extent, target_extent, Gather, 2);

    Portage::SearchSimplePoints<2, Swarm<2>, Swarm<2>>
      simplesearch(source_swarm, target_swarm, source_extent, target_extent, Gather);

    Wonton::vector<std::vector<int>> candidates(num_target_points);
    Wonton::transform(target_swarm.begin(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                      target_swarm.end(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                      candidates.begin(), cellsearch);

    for (int tp = 0; tp < num_target_points; tp++) {
      std::vector<int> xnbr = candidates[tp];
      auto cnbr = cellsearch(tp);
      auto snbr = simplesearch(tp);

      ASSERT_EQ(snbr.size(), cnbr.size());
      ASSERT_EQ(snbr.size(), xnbr.size());

      std::sort(cnbr.begin(), cnbr.end());
      std::sort(xnbr.begin(), xnbr.end());
      int const num_snbr = snbr.size();
      for (int j =  0; j < num_snbr; j++) {
        ASSERT_EQ(snbr[j], cnbr[j]);
        ASSERT_EQ(snbr[j], xnbr[j]);
      }
    }

  }

  static void test_random_3d(int num_source_points, int num_target_points) {

    // random point sets and test against SearchSimplePoints
    Wonton::vector<Wonton::Point<3>> source_points(num_source_points);
    Wonton::vector<Wonton::Point<3>> source_extent(num_source_points);
    Wonton::vector<Wonton::Point<3>> target_points(num_target_points);
    Wonton::vector<Wonton::Point<3>> target_extent(num_target_points);

    for (int j = 0; j < num_source_points; ++j) {
      double x = 1.2 * rand()/RAND_MAX-.1;
      double y = 1.2 * rand()/RAND_MAX-.1;
      double z = 1.2 * rand()/RAND_MAX-.1;
      double ext = 1./pow(num_source_points * 1.0, 1. / 3.);
      source_points[j] = Wonton::Point<3>(x, y, z);
      source_extent[j] = Wonton::Point<3>(ext, ext, ext);
    }

    for (int j = 0; j < num_target_points; ++j) {
      double x = 1.0 * rand()/RAND_MAX;
      double y = 1.0 * rand()/RAND_MAX;
      double z = 1.0 * rand()/RAND_MAX;
      double ext = 1./pow(num_target_points * 1.0, 1. / 3.);
      target_points[j] = Wonton::Point<3>(x, y, z);
      target_extent[j] = Wonton::Point<3>(ext, ext, ext);
    }

    Swarm<3> source_swarm(source_points);
    Swarm<3> target_swarm(target_points);

    Portage::SearchPointsBins<3, Swarm<3>, Swarm<3>>
      cellsearch(source_swarm, target_swarm, source_extent, target_extent, Gather, 2);

    Wonton::vector<std::vector<int>> candidates(num_target_points);
    Wonton::transform(target_swarm.begin(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                      target_swarm.end(Wonton::PARTICLE, Wonton::PARALLEL_OWNED),
                      candidates.begin(), cellsearch);

    Portage::SearchSimplePoints<3, Wonton::Swarm<3>, Wonton::Swarm<3>>
      simplesearch(source_swarm, target_swarm, source_extent, target_extent, Portage::Meshfree::Gather);

    for (int tp = 0; tp < num_target_points; tp++) {
      std::vector<int> xnbr = candidates[tp];
      auto cnbr = cellsearch(tp);
      auto snbr = simplesearch(tp);

      ASSERT_EQ(snbr.size(), cnbr.size());
      ASSERT_EQ(snbr.size(), xnbr.size());

      std::sort(cnbr.begin(), cnbr.end());
      std::sort(xnbr.begin(), xnbr.end());
      int const num_snbr = snbr.size();
      for (int j =  0; j < num_snbr; j++) {
        ASSERT_EQ(snbr[j], cnbr[j]);
        ASSERT_EQ(snbr[j], xnbr[j]);
      }
    }
  }
};

TEST_F(SearchBins, regular_2d) { test_regular_2d(); }

TEST_F(SearchBins, regular_3d) { test_regular_3d(); }

TEST_F(SearchBins, random_2d_case1) { test_random_2d(128, 128); }

TEST_F(SearchBins, random_2d_case2) { test_random_2d(256, 128); }

TEST_F(SearchBins, random_2d_case3) { test_random_2d(128, 256); }

TEST_F(SearchBins, random_3d_case1) { test_random_3d(1000, 1000); }

TEST_F(SearchBins, random_3d_case2) { test_random_3d(1000, 729); }

TEST_F(SearchBins, random_3d_case3) { test_random_3d(729, 1000); }
