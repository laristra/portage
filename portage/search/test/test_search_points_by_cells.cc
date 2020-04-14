/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <memory>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "gtest/gtest.h"

#include "portage/support/portage.h"
#include "portage/swarm/swarm.h"
#include "portage/search/search_simple_points.h"
#include "portage/search/search_points_by_cells.h"

#include "wonton/support/Point.h"

TEST(search_by_cells, scatter_2d) {

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

  Portage::SearchPointsByCells<2, Swarm<2>, Swarm<2>>
    search(source_swarm, target_swarm,
           source_extent, target_extent);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int tp = i + j * 3;
      auto candidates = search(tp);

      // there should be four candidate source points, in a square
      ASSERT_EQ(unsigned(4), candidates.size());
      // compute spbase = index of lower left source point
      const int spbase = i + j * 4;
      ASSERT_EQ(spbase,     candidates[0]);
      ASSERT_EQ(spbase + 1, candidates[1]);
      ASSERT_EQ(spbase + 4, candidates[2]);
      ASSERT_EQ(spbase + 5, candidates[3]);
    }
  }

} // TEST(search_by_cells, scatter_2d)


TEST(search_by_cells, gather_2d) {

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

  Portage::SearchPointsByCells<2, Swarm<2>, Swarm<2>>
    search(source_swarm, target_swarm,
           source_extent, target_extent, Portage::Meshfree::Gather);

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

} // TEST(search_by_cells, gather_2d)


TEST(search_by_cells, scatter_3d) {

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

  Portage::SearchPointsByCells<3, Swarm<3>, Swarm<3>>
    search(source_swarm, target_swarm,
           source_extent, target_extent);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        int const target_point_index = i + j * 2 + k * 4;
        auto candidates = search(target_point_index);
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
} // TEST(search_by_cells, scatter_3d)


TEST(search_by_cells, gather_3d) {

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

  Portage::SearchPointsByCells<3, Swarm<3>, Swarm<3>>
    search(source_swarm, target_swarm,
           source_extent, target_extent, Portage::Meshfree::Gather);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        int const target_point_index = i + j * 2 + k * 4;
        auto candidates = search(target_point_index);
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
} // TEST(search_by_cells, gather_3d)


void test_scatter_2d_random(int nsrc, int ntgt) {

  using Portage::Meshfree::Swarm;

  // random point sets and test against SearchSimplePoints
  Portage::vector<Wonton::Point<2>> source_points(nsrc);
  Portage::vector<Wonton::Point<2>> source_extent(nsrc);
  Portage::vector<Wonton::Point<2>> target_points(ntgt);
  Portage::vector<Wonton::Point<2>> target_extent(ntgt);

  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2 * rand()/RAND_MAX-.1;
    double y = 1.2 * rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    source_points[j] = Wonton::Point<2>(x, y);
    source_extent[j] = Wonton::Point<2>(ext, ext);
  }


  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0 * rand()/RAND_MAX;
    double y = 1.0 * rand()/RAND_MAX;
    double ext = 1./sqrt(ntgt*1.0);
    target_points[j] = Wonton::Point<2>(x, y);
    target_extent[j] = Wonton::Point<2>(ext, ext);
  }

  Swarm<2> source_swarm(source_points);
  Swarm<2> target_swarm(target_points);

  Portage::SearchPointsByCells<2, Swarm<2>, Swarm<2>>
    cellsearch(source_swarm, target_swarm, source_extent, target_extent);

  Portage::SearchSimplePoints<2, Swarm<2>, Swarm<2>>
    simplesearch(source_swarm, target_swarm, source_extent, target_extent);

  Portage::vector<std::vector<int>> candidates(ntgt);
  Portage::transform(target_swarm.begin(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     target_swarm.end(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     candidates.begin(), cellsearch);

  // use direct application
  for (int tp = 0; tp < ntgt; tp++) {
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

} // test_scatter_2d_random

TEST(search_by_cells, scatter_2d_random_case1)
{
  test_scatter_2d_random(128, 128);
}

TEST(search_by_cells, scatter_2d_random_case2)
{
  test_scatter_2d_random(256, 128);
}

TEST(search_by_cells, scatter_2d_random_case3)
{
  test_scatter_2d_random(128, 256);
}


void test_scatter_3d_random(int nsrc, int ntgt, bool check=true) {

  using Portage::Meshfree::Swarm;

  // random point sets and test against SearchSimplePoints
  Portage::vector<Wonton::Point<3>> source_points(nsrc);
  Portage::vector<Wonton::Point<3>> source_extent(nsrc);
  Portage::vector<Wonton::Point<3>> target_points(ntgt);
  Portage::vector<Wonton::Point<3>> target_extent(ntgt);

  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2 * rand()/RAND_MAX-.1;
    double y = 1.2 * rand()/RAND_MAX-.1;
    double z = 1.2 * rand()/RAND_MAX-.1;
    double ext = 1./pow(nsrc*1.0,1./3.);
    source_points[j] = Wonton::Point<3>(x, y, z);
    source_extent[j] = Wonton::Point<3>(ext, ext, ext);
  }

  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0 * rand()/RAND_MAX;
    double y = 1.0 * rand()/RAND_MAX;
    double z = 1.0 * rand()/RAND_MAX;
    double ext = 1./pow(ntgt*1.0,1./3.);
    target_points[j] = Wonton::Point<3>(x, y, z);
    target_extent[j] = Wonton::Point<3>(ext, ext, ext);
  }

  Swarm<3> source_swarm(source_points);
  Swarm<3> target_swarm(target_points);

  Portage::SearchPointsByCells<3, Swarm<3>, Swarm<3>>
    cellsearch(source_swarm, target_swarm, source_extent, target_extent);

  Portage::vector<std::vector<int>> candidates(ntgt);
  Portage::transform(target_swarm.begin(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     target_swarm.end(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     candidates.begin(), cellsearch);

  if (check) {
    Portage::SearchSimplePoints<3, Swarm<3>, Swarm<3>>
      simplesearch(source_swarm, target_swarm, source_extent, target_extent);

    for (int tp = 0; tp < ntgt; tp++) {
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
} // test_scatter_3d_random

TEST(search_by_cells, scatter_3d_random_case1)
{
  test_scatter_3d_random(1000, 1000);
}

TEST(search_by_cells, scatter_3d_random_case2)
{
  test_scatter_3d_random(1000, 729);
}

TEST(search_by_cells, scatter_3d_random_case3)
{
  test_scatter_3d_random(729, 1000);
}

TEST(search_by_cells, scatter_3d_random_case3_nocheck)
{
  test_scatter_3d_random(2744, 2744, false);
}


void test_gather_2d_random(const int nsrc, const int ntgt) {
  // random point sets and test against SearchSimplePoints

  using Portage::Meshfree::Swarm;
  using Portage::Meshfree::Gather;

  // random point sets and test against SearchSimplePoints
  Portage::vector<Wonton::Point<2>> source_points(nsrc);
  Portage::vector<Wonton::Point<2>> source_extent(nsrc);
  Portage::vector<Wonton::Point<2>> target_points(ntgt);
  Portage::vector<Wonton::Point<2>> target_extent(ntgt);

  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2 * rand()/RAND_MAX-.1;
    double y = 1.2 * rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    source_points[j] = Wonton::Point<2>(x, y);
    source_extent[j] = Wonton::Point<2>(ext, ext);
  }

  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0 * rand()/RAND_MAX;
    double y = 1.0 * rand()/RAND_MAX;
    double ext = 1./sqrt(ntgt*1.0);
    target_points[j] = Wonton::Point<2>(x, y);
    target_extent[j] = Wonton::Point<2>(ext, ext);
  }

  Swarm<2> source_swarm(source_points);
  Swarm<2> target_swarm(target_points);

  Portage::SearchPointsByCells<2, Swarm<2>, Swarm<2>>
    cellsearch(source_swarm, target_swarm, source_extent, target_extent, Gather);

  Portage::SearchSimplePoints<2, Swarm<2>, Swarm<2>>
    simplesearch(source_swarm, target_swarm, source_extent, target_extent, Gather);

  Portage::vector<std::vector<int>> candidates(ntgt);
  Portage::transform(target_swarm.begin(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     target_swarm.end(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     candidates.begin(), cellsearch);

  for (int tp = 0; tp < ntgt; tp++) {
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

} // test_gather_2d_random

TEST(search_by_cells, gather_2d_random_case1)
{
  test_gather_2d_random(128, 128);
}

TEST(search_by_cells, gather_2d_random_case2)
{
  test_gather_2d_random(256, 128);
}

TEST(search_by_cells, gather_2d_random_case3)
{
  test_gather_2d_random(128, 256);
}


void test_gather_3d_random(const int nsrc, const int ntgt, bool check = true) {

  // random point sets and test against SearchSimplePoints
  using Portage::Meshfree::Swarm;
  using Portage::Meshfree::Gather;

  // random point sets and test against SearchSimplePoints
  Portage::vector<Wonton::Point<3>> source_points(nsrc);
  Portage::vector<Wonton::Point<3>> source_extent(nsrc);
  Portage::vector<Wonton::Point<3>> target_points(ntgt);
  Portage::vector<Wonton::Point<3>> target_extent(ntgt);

  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2 * rand()/RAND_MAX-.1;
    double y = 1.2 * rand()/RAND_MAX-.1;
    double z = 1.2 * rand()/RAND_MAX-.1;
    double ext = 1./pow(nsrc*1.0,1./3.);
    source_points[j] = Wonton::Point<3>(x, y, z);
    source_extent[j] = Wonton::Point<3>(ext, ext, ext);
  }

  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0 * rand()/RAND_MAX;
    double y = 1.0 * rand()/RAND_MAX;
    double z = 1.0 * rand()/RAND_MAX;
    double ext = 1./pow(ntgt*1.0,1./3.);
    target_points[j] = Wonton::Point<3>(x, y, z);
    target_extent[j] = Wonton::Point<3>(ext, ext, ext);
  }

  Swarm<3> source_swarm(source_points);
  Swarm<3> target_swarm(target_points);

  Portage::SearchPointsByCells<3, Swarm<3>, Swarm<3>>
    cellsearch(source_swarm, target_swarm, source_extent, target_extent, Gather);

  Portage::vector<std::vector<int>> candidates(ntgt);
  Portage::transform(target_swarm.begin(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     target_swarm.end(Portage::Entity_kind::PARTICLE, Portage::Entity_type::PARALLEL_OWNED),
                     candidates.begin(), cellsearch);

  if (check) {
    Portage::SearchSimplePoints<
      3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
      simplesearch(source_swarm, target_swarm, source_extent, target_extent, Portage::Meshfree::Gather);

    for (int tp = 0; tp < ntgt; tp++) {
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
} // test_gather_3d_random

TEST(search_by_cells, gather_3d_random_case1)
{
  test_gather_3d_random(1000, 1000);
}

TEST(search_by_cells, gather_3d_random_case2)
{
  test_gather_3d_random(1000, 729);
}

TEST(search_by_cells, gather_3d_random_case3)
{
  test_gather_3d_random(729, 1000);
}

TEST(search_by_cells, gather_3d_random_case1_nocheck)
{
  test_gather_3d_random(2744, 2744, false);
}


TEST(search_by_cells, scatter_2d_random_disjoint) {

  using Portage::Meshfree::Swarm;

  int const nsrc = 256;
  int const ntgt = 128;
  // random point sets and test against SearchSimplePoints
  Portage::vector<Wonton::Point<2>> source_points(nsrc);
  Portage::vector<Wonton::Point<2>> source_extent(nsrc);
  Portage::vector<Wonton::Point<2>> target_points(ntgt);
  Portage::vector<Wonton::Point<2>> target_extent(ntgt);

  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2 * rand()/RAND_MAX-.1;
    double y = 1.2 * rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    source_points[j] = Wonton::Point<2>(x, y);
    source_extent[j] = Wonton::Point<2>(ext, ext);
  }

  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0 * rand()/RAND_MAX + 10.;
    double y = 1.0 * rand()/RAND_MAX + 10.;
    double ext = 1./sqrt(ntgt*1.0);
    target_points[j] = Wonton::Point<2>(x, y);
    target_extent[j] = Wonton::Point<2>(ext, ext);
  }

  Swarm<2> source_swarm(source_points);
  Swarm<2> target_swarm(target_points);

  Portage::SearchPointsByCells<2, Swarm<2>, Swarm<2>>
    cellsearch(source_swarm, target_swarm, source_extent, target_extent);

  Portage::SearchSimplePoints<2, Swarm<2>, Swarm<2>>
    simplesearch(source_swarm, target_swarm, source_extent, target_extent);

  for (int tp = 0; tp < ntgt; tp++) {
    auto cnbr = cellsearch(tp);
    auto snbr = simplesearch(tp);
    ASSERT_EQ(snbr.size(), unsigned(0));
    ASSERT_EQ(unsigned(0), cnbr.size());
  }

} // TEST(search_by_cells, scatter_2d_random_disjoint)


TEST(search_by_cells, scatter_2d_random_edge) {

  using Portage::Meshfree::Swarm;

  int const nsrc = 256;
  int const ntgt = 128;

  // random point sets and test against SearchSimplePoints
  Portage::vector<Wonton::Point<2>> source_points(nsrc);
  Portage::vector<Wonton::Point<2>> source_extent(nsrc);
  Portage::vector<Wonton::Point<2>> target_points(ntgt);
  Portage::vector<Wonton::Point<2>> target_extent(ntgt);

  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2 * rand()/RAND_MAX-.1;
    double y = 1.2 * rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    source_points[j] = Wonton::Point<2>(x, y);
    source_extent[j] = Wonton::Point<2>(ext, ext);
  }

  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0 * rand()/RAND_MAX + .5;
    double y = 1.0 * rand()/RAND_MAX + .5;
    double ext = 1./sqrt(ntgt*1.0);
    target_points[j] = Wonton::Point<2>(x, y);
    target_extent[j] = Wonton::Point<2>(ext, ext);
  }

  Swarm<2> source_swarm(source_points);
  Swarm<2> target_swarm(target_points);

  Portage::SearchPointsByCells<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  cellsearch(source_swarm, target_swarm, source_extent, target_extent);

  Portage::SearchSimplePoints<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  simplesearch(source_swarm, target_swarm, source_extent, target_extent);

  for (int tp = 0; tp < ntgt; tp++) {
    auto cnbr = cellsearch(tp);
    auto snbr = simplesearch(tp);

    ASSERT_EQ(snbr.size(), cnbr.size());

    std::sort(cnbr.begin(), cnbr.end());
    int const num_snbr = snbr.size();
    for (int j =  0; j < num_snbr; j++) {
      ASSERT_EQ(snbr[j], cnbr[j]);
    }
  }

} // TEST(search_by_cells, scatter_2d_random_edge)


