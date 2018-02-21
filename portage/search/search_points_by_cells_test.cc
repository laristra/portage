/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "search_points_by_cells.h"

#include <memory>
#include <vector>
#include <cstdlib>
#include <algorithm>

#include "gtest/gtest.h"

#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "portage/search/search_simple_points.h"


TEST(search_by_cells, scatter_2d)
{
  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points

  std::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<2>>>(srce);
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      double x = (i + 0.5);
      double y = (j + 0.5);
      double ext = 0.375;
      srcpts->push_back(Portage::Point<2>{x, y});
      srcexts->push_back(Portage::Point<2>{ext, ext});
    }
  }
  Portage::Meshfree::Swarm<2> srcswarm(srcpts);

  std::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<2>>>(tgte);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      double x = (i+1.);
      double y = (j+1.);
      double ext = 0.375;
      tgtpts->push_back(Portage::Point<2>{x, y});
      tgtexts->push_back(Portage::Point<2>{ext, ext});
    }
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
    2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
    search(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int tp = i + j * 3;
      std::vector<unsigned int> candidates;
      candidates = search(tp);

      // there should be four candidate source points, in a square
      ASSERT_EQ(4, candidates.size());
      // compute spbase = index of lower left source point
      const int spbase = i + j * 4;
      ASSERT_EQ(spbase,     candidates[0]);
      ASSERT_EQ(spbase + 1, candidates[1]);
      ASSERT_EQ(spbase + 4, candidates[2]);
      ASSERT_EQ(spbase + 5, candidates[3]);
    }
  }

} // TEST(search_by_cells, scatter_2d)


TEST(search_by_cells, gather_2d)
{
  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points

  std::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<2>>>(srce);
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      double x = (i + 0.5);
      double y = (j + 0.5);
      double ext = 0.375;
      srcpts->push_back(Portage::Point<2>{x, y});
      srcexts->push_back(Portage::Point<2>{ext, ext});
    }
  }
  Portage::Meshfree::Swarm<2> srcswarm(srcpts);

  std::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<2>>>(tgte);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      double x = (i+1.);
      double y = (j+1.);
      double ext = 0.375;
      tgtpts->push_back(Portage::Point<2>{x, y});
      tgtexts->push_back(Portage::Point<2>{ext, ext});
    }
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
    2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
    search(srcswarm, tgtswarm, srcexts, tgtexts, Portage::Meshfree::Gather);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int tp = i + j * 3;
      std::vector<unsigned int> candidates;
      candidates = search(tp);

      // there should be four candidate source points, in a square
      ASSERT_EQ(4, candidates.size());
      // compute spbase = index of lower left source point
      const int spbase = i + j * 4;
      ASSERT_EQ(spbase,     candidates[0]);
      ASSERT_EQ(spbase + 1, candidates[1]);
      ASSERT_EQ(spbase + 4, candidates[2]);
      ASSERT_EQ(spbase + 5, candidates[3]);
    }
  }

} // TEST(search_by_cells, gather_2d)


TEST(search_by_cells, scatter_3d)
{
  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
  // each target point should have eight candidate source points

  std::vector<Portage::Point<3>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<3>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<3>>>(srce);
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double x = (i + 0.5);
        double y = (j + 0.5);
        double z = (k + 0.5);
        double ext = .375;
        srcpts->push_back(Portage::Point<3>{x, y, z});
        srcexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> srcswarm(srcpts);

  std::vector<Portage::Point<3>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<3>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<3>>>(tgte);
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        double x = (i + 1.0);
        double y = (j + 1.0);
        double z = (k + 1.0);
        double ext = .375;
        tgtpts->push_back(Portage::Point<3>{x, y, z});
        tgtexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
    3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
    search(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        const int tp = i + j * 2 + k * 4;
        std::vector<unsigned int> candidates;
        candidates = search(tp);

        // there should be eight candidate source points, in a cube
        ASSERT_EQ(8, candidates.size());
        // compute spbase = index of lower left source point
        const int spbase = i + j * 3 + k * 9;
        ASSERT_EQ(spbase,      candidates[0]);
        ASSERT_EQ(spbase + 1,  candidates[1]);
        ASSERT_EQ(spbase + 3,  candidates[2]);
        ASSERT_EQ(spbase + 4,  candidates[3]);
        ASSERT_EQ(spbase + 9,  candidates[4]);
        ASSERT_EQ(spbase + 10, candidates[5]);
        ASSERT_EQ(spbase + 12, candidates[6]);
        ASSERT_EQ(spbase + 13, candidates[7]);
      }
    }
  }

} // TEST(search_by_cells, scatter_3d)


TEST(search_by_cells, gather_3d)
{
  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
  // each target point should have eight candidate source points

  std::vector<Portage::Point<3>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<3>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<3>>>(srce);
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double x = (i + 0.5);
        double y = (j + 0.5);
        double z = (k + 0.5);
        double ext = .375;
        srcpts->push_back(Portage::Point<3>{x, y, z});
        srcexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> srcswarm(srcpts);

  std::vector<Portage::Point<3>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<3>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<3>>>(tgte);
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        double x = (i + 1.0);
        double y = (j + 1.0);
        double z = (k + 1.0);
        double ext = .375;
        tgtpts->push_back(Portage::Point<3>{x, y, z});
        tgtexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
    3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
    search(srcswarm, tgtswarm, srcexts, tgtexts, Portage::Meshfree::Gather);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        const int tp = i + j * 2 + k * 4;
        std::vector<unsigned int> candidates;
        candidates = search(tp);

        // there should be eight candidate source points, in a cube
        ASSERT_EQ(8, candidates.size());
        // compute spbase = index of lower left source point
        const int spbase = i + j * 3 + k * 9;
        ASSERT_EQ(spbase,      candidates[0]);
        ASSERT_EQ(spbase + 1,  candidates[1]);
        ASSERT_EQ(spbase + 3,  candidates[2]);
        ASSERT_EQ(spbase + 4,  candidates[3]);
        ASSERT_EQ(spbase + 9,  candidates[4]);
        ASSERT_EQ(spbase + 10, candidates[5]);
        ASSERT_EQ(spbase + 12, candidates[6]);
        ASSERT_EQ(spbase + 13, candidates[7]);
      }
    }
  }

} // TEST(search_by_cells, gather_3d)


void test_scatter_2d_random(const size_t nsrc, const size_t ntgt)
{ // random point sets and test against SearchSimplePoints

  std::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<2>>>(srce);
  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2*rand()/RAND_MAX-.1;
    double y = 1.2*rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    srcpts->push_back(Portage::Point<2>{x, y});
    srcexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> srcswarm(srcpts);

  std::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<2>>>(tgte);
  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0*rand()/RAND_MAX;
    double y = 1.0*rand()/RAND_MAX;
    double ext = 1./sqrt(ntgt*1.0);
    tgtpts->push_back(Portage::Point<2>{x, y});
    tgtexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  cellsearch(srcswarm, tgtswarm, srcexts, tgtexts);

  Portage::SearchSimplePoints<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  simplesearch(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int tp = 0; tp < ntgt; tp++) {
    std::vector<unsigned int> cnbr, snbr;
    cnbr = cellsearch(tp);
    snbr = simplesearch(tp);

    ASSERT_EQ(snbr.size(), cnbr.size());

    std::sort(cnbr.begin(), cnbr.end());
    for (int j=0; j<snbr.size(); j++) {
      ASSERT_EQ(snbr[j], cnbr[j]);
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


void test_scatter_3d_random(const size_t nsrc, const size_t ntgt)
{ // random point sets and test against SearchSimplePoints

  std::vector<Portage::Point<3>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<3>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<3>>>(srce);
  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2*rand()/RAND_MAX-.1;
    double y = 1.2*rand()/RAND_MAX-.1;
    double z = 1.2*rand()/RAND_MAX-.1;
    double ext = 4./pow(nsrc*1.0,1./3.);
    srcpts->push_back(Portage::Point<3>{x, y, z});
    srcexts->push_back(Portage::Point<3>{ext, ext, ext});
  }
  Portage::Meshfree::Swarm<3> srcswarm(srcpts);

  std::vector<Portage::Point<3>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<3>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<3>>>(tgte);
  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0*rand()/RAND_MAX;
    double y = 1.0*rand()/RAND_MAX;
    double z = 1.0*rand()/RAND_MAX;
    double ext = 4./pow(ntgt*1.0,1./3.);
    tgtpts->push_back(Portage::Point<3>{x, y, z});
    tgtexts->push_back(Portage::Point<3>{ext, ext, ext});
  }
  Portage::Meshfree::Swarm<3> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
  3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
  cellsearch(srcswarm, tgtswarm, srcexts, tgtexts);

  Portage::SearchSimplePoints<
  3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
  simplesearch(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int tp = 0; tp < ntgt; tp++) {
    std::vector<unsigned int> cnbr, snbr;
    cnbr = cellsearch(tp);
    snbr = simplesearch(tp);

    ASSERT_EQ(snbr.size(), cnbr.size());

    std::sort(cnbr.begin(), cnbr.end());
    for (int j=0; j<snbr.size(); j++) {
      ASSERT_EQ(snbr[j], cnbr[j]);
    }
  }

} // test_scatter_3d_random

TEST(search_by_cells, scatter_3d_random_case1)
{
  test_scatter_3d_random(500, 500);
}

TEST(search_by_cells, scatter_3d_random_case2)
{
  test_scatter_3d_random(1000, 500);
}

TEST(search_by_cells, scatter_3d_random_case3)
{
  test_scatter_3d_random(500, 1000);
}


void test_gather_2d_random(const size_t nsrc, const size_t ntgt)
{ // random point sets and test against SearchSimplePoints

  std::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<2>>>(srce);
  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2*rand()/RAND_MAX-.1;
    double y = 1.2*rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    srcpts->push_back(Portage::Point<2>{x, y});
    srcexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> srcswarm(srcpts);

  std::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<2>>>(tgte);
  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0*rand()/RAND_MAX;
    double y = 1.0*rand()/RAND_MAX;
    double ext = 1./sqrt(ntgt*1.0);
    tgtpts->push_back(Portage::Point<2>{x, y});
    tgtexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  cellsearch(srcswarm, tgtswarm, srcexts, tgtexts, Portage::Meshfree::Gather);

  Portage::SearchSimplePoints<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  simplesearch(srcswarm, tgtswarm, srcexts, tgtexts, Portage::Meshfree::Gather);

  for (int tp = 0; tp < ntgt; tp++) {
    std::vector<unsigned int> cnbr, snbr;
    cnbr = cellsearch(tp);
    snbr = simplesearch(tp);

    ASSERT_EQ(snbr.size(), cnbr.size());

    std::sort(cnbr.begin(), cnbr.end());
    for (int j=0; j<snbr.size(); j++) {
      ASSERT_EQ(snbr[j], cnbr[j]);
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


void test_gather_3d_random(const size_t nsrc, const size_t ntgt)
{ // random point sets and test against SearchSimplePoints

  std::vector<Portage::Point<3>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<3>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<3>>>(srce);
  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2*rand()/RAND_MAX-.1;
    double y = 1.2*rand()/RAND_MAX-.1;
    double z = 1.2*rand()/RAND_MAX-.1;
    double ext = 4./pow(nsrc*1.0,1./3.);
    srcpts->push_back(Portage::Point<3>{x, y, z});
    srcexts->push_back(Portage::Point<3>{ext, ext, ext});
  }
  Portage::Meshfree::Swarm<3> srcswarm(srcpts);

  std::vector<Portage::Point<3>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<3>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<3>>>(tgte);
  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0*rand()/RAND_MAX;
    double y = 1.0*rand()/RAND_MAX;
    double z = 1.0*rand()/RAND_MAX;
    double ext = 4./pow(ntgt*1.0,1./3.);
    tgtpts->push_back(Portage::Point<3>{x, y, z});
    tgtexts->push_back(Portage::Point<3>{ext, ext, ext});
  }
  Portage::Meshfree::Swarm<3> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
  3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
  cellsearch(srcswarm, tgtswarm, srcexts, tgtexts, Portage::Meshfree::Gather);

  Portage::SearchSimplePoints<
  3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
  simplesearch(srcswarm, tgtswarm, srcexts, tgtexts, Portage::Meshfree::Gather);

  for (int tp = 0; tp < ntgt; tp++) {
    std::vector<unsigned int> cnbr, snbr;
    cnbr = cellsearch(tp);
    snbr = simplesearch(tp);

    ASSERT_EQ(snbr.size(), cnbr.size());

    std::sort(cnbr.begin(), cnbr.end());
    for (int j=0; j<snbr.size(); j++) {
      ASSERT_EQ(snbr[j], cnbr[j]);
    }
  }

} // test_gather_3d_random

TEST(search_by_cells, gather_3d_random_case1)
{
  test_gather_3d_random(500, 500);
}

TEST(search_by_cells, gather_3d_random_case2)
{
  test_gather_3d_random(1000, 500);
}

TEST(search_by_cells, gather_3d_random_case3)
{
  test_gather_3d_random(500, 1000);
}


TEST(search_by_cells, scatter_2d_random_disjoint)
{ // random point sets and test against SearchSimplePoints

  std::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<2>>>(srce);
  const size_t nsrc = 256;
  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2*rand()/RAND_MAX-.1;
    double y = 1.2*rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    srcpts->push_back(Portage::Point<2>{x, y});
    srcexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> srcswarm(srcpts);

  std::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<2>>>(tgte);
  const size_t ntgt = 128;
  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0*rand()/RAND_MAX + 10.;
    double y = 1.0*rand()/RAND_MAX + 10.;
    double ext = 1./sqrt(ntgt*1.0);
    tgtpts->push_back(Portage::Point<2>{x, y});
    tgtexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  cellsearch(srcswarm, tgtswarm, srcexts, tgtexts);

  Portage::SearchSimplePoints<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  simplesearch(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int tp = 0; tp < ntgt; tp++) {
    std::vector<unsigned int> cnbr, snbr;
    cnbr = cellsearch(tp);
    snbr = simplesearch(tp);

    ASSERT_EQ(snbr.size(), 0);
    ASSERT_EQ(0, cnbr.size());
  }

} // TEST(search_by_cells, scatter_2d_random_disjoint)


TEST(search_by_cells, scatter_2d_random_edge)
{ // random point sets and test against SearchSimplePoints

  std::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<2>>>(srce);
  const size_t nsrc = 256;
  for (int j = 0; j < nsrc; ++j) {
    double x = 1.2*rand()/RAND_MAX-.1;
    double y = 1.2*rand()/RAND_MAX-.1;
    double ext = 1./sqrt(nsrc*1.0);
    srcpts->push_back(Portage::Point<2>{x, y});
    srcexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> srcswarm(srcpts);

  std::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<2>>>(tgte);
  const size_t ntgt = 128;
  for (int j = 0; j < ntgt; ++j) {
    double x = 1.0*rand()/RAND_MAX + .5;
    double y = 1.0*rand()/RAND_MAX + .5;
    double ext = 1./sqrt(ntgt*1.0);
    tgtpts->push_back(Portage::Point<2>{x, y});
    tgtexts->push_back(Portage::Point<2>{ext, ext});
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchPointsByCells<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  cellsearch(srcswarm, tgtswarm, srcexts, tgtexts);

  Portage::SearchSimplePoints<
  2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
  simplesearch(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int tp = 0; tp < ntgt; tp++) {
    std::vector<unsigned int> cnbr, snbr;
    cnbr = cellsearch(tp);
    snbr = simplesearch(tp);

    ASSERT_EQ(snbr.size(), cnbr.size());

    std::sort(cnbr.begin(), cnbr.end());
    for (int j=0; j<snbr.size(); j++) {
      ASSERT_EQ(snbr[j], cnbr[j]);
    }
  }

} // TEST(search_by_cells, scatter_2d_random_edge)


