/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "search_simple_points.h"

#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"
#include "portage/accumulate/accumulate.h"


TEST(search_simple_points, scatter_2d)
{
  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points

  Portage::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<Portage::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<Portage::vector<Portage::Point<2>>>(srce);
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

  Portage::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<Portage::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<Portage::vector<Portage::Point<2>>>(tgte);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      double x = (i + 1.0);
      double y = (j + 1.0);
      double ext = 0.375;
      tgtpts->push_back(Portage::Point<2>{x, y});
      tgtexts->push_back(Portage::Point<2>{ext, ext});
    }
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchSimplePoints<
    2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
    search(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int tp = i + j * 3;
      std::vector<unsigned int> candidates = search(tp);

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

} // TEST(search_simple_points, scatter_2d)


TEST(search_simple_points, scatter_3d)
{
  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
  // each target point should have eight candidate source points

  Portage::vector<Portage::Point<3>> srcp, srce;
  auto srcpts = std::make_shared<Portage::vector<Portage::Point<3>>>(srcp);
  auto srcexts = std::make_shared<Portage::vector<Portage::Point<3>>>(srce);
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double x = (i + 0.5);
        double y = (j + 0.5);
        double z = (k + 0.5);
        double ext = 0.375;
        srcpts->push_back(Portage::Point<3>{x, y, z});
        srcexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> srcswarm(srcpts);

  Portage::vector<Portage::Point<3>> tgtp, tgte;
  auto tgtpts = std::make_shared<Portage::vector<Portage::Point<3>>>(tgtp);
  auto tgtexts = std::make_shared<Portage::vector<Portage::Point<3>>>(tgte);
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        double x = (i + 1.0);
        double y = (j + 1.0);
        double z = (k + 1.0);
        double ext = 0.375;
        tgtpts->push_back(Portage::Point<3>{x, y, z});
        tgtexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> tgtswarm(tgtpts);

  Portage::SearchSimplePoints<
    3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
    search(srcswarm, tgtswarm, srcexts, tgtexts);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        const int tp = i + j * 2 + k * 4;
        std::vector<unsigned int> candidates = search(tp);

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

} // TEST(search_simple_points, scatter_3d)


TEST(search_simple_points, gather_2d)
{
  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points

  Portage::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<Portage::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<Portage::vector<Portage::Point<2>>>(srce);
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

  Portage::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<Portage::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<Portage::vector<Portage::Point<2>>>(tgte);
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      double x = (i + 1.0);
      double y = (j + 1.0);
      double ext = 0.375;
      tgtpts->push_back(Portage::Point<2>{x, y});
      tgtexts->push_back(Portage::Point<2>{ext, ext});
    }
  }
  Portage::Meshfree::Swarm<2> tgtswarm(tgtpts);

  Portage::SearchSimplePoints<
    2, Portage::Meshfree::Swarm<2>, Portage::Meshfree::Swarm<2>>
    search(srcswarm, tgtswarm, srcexts, tgtexts,
           Portage::Meshfree::Gather);

  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      const int tp = i + j * 3;
      std::vector<unsigned int> candidates = search(tp);

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

} // TEST(search_simple_points, gather_2d)


TEST(search_simple_points, gather_3d)
{
  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
  // each target point should have eight candidate source points

  Portage::vector<Portage::Point<3>> srcp, srce;
  auto srcpts = std::make_shared<Portage::vector<Portage::Point<3>>>(srcp);
  auto srcexts = std::make_shared<Portage::vector<Portage::Point<3>>>(srce);
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double x = (i + 0.5);
        double y = (j + 0.5);
        double z = (k + 0.5);
        double ext = 0.375;
        srcpts->push_back(Portage::Point<3>{x, y, z});
        srcexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> srcswarm(srcpts);

  Portage::vector<Portage::Point<3>> tgtp, tgte;
  auto tgtpts = std::make_shared<Portage::vector<Portage::Point<3>>>(tgtp);
  auto tgtexts = std::make_shared<Portage::vector<Portage::Point<3>>>(tgte);
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        double x = (i + 1.0);
        double y = (j + 1.0);
        double z = (k + 1.0);
        double ext = 0.375;
        tgtpts->push_back(Portage::Point<3>{x, y, z});
        tgtexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> tgtswarm(tgtpts);

  Portage::SearchSimplePoints<
    3, Portage::Meshfree::Swarm<3>, Portage::Meshfree::Swarm<3>>
    search(srcswarm, tgtswarm, srcexts, tgtexts,
           Portage::Meshfree::Gather);

  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        const int tp = i + j * 2 + k * 4;
        std::vector<unsigned int> candidates = search(tp);

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

} // TEST(search_simple_points, gather_3d)


