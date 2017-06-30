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


#include "search_simple_points.h"

#include <memory>
#include <vector>

#include "gtest/gtest.h"

#include "portage/support/Point.h"
#include "portage/swarm/swarm.h"


TEST(search_simple_points, case2d)
{
  // overlay a 3x3 target swarm on a 4x4 source swarm
  // each target point should have four candidate source points

  std::vector<Portage::Point<2>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<2>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<2>>>(srce);
  const double srclen = 1./4.;
  for (int j = 0; j < 4; ++j) {
    for (int i = 0; i < 4; ++i) {
      double x = (i + 0.5) * srclen;
      double y = (j + 0.5) * srclen;
      double ext = srclen / 4.;
      srcpts->push_back(Portage::Point<2>{x, y});
      srcexts->push_back(Portage::Point<2>{ext, ext});
    }
  }
  Portage::Meshfree::Swarm<2> srcswarm(srcpts);

  std::vector<Portage::Point<2>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<2>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<2>>>(tgte);
  const double tgtlen = 1./3.;
  for (int j = 0; j < 3; ++j) {
    for (int i = 0; i < 3; ++i) {
      double x = (i + 0.5) * tgtlen;
      double y = (j + 0.5) * tgtlen;
      double ext = tgtlen / 4.;
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

} // TEST(search_simple_points, case2d)


TEST(search_simple_points, case3d)
{
  // overlay a 2x2x2 target swarm on a 3x3x3 source swarm
  // each target point should have eight candidate source points

  std::vector<Portage::Point<3>> srcp, srce;
  auto srcpts = std::make_shared<std::vector<Portage::Point<3>>>(srcp);
  auto srcexts = std::make_shared<std::vector<Portage::Point<3>>>(srce);
  const double srclen = 1./3.;
  for (int k = 0; k < 3; ++k) {
    for (int j = 0; j < 3; ++j) {
      for (int i = 0; i < 3; ++i) {
        double x = (i + 0.5) * srclen;
        double y = (j + 0.5) * srclen;
        double z = (k + 0.5) * srclen;
        double ext = srclen / 4.;
        srcpts->push_back(Portage::Point<3>{x, y, z});
        srcexts->push_back(Portage::Point<3>{ext, ext, ext});
      }
    }
  }
  Portage::Meshfree::Swarm<3> srcswarm(srcpts);

  std::vector<Portage::Point<3>> tgtp, tgte;
  auto tgtpts = std::make_shared<std::vector<Portage::Point<3>>>(tgtp);
  auto tgtexts = std::make_shared<std::vector<Portage::Point<3>>>(tgte);
  const double tgtlen = 1./2.;
  for (int k = 0; k < 2; ++k) {
    for (int j = 0; j < 2; ++j) {
      for (int i = 0; i < 2; ++i) {
        double x = (i + 0.5) * tgtlen;
        double y = (j + 0.5) * tgtlen;
        double z = (k + 0.5) * tgtlen;
        double ext = tgtlen / 4.;
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

} // TEST(search_simple_points, case3d)


