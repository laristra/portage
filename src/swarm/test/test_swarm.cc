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

/*
 * test_swarm.cc
 *
 *  Created on: Apr 19, 2017
 *      Author: gad
 *      Updated: rao
 */

#include <cstdlib>
#include <ctime>
#include <algorithm>

#include "portage/swarm/swarm.h"

#include "gtest/gtest.h"

TEST(Swarm, Sanity_Check) {
  std::vector<Portage::Point<3>> points(10), extents(10);

  double h = 0.01;

  srand(time(NULL));
  for (int i = 0; i < 10; i++) {
    points[i] = Portage::Point<3>(((double)rand()/RAND_MAX),
                                  ((double)rand()/RAND_MAX),
                                  ((double)rand()/RAND_MAX));
    extents[i] = Portage::Point<3>(h, h, h);
  }

  auto p_ptr = std::make_shared<std::vector<Portage::Point<3>>>(points);
  auto e_ptr = std::make_shared<std::vector<Portage::Point<3>>>(extents);
                               
  Portage::Meshfree::Swarm<3> swarm(p_ptr, e_ptr);

  // Did we get back 10 points?
  ASSERT_EQ(10, swarm.num_owned_cells());

  // Check the bounding box and center of the hexahedral "cell"
  // surrounding each point

  for (int i = 0; i < 10; i++) {
    Portage::Point<3> exp_bbox[2];
    exp_bbox[0] = Portage::Point<3>(points[i][0]-0.5*h,
                                    points[i][1]-0.5*h,
                                    points[i][2]-0.5*h);
    exp_bbox[1] = Portage::Point<3>(points[i][0]+0.5*h,
                                    points[i][1]+0.5*h,
                                    points[i][2]+0.5*h);

    std::vector<Portage::Point<3>> cpoints;
    swarm.cell_get_coordinates(i, &cpoints);

    // Did we get all corners of the hex
    int np = cpoints.size();
    ASSERT_EQ(8, np);

    // What is the bounding box and center of this hex
    Portage::Point<3> bbox[2], cen;
    bbox[0] = Portage::Point<3>( 1.e99,  1.e99,  1.e99);
    bbox[1] = Portage::Point<3>(-1.e99, -1.e99, -1.e99);
    cen = Portage::Point<3>(0.0, 0.0, 0.0);

    for (int p = 0; p < np; p++) {
      cen += cpoints[p];
      for (int d = 0; d < 3; d++) {
        bbox[0][d] = std::min(cpoints[p][d], bbox[0][d]);
        bbox[1][d] = std::max(cpoints[p][d], bbox[1][d]);
      }
    }
    cen /= np;  // average of all coordinates

    for (int d = 0; d < 3; d++) {
      ASSERT_NEAR(points[i][d], cen[d], 1.0e-8);
      ASSERT_NEAR(h, bbox[1][d]-bbox[0][d], 1.0e-8);
    }  
  }  // for each point used to construct the swarm

}  // TEST


