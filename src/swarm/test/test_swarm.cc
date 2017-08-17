/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
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
  std::vector<Portage::Point<3>> points(10);

  srand(time(NULL));
  for (int i = 0; i < 10; i++)
    points[i] = Portage::Point<3>(((double)rand()/RAND_MAX),
                                  ((double)rand()/RAND_MAX),
                                  ((double)rand()/RAND_MAX));

  auto p_ptr = std::make_shared<std::vector<Portage::Point<3>>>(points);
                               
  Portage::Meshfree::Swarm<3> swarm(p_ptr);

  // Did we get back 10 points?
  ASSERT_EQ(10, swarm.num_owned_particles());

  // Check the bounding box and center of the hexahedral "cell"
  // surrounding each point

  for (int i = 0; i < 10; i++) {
    // Are the point coordinates correct?
    auto pt = swarm.get_particle_coordinates(i);
    for (size_t j=0; j<3; j++) ASSERT_EQ(pt[j], points[i][j]);
  }

}  // TEST


