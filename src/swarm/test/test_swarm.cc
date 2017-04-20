/*
 * test_swarm.cc
 *
 *  Created on: Apr 19, 2017
 *      Author: gad
 */

#include "gtest/gtest.h"
#include "swarm.h"

TEST(Swarm, placeholder) {
  std::vector<Portage::Point<3>> points(100), extents(100);
  auto p_ptr=std::make_shared<std::vector<Portage::Point<3>>>(points);
  auto e_ptr=std::make_shared<std::vector<Portage::Point<3>>>(extents);
  Portage::Meshfree::Swarm<3> swarm(p_ptr, e_ptr);
}


