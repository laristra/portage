/*----------------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include <memory>
#include <iterator>

#include "portage/simple_mesh/simple_state.h"
#include "portage/simple_mesh/simple_mesh.h"

#include "portage/support/portage.h"
#include "portage/support/Point.h"

#include "gtest/gtest.h"


TEST(Simple_State, MultiCell) {
  Portage::Simple_Mesh mymesh(0.0, 0.0, 0.0,
                              1.0, 1.0, 1.0,
                              2, 2, 2);
  Portage::Simple_State mystate(std::make_shared<Portage::Simple_Mesh>(mymesh));

  auto numcells = mymesh.num_entities(Portage::Entity_kind::CELL,
                                      Portage::Entity_type::PARALLEL_OWNED);
  auto numnodes = mymesh.num_entities(Portage::Entity_kind::NODE,
                                      Portage::Entity_type::PARALLEL_OWNED);
  std::cout << "numcells " << numcells << std::endl;
  // Add a cell-centered variable, which should initialize to zeros
  auto& cellvar1 = mystate.add("cellvar1", Portage::Entity_kind::CELL);
  ASSERT_EQ(numcells, cellvar1.size());
  for (auto const cv : cellvar1) {
    ASSERT_EQ(0.0, cv);
    std::cout << " cv " << cv << std::endl;
  }

  // Add a node centered variable with an initialized array
  std::vector<double> nodevar1(numnodes);
  for (int i(0); i < numnodes; ++i)
    nodevar1[i] = i;
  auto retnodevar1 = mystate.add("nodevar1", Portage::Entity_kind::NODE,
                                 &(nodevar1[0]));
  ASSERT_EQ(numnodes, retnodevar1.size());
  for (int i(0); i < numnodes; ++i)
    ASSERT_EQ(nodevar1[i], retnodevar1[i]);

  auto numVars = mystate.numVars();
  ASSERT_EQ(2, numVars);

  // Now try modifying the data and then getting a new copy to make sure
  // it was actually modified.
  std::cout << "setting var" << std::endl;
  std::copy(cellvar1.begin(), cellvar1.end(),
            std::ostream_iterator<double>(std::cout , " "));
  std::cout << std::endl;
  for (int i(0); i < numcells; ++i)
    std::cout << " " << cellvar1[i];
  std::cout << std::endl;

  
  cellvar1[4] = 10.0;
  std::cout << "getting it again" << std::endl;
  auto cellvar2 = mystate.get_data("cellvar1", Portage::Entity_kind::CELL);
  ASSERT_EQ(10.0, cellvar2[4]);
}
