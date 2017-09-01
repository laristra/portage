/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <vector>
#include <memory>
#include <iterator>
#include <stdexcept>

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

  // Add a cell-centered variable, which should initialize to zeros
  auto& cellvar1 = mystate.add("cellvar1", Portage::Entity_kind::CELL);
  ASSERT_EQ(numcells, cellvar1.size());
  for (auto const cv : cellvar1) {
    ASSERT_EQ(0.0, cv);
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
  cellvar1[4] = 10.0;

  auto & cellvar2 = mystate.get("cellvar1", Portage::Entity_kind::CELL);
  ASSERT_EQ(10.0, cellvar2[4]);

  // Make sure this throws an exception
  ASSERT_THROW(mystate.get("missing variable",
                           Portage::Entity_kind::CELL),
               std::runtime_error);

  // Add a variable using just a constant value initializer
  auto& cellvar3 = mystate.add("cellvar3", Portage::Entity_kind::CELL, 7.0);
  for (auto const cv : cellvar3)
    ASSERT_EQ(7.0, cv);
}
