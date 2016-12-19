/*---------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include <memory>
#include <vector>
#include <string>

#include "portage/wrappers/state/simple_state/simple_state_wrapper.h"

#include "gtest/gtest.h"

#include "portage/support/portage.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/simple_mesh/simple_state.h"
#include "portage/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"


TEST(Simple_State_Wrapper, WrapperTest) {
  Portage::Simple_Mesh mymesh(0.0, 0.0, 0.0,
                              1.0, 1.0, 1.0,
                              2, 2, 2);
  Portage::Simple_State mystate(std::make_shared<Portage::Simple_Mesh>(mymesh));

  // Wrappers
  Portage::Simple_Mesh_Wrapper mymesh_wrapper(mymesh);
  Portage::Simple_State_Wrapper mystate_wrapper(mystate);

  auto numcells = mymesh_wrapper.num_owned_cells();
  auto numnodes = mymesh_wrapper.num_owned_nodes();

  std::vector<double> cellvar(numcells);
  for (int i(0); i < numcells; ++i)
    cellvar[i] = i;
  mystate.add("cellvar1", Portage::Entity_kind::CELL, &(cellvar[0]));

  double* celldata;
  mystate_wrapper.get_data(Portage::Entity_kind::CELL, "cellvar1",
                           &celldata);
  for (int i(0); i < numcells; ++i) {
    ASSERT_EQ(cellvar[i], celldata[i]);
  }


  std::vector<double> nodevar(numnodes);
  for (int i(0); i < numnodes; ++i)
    nodevar[i] = i;
  mystate.add("nodevar1", Portage::Entity_kind::NODE, &(nodevar[0]));

  double* nodedata;
  mystate_wrapper.get_data(Portage::Entity_kind::NODE, "nodevar1",
                           &nodedata);
  for (int i(0); i < numnodes; ++i) {
    ASSERT_EQ(nodevar[i], nodedata[i]);
  }

  // Check get_entity
  std::vector<std::string> names = {"cellvar1", "nodevar1"};
  std::vector<Portage::Entity_kind>
      expected_kinds = {Portage::Entity_kind::CELL,
                        Portage::Entity_kind::NODE};
  Portage::Entity_kind kind;
  for (int i(0); i < names.size(); ++i) {
    kind = mystate_wrapper.get_entity(names[i]);
    ASSERT_EQ(expected_kinds[i], kind);
  }
}
