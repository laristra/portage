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
