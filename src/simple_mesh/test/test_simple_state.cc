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
