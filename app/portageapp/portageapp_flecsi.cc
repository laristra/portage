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



#include <vector>
#include <iostream>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <iterator>

#define PORTAGE_SERIAL_ONLY
#include "portage/support/Point.h"
#include "portage/wrappers/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "portage/wrappers/state/flecsi/flecsi_state_wrapper.h"
#include "portage/driver/driver.h"

#include "flecsi/specializations/burton/burton.h"
#include "flecsi/io/io.h"
#include "flecsi/specializations/burton/burton_io_exodus.h"

using mesh_t = flecsi::burton_mesh_t;

using real_t = flecsi::mesh_t::real_t;

/*!
  @file main_flecsi.cc
  @brief A simple application that drives our remap routines on a FleCSI mesh.

  The program is used to showcase our capabilities with various types of remap
  operations on a 2d FleCSI mesh with simple linear.  For the
  cases of remapping linear data with a second-order interpolator, the L2 norm
  output at the end should be identically zero.
 */


void print_usage() {
  std::printf("usage: portage_flecsiapp ncellsx ncellsy [order=1]\n");  // [y];
  //  std::printf("    if 'y' is specified, then dump the meshes to Exodus\n");
}

int main(int argc, char** argv) {
  if (argc <= 2) {
    print_usage();
    return 0;
  }

  // number of CELLS in x and y for input mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  auto order = (argc > 3) ? atoi(argv[3]) : 1;
  // bool dump_output = (argc > 4) ?
  //     (std::string(argv[4]) == "y") ? true : false
  //     : false;

  // dimensions of meshes
  real_t xmin = 0.0, xmax = 1.0;
  real_t ymin = 0.0, ymax = 1.0;

  std::printf("starting portage_flecsiapp...\n");

  // Setup the meshes
  mesh_t inputMesh, targetMesh;
  Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx, ny, inputMesh);
  Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx+1, ny+1, targetMesh);

  // Setup the wrappers
  Portage::Flecsi_Mesh_Wrapper inputMeshWrapper(inputMesh);
  Portage::Flecsi_Mesh_Wrapper targetMeshWrapper(targetMesh);

  // Register the state with FleCSI
  // register_state is a macro from burton.h
  // we keep the inputState accessor for filling of data
  auto inputState = register_state(inputMesh, "celldata", cells, real_t,
                                   flecsi::persistent);
  register_state(targetMesh, "celldata", cells, real_t, flecsi::persistent);

  // Fill with linear function
  for (auto c : inputMesh.cells()) {
    auto cen = c->centroid();
    inputState[c] = cen[0] + cen[1];
  }

  // Build the state wrappers
  Portage::Flecsi_State_Wrapper inputStateWrapper(inputMesh);
  Portage::Flecsi_State_Wrapper targetStateWrapper(targetMesh);

  // Setup the main driver for this mesh type
  if(order==2){

    Portage::Driver<
      Portage::SearchKDTree,
          Portage::IntersectR2D,
      Portage::Interpolate_2ndOrder,
      2,
      Portage::Flecsi_Mesh_Wrapper,
          Portage::Flecsi_State_Wrapper>
      d(inputMeshWrapper, inputStateWrapper, targetMeshWrapper, targetStateWrapper);
  
  // Declare which variables are remapped
  std::vector<std::string> varnames(1, "celldata");
  d.set_remap_var_names(varnames);

  // Do the remap
  d.run(false);
  }

  // Setup the main driver for this mesh type
  if(order==1){

    Portage::Driver<
      Portage::SearchKDTree,
          Portage::IntersectR2D,
      Portage::Interpolate_1stOrder,
      2,
      Portage::Flecsi_Mesh_Wrapper,
          Portage::Flecsi_State_Wrapper>
      d(inputMeshWrapper, inputStateWrapper, targetMeshWrapper, targetStateWrapper);
  
  // Declare which variables are remapped
  std::vector<std::string> varnames(1, "celldata");
  d.set_remap_var_names(varnames);

  // Do the remap
  d.run(false);
  }


  // Get the new data and calculate the error
  // use flecsi intrinsics - access_state is a macro in burton.h
  auto outData = access_state(targetMesh, "celldata", real_t);

  double toterr = 0.0;
  for (auto c : targetMesh.cells()) {
    auto cen = c->centroid();
    double error = cen[0] + cen[1] - outData[c];
    std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c.id(),
                cen[0], cen[1]);
    std::printf("  Value = % 10.6lf  Err = % lf\n",
                outData[c], error);
    toterr += error*error;
  }
  std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));

  /// @TODO Dumping the targetMesh is still broken with FleCSI
  // if (dump_output) {
  //   flecsi::write_mesh("input.exo", inputMesh);
  //   flecsi::write_mesh("output.exo", targetMesh);
  // }

  return 0;
}
