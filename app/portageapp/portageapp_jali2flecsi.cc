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

#include <mpi.h>

#define PORTAGE_SERIAL_ONLY
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "portage/wrappers/state/flecsi/flecsi_state_wrapper.h"
#include "portage/driver/driver.h"

// FleCSI includes
#include "flecsi/specializations/burton/burton.h"

// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

using mesh_t = flecsi::burton_mesh_t;

using real_t = flecsi::mesh_t::real_t;

/*!
  @file main_jali_to_flecsi.cc
  @brief A simple application that drives our remap routines from Jali mesh to
  a FleCSI mesh.

  The program is used to showcase our capabilities with various types of remap
  operations on 2d meshes with simple linear or quadratic data.For the
  cases of remapping linear data with a second-order interpolator, the L2 norm
  output at the end should be identically zero.
 */

void print_usage() {
  std::printf("usage: portage_jali2flecsiapp ncellsx ncellsy order=1\n");
}

int main(int argc, char** argv) {
  if (argc <= 2) {
    print_usage();
    return 0;
  }

  // number of CELLS in x and y for input mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  auto order = argc == 4 ? atoi(argv[3]) : 1;
  // dimensions of meshes
  real_t xmin = 0.0, xmax = 1.0;
  real_t ymin = 0.0, ymax = 1.0;

  // Jali needs a COMM
  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  if (numpe > 1) {
    std::printf("error - only 1 mpi rank is allowed\n");
    std::exit(1);
  }

  std::printf("starting portageapp_jali_to_flecsi...\n");

  // Construct the input Jali mesh
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> inputMesh;
  inputMesh = mf(xmin, ymin, xmax, ymax, nx, ny);

  // Construct the output FleCSI mesh
  mesh_t targetMesh;
  Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx+1, ny+1, targetMesh);

  // Create the mesh wrappers
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Portage::Flecsi_Mesh_Wrapper targetMeshWrapper(targetMesh);

  // Fill the input state with linear func
  const int nsrccells = inputMeshWrapper.num_owned_cells();
  Jali::State inputState(inputMesh);
  std::vector<double> inputData(nsrccells);
  JaliGeometry::Point cen;
  for (auto c = 0; c < nsrccells; ++c) {
    cen = inputMesh.cell_centroid(c);
    inputData[c] = cen[0] + cen[1];
  }
  inputState.add("celldata", inputMesh, Jali::Entity_kind::CELL,
          Jali::Entity_type::ALL, &(inputData[0]));
  Portage::Jali_State_Wrapper inputStateWrapper(inputState);

  // Declare the target storage
  register_state(targetMesh, "celldata", cells, real_t, flecsi::persistent);
  Portage::Flecsi_State_Wrapper targetStateWrapper(targetMesh);

  // Setup the main driver for this mesh type
  Portage::Driver<Portage::Jali_Mesh_Wrapper,
                  Portage::Jali_State_Wrapper,
                  Portage::Flecsi_Mesh_Wrapper,
                  Portage::Flecsi_State_Wrapper> d(inputMeshWrapper,
                                                   inputStateWrapper,
                                                   targetMeshWrapper,
                                                   targetStateWrapper);
  // Register the variable to be remapped - this is the same on both meshes
  std::vector<std::string> remap_fields;
  remap_fields.push_back("celldata");
  d.set_remap_var_names(remap_fields);

  // Set interpolation order
  d.set_interpolation_order(order);

  // Do the remap
  d.run();

  // Get the new data - CELL is ignored here
  double * outData;
  targetStateWrapper.get_data(Portage::CELL, "celldata", &outData);

  const auto ntarcells = targetMeshWrapper.num_owned_cells();

  double toterr = 0.0;
  for (auto c = 0; c < ntarcells; ++c) {
    Portage::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double error = cen[0] + cen[1] - outData[c];
    std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                cen[0], cen[1]);
    std::printf("  Value = % 10.6lf  Err = % lf\n",
                outData[c], error);
    toterr += error*error;
  }
  std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));

  MPI_Finalize();

  return 0;
}
