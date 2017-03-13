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
#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <memory>
#include <stdexcept>

#ifdef ENABLE_MPI
#include <mpi.h>
#else
#define PORTAGE_SERIAL_ONLY
#endif

#include "portage/support/portage.h"
#include "portage/driver/driver.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/simple_mesh/simple_state.h"
#include "portage/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/wrappers/state/simple_state/simple_state_wrapper.h"

using Portage::Simple_Mesh;
using Portage::Simple_State;
using Portage::Simple_Mesh_Wrapper;
using Portage::Simple_State_Wrapper;



//////////////////////////////////////////////////////////////////////
// Helper routines and data structures

struct example_properties {
  example_properties(const int order, const bool cell_centered,
                     const int field_order, const bool conformal = true)
      : order(order), cell_centered(cell_centered), field_order(field_order),
        conformal(conformal) { }
  int order;           // interpolation order in example
  bool cell_centered;  // is this example a cell-centered remap?
  int field_order;     // order of the field to map
  bool conformal;      // are the two meshes boundary-conformal?
};

// Use this to add new problems.  If needed, we can extend the
// example_properties struct to contain more information.
std::vector<example_properties> setup_examples() {
  std::vector<example_properties> examples;

  // Cell-centered remaps:

  // 1st order cell-centered remap of a const func
  examples.emplace_back(1, true, 0);

  // 1st order cell-centered remap of linear func
  examples.emplace_back(1, true, 1);

  // 1st order cell-centered remap of quad func
  examples.emplace_back(1, true, 2);

  // 2nd order cell-centered remap of linear func
  examples.emplace_back(2, true, 1);

  // 2nd order cell-centered remap of quad func
  examples.emplace_back(2, true, 2);

  // 2nd order cell-centered remap of linear func on non-conformal mesh
  examples.emplace_back(2, true, 1, false);

  // Node-centered remaps

  // 1st order node-centered remap of a const func
  examples.emplace_back(1, false, 0);

  // 1st order node-centered remap of linear func
  examples.emplace_back(1, false, 1);

  // 1st order node-centered remap of quad func
  examples.emplace_back(1, false, 2);

  // 2nd order node-centered remap of linear func
  examples.emplace_back(2, false, 1);

  // 2nd order node-centered remap of quad func
  examples.emplace_back(2, false, 2);

  // 2nd order node-centered remap of linear func on non-conformal mesh
  examples.emplace_back(2, false, 1, false);
  return examples;
}

void usage() {
  auto examples = setup_examples();
  std::cout << "Usage: simple_mesh example-number nsourcecells ntargetcells"
            << std::endl;
  std::cout << "List of example numbers:" << std::endl;
  int i = 0;
  bool separated = false;
  std::cout << "CELL-CENTERED EXAMPLES" << std::endl;
  for (const auto &example : examples) {
    if (!separated && !example.cell_centered) {
      std::cout << std::endl;
      std::cout << "NODE-CENTERED EXAMPLES:" << std::endl;
      separated = true;
    }
    std::printf("  %d: %s order %s-centered remap of order %d func %s\n",
                i, (example.order == 1) ? "1st" : "2nd",
                example.cell_centered ? "cell" : "node",
                example.field_order,
                !example.conformal ? "on non-conformal mesh" : "");
    ++i;
  }
}

int main(int argc, char** argv) {
  int example_num, n_source, n_target;
  if (argc <= 3) {
    usage();
    return 0;
  }

  example_num = atoi(argv[1]);
  n_source = atoi(argv[2]);
  n_target = atoi(argv[3]);

  // Even though Simple_Mesh is serial only, we still need to initialize MPI
  // for other Portage code.

#ifdef ENABLE_MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  if (numpe > 1) {
    std::cerr << "Simple Mesh is only designed for a single process!"
              << std::endl;
    return 1;
  }
#endif

  std::cout << "starting simple_mesh app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  example_properties example = setup_examples()[example_num];

  std::shared_ptr<Simple_Mesh> inputMesh, targetMesh;

  // Value for the constant function
  const double fval = 42.0;

  if (example.cell_centered) {
    inputMesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                              1.0, 1.0, 1.0,
                                              n_source, n_source, n_source);
    if (example.conformal) {
      targetMesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                                 1.0, 1.0, 1.0,
                                                 n_target, n_target, n_target);
    } else {
      double dx = 1.0/static_cast<double>(n_target);
      targetMesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                                 1.0+1.5*dx, 1.0, 1.0,
                                                 n_target, n_target, n_target);
    }

    Simple_Mesh_Wrapper inputMeshWrapper(*inputMesh);
    Simple_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    const int ninpcells = inputMeshWrapper.num_owned_cells();
    const int ntarcells = targetMeshWrapper.num_owned_cells();

    Simple_State inputState(inputMesh);
    std::vector<double> inputData(ninpcells);

    Portage::Point<3> cen;

    switch (example.field_order) {
      case 0:
        for (int c(0); c < ninpcells; ++c) {
          inputData[c] = fval;
        }
        break;
      case 1:
        for (int c(0); c < ninpcells; ++c) {
          inputMeshWrapper.cell_centroid(c, &cen);
          inputData[c] = cen[0] + cen[1] + cen[2];
        }
        break;
      case 2:
        for (int c(0); c < ninpcells; ++c) {
          inputMeshWrapper.cell_centroid(c, &cen);
          inputData[c] = cen[0]*cen[0] + cen[1]*cen[1] + cen[2]*cen[2];
        }
        break;
      default:
        throw std::runtime_error("Unknown field order!");
    }

    inputState.add("celldata", Portage::Entity_kind::CELL, &(inputData[0]));
    Simple_State_Wrapper inputStateWrapper(inputState);

    Simple_State targetState(targetMesh);
    std::vector<double> targetData(ntarcells, 0.0);
    auto& cellvecout = targetState.add("celldata",
                                       Portage::Entity_kind::CELL,
                                       &(targetData[0]));
    Simple_State_Wrapper targetStateWrapper(targetState);

    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");

    if (example.order == 1) {
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_1stOrder,
        3,
        Simple_Mesh_Wrapper,
        Simple_State_Wrapper>
          d(inputMeshWrapper, inputStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(false);
    } else {  // 2nd order
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_2ndOrder,
        3,
        Simple_Mesh_Wrapper,
        Simple_State_Wrapper>
          d(inputMeshWrapper, inputStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(false);
    }

    double toterr = 0.0;
    Portage::Point<3> ccen;
    for (int c(0); c < ntarcells; ++c) {
      targetMeshWrapper.cell_centroid(c, &ccen);

      double error;
      switch (example.field_order) {
        case 0:
          error = fval - cellvecout[c];
          break;
        case 1:
          error = ccen[0] + ccen[1] + ccen[2] - cellvecout[c];
          break;
        case 2:
          error = ccen[0]*ccen[0] + ccen[1]*ccen[1] + ccen[2]*ccen[2]
              - cellvecout[c];
          break;
        default:
          throw std::runtime_error("Unknown field_order!");
      }

      std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf,% 5.3lf)", c,
                  ccen[0], ccen[1], ccen[2]);
      std::printf("  Value = % 10.6lf  Err = % lf\n",
                  cellvecout[c], error);

      toterr += error*error;
    }

    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));

    // Create the output file for comparison.
    // The `static_cast` is a workaround for Intel compiler, which are mising a
    // `to_string` function for ints.
    std::ofstream fout("field"
                       + std::to_string(static_cast<long long>(example_num))
                       + ".txt");
    fout << std::scientific;
    fout.precision(17);
    for (int c(0); c < ntarcells; ++c)
      fout << c << " " << cellvecout[c] << std::endl;
  } else {
    // Node-centered remaps.
    inputMesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                              1.0, 1.0, 1.0,
                                              n_source, n_source, n_source);
    if (example.conformal) {
      targetMesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                                 1.0, 1.0, 1.0,
                                                 n_target, n_target, n_target);
    } else {
      double dx = 1.0/static_cast<double>(n_target);
      targetMesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                                 1.0+1.5*dx, 1.0, 1.0,
                                                 n_target, n_target, n_target);
    }

    Simple_Mesh_Wrapper inputMeshWrapper(*inputMesh);
    Simple_Mesh_Wrapper targetMeshWrapper(*targetMesh);

    const int ninpnodes = inputMeshWrapper.num_owned_nodes();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();

    Simple_State inputState(inputMesh);
    std::vector<double> inputData(ninpnodes);

    Portage::Point<3> nodexyz;
    switch (example.field_order) {
      case 0:
        for (int i(0); i < ninpnodes; ++i) {
          inputData[i] = fval;
        }
        break;
      case 1:
        for (int i(0); i < ninpnodes; ++i) {
          inputMeshWrapper.node_get_coordinates(i, &nodexyz);
          inputData[i] = nodexyz[0] + nodexyz[1] + nodexyz[2];
        }
        break;
      case 2:
        for (int i(0); i < ninpnodes; ++i) {
          inputMeshWrapper.node_get_coordinates(i, &nodexyz);
          inputData[i] = nodexyz[0]*nodexyz[0] + nodexyz[1]*nodexyz[1]
              + nodexyz[2]*nodexyz[2];
        }
        break;
      default:
        throw std::runtime_error("Unknown field_order!");
    }

    inputState.add("nodedata", Portage::Entity_kind::NODE, &(inputData[0]));
    Simple_State_Wrapper inputStateWrapper(inputState);

    Simple_State targetState(targetMesh);
    auto& nodevecout = targetState.add("nodedata",
                                       Portage::Entity_kind::NODE,
                                       0.0);
    Simple_State_Wrapper targetStateWrapper(targetState);

    std::vector<std::string> remap_fields;
    remap_fields.push_back("nodedata");

    if (example.order == 1) {
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_1stOrder,
        3,
        Simple_Mesh_Wrapper,
        Simple_State_Wrapper>
          d(inputMeshWrapper, inputStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(false);
    } else {
      Portage::Driver<
        Portage::SearchKDTree,
        Portage::IntersectR3D,
        Portage::Interpolate_2ndOrder,
        3,
        Simple_Mesh_Wrapper,
        Simple_State_Wrapper>
          d(inputMeshWrapper, inputStateWrapper,
            targetMeshWrapper, targetStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run(false);
    }

    double toterr(0.0);
    double stdval, err;
    Portage::Point<3> nnodexyz;
    for (int i(0); i < ntarnodes; ++i) {
      targetMeshWrapper.node_get_coordinates(i, &nnodexyz);

      double err;
      switch (example.field_order) {
        case 0:
          err = fval - nodevecout[i];
          break;
        case 1:
          err = nnodexyz[0] + nnodexyz[1] + nnodexyz[2] - nodevecout[i];
          break;
        case 2:
          err = nnodexyz[0]*nnodexyz[0] + nnodexyz[1]*nnodexyz[1]
              + nnodexyz[2]*nnodexyz[2] - nodevecout[i];
          break;
        default:
          throw std::runtime_error("Unknown field_order!");
      }

      std::printf("Node=% 4d Coords = (% 5.3lf, % 5.3lf, % 5.3lf) ", i,
                  nnodexyz[0], nnodexyz[1], nnodexyz[2]);
      std::printf("  Value = %10.6lf Err = % lf\n", nodevecout[i], err);
      toterr += err*err;
    }
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));

    // Create the output file for comparison.
    // The `static_cast` is a workaround for Intel compiler, which are missing a
    // `to_string` function for ints.
    std::ofstream fout("field"
                       + std::to_string(static_cast<long long>(example_num))
                       + ".txt");
    fout << std::scientific;
    fout.precision(17);
    for (int n(0); n < ntarnodes; ++n)
      fout << n << " " << nodevecout[n] << std::endl;
  }
  std::cout << "finishing simple_mesh app..." << std::endl;

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif
}
