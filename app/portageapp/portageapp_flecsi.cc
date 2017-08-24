/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <vector>
#include <iostream>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <iterator>

#include "portage/support/Point.h"
#include "portage/driver/driver.h"
#include "portage/wonton/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "portage/wonton/state/flecsi/flecsi_state_wrapper.h"

#include "flecsi-sp.h"
#include "flecsi/io/io.h"
#include "flecsi-sp/burton/burton.h"
#include "flecsi-sp/burton/factory.h"
#include "flecsi-sp/burton/burton_io_exodus.h"

namespace math = flecsi::sp::math;
namespace mesh = flecsi::sp::burton;

using mesh_t = mesh::burton_mesh_2d_t;
using real_t = typename mesh_t::real_t;

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

  // length of meshes
  constexpr size_t lenx = 1.;
  constexpr size_t leny = 1.;

  std::printf("starting portage_flecsiapp...\n");

  // Setup::flecsi meshes
  mesh_t inputMesh, targetMesh;
  auto inputMesh  = mesh::box<mesh_t>(nx, ny, 0, 0, lenx, leny);
  auto targetMesh = mesh::box<mesh_t>(nx+1, ny+1, 0, 0, lenx, leny);

  // Setup::portage-flecsi mesh wrappers
  wonton::flecsi_mesh_t inputMeshWrapper(inputMesh);
  wonton::flecsi_mesh_t targetMeshWrapper(targetMesh);
  
  // Setup::portage-flecsi  state wrappers
  wonton::flecsi_state_t inputStateWrapper(inputMesh);
  wonton::flecsi_state_t targetStateWrapper(targetMesh);

  //Register data on input and target mesh
  flecsi_register_data(inputMesh, hydro, cell_data, real_t, dense, 1, cells);
  flecsi_register_data(targetMesh, hydro, cell_data, real_t, dense, 1, cells);

  //Get accessors to data on input and target mesh
  auto inputMeshAccessor  = flecsi_get_accessor(intputMesh, hydro, cell_data, real_t, dense, 0);
  auto targetMeshAccessor = flecsi_get_accesor(targetMesh, hydro, cell_data, real_t, dense, 0); 

  inputMeshAccessor.attributes().set(persistent);
  targetMeshAccessor.attributes().set(persistent);

  // Fill data on input mesh with linear function
  for (auto c : inputMesh.cells()) {
    auto cen = c->centroid();
    inputMeshAccessor[c] = cen[0] + cen[1];
  }
  
  // Setup the main driver for this mesh type 2
  if(order==2){
    Portage::Driver<
      Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_2ndOrder,
      2,
      wonton::flecsi_mesh_t,
      wonton::flecsi_state_t>
      d(inputMeshWrapper, inputStateWrapper, targetMeshWrapper, targetStateWrapper);
  
    // Declare which variables are remapped
    std::vector<std::string> varnames(1, "cell_data");
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
      wonton::flecsi_mesh_t,
      wonton::flecsi_state_t>
      d(inputMeshWrapper, inputStateWrapper, targetMeshWrapper, targetStateWrapper);
  
     // Declare which variables are remapped
     std::vector<std::string> varnames(1, "cell_data");
     d.set_remap_var_names(varnames);

     // Do the remap
     d.run(false);
  }


  // Get the new data and calculate the error
  double toterr = 0.0;
  for (auto c : targetMesh.cells()) {
    auto cen = c->centroid();
    double error = cen[0] + cen[1] - targetMeshAccessor[c];
    std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c.id(),cen[0], cen[1]);
    std::printf("  Value = % 10.6lf  Err = % lf\n",targetMeshAccessor[c], error);
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
