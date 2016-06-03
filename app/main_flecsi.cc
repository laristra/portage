#include <vector>
#include <iostream>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <iterator>

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
  Portage::Driver<Portage::Flecsi_Mesh_Wrapper,
                  Portage::Flecsi_State_Wrapper> d(inputMeshWrapper,
                                                   inputStateWrapper,
                                                   targetMeshWrapper,
                                                   targetStateWrapper);

  // Declare which variables are remapped
  std::vector<std::string> varnames(1, "celldata");
  d.set_remap_var_names(varnames);

  // Declare interpolation order
  d.set_interpolation_order(order);

  // Do the remap
  d.run();

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
