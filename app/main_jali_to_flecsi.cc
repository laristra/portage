#include <vector>
// note that flecsi is missing this in some files - e.g. data.h
#include <iostream>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <iterator>

#include <mpi.h>

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "portage/wrappers/state/flecsi/flecsi_state_wrapper.h"
#include "portage/driver/driver.h"

// FleCSI includes
#include "flecsi/specializations/burton/burton.h"
#include "flecsi/specializations/burton/burton_io_exodus.h"
// IO needs fixed
//#include "flecsi/io/io.h"
//#include "flecsi/io/io_exodus.h"
// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"
// Needed for IO
// #define MSTK_HAVE_MPI
// #include "Mesh_MSTK.hh"

using mesh_t = flecsi::burton_mesh_t;

using vertex_t = flecsi::mesh_t::vertex_t;
using real_t = flecsi::mesh_t::real_t;
//using vector_t = flecsi::mesh_t::vector_t;

/*!
  @file main_jali_to_flecsi.cc
  @brief A simple application that drives our remap routines from Jali mesh to
  a FleCSI mesh.

  The program is used to showcase our capabilities with various types of remap
  operations on 2d meshes with simple linear or quadratic data.For the
  cases of remapping linear data with a second-order interpolator, the L2 norm
  output at the end should be identically zero.
 */

/// @TODO add multiple examples

void print_usage() {
  std::printf("usage: portageapp_flecsi ncellsx ncellsy order=1\n");
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
  std::vector<double> cen;
  for (auto c = 0; c < nsrccells; ++c) {
    inputMeshWrapper.cell_centroid(c, &cen);
    inputData[c] = cen[0] + cen[1];
  }
  inputState.add("celldata", inputMesh, Jali::Entity_kind::CELL,
          Jali::Parallel_type::ALL, &(inputData[0]));
  Portage::Jali_State_Wrapper inputStateWrapper(inputState);

  // Declare the target storage
  register_state(targetMesh, "celldata", cells, real_t, flecsi::persistent);
  Portage::Flecsi_State_Wrapper targetStateWrapper(targetMesh);

  // Setup the main driver for this mesh type
  Portage::Driver<Portage::Jali_Mesh_Wrapper,
                  Portage::Jali_State_Wrapper,
                  Portage::Flecsi_Mesh_Wrapper,
                  Portage::Flecsi_State_Wrapper> d(Portage::CELL,
                                                   inputMeshWrapper,
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
    std::vector<double> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double error = cen[0] + cen[1] - outData[c];
    std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                cen[0], cen[1]);
    std::printf("  Value = % 10.6lf  Err = % lf\n",
                outData[c], error);
    toterr += error*error;
  }
  std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));


  // output -- flecsi output not working at the moment
  // inputState.export_to_mesh();
  // dynamic_cast<Jali::Mesh_MSTK*>(inputMesh)->write_to_exodus_file("input.exo");
  // std::cout << "done inputState mesh" << std::endl;
  // std::cout << flecsi::burton_exodus_exo_registered << std::endl;
  // // FlecSI mesh
  // auto writer = flecsi::io_exodus_t<mesh_t>();
  // writer.write("output.exo", targetMesh);
  //  execute(flecsi::write<mesh_t>, "output.exo", targetMesh);
  // //  flecsi::io_exodus_t<mesh_t>::write("output.exo", targetMesh);

  MPI_Finalize();

  return 0;
}
