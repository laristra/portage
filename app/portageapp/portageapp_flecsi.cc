/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <vector>
#include <fstream>
#include <iostream>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <iterator>
#include <cstring>

#include "flecsi-sp.h"
#include "flecsi/io/io.h"
#include "flecsi-sp/burton/burton.h"
#include "flecsi-sp/burton/factory.h"
#include "flecsi-sp/burton/burton_io_exodus.h"

#include "wonton/support/wonton.h"
#include "wonton/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "wonton/state/flecsi/flecsi_state_wrapper.h"

#include "portage/driver/mmdriver.h"
#include "portage/support/mpi_collate.h"


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
  std::printf("usage: portage_flecsiapp ncellsx ncellsy [order=1] [dump_output=y]\n");
}

int main(int argc, char** argv) {
  if (argc <= 2) {
    print_usage();
    return 0;
  }

  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
     MPI_Init(&argc, &argv);
  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // number of CELLS in x and y for input mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  auto order = (argc > 3) ? atoi(argv[3]) : 1;
  auto dump_output = (argc > 4) ? ((strcmp(argv[4],"n") == 0) || (strcmp(argv[4],"N") == 0)) : true;

  // length of meshes
  constexpr size_t lenx = 1.;
  constexpr size_t leny = 1.;

  std::printf("starting portage_flecsiapp...\n");
  std::cout << "nx: " << nx << std::endl;
  std::cout << "ny: " << ny << std::endl;
  std::cout << "order: " << order << std::endl;
  std::cout << "dump_output: " << dump_output << std::endl;

  // Setup::flecsi meshes
  mesh_t inputMesh, targetMesh;
  inputMesh  = mesh::box<mesh_t>(nx, ny, 0, 0, lenx, leny);
  targetMesh = mesh::box<mesh_t>(nx+1, ny+1, 0, 0, lenx, leny);

  // Setup::portage-flecsi mesh wrappers
  Wonton::flecsi_mesh_t<mesh_t> inputMeshWrapper(inputMesh);
  Wonton::flecsi_mesh_t<mesh_t> targetMeshWrapper(targetMesh);

  // Setup::portage-flecsi  state wrappers
  Wonton::flecsi_state_t<mesh_t> inputStateWrapper(inputMesh);
  Wonton::flecsi_state_t<mesh_t> targetStateWrapper(targetMesh);

  // Register data on input and target mesh
  flecsi_register_data(inputMesh, hydro, cell_data, real_t, dense, 1, cells);
  flecsi_register_data(targetMesh, hydro, cell_data, real_t, dense, 1, cells);

  // Get accessors to data on input and target mesh
  auto inputMeshAccessor = flecsi_get_accessor(inputMesh, hydro, cell_data,
                                               real_t, dense, 0);
  auto targetMeshAccessor = flecsi_get_accessor(targetMesh, hydro, cell_data,
                                                real_t, dense, 0);

  inputMeshAccessor.attributes().set(persistent);
  targetMeshAccessor.attributes().set(persistent);

  // Fill data on input mesh with linear function
  for (auto c : inputMesh.cells()) {
    auto cen = c->centroid();

    if (order == 1)
      inputMeshAccessor[c] = cen[0] + cen[1];
    else
      inputMeshAccessor[c] = cen[0]*cen[0] + cen[1]*cen[1];
  }

  // Setup the main driver for this mesh type 2
  if (order == 2) {
    Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_2ndOrder,
      2,
      Wonton::flecsi_mesh_t<mesh_t>,
      Wonton::flecsi_state_t<mesh_t> >
      d(inputMeshWrapper, inputStateWrapper,
        targetMeshWrapper, targetStateWrapper);

    // Declare which variables are remapped
    std::vector<std::string> varnames(1, "cell_data");
    d.set_remap_var_names(varnames);

    // Do the remap
    d.run();  // executor argument defaults to false -> serial
  }

  // Setup the main driver for this mesh type
  if (order == 1) {
     Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_1stOrder,
      2,
      Wonton::flecsi_mesh_t<mesh_t>,
      Wonton::flecsi_state_t<mesh_t> >
      d(inputMeshWrapper, inputStateWrapper,
        targetMeshWrapper, targetStateWrapper);

     // Declare which variables are remapped
     std::vector<std::string> varnames(1, "cell_data");
     d.set_remap_var_names(varnames);

     // Do the remap
     d.run();  // executor arguments defaults to false -> serial
  }


  // Get the new data and calculate the error
  double toterr = 0.0;
  for (auto c : targetMesh.cells()) {
    auto cen = c->centroid();
    double error;
    if (order ==1)
       error = cen[0] + cen[1] - targetMeshAccessor[c];
    else
       error = cen[0]*cen[0] + cen[1]*cen[1] - targetMeshAccessor[c];
    std::printf("Cell=%zu Centroid = (% 5.3lf,% 5.3lf)",
                c.id(), cen[0], cen[1]);
    std::printf("  Value = % 10.6lf  Err = % lf\n",
                targetMeshAccessor[c], error);
    toterr += error*error;
  }
  std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));

  if (dump_output) {
    // Dump the output for comparison
    std::string entstr("cell");
    // we only allow 2d
    auto lorder = static_cast<long long>(order);
    std::string fieldfilename = "flecsi_field_2d_" +
        entstr + "_f" + std::to_string(lorder) + "_r" +
        std::to_string(lorder) + ".txt";

    std::vector<int> idx, lgid;
    std::vector<real_t> lvalues;
    for (auto c : targetMesh.cells()) {
      lgid.push_back(c.id());
      lvalues.push_back(targetMeshAccessor[c]);
    }
    // Sort things, even though we should only be on a single proc
    Portage::argsort(lgid, idx);
    Portage::reorder(lgid, idx);
    Portage::reorder(lvalues, idx);

    std::ofstream fout(fieldfilename);
    fout << std::scientific;
    fout.precision(17);

    for (unsigned int i = 0; i < lgid.size(); ++i)
      fout << lgid[i] << " " << lvalues[i] << std::endl;
  }  // if (dump_output)

  MPI_Finalize();
}
