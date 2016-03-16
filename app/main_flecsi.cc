#include <vector>
// note that flecsi is missing this in some files - e.g. data.h
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

using mesh_t = flecsi::burton_mesh_t;

using real_t = flecsi::mesh_t::real_t;

/*!
  @file main_flecsi.cc
  @brief A simple application that drives our remap routines on a FleCSI mesh.

  The program is used to showcase our capabilities with various types of remap
  operations on a 2d FleCSI mesh with simple linear or quadratic data.For the
  cases of remapping linear data with a second-order interpolator, the L2 norm
  output at the end should be identically zero.
 */

/// @TODO add multiple examples

void print_usage() {
  std::printf("usage: portageapp_flecsi ncellsx ncellsy\n");
}

int main(int argc, char** argv) {

  if (argc <= 2) {
    print_usage();
    return 0;
  }

  // number of CELLS in x and y for input mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  // dimensions of meshes
  real_t xmin = 0.0, xmax = 1.0;
  real_t ymin = 0.0, ymax = 1.0;

  std::printf("starting portageapp_flecsi...\n");

  ///@TODO swap to a smart pointer

  // Setup the meshes
  mesh_t inputMesh, targetMesh;
  Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx, ny, inputMesh);
  // Example of how to get information directly from mesh object
  // std::cout << "input mesh before target mesh creattion" << std::endl;
  // for (auto v : inputMesh.vertices()) {
  //   std::cout << "coords: " << v->coordinates() << std::endl;
  // }

  Portage::Flecsi_Mesh_Wrapper inputMeshWrapper(inputMesh);

  // Compare data from wrapper to Mesh
  auto someCellIndex = 5;
  auto thisCell = inputMesh.cells()[someCellIndex];

  // Examples of getting data from mesh object or wrapper
  std::vector<int> nodeIDs;
  inputMeshWrapper.cell_get_nodes(someCellIndex, &nodeIDs);
  std::cout << "node IDs from wrapper: ";
  std::copy(nodeIDs.begin(), nodeIDs.end(),
            std::ostream_iterator<int>(std::cout, " "));
  std::cout << std::endl;
  std::cout << "node IDs from mesh:    ";
  for (auto n : inputMesh.vertices(thisCell))
    std::cout << n.global_id() << " ";
  std::cout << std::endl;

  std::vector<double> centroid;
  inputMeshWrapper.cell_centroid(someCellIndex, &centroid);
  std::cout << "centroid from wrapper: [ ";
  std::copy(centroid.begin(), centroid.end(),
            std::ostream_iterator<double>(std::cout, " "));
  std::cout << " ]" << std::endl;
  std::cout << "centroid from mesh:    [ ";
  std::cout << thisCell->centroid() << std::endl;

  //  std::vector<std::pair<double,double>> nodeCoords;
  std::vector<Portage::Point2> nodeCoords;
  inputMeshWrapper.cell_get_coordinates(someCellIndex, &nodeCoords);
  std::cout << "node coordinates from wrapper:" << std::endl;
  for (auto nc : nodeCoords)
    std::cout << nc << std::endl;
    // std::cout << "[ " << nc.first << " " << nc.second << " ]" << std::endl;
  std::cout << "node coordinates from mesh:" << std::endl;
  for (auto n : inputMesh.vertices(thisCell))
    std::cout << n->coordinates() << std::endl;

  return 0;

  // This is broken because of the singleton issue in FleCSI state manager.
  // Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx+1, ny+1, targetMesh);
  // std::cout << "input mesh after target mesh creattion" << std::endl;
  // for (auto v : inputMesh.vertices()) {
  //   std::cout << "coords: " << v->coordinates() << std::endl;
  // }

  return 0;
}
