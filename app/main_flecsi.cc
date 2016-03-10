#include <vector>
// note that flecsi is missing this in some files - e.g. data.h
#include <iostream>
#include <cstdio>
#include <memory>
#include <algorithm>
#include <iterator>

#include "portage/wrappers/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "portage/wrappers/state/flecsi/flecsi_state_wrapper.h"
#include "portage/driver/driver.h"

#include "flecsi/specializations/burton/burton.h"

using mesh_t = flecsi::burton_mesh_t;

using vertex_t = flecsi::mesh_t::vertex_t;
using real_t = flecsi::mesh_t::real_t;
//using vector_t = flecsi::mesh_t::vector_t;

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
  // This isn't working because it resets inputMesh to the same as targetMesh...
  Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx, ny, inputMesh);
  // std::cout << "input mesh before target mesh creattion" << std::endl;
  // for (auto v : inputMesh.vertices()) {
  //   std::cout << "coords: " << v->coordinates() << std::endl;
  // }

  Portage::Flecsi_Mesh_Wrapper inputMeshWrapper(inputMesh);

  // Compare data from wrapper to Mesh
  auto someCellIndex = 5;
  auto thisCell = inputMesh.cells()[someCellIndex];

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

  std::vector<std::pair<double,double>> nodeCoords;
  inputMeshWrapper.cell_get_coordinates(someCellIndex, &nodeCoords);
  std::cout << "node coordinates from wrapper:" << std::endl;
  for (auto nc : nodeCoords)
    std::cout << "[ " << nc.first << " " << nc.second << " ]" << std::endl;
  std::cout << "node coordinates from mesh:" << std::endl;
  for (auto n : inputMesh.vertices(thisCell))
    std::cout << n->coordinates() << std::endl;

  return 0;

  // Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx+1, ny+1, targetMesh);
  // std::cout << "input mesh after target mesh creattion" << std::endl;
  // for (auto v : inputMesh.vertices()) {
  //   std::cout << "coords: " << v->coordinates() << std::endl;
  // }








  // Do this by hand...
  // auto dx = (xmax - xmin) / real_t(nx);
  // auto dy = (ymax - ymin) / real_t(ny);
  // std::cout << "dx dy " << dx << " " << dy << std::endl;
  // auto num_verts = (nx+1)*(ny+1);
  // inputMesh.init_parameters(num_verts);
  // std::vector<vertex_t*> verts;
  // for (auto j = 0; j < ny+1; ++j) {
  //   for (auto i = 0; i < nx+1; ++i) {
  //     auto vert = inputMesh.create_vertex(
  //         {xmin + dx*real_t(i), ymin + dy*real_t(j)});
  //     verts.push_back(vert);
  //   }
  // }
  // auto nx1 = nx + 1;
  // for (auto j = 0; j < ny; ++j) {
  //   for (auto i = 0; i < nx; ++i) {
  //     auto c = inputMesh.create_cell({verts[i + j*nx1],
  //             verts[i + 1 + j*nx1],
  //             verts[i + 1 + (j + 1)*nx1],
  //             verts[i + (j + 1)*nx1]});
  //   }
  // }
  // inputMesh.init();
  // std::cout << "input mesh " << std::endl;
  // for (auto v : inputMesh.vertices()) {
  //   std::cout << "coords: " << v->coordinates() << std::endl;
  // }

  // // target mesh
  // ++nx;
  // ++ny;
  // dx = (xmax - xmin) / real_t(nx);
  // dy = (ymax - ymin) / real_t(ny);
  // std::cout << "dx dy " << dx << " " << dy << std::endl;
  // num_verts = (nx+1)*(ny+1);
  // targetMesh.init_parameters(num_verts);
  // verts.clear();
  // for (auto j = 0; j < ny+1; ++j) {
  //   for (auto i = 0; i < nx+1; ++i) {
  //     auto vert = targetMesh.create_vertex(
  //         {xmin + dx*real_t(i), ymin + dy*real_t(j)});
  //     verts.push_back(vert);
  //   }
  // }
  // nx1 = nx + 1;
  // for (auto j = 0; j < ny; ++j) {
  //   for (auto i = 0; i < nx; ++i) {
  //     auto c = targetMesh.create_cell({verts[i + j*nx1],
  //             verts[i + 1 + j*nx1],
  //             verts[i + 1 + (j + 1)*nx1],
  //             verts[i + (j + 1)*nx1]});
  //   }
  // }
  // targetMesh.init();

  // std::cout << "input mesh " << std::endl;
  // for (auto v : inputMesh.vertices()) {
  //   std::cout << "coords: " << v->coordinates() << std::endl;
  // }
  // // mesh_t targetMesh;
  // // Portage::make_mesh_cart2d(xmin, xmax, ymin, ymax, nx+1, ny+1, targetMesh);
  // // std::cout << "input mesh " << std::endl;
  // // for (auto v : inputMesh.vertices()) {
  // //   std::cout << "coords: " << v->coordinates() << std::endl;
  // // }

  // std::cout << "target mesh " << std::endl;
  // for (auto v : targetMesh.vertices()) {
  //   std::cout << "coords: " << v->coordinates() << std::endl;
  // }



  // // Setup the wrappers
  // Portage::Flecsi_Mesh_Wrapper inputMeshWrapper(inputMesh);
  // Portage::Flecsi_Mesh_Wrapper targetMeshWrapper(targetMesh);

  // // Register the state data with the meshes; this just lets the statemanager
  // // know that there will be data
  // // NOTES: register_state is a macro defined in burton.h
  // //        cells is an enum in burton_mesh_traits.h
  // register_state(inputMesh, "celldata", cells, real_t, flecsi::persistent);
  // // NOTE:
  // // There is a naming conflict here for some reason - must change variable name
  // //  register_state(targetMesh, "celldata", cells, real_t, flecsi::persistent);
  // register_state(targetMesh, "celldata_tar", cells, real_t, flecsi::persistent);

  // // Setup the iniital data
  // // NOTE: access_state is a macro defined in burton.h
  // auto cd = access_state(inputMesh, "celldata", real_t);
  // // linear function
  // for (auto c : inputMesh.cells()) {
  //   auto cen = c->centroid();
  //   cd[c] = cen[0] + cen[1];
  // }

  // // Make the state wrappers
  // Portage::Flecsi_State_Wrapper sourceStateWrapper(inputMesh);
  // Portage::Flecsi_State_Wrapper targetStateWrapper(targetMesh);

  // // Build the main driver data for this remap
  // Portage::Driver<Portage::Flecsi_Mesh_Wrapper,
  //                 Portage::Flecsi_State_Wrapper> d(Portage::CELL,
  //                                                  inputMeshWrapper,
  //                                                  sourceStateWrapper,
  //                                                  targetMeshWrapper,
  //                                                  targetStateWrapper);

  // // Register the variables to be remapped
  // std::vector<std::string> src_remap_fields, tar_remap_fields;
  // src_remap_fields.push_back("celldata");
  // tar_remap_fields.push_back("celldata_tar");
  // d.set_remap_var_names(src_remap_fields, tar_remap_fields);

  // d.set_interpolation_order(1);

  // d.run();

  // for (auto c: inputMesh.cells()) {
  //   auto vol = c->area();
  //   //    auto coords = c->coordinates();
  //   //    std::cout << "area: " << vol << " " << coords << std::endl;
  //   std::cout << "area: " << vol << std::endl;
  // }

  // auto cell5 = inputMesh.cells()[4];
  // auto vert5 = inputMesh.vertices()[4];

  // // std::cout << "cell 5 " << &cell5 << std::endl;
  // // for (auto v : inputMesh.vertices(cell5)) {
  // //     std::cout << "coords: " << v->coordinates() << std::endl;
  // // }

  // std::cout << "testing..." << std::endl;
  // // inputMesh.cells(cell5);   Fails
  // inputMesh.vertices(cell5);
  // inputMesh.cells(vert5);
  // inputMesh.wedges(cell5);
  // std::cout << "done testing..." << std::endl;

  // for (auto vi : inputMesh.cells(cell5))
  //   std::cout << "ids: " << vi << std::endl;



  // std::cout << "meshwrapper dim: " << inputMeshWrapper.space_dimension()
  //           << std::endl;

  // std::vector<std::pair<double,double>> vert_coords;
  // inputMeshWrapper.cell_get_coordinates(4, &vert_coords);
  // for (auto c : vert_coords)
  //   std::cout << c.first << " " << c.second << std::endl;
  // std::vector<int> vertIds;
  // inputMeshWrapper.cell_get_nodes(4, &vertIds);
  // for (auto vi : vertIds)
  //   std::cout << "id: " << vi << std::endl;


  return 0;
}
