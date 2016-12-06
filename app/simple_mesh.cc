/*----------------------------------------------------------------------------*
 * Copyright (c) 2016 Los Alamos National Security, LLC
 * All rights reserved.
 *----------------------------------------------------------------------------*/

#include <vector>
#include <string>
#include <iostream>
#include <cstdio>
#include <memory>

#include "portage/support/portage.h"
#include "portage/driver/driver.h"
#include "portage/wrappers/mesh/simple_mesh/simple_mesh.h"
#include "portage/wrappers/state/simple_state/simple_state.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/simple_mesh/simple_state.h"

using Portage::Simple_Mesh;
using Portage::Simple_State;
using Portage::Simple_Mesh_Wrapper;
using Portage::Simple_State_Wrapper;


//////////////////////////////////////////////////////////////////////
// Helper routines and data structures

struct example_properties {
  example_properties(const int order, const bool cell_centered,
		     const bool linear, const bool conformal = true)
    : order(order), cell_centered(cell_centered), linear(linear),
      conformal(conformal) { }
  int order;           // interpolation order in example
  bool cell_centered;  // is this example a cell-centered remap?
  bool linear;         // is this example a remap of linear data?
  bool conformal;      // are the two meshes boundary-conformal?
};

// Use this to add new problems.  If needed, we can extend the
// example_properties struct to contain more information.
std::vector<example_properties> setup_examples() {
  std::vector<example_properties> examples;

  // Cell-centered remaps:

  // 1st order cell-centered remap of linear func
  examples.emplace_back(1, true, true);

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
    std::printf("  %d: %s order %s-centered remap of %s func %s\n",
		i, (example.order=1) ? "1st" : "2nd",
		example.cell_centered ? "cell" : "node",
		example.linear ? "linear" : "quadratic",
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

  std::cout << "starting simple_mesh app..." << std::endl;
  std::cout << "running example " << example_num << std::endl;

  example_properties example = setup_examples()[example_num];

  std::shared_ptr<Simple_mesh> inputMesh, targetMesh;

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

    const int nsrccells = inputMeshWrapper.num_owned_cells();
    const int ntarcells = targetMeshWRapper.num_owned_cells();

    Simple_State srcState(inputMesh);
    std::vector<double> srcData(nsrccells);

    if (example.linear) {
      for (int c(0); c < nsrccells; ++c) {
	Portage::Point<3> cen = inputMesh->cell_centroid(c);
	srcData[c] = cen[0] + cen[1] + cen[2];
      }
    } else {
      for (int c(0); c < nsrccells; ++c) {
	Portage::Point<3> cen = inputMesh->cell_centroid(c);
	srcData[c] = cen[0]*cen[0] + cen[1]*cen[1] + cen[2]*cen[2];
      }
    }

    srcState.add("celldata", inputMesh, Portage::Entity_kind::CELL,
		 Portage::Entity_type::ALL, &(srcData[0]));
    Simple_State_Wrapper srcStateWrapper(srcState);

    Simple_State tarState(targetMesh);
    std::vector<double> tarData(ntarcells, 0.0);
    auto& cellvecout = tarState.add("celldata", targetMesh,
				    Portage::Entity_kind::CELL,
				    Portage::Entity_type::ALL,
				    &(tarData[0]));
    Simple_State_Wrapper tarStateWrapper(tarState);

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
	d(inputMeshWrapper, srcStateWrapper, targetMeshWraper, tarStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run();
    } else {  // 2nd order
      Portage::Driver<
	Portage::SearchKDTree,
	Portage::IntersectR3D,
	Portage::Interpolate_2ndOrder,
	3,
	Simple_Mesh_Wrapper,
	Simple_State_Wrapper>
	d(inputMeshWrapper, srcStateWrapper, targetMeshWrapper, tarStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run();
    }

    double toterr = 0.0;
    for (int c(0); c < ntarcells; ++c) {
      Portage::Point<3> ccen = targetMesh->cell_centroid(c);

      double error;
      if (example.linear) {
	error = ccen[0] + ccen[1] + ccen[2] - cellvecout[c];
	std::cout << "error is " << error << std::endl;
      } else {
	error = ccen[0]*ccen[0] + ccen[1]*ccen[1] + ccen[2]*ccen[2]
	  - cellvecout[c];
	std::cout << "error is " << error << std::endl;
      }

      std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf,% 5.3lf)", c,
		  ccen[0], ccen[1], ccen[2]);
      std::printf("  Value = % 10.6lf  Err = % lf\n",
		  cellvecout[c], error);

      toterr += error*error;
    }

    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
  } else {
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
    Simple_Mesh_Wrapper targetMeshWRapper(*targetMesh);

    const int nsrcnodes = inputMeshWrapper.num_owned_nodes();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();

    Simple_State srcState(inputMesh);
    std::vector<double> srcData(nsrcnodes);

    Portage::Point<3> nodxyz;
    for (int i(0); i < ntarnodes; ++i) {
      inputMeshWrapper.node_get_coordinates(i, &nodexyz);
      srcData[i] = nodexyz[o]*nodexyz[0] + nodexyz[1]*nodexyz[1]
	+ nodexyz[2]*nodexyz[2];
    }

    srcState.add("nodedata", inputMesh, Portage::Entity_kind::NODE,
		 Portage::Entity_type::ALL, &(srcData[0]));
    Simple_State_Wrapper srcStateWrapper(srcState);

    Simple_State tarState(targetMesh);
    auto& nodevecout = tarState.add("nodedata", targetMesh,
				    Portage::Entity_kind::NODE,
				    Portage::Entity_type::ALL,
				    0.0);
    Simple_State_Wrapper tarStateWrapper(tarState);

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
	d(inputMeshWrapper, srcStateWrapper,
	  targetMeshWrapper, tarStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run();
    } else {
      Portage::Driver<
	Portage::SearchKDTree,
	Portage::IntersectR3D,
	Portage::Interpolate_2ndOrder,
	3,
	Simple_Mesh_Wrapper,
	Simple_State_Wrapper>
	d(inputMeshWrapper, srcStateWrapper,
	  targetMeshWrapper, tarStateWrapper);
      d.set_remap_var_names(remap_fields);
      d.run();
    }

    double toterr(0.0);
    double stdval, err;
    Portage::Point<3> nodexyz;
    for (int i(0); i < ntarnodes; ++i) {
      targetMeshWrapper.node_get_coordinates(i, &nodexyz);
      stdval = nodexyz[0]*nodexyz[0] + nodexyz[1]*nodexyz[1]
	+ nodexyz[2]*nodexyz[2];
      err = fabs(stdval-nodevecout[i]);
      std::printf("Node=% 4d Coords = (% 5.3lf, % 5.3lf, % 5.3lf) ", i,
		  nodexyz[0], nodexyz[1], nodexyz[2]);
      std::printf("  Value = %10.6lf Err = % lf\n", nodevecout[i], err);
      toterr += err;
    }
    std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
  }
  std::cout << "finishing simple_mesh app..." << std::endl;
}
