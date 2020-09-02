/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <vector>

#include "gtest/gtest.h"

// wonton includes
#include "wonton/support/CoordinateSystem.h"
#include "wonton/support/Point.h"
#include "wonton/support/wonton.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

// portage includes
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/support/portage.h"
#include "portage/driver/coredriver.h"

class Order2Test : public ::testing::TestWithParam<int> {
 protected:
  int itest;
};

/// Second order interpolation of cell-centered field in 2D
using CoordSys = Wonton::CylindricalAxisymmetricCoordinates;

TEST_P(Order2Test, SimpleMesh) {
  int itest = GetParam();

  int mx(4), nx(5);
  double TOL_L8(1e-12), TOL_L2(1e-12);
  std::vector<double> coefs(4);
  if (itest == 1) {
    coefs = { 1.25, 0.0, 0.0, 0.0 };  // constant function
  } else if (itest == 2) {
    coefs = { 0.0, 1.0, 0.0, 0.0 }; 
    TOL_L8 = 0.01;
    TOL_L2 = 0.01;
  } else { 
    mx *= std::pow(2, itest - 3);
    nx *= std::pow(2, itest - 3);
    coefs = { 0.0, 0.0, 0.0, 1.0 };  // cubic function
    TOL_L8 = 0.8 / nx / nx;
    TOL_L2 = 0.4 / nx / nx;
  }

  // Create simple meshes
  std::shared_ptr<Wonton::Simple_Mesh> source_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, mx, mx);
  std::shared_ptr<Wonton::Simple_Mesh> target_mesh =
    std::make_shared<Wonton::Simple_Mesh>(0.0, 0.0, 1.0, 1.0, nx, nx);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper(*source_mesh, true, true, true, Wonton::CoordSysType::CylindricalAxisymmetric);
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(*target_mesh, true, true, true, Wonton::CoordSysType::CylindricalAxisymmetric);

  const int ncells_source = sourceMeshWrapper.num_owned_cells();
  const int ncells_target = targetMeshWrapper.num_owned_cells();

  // Create a state object
  Wonton::Simple_State source_state(source_mesh);
  Wonton::Simple_State target_state(target_mesh);

  // Define a state vector with constant value and add it to the source state
  double vol, mass0(0.0);
  Wonton::Point<2> xc, xmin, xmax;
  std::vector<double> data(ncells_source);

  for (int c = 0; c < ncells_source; ++c) {
    sourceMeshWrapper.cell_centroid(c, &xc);
    vol = sourceMeshWrapper.cell_volume(c);

    data[c] = coefs[0] 
            + coefs[1] * xc[0] + coefs[2] * xc[1] 
            + coefs[3] * xc[0] * xc[1] * xc[1];

    mass0 += data[c] * vol;
  }
  source_state.add("cellvars", Wonton::Entity_kind::CELL, &(data[0]));

  // Create state wrapper
  Wonton::Simple_State_Wrapper sourceStateWrapper(source_state);
  Wonton::Simple_State_Wrapper targetStateWrapper(target_state);

  // Gather the cell coordinates as Portage Points for source and target meshes
  // for intersection. The outer vector is the cells, the inner vector is the
  // points of the vertices of that cell.
  std::vector<std::vector<Wonton::Point<2>>> source_cell_coords(ncells_source);
  std::vector<std::vector<Wonton::Point<2>>> target_cell_coords(ncells_target);

  for (int c = 0; c < ncells_source; ++c)
    sourceMeshWrapper.cell_get_coordinates(c, &(source_cell_coords[c]));
  for (int c = 0; c < ncells_target; ++c)
    targetMeshWrapper.cell_get_coordinates(c, &(target_cell_coords[c]));

  // Interpolate from source to target mesh using the independent calculation
  // in simple_intersect_for_tests.h
  std::vector<double> outvals(ncells_target);
  std::vector<std::vector<Portage::Weights_t>> sources_and_weights(ncells_target);

  // Loop over target cells
  for (int c = 0; c < ncells_target; ++c) {
    std::vector<int> xcells;
    std::vector<std::vector<double>> xwts;

    // Compute the moments
    // xcells is the source cell indices that intersect
    // xwts is the moments vector for each cell that intersects
    BOX_INTERSECT::intersection_moments<2,Wonton::CylindricalAxisymmetricCoordinates>(
        target_cell_coords[c],
        source_cell_coords,
        &xcells, &xwts);

    // Pack the results into a vector of true Portage::Weights_t
    int const num_intersect_cells = xcells.size();
    std::vector<Portage::Weights_t> wtsvec(num_intersect_cells);
    for (int i = 0; i < num_intersect_cells; ++i) {
      wtsvec[i].entityID = xcells[i];
      wtsvec[i].weights = xwts[i];
    }

    // Put the weights in final form
    sources_and_weights[c] = wtsvec;
  }

  // compute gradient field to pass to the interpolator
  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Simple_Mesh_Wrapper,
                      Wonton::Simple_State_Wrapper> driver(
     sourceMeshWrapper, sourceStateWrapper,
     targetMeshWrapper, targetStateWrapper);

  auto gradients = driver.compute_source_gradient("cellvars");

  // create interpolation object
  Portage::NumericTolerances_t num_tols = Portage::DEFAULT_NUMERIC_TOLERANCES<2>;
 
  Portage::Interpolate_2ndOrder<2, Wonton::Entity_kind::CELL,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_Mesh_Wrapper,
                                Wonton::Simple_State_Wrapper,
                                Wonton::Simple_State_Wrapper,
                                double,
                                Portage::DummyInterfaceReconstructor,
                                void,
                                void,
                                Wonton::CylindricalAxisymmetricCoordinates>
      interpolator(sourceMeshWrapper, targetMeshWrapper, sourceStateWrapper,
                   num_tols);

  interpolator.set_interpolation_variable("cellvars", &gradients);

  Wonton::transform(targetMeshWrapper.begin(Wonton::Entity_kind::CELL),
                     targetMeshWrapper.end(Wonton::Entity_kind::CELL),
                     sources_and_weights.begin(),
                     outvals.begin(), interpolator);

  // error analysis
  double stdval, mass1(0.0), errl2(0.0);
  for (int c = 0; c < ncells_target; ++c) {
    targetMeshWrapper.cell_centroid(c, &xc);
    vol = targetMeshWrapper.cell_volume(c);

    stdval = coefs[0]
           + coefs[1] * xc[0] + coefs[2] * xc[1] 
           + coefs[3] * xc[0] * xc[1] * xc[1];
    ASSERT_NEAR(stdval, outvals[c], TOL_L8);

    mass1 += outvals[c] * vol;
    errl2 += std::pow(stdval - outvals[c], 2);
  }
  errl2 = sqrt(errl2 / ncells_target);

  std::cout << "masses: " << mass0 << " " << mass1
            << " tols: " << TOL_L2 << " " << TOL_L8 << std::endl;
  std::cout << "mass error:     " << mass0 - mass1 << std::endl;
  std::cout << "solution error: " << errl2 << std::endl;

  ASSERT_NEAR(mass0, mass1, 1.0e-14);
  ASSERT_NEAR(0.0, errl2, TOL_L2);

}


INSTANTIATE_TEST_SUITE_P(
  Order2TestAll,
  Order2Test,
  ::testing::Values(1, 2, 3, 4, 5));

