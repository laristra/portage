/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

// this should be included prior to the use of Portage macros
#include "portage/support/portage.h"

#include <iostream>
#include <memory>

#include "gtest/gtest.h"

#include "wonton/support/wonton.h"

#ifdef WONTON_ENABLE_MPI
#include "mpi.h"
#endif

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"

#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "portage/driver/coredriver.h"

double TOL = 1e-6;


// Tests for single material remap on cells with 2nd Order Accurate Remap

TEST(CellDriver, 2D_2ndOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  MPI_Comm comm = MPI_COMM_WORLD;

  // This is bypassing the redistribution step in parallel runs So the
  // two meshes have to be nearly identical and the partitioner should
  // give a nearly identical partitioning (best to use BLOCK
  // partitioner in Jali)
  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.partitioner(Jali::Partitioner_type::BLOCK);
  
  sourceMesh = mesh_factory(0.0, 0.0, 1.0, 1.0, 4, 2);
  targetMesh = mesh_factory(0.0, 0.0, 1.0, 1.0, 5, 3);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  //-------------------------------------------------------------------
  // Now add density field to the mesh
  //-------------------------------------------------------------------

  double source_mass_local = 0.0, source_mass_global = 0.0;
  std::vector<double> source_density(nsrccells);
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    source_density[c] = cen[0] + 2 * cen[1];

    if (sourceMeshWrapper.cell_get_type(c) == Wonton::PARALLEL_OWNED) {
      source_mass_local += (sourceMeshWrapper.cell_volume(c) * source_density[c]);
    }
  }

  MPI_Allreduce(&source_mass_local, &source_mass_global, 1, MPI_DOUBLE, MPI_SUM, comm);
  
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL, "density", source_density.data());
  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,"density", 0.0);

  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, default mismatch fixup options

  Wonton::MPIExecutor_type executor(comm);
  
  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper, &executor);

  auto candidates = d.search<Portage::SearchKDTree>();
  auto srcwts = d.intersect_meshes<Portage::IntersectRnD>(candidates);

  auto gradients = d.compute_source_gradient("density");

  // field gradient on ghost cells should be already updated at this point
  double exact_gradient[2] = {1.0, 2.0};
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Vector<2> const& approx_gradient = gradients[c];
    for (int dim = 0; dim < 2; dim++) {
      ASSERT_NEAR(exact_gradient[dim], approx_gradient[dim], 1.0e-12);
    }
  }

  
  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
      "density", "density", srcwts, &gradients
  );


  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------


  // Finally check that we got the right target density values
  double* target_density;
  targetStateWrapper.mesh_get_data(Wonton::CELL, "density", &target_density);

  int const num_owned_target_cells = targetMeshWrapper.num_owned_cells();

  double target_mass_local = 0.0, target_mass_global = 0.0;
  for (int c = 0; c < num_owned_target_cells; c++) {
    double cvol = targetMeshWrapper.cell_volume(c);
    double cmass = target_density[c] * cvol;
    target_mass_local += cmass;
  }

  MPI_Allreduce(&target_mass_local, &target_mass_global, 1, MPI_DOUBLE, MPI_SUM, comm);

  ASSERT_NEAR(source_mass_global, target_mass_global, 1.0e-10);

  for (int c = 0; c < num_owned_target_cells; c++) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double const expected_density = cen[0] + 2 * cen[1];
    ASSERT_NEAR(target_density[c], expected_density, 1.0e-10);
  }

}  // CellDriver_2D_2ndOrder





TEST(CellDriver, 3D_2ndOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  MPI_Comm comm = MPI_COMM_WORLD;

  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.partitioner(Jali::Partitioner_type::BLOCK);
  
  sourceMesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 8, 8);
  targetMesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);


  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  double source_mass_local = 0.0, source_mass_global = 0.0;
  std::vector<double> srcdensity(nsrccells);
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Point<3> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    srcdensity[c] = cen[0] + 2*cen[1] + 3*cen[2];

    if (sourceMeshWrapper.cell_get_type(c) == Wonton::PARALLEL_OWNED) {
      source_mass_local += (sourceMeshWrapper.cell_volume(c) * srcdensity[c]);
    }
  }

  MPI_Allreduce(&source_mass_local, &source_mass_global, 1, MPI_DOUBLE, MPI_SUM, comm);
  
  //-------------------------------------------------------------------
  // Now add density field to the mesh
  //-------------------------------------------------------------------

  sourceStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL, "density", srcdensity.data());
  

  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL, "density", 0.0);

  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, default mismatch fixup options
  
  Wonton::MPIExecutor_type executor(comm);

  Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper, &executor);

  auto candidates = d.search<Portage::SearchKDTree>();
  auto srcwts = d.intersect_meshes<Portage::IntersectRnD>(candidates);

  auto gradients = d.compute_source_gradient("density");

  // field gradient on ghost cells should be already updated at this point
  // we should get a gradient of 1.0, 2.0, 3.0
  double exact_gradient[3] = {1.0, 2.0, 3.0};
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Vector<3> const& approx_gradient = gradients[c];
    for (int dim = 0; dim < 3; dim++) {
      ASSERT_NEAR(exact_gradient[dim], approx_gradient[dim], 1.0e-12);
    }
  }

  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "density", "density", srcwts, &gradients
  );
  
  // Finally check that we got the right target density values
  double* target_density;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "density", &target_density);
  
  int num_owned_target_cells = targetMeshWrapper.num_owned_cells();

  double target_mass_local = 0.0, target_mass_global = 0.0;
  for (int c = 0; c < num_owned_target_cells; c++) {
    target_mass_local += (targetMeshWrapper.cell_volume(c) * target_density[c]);
  }

  MPI_Allreduce(&target_mass_local, &target_mass_global, 1, MPI_DOUBLE, MPI_SUM, comm);

  ASSERT_NEAR(source_mass_global, target_mass_global, 1.0e-10);

  for (int c = 0; c < num_owned_target_cells; c++) {
    Wonton::Point<3> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double expected_density = cen[0] + 2 * cen[1] + 3 * cen[2];
    ASSERT_NEAR(target_density[c], expected_density, 1.0e-10);
  }

}  // CellDriver_3D_2ndOrder


