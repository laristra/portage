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

  // This is bypassing the redistribution step in parallel runs So the
  // two meshes have to be nearly identical and the partitioner should
  // give a nearly identical partitioning (best to use BLOCK
  // partitioner in Jali)

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  
  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 4, 2);
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 5, 3);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  int nsrccells_all = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  //-------------------------------------------------------------------
  // Now add density field to the mesh
  //-------------------------------------------------------------------

  double source_mass = 0.0, source_mass_global = 0.0;
  std::vector<double> srcdensity(nsrccells_all);
  for (int c = 0; c < nsrccells_all; c++) {
    Wonton::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    srcdensity[c] = cen[0] + 2*cen[1];

    if (sourceMeshWrapper.cell_get_type(c) == Wonton::PARALLEL_OWNED) {
      double cvol = sourceMeshWrapper.cell_volume(c);
      double cmass = cvol*srcdensity[c];
      source_mass += cmass;
    }
  }

#ifdef WONTON_ENABLE_MPI
  MPI_Allreduce(&source_mass, &source_mass_global, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
#else
  source_mass_global = source_mass;
#endif

  
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "density", srcdensity.data());

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "density", 0.0);

  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, default mismatch fixup options

  Wonton::MPIExecutor_type executor(MPI_COMM_WORLD);
  
  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper, &executor);

  auto candidates = d.search<Portage::SearchKDTree>();
  auto srcwts = d.intersect_meshes<Portage::IntersectRnD>(candidates);

  auto gradients = d.compute_source_gradient("density");

  double exact_gradient[2] = {1.0, 2.0};
  for (int c = 0; c < nsrccells_all; c++) {
    Wonton::Vector<2> const& nabla = gradients[c];
    for (int dim = 0; dim < 2; dim++) {
      ASSERT_NEAR(exact_gradient[dim], nabla[dim], 1.0e-12);
    }
  }

  
  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
      "density", "density", srcwts, &gradients
  );


  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------


  // Finally check that we got the right target density values
  double *targetdensity;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "density",
                                   &targetdensity);

  int ntrgcells_owned =
      targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                     Wonton::Entity_type::PARALLEL_OWNED);

  double target_mass = 0.0, target_mass_global = 0.0;
  for (int c = 0; c < ntrgcells_owned; c++) {
    double cvol = targetMeshWrapper.cell_volume(c);
    double cmass = targetdensity[c]*cvol;
    target_mass += cmass;
  }

#ifdef WONTON_ENABLE_MPI
  MPI_Allreduce(&target_mass, &target_mass_global, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
#else
  target_mass_global = target_mass;
#endif

  ASSERT_NEAR(source_mass_global, target_mass_global, 1.0e-10);

  for (int c = 0; c < ntrgcells_owned; c++) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double trgdensity = cen[0] + 2*cen[1];
    ASSERT_NEAR(targetdensity[c], trgdensity, 1.0e-10);
  }

}  // CellDriver_2D_2ndOrder





TEST(CellDriver, 3D_2ndOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  
  sourceMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 8, 8, 8);
  targetMesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 10, 10, 10);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);


  int nsrccells_all = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                     Wonton::Entity_type::ALL);

  double source_mass = 0.0, source_mass_global = 0.0;
  std::vector<double> srcdensity(nsrccells_all);
  for (int c = 0; c < nsrccells_all; c++) {
    Wonton::Point<3> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    srcdensity[c] = cen[0] + 2*cen[1] + 3*cen[2];

    if (sourceMeshWrapper.cell_get_type(c) == Wonton::PARALLEL_OWNED) {
      double cvol = sourceMeshWrapper.cell_volume(c);
      double cmass = cvol*srcdensity[c];
      source_mass += cmass;
    }
  }

#ifdef WONTON_ENABLE_MPI
  MPI_Allreduce(&source_mass, &source_mass_global, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
#else
  source_mass_global = source_mass;
#endif
  
  
  //-------------------------------------------------------------------
  // Now add density field to the mesh
  //-------------------------------------------------------------------

  sourceStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "density", srcdensity.data());
  

  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "DENS", 0.0);

  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, default mismatch fixup options
  
  Wonton::MPIExecutor_type executor(MPI_COMM_WORLD);

  Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper, &executor);

  auto candidates = d.search<Portage::SearchKDTree>();
  auto srcwts = d.intersect_meshes<Portage::IntersectRnD>(candidates);

  auto gradients = d.compute_source_gradient("density");

  // Should get a gradient of 1.0, 2.0, 3.0
  double exact_gradient[3] = {1.0, 2.0, 3.0};
  for (int c = 0; c < nsrccells_all; c++) {
    Wonton::Vector<3> const& nabla = gradients[c];
    for (int dim = 0; dim < 3; dim++) {
      ASSERT_NEAR(exact_gradient[dim], nabla[dim], 1.0e-12);
    }
  }

  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "density", "DENS", srcwts, &gradients
  );
  
  // Finally check that we got the right target density values
  double *targetdensity;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "DENS",
                                   &targetdensity);
  
  int ntrgcells_owned =
      targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                     Wonton::Entity_type::PARALLEL_OWNED);
  double target_mass = 0.0, target_mass_global = 0.0;
  for (int c = 0; c < ntrgcells_owned; c++) {
    double cvol = targetMeshWrapper.cell_volume(c);
    double cmass = targetdensity[c]*cvol;
    target_mass += cmass;
  }

#ifdef WONTON_ENABLE_MPI
  MPI_Allreduce(&target_mass, &target_mass_global, 1, MPI_DOUBLE,
                MPI_SUM, MPI_COMM_WORLD);
#else
  target_mass_global = target_mass;
#endif

  ASSERT_NEAR(source_mass_global, target_mass_global, 1.0e-10);

  for (int c = 0; c < ntrgcells_owned; c++) {
    Wonton::Point<3> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double trgdensity = cen[0] + 2*cen[1] + 3*cen[2];
    ASSERT_NEAR(targetdensity[c], trgdensity, 1.0e-10);
  }

}  // CellDriver_3D_2ndOrder


