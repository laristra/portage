/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifdef HAVE_TANGRAM

#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
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

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 7, 6);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  //-------------------------------------------------------------------
  // Now add temperature field to the mesh
  //-------------------------------------------------------------------

  std::vector<double> srctemp(nsrccells);
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    srctemp[c] = cen[0] + 2*cen[1];
  }
  
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "temperature", srctemp.data());

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "temperature", 0.0);

  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, no mismatch fixup

  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  auto candidates = d.search<Portage::SearchKDTree>();
  auto srcwts = d.intersect_meshes<Portage::IntersectR2D>(candidates);

  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>("temperature",
                                                                "temperature",
                                                                srcwts,
                                                                dblmin, dblmax);


  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------


  // Finally check that we got the right target temperature values
  double *targettemp;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "temperature",
                                   &targettemp);

  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < ntrgcells; c++) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double trgtemp = cen[0] + 2*cen[1];
    ASSERT_NEAR(targettemp[c], trgtemp, 1.0e-10);
  }

}  // CellDriver_2D_2ndOrder





TEST(CellDriver, 3D_2ndOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);

  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);


  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  std::vector<double> srctemp(nsrccells);
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Point<3> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    srctemp[c] = cen[0] + 2*cen[1] + 3*cen[2];
  }
  
  
  //-------------------------------------------------------------------
  // Now add temperature field to the mesh
  //-------------------------------------------------------------------

  sourceStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "temperature", srctemp.data());
  

  //-------------------------------------------------------------------
  // Field(s) we have to remap
  //-------------------------------------------------------------------

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "TEMP", 0.0);

  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, no mismatch fixup
  
  Wonton::SerialExecutor_type executor;

  Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper, &executor);

  auto candidates = d.search<Portage::SearchKDTree>();

  auto srcwts = d.intersect_meshes<Portage::IntersectR3D>(candidates);

  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  d.interpolate_mesh_var<double,
                         Portage::Interpolate_2ndOrder>("temperature", "TEMP",
                                                        srcwts, dblmin, dblmax);
  
  // Finally check that we got the right target temperature values
  double *targettemp;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "TEMP",
                                   &targettemp);
  
  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  for (int c = 0; c < ntrgcells; c++) {
    Wonton::Point<3> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double trgtemp = cen[0] + 2*cen[1] + 3*cen[2];
    ASSERT_NEAR(targettemp[c], trgtemp, 1.0e-10);
  }
}  // CellDriver_3D_2ndOrder


#endif  // ifdef HAVE_TANGRAM
