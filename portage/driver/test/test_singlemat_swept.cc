/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/search/search_swept_face.h"
#include "portage/intersect/intersect_swept_face.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "portage/driver/coredriver.h"
#include "portage/support/portage.h"

double TOL = 1e-6;


// Integrated tests for single material swept-face remap


// 1st order remap of constant field preserves both the field and the
// integral value but only integral value of general field

TEST(SweptFaceRemap, 2D_1stOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  // For sweptface remap, target and source mesh must be identical with
  // identical numbering of all entities

  double minxyz = 0.0, maxxyz = 1.0; 
  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);

  //-------------------------------------------------------------------
  // Shift internal nodes of the targetmesh
  //-------------------------------------------------------------------

  int ntrgnodes = targetMesh->num_entities(Jali::Entity_kind::NODE,
                                          Jali::Entity_type::ALL);

  // node coordinate modification has to be done in Jali (the mesh wrappers
  // are query only)
  for (int i = 0; i < ntrgnodes; i++) {
    std::array<double, 2> pnt;
    targetMesh->node_get_coordinates(i, &pnt);

    // Move only the internal nodes because we don't want to mess with
    // boundary conditions
    
    if (fabs(pnt[0]-minxyz) < 1.0e-16 || fabs(pnt[0]-maxxyz) < 1.0e-16 ||
        fabs(pnt[1]-minxyz) < 1.0e-16 || fabs(pnt[1]-maxxyz) < 1.0e-16)
      continue;  // boundary point - don't move

    pnt[0] += 0.1*sin(2*M_PI*pnt[0]);
    pnt[1] += 0.1*sin(2*M_PI*pnt[1]);
    targetMesh->node_set_coordinates(i, pnt.data());
  }

  // Create state managers
  
  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  // Create the Portage mesh/state wrappers - MUST BE AFTER WE SHIFT
  // THE NODES OTHERWISE THE CENTROIDS GIVEN BY THE WRAPPER FOR THE
  // TARGET MESH WILL BE WRONG (BECAUSE THEY ARE CACHED AND NOT UPDATED)

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);


  
  //-------------------------------------------------------------------
  // Now add a constant temperature field and a general density field
  // to the mesh
  // -------------------------------------------------------------------

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  double source_integral_temp = 0.0, source_integral_density = 0.0;
  std::vector<double> srctemp(nsrccells), srcdensity(nsrccells);
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    double cellvol = sourceMeshWrapper.cell_volume(c);
    srctemp[c] = 42;
    source_integral_temp += srctemp[c]*cellvol;
    srcdensity[c] = cen[0]*cen[0] + cen[1]*cen[1]*cen[1];
    source_integral_density += srcdensity[c]*cellvol;
  }
  
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "temperature", srctemp.data());
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "density", srcdensity.data());
  


  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "temperature", 0.0);
  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "density", 0.0);

  
  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, default mismatch fixup options

  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  Portage::NumericTolerances_t default_num_tols;
  default_num_tols.use_default();
  d.set_num_tols(default_num_tols);

  auto candidates = d.search<Portage::SearchSweptFace>();
  auto srcwts = d.intersect_meshes<Portage::IntersectSweptFace2D>(candidates);

  // Remap temperature and check that we got the right target
  // temperatures - both cell wise values and integral should match
  auto temp_gradients = d.compute_source_gradient("temperature");

  d.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
    "temperature", "temperature", srcwts, &temp_gradients);

  double *targettemp;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "temperature",
                                   &targettemp);

  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  double target_integral_temp = 0.0;
  for (int c = 0; c < ntrgcells; c++) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double trgtemp = 42;
    ASSERT_NEAR(targettemp[c], trgtemp, 1.0e-12);
    target_integral_temp += targettemp[c]*targetMeshWrapper.cell_volume(c);
  }

  ASSERT_NEAR(source_integral_temp, target_integral_temp, 1.0e-12);

  // Remap density and check that we got the right target densities -
  // cell wise values won't necessarily match as it is a general
  // function but integral values should match
  
  auto density_gradients = d.compute_source_gradient("density");

  d.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>("density",
                                                                "density",
                                                                srcwts,
                                                                &density_gradients);

  double *targetdensity;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "density",
                                   &targetdensity);

  double target_integral_density = 0.0;
  for (int c = 0; c < ntrgcells; c++)
    target_integral_density += targetdensity[c]*targetMeshWrapper.cell_volume(c);

  ASSERT_NEAR(source_integral_density, target_integral_density, 1.0e-12);

}  // SweptFace_2D_1stOrder



// 2nd order remap of linear field preserves both the field and the
// integral value but only integral value of general field

TEST(SweptFaceRemap, 2D_2ndOrder) {
  // Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  // Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;

  // For sweptface remap, target and source mesh must be identical with
  // identical numbering of all entities

  double minxyz = 0.0, maxxyz = 1.0; 
  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);
  targetMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);

  //-------------------------------------------------------------------
  // Shift internal nodes of the targetmesh
  //-------------------------------------------------------------------

  int ntrgnodes = targetMesh->num_entities(Jali::Entity_kind::NODE,
                                          Jali::Entity_type::ALL);

  // node coordinate modification has to be done in Jali (the mesh wrappers
  // are query only)
  for (int i = 0; i < ntrgnodes; i++) {
    std::array<double, 2> pnt;
    targetMesh->node_get_coordinates(i, &pnt);

    // Move only the internal nodes because we don't want to mess with
    // boundary conditions
    
    if (fabs(pnt[0]-minxyz) < 1.0e-16 || fabs(pnt[0]-maxxyz) < 1.0e-16 ||
        fabs(pnt[1]-minxyz) < 1.0e-16 || fabs(pnt[1]-maxxyz) < 1.0e-16)
      continue;  // boundary point - don't move

    pnt[0] += 0.1*sin(2*M_PI*pnt[0]);
    pnt[1] += 0.1*sin(2*M_PI*pnt[1]);
    targetMesh->node_set_coordinates(i, pnt.data());
  }

  // Create state managers
  
  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  // Create the Portage mesh/state wrappers - MUST BE AFTER WE SHIFT
  // THE NODES OTHERWISE THE CENTROIDS GIVEN BY THE WRAPPER FOR THE
  // TARGET MESH WILL BE WRONG (BECAUSE THEY ARE CACHED AND NOT UPDATED)

  Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);


  
  //-------------------------------------------------------------------
  // Now add a constant temperature field and a general density field
  // to the mesh
  // -------------------------------------------------------------------

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  double source_integral_temp = 0.0, source_integral_density = 0.0;
  std::vector<double> srctemp(nsrccells), srcdensity(nsrccells);
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    double cellvol = sourceMeshWrapper.cell_volume(c);
    srctemp[c] = cen[0] + 2*cen[1];
    source_integral_temp += srctemp[c]*cellvol;
    srcdensity[c] = cen[0]*cen[0] + cen[1]*cen[1]*cen[1];
    source_integral_density += srcdensity[c]*cellvol;
  }
  
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "temperature", srctemp.data());
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "density", srcdensity.data());
  


  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "temperature", 0.0);
  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "density", 0.0);

  
  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, default mismatch fixup options

  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  Portage::NumericTolerances_t default_num_tols;
  default_num_tols.use_default();
  d.set_num_tols(default_num_tols);

  auto candidates = d.search<Portage::SearchSweptFace>();
  auto srcwts = d.intersect_meshes<Portage::IntersectSweptFace2D>(candidates);

  bool has_mismatch = d.check_mesh_mismatch(srcwts);

  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  
  // Remap temperature and check that we got the right target
  // temperatures - both cell wise values and integral should match
  auto temp_gradients = d.compute_source_gradient("temperature");

  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "temperature", "temperature", srcwts, &temp_gradients);

  double *targettemp;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "temperature",
                                   &targettemp);

  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  double target_integral_temp = 0.0;
  for (int c = 0; c < ntrgcells; c++) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double trgtemp = cen[0] + 2*cen[1];
    ASSERT_NEAR(targettemp[c], trgtemp, 1.0e-12);
    target_integral_temp += targettemp[c]*targetMeshWrapper.cell_volume(c);
  }

  ASSERT_NEAR(source_integral_temp, target_integral_temp, 1.0e-12);



  // Remap density and check that we got the right target densities -
  // cell wise values won't necessarily match as it is a general
  // function but integral values should match
  auto density_gradients = d.compute_source_gradient("density");

  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>("density",
                                                                "density",
                                                                srcwts,
                                                                &density_gradients);

  double *targetdensity;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "density",
                                   &targetdensity);

  double target_integral_density = 0.0;
  for (int c = 0; c < ntrgcells; c++) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    target_integral_density += targetdensity[c]*targetMeshWrapper.cell_volume(c);
  }

  ASSERT_NEAR(source_integral_density, target_integral_density, 1.0e-12);

}  // SweptFace_2D_1stOrder
