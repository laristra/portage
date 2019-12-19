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
#include "portage/interpolate/interpolate_2nd_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "portage/driver/coredriver.h"
#include "portage/support/portage.h"

double TOL = 1e-6;


// Integrated test for single material swept-face remap

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
  // Now add a linear temperature field to the mesh
  //-------------------------------------------------------------------

  int nsrccells = sourceMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);

  double source_integral = 0.0;
  std::vector<double> srctemp(nsrccells);
  for (int c = 0; c < nsrccells; c++) {
    Wonton::Point<2> cen;
    sourceMeshWrapper.cell_centroid(c, &cen);
    srctemp[c] = cen[0] + 2*cen[1];
    source_integral += srctemp[c]*sourceMeshWrapper.cell_volume(c);
  }
  
  sourceStateWrapper.mesh_add_data(Wonton::Entity_kind::CELL,
                                   "temperature", srctemp.data());

  targetStateWrapper.mesh_add_data<double>(Wonton::Entity_kind::CELL,
                                           "temperature", 0.0);

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

  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  auto gradients = d.compute_source_gradient("temperature");

  d.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "temperature", "temperature", srcwts, &gradients
  );

  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------


  // Finally check that we got the right target temperature values
  double *targettemp;
  targetStateWrapper.mesh_get_data(Wonton::Entity_kind::CELL, "temperature",
                                   &targettemp);

  int ntrgcells = targetMeshWrapper.num_entities(Wonton::Entity_kind::CELL,
                                                 Wonton::Entity_type::ALL);
  double target_integral = 0.0;
  for (int c = 0; c < ntrgcells; c++) {
    Wonton::Point<2> cen;
    targetMeshWrapper.cell_centroid(c, &cen);
    double trgtemp = cen[0] + 2*cen[1];
    ASSERT_NEAR(targettemp[c], trgtemp, 1.0e-12);
    target_integral += targettemp[c]*targetMeshWrapper.cell_volume(c);
  }

  ASSERT_NEAR(source_integral, target_integral, 1.0e-12);

}  // CellDriver_2D_2ndOrder


TEST(SweptFaceRemap, 3D_2ndOrder) {

  // useful constants
  static constexpr double const min = 0.0;
  static constexpr double const max = 1.0;
  static constexpr auto const CELL = Wonton::Entity_kind::CELL;

  // useful shortcuts
  using Remapper = Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
                                       Wonton::Jali_Mesh_Wrapper,
                                       Wonton::Jali_State_Wrapper>;

  // create meshes and related states
  auto source_mesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0,0.0,0.0,1.0,1.0,1.0,5,5,5);
  auto target_mesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0,0.0,0.0,1.0,1.0,1.0,5,5,5);

  int const nb_target_nodes = target_mesh->num_nodes<Jali::Entity_type::ALL>();
  int const nb_source_cells = source_mesh->num_cells<Jali::Entity_type::ALL>();
  int const nb_target_cells = target_mesh->num_cells<Jali::Entity_type::ALL>();

  // check if a point is on boundary
  auto is_boundary = [&](auto const& point) -> bool {
    for (int d = 0; d < 3; ++d) {
      if (std::abs(point[d] - min) < 1.E-16
       or std::abs(max - point[d]) < 1.E-16)
        return false;
    }
    return true;
  };

  // move target points
  for (int i = 0; i < nb_target_nodes; i++) {
    std::array<double, 3> point = {0.,0.,0.};
    target_mesh->node_get_coordinates(i, &point);

    // move only internal nodes to avoid dealing with boundary conditions.
    if (not is_boundary(point)) {
      point[0] += 0.1 * sin(2 * M_PI * point[0]);
      point[1] += 0.1 * sin(2 * M_PI * point[1]);
      point[2] += 0.1 * sin(2 * M_PI * point[2]);
      target_mesh->node_set_coordinates(i, point.data());
    }
  }

  // create states then
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  // create related wrappers
  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  // now add a linear temperature field on the source mesh
  double field[nb_source_cells];
  double *remapped = nullptr;

  double source_integral = 0.0;
  for (int c = 0; c < nb_source_cells; c++) {
    auto const centroid = source_mesh->cell_centroid(c);
    auto const volume = source_mesh_wrapper.cell_volume(c);
    field[c] = centroid[0] + 2 * centroid[1] + 3 * centroid[2];
    source_integral += field[c] * volume;
  }

  source_state_wrapper.mesh_add_data<double>(CELL, "temperature", field);
  target_state_wrapper.mesh_add_data<double>(CELL, "temperature", 0.0);

  // now instantiate the driver
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  Portage::NumericTolerances_t tolerances;
  tolerances.use_default();
  remapper.set_num_tols(tolerances);

  auto candidates = remapper.search<Portage::SearchSweptFace>();
  auto gradients  = remapper.compute_source_gradient("temperature");
  auto weights    = remapper.intersect_meshes<Portage::IntersectSweptFace3D>(candidates);
  auto has_mismatch = remapper.check_mesh_mismatch(weights);

  remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "temperature", "temperature", weights, 0.0, 1.0,
    Portage::DEFAULT_PARTIAL_FIXUP_TYPE, Portage::DEFAULT_EMPTY_FIXUP_TYPE,
    Portage::DEFAULT_CONSERVATION_TOL,
    Portage::DEFAULT_MAX_FIXUP_ITER, nullptr, &gradients
  );

  // check remapped values on target mesh
  target_state_wrapper.mesh_get_data<double>(CELL, "temperature", &remapped);

  double target_integral = 0.0;
  // compare remapped values for each cell.
  for (int c = 0; c < nb_target_cells; c++) {
    Wonton::Point<3> centroid;
    target_mesh_wrapper.cell_centroid(c, &centroid);
    auto const volume = target_mesh_wrapper.cell_volume(c);
    auto const expected = centroid[0] + 2 * centroid[1] + 3 * centroid[2];
    ASSERT_NEAR(remapped[c], expected, 1.E-12);
    target_integral += remapped[c] * volume;
  }
  // check for integral quantities conservation.
  ASSERT_DOUBLE_EQ(source_integral, target_integral);
}
