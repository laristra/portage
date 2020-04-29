/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

// this should be included prior to the use of Portage macros
#include "portage/support/portage.h"

#ifdef HAVE_TANGRAM
// system
#include "gtest/gtest.h"
#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

// wonton
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// jali
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

// portage
#include "portage/driver/coredriver.h"
#include "portage/search/search_swept_face.h"
#include "portage/intersect/intersect_swept_face.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"

// tangram
#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/reconstruct/MOF.h"

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
  double dst_tol = Portage::DEFAULT_NUMERIC_TOLERANCES<2>.min_absolute_distance;

  //-------------------------------------------------------------------
  // Shift internal nodes of the targetmesh
  //-------------------------------------------------------------------

  int ntrgnodes = targetMesh->num_entities(Jali::Entity_kind::NODE,
                                          Jali::Entity_type::ALL);

  // node coordinate modification has to be done in Jali (the mesh wrappers
  // are query only)
  for (int i = 0; i < ntrgnodes; i++) {
    std::array<double, 2> pnt {};
    targetMesh->node_get_coordinates(i, &pnt);

    // Move only the internal nodes because we don't want to mess with
    // boundary conditions
    
    if (fabs(pnt[0]-minxyz) <= dst_tol || fabs(pnt[0]-maxxyz) <= dst_tol ||
        fabs(pnt[1]-minxyz) <= dst_tol || fabs(pnt[1]-maxxyz) <= dst_tol)
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

  // IMPORTANT: we need to set the right tangram options even in single material
  // since we explicitly use MatPoly_clipper in face_group_moments.
  // It will trigger a compilation error otherwise - even in single material setting.
  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                      Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>

      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  auto candidates = d.search<Portage::SearchSweptFace>();
  auto srcwts = d.intersect_meshes<Portage::IntersectSweptFace2D>(candidates);

  // Remap temperature and check that we got the right target
  // temperatures - both cell wise values and integral should match
  d.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
    "temperature", "temperature", srcwts);

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
  d.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>("density",
                                                                "density",
                                                                srcwts); 

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
    std::array<double, 2> pnt {};
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

  // IMPORTANT: we need to set the right tangram options even in single material
  // since we explicitly use MatPoly_clipper in face_group_moments.
  // It will trigger a compilation error otherwise - even in single material setting.
  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                      Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>

      d(sourceMeshWrapper, sourceStateWrapper,
        targetMeshWrapper, targetStateWrapper);

  auto candidates = d.search<Portage::SearchSweptFace>();
  auto srcwts = d.intersect_meshes<Portage::IntersectSweptFace2D>(candidates);

  //bool has_mismatch = d.check_mesh_mismatch(srcwts);

  //double dblmin = -std::numeric_limits<double>::max();
  //double dblmax =  std::numeric_limits<double>::max();

  
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

TEST(SweptFaceRemap, 3D_2ndOrder) {

  // useful constants
  double const p_min = 0.0;
  double const p_max = 1.0;
  using Wonton::Entity_kind::CELL;

  // IMPORTANT: we need to set the right tangram options even in single material
  // since we explicitly use MatPoly_clipper in face_group_moments.
  // It will trigger a compilation error otherwise - even in single material setting.
  using Remapper = Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
                                        Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                        Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                        Tangram::MOF, Tangram::SplitR3D, Tangram::ClipR3D>;

  // create meshes and related states
  auto source_mesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0,0.0,0.0,1.0,1.0,1.0,5,5,5);
  auto target_mesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0,0.0,0.0,1.0,1.0,1.0,5,5,5);

  int const nb_target_nodes = target_mesh->num_nodes<Jali::Entity_type::ALL>();
  int const nb_source_cells = source_mesh->num_cells<Jali::Entity_type::ALL>();
  int const nb_target_cells = target_mesh->num_cells<Jali::Entity_type::ALL>();

  // check if a point is on boundary
  auto is_boundary = [&](auto const& point) -> bool {
    for (int d = 0; d < 3; ++d) {
      if (std::abs(point[d] - p_min) < 1.E-16
       or std::abs(p_max - point[d]) < 1.E-16)
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

  auto candidates = remapper.search<Portage::SearchSweptFace>();
  auto gradients  = remapper.compute_source_gradient("temperature");
  auto weights    = remapper.intersect_meshes<Portage::IntersectSweptFace3D>(candidates);

  remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "temperature", "temperature", weights, &gradients
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

TEST(SweptFaceRemap, MassConservationConstantField2D) {

  using Wonton::Entity_kind::CELL;
  using Wonton::Entity_type::ALL;

  // IMPORTANT: we need to set the right tangram options even in single material
  // since we explicitly use MatPoly_clipper in face_group_moments.
  // It will trigger a compilation error otherwise - even in single material setting.
  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Tangram::MOF, Tangram::SplitR2D, Tangram::ClipR2D>;

  MPI_Comm comm = MPI_COMM_WORLD;

  auto source_mesh  = Jali::MeshFactory(comm)(0.0,0.0,1.0,1.0,5,5);
  auto target_mesh  = Jali::MeshFactory(comm)(0.5,0.5,1.5,1.5,5,5);
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  int const nb_source_cells = source_mesh_wrapper.num_entities(CELL, ALL);
  int const nb_target_cells = target_mesh_wrapper.num_entities(CELL, ALL);

  // now add a linear temperature field on the source mesh
  double source_density[nb_source_cells];
  double* target_density = nullptr;

  double source_mass = 0.0;
  for (int c = 0; c < nb_source_cells; c++) {
    source_density[c] = 42;
    double const area = source_mesh_wrapper.cell_volume(c);
    source_mass += source_density[c] * area;
  }

  source_state_wrapper.mesh_add_data<double>(CELL, "density", source_density);
  target_state_wrapper.mesh_add_data<double>(CELL, "density", 0.0);
  target_state_wrapper.mesh_get_data<double>(CELL, "density", &target_density);

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchSweptFace>();
  auto gradients  = remapper.compute_source_gradient("density");
  auto weights    = remapper.intersect_meshes<Portage::IntersectSweptFace2D>(candidates);

  remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "density", "density", weights, &gradients
  );

  double target_mass = 0.0;

  for (int c = 0; c < nb_target_cells; c++) {
    double const area = target_mesh_wrapper.cell_volume(c);
    target_mass += target_density[c] * area;
  }

  ASSERT_DOUBLE_EQ(source_mass, target_mass);
}

TEST(SweptFaceRemap, MassConservationConstantField3D) {

  using Wonton::Entity_kind::CELL;
  using Wonton::Entity_type::ALL;

  // IMPORTANT: we need to set the right tangram options even in single material
  // since we explicitly use MatPoly_clipper in face_group_moments.
  // It will trigger a compilation error otherwise - even in single material setting.
  using Remapper = Portage::CoreDriver<3, Wonton::Entity_kind::CELL,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper,
                                       Tangram::MOF, Tangram::SplitR3D, Tangram::ClipR3D>;

  MPI_Comm comm = MPI_COMM_WORLD;

  auto source_mesh  = Jali::MeshFactory(comm)(0.0,0.0,0.0,1.0,1.0,1.0,5,5,5);
  auto target_mesh  = Jali::MeshFactory(comm)(0.5,0.5,0.5,1.5,1.5,1.5,5,5,5);
  auto source_state = Jali::State::create(source_mesh);
  auto target_state = Jali::State::create(target_mesh);

  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Wonton::Jali_State_Wrapper source_state_wrapper(*source_state);
  Wonton::Jali_State_Wrapper target_state_wrapper(*target_state);

  int const nb_source_cells = source_mesh_wrapper.num_entities(CELL, ALL);
  int const nb_target_cells = target_mesh_wrapper.num_entities(CELL, ALL);

  // now add a linear temperature field on the source mesh
  double source_density[nb_source_cells];
  double* target_density = nullptr;

  double source_mass = 0.0;
  for (int c = 0; c < nb_source_cells; c++) {
    source_density[c] = 42;
    double const volume = source_mesh_wrapper.cell_volume(c);
    source_mass += source_density[c] * volume;
  }

  source_state_wrapper.mesh_add_data<double>(CELL, "density", source_density);
  target_state_wrapper.mesh_add_data<double>(CELL, "density", 0.0);
  target_state_wrapper.mesh_get_data<double>(CELL, "density", &target_density);

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  auto candidates = remapper.search<Portage::SearchSweptFace>();
  auto gradients  = remapper.compute_source_gradient("density");
  auto weights    = remapper.intersect_meshes<Portage::IntersectSweptFace3D>(candidates);

  remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
    "density", "density", weights, &gradients
  );

  double target_mass = 0.0;

  for (int c = 0; c < nb_target_cells; c++) {
    auto const volume = target_mesh_wrapper.cell_volume(c);
    target_mass += target_density[c] * volume;
  }

  ASSERT_DOUBLE_EQ(source_mass, target_mass);
}
#endif