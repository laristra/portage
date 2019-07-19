//
// Created by Hoby Rakotarivelo on 2019-07-17.
//

#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#ifdef PORTAGE_ENABLE_MPI
  #include "mpi.h"
#endif

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "portage/driver/coredriver.h"

TEST(PartDriver, SingleMat_NoMismatch) {

  // useful constants
  double const dblmax  = std::numeric_limits<double>::max();
  double const dblmin  = -dblmax;
  double const epsilon = 1.E-10;

  // Source and target meshes and states
  std::shared_ptr<Jali::Mesh>  sourceMesh;
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::Mesh>  targetMesh;
  std::shared_ptr<Jali::State> targetState;

  sourceMesh  = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);
  targetMesh  = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 10, 10);
  sourceState = Jali::State::create(sourceMesh);
  targetState = Jali::State::create(targetMesh);

  Wonton::Jali_Mesh_Wrapper  sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_Mesh_Wrapper  targetMeshWrapper(*targetMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Wonton::Jali_State_Wrapper targetStateWrapper(*targetState);

  int const nb_source_cells = sourceMeshWrapper.num_entities(
    Wonton::Entity_kind::CELL, Wonton::Entity_type::ALL
  );

  //-------------------------------------------------------------------
  // Now add temperature field to the mesh
  //-------------------------------------------------------------------
  auto compute_temperature = [&](int c, Wonton::Jali_Mesh_Wrapper const& mesh) {
    double const x_gap = 0.4;
    double const t_min = 30.;
    double const t_max = 100.;

    Wonton::Point<2> centroid;
    mesh.cell_centroid(c, &centroid);
    return (centroid[0] < x_gap ? t_min : t_max);
  };

  double source_temperature[nb_source_cells];
  for (int c = 0; c < nb_source_cells; c++) {
    source_temperature[c] = compute_temperature(c, sourceMeshWrapper);
  }

  sourceStateWrapper.mesh_add_data(
    Wonton::Entity_kind::CELL, "temperature", source_temperature
  );

  targetStateWrapper.mesh_add_data<double>(
    Wonton::Entity_kind::CELL, "temperature", 0.0
  );

  // assuming that entities ordering are column-wise,
  // let's pick a couple of aligned subsets of entites
  // from source and target meshes (5x5, 10x10).
  std::vector<int> source_entities {
    6, 11,
    7, 12,
    8, 13
  };

  std::vector<int> target_entities {
    22, 32, 42,	52,
    23, 33, 43,	53,
    24, 34, 44,	54,
    25, 35, 45,	55,
    26, 36, 46,	56,
    27, 37, 47,	57
  };

  // Do the basic remap algorithm (search, intersect, interpolate) -
  // no redistribution, no mismatch fixup
  Portage::CoreDriver<
    2, Wonton::Entity_kind::CELL,
    Wonton::Jali_Mesh_Wrapper,
    Wonton::Jali_State_Wrapper
  > d(sourceMeshWrapper, sourceStateWrapper,
      targetMeshWrapper, targetStateWrapper);

  auto candidates = d.search<Portage::SearchKDTree>();
  auto source_weights = d.intersect_meshes<Portage::IntersectR2D>(candidates);

  d.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
    "temperature", "temperature", source_weights,
    source_entities, target_entities, dblmin, dblmax
  );

  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------
  // Finally check that we got the right target temperature values
  double* remapped;
  targetStateWrapper.mesh_get_data(
    Wonton::Entity_kind::CELL, "temperature", &remapped
  );

  for (auto&& c : target_entities) {
    double const expected = compute_temperature(c, targetMeshWrapper);
#ifdef DEBUG
    std::printf("remap[%d]: %.3f, exact[%d]: %.3f\n",
      c, remapped[c], c, expected
    );
#endif
    ASSERT_NEAR(remapped[c], expected, epsilon);
  }
}