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

TEST(PartDriverTest, Basic) {

  using Wonton::Entity_kind;
  using Wonton::Entity_type;
  using Wonton::Jali_Mesh_Wrapper;
  using Wonton::Jali_State_Wrapper;

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

  Jali_Mesh_Wrapper  sourceMeshWrapper(*sourceMesh);
  Jali_Mesh_Wrapper  targetMeshWrapper(*targetMesh);
  Jali_State_Wrapper sourceStateWrapper(*sourceState);
  Jali_State_Wrapper targetStateWrapper(*targetState);

  int const nb_source_cells =
    sourceMeshWrapper.num_entities(Entity_kind::CELL, Entity_type::ALL);

  //-------------------------------------------------------------------
  // Now add density field to the mesh
  //-------------------------------------------------------------------
  auto compute_density = [&](int c, Wonton::Jali_Mesh_Wrapper const& mesh) {
    double const x_gap = 0.4;
    double const t_min = 30.;
    double const t_max = 100.;

    Wonton::Point<2> centroid;
    mesh.cell_centroid(c, &centroid);
    return (centroid[0] < x_gap ? t_min : t_max);
  };

  double source_density[nb_source_cells];
  for (int c = 0; c < nb_source_cells; c++) {
    source_density[c] = compute_density(c, sourceMeshWrapper);
  }

  sourceStateWrapper.mesh_add_data(Entity_kind::CELL, "density", source_density);
  targetStateWrapper.mesh_add_data<double>(Entity_kind::CELL, "density", 0.);

  // assuming that entities ordering are column-wise,
  // let's pick a couple of aligned subsets of entites
  // from source and target meshes (5x5, 10x10).
  // FIXME update according to coordinates
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
    2, Entity_kind::CELL, Jali_Mesh_Wrapper, Jali_State_Wrapper
  > d(sourceMeshWrapper, sourceStateWrapper,
      targetMeshWrapper, targetStateWrapper);

  auto candidates = d.search<Portage::SearchKDTree>();
  auto source_weights = d.intersect_meshes<Portage::IntersectR2D>(candidates);

  d.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
    "density", "density", source_weights,
    source_entities, target_entities, dblmin, dblmax
  );

  //-------------------------------------------------------------------
  // CHECK REMAPPING RESULTS ON TARGET MESH SIDE
  //-------------------------------------------------------------------
  // Finally check that we got the right target density values
  double* remapped;
  targetStateWrapper.mesh_get_data(Entity_kind::CELL, "density", &remapped);

  for (auto&& c : target_entities) {
    auto const& obtained = remapped[c];
    auto const  expected = compute_density(c, targetMeshWrapper);
#ifdef DEBUG
std::printf("cell[%d]: remapped: %.3f, expected: %.3f\n", c, obtained, expected);
#endif
    ASSERT_NEAR(obtained, expected, epsilon);
  }
}
