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

  // useful shortcuts
  using Wonton::Entity_kind;
  using Wonton::Entity_type;
  using Wonton::Jali_Mesh_Wrapper;
  using Wonton::Jali_State_Wrapper;

  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                                          Wonton::Jali_Mesh_Wrapper,
                                          Wonton::Jali_State_Wrapper>;
  using Partition = Portage::Parts<2, Wonton::Entity_kind::CELL,
                                      Wonton::Jali_Mesh_Wrapper,
                                      Wonton::Jali_State_Wrapper>;

  // useful constants
  double const upper_bound  = std::numeric_limits<double>::max();
  double const lower_bound  = -upper_bound;
  double const epsilon = 1.E-10;

  // Source and target meshes and states
  std::shared_ptr<Jali::Mesh>  source_mesh;
  std::shared_ptr<Jali::State> source_state;
  std::shared_ptr<Jali::Mesh>  target_mesh;
  std::shared_ptr<Jali::State> target_state;

  source_mesh  = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5);
  target_mesh  = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 10, 10);
  source_state = Jali::State::create(source_mesh);
  target_state = Jali::State::create(target_mesh);

  Jali_Mesh_Wrapper  source_mesh_wrapper(*source_mesh);
  Jali_Mesh_Wrapper  target_mesh_wrapper(*target_mesh);
  Jali_State_Wrapper source_state_wrapper(*source_state);
  Jali_State_Wrapper target_state_wrapper(*target_state);

  int const nb_source_cells =
    source_mesh_wrapper.num_entities(Entity_kind::CELL, Entity_type::ALL);

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
    source_density[c] = compute_density(c, source_mesh_wrapper);
  }

  source_state_wrapper.mesh_add_data(Entity_kind::CELL, "density", source_density);
  target_state_wrapper.mesh_add_data<double>(Entity_kind::CELL, "density", 0.);

  // assuming that entities ordering are column-wise,
  // let's pick a couple of aligned subsets of entites
  // from source and target meshes (5x5, 10x10).
  std::vector<int> source_cells {
    6, 11,
    7, 12,
    8, 13
  };

  std::vector<int> target_cells {
    22, 32, 42,	52,
    23, 33, 43,	53,
    24, 34, 44,	54,
    25, 35, 45,	55,
    26, 36, 46,	56,
    27, 37, 47,	57
  };

  // perform remap without redistribution, nor mismatch fixup
  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  Partition partition(source_mesh_wrapper, target_mesh_wrapper,
                      source_state_wrapper,target_state_wrapper,
                      source_cells, target_cells, nullptr);

  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto source_weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);

  // check mismatch and compute cell volumes before
  partition.test_mismatch(source_weights);
  assert(not partition.has_mismatch());

  remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
    "density", "density", source_weights, lower_bound, upper_bound,
    Portage::DEFAULT_LIMITER, Portage::DEFAULT_PARTIAL_FIXUP_TYPE,
    Portage::DEFAULT_EMPTY_FIXUP_TYPE, Portage::DEFAULT_CONSERVATION_TOL,
    Portage::DEFAULT_MAX_FIXUP_ITER, &partition
  );

  // Finally check that we got the right target density values
  double* remapped;
  target_state_wrapper.mesh_get_data(Entity_kind::CELL, "density", &remapped);

  for (auto&& c : target_cells) {
    auto const& obtained = remapped[c];
    auto const  expected = compute_density(c, target_mesh_wrapper);
    #ifdef DEBUG
      std::printf("target[%02d]: remapped: %7.3f, expected: %7.3f\n",
                  c, obtained, expected);
    #endif
    ASSERT_NEAR(obtained, expected, epsilon);
  }
}
