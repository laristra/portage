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

/**
 * @brief Basic test for part-by-part interpolation.
 *
 * Here, source and target meshes are partitioned into two parts,
 * and each part is remapped independently. Remapped values are then
 * compared to the exact values given by an analytical function.
 * Here source and target parts are perfectly aligned (no mismatch),
 * but target mesh resolution is twice that of source mesh.
 * The generated parts looks like below:
 *
 *  0,1                  1,1
 *    ---------:---------
 *   |   s1    :    t1   |
 *   |    _____:__       |
 *   |   |     :  |      |
 *   |   |     :  |      |
 *   |   | s0  :t0|      |
 *   |   |     :  |      |
 *   |   |     :  |      |
 *   |    -----:--       |
 *   | r=30    : r=100   |
 *    ---------:---------
 *  0,0                 1,0
 */
class PartDriverTest : public testing::Test {

protected:
  // useful shortcuts
  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                                          Wonton::Jali_Mesh_Wrapper,
                                          Wonton::Jali_State_Wrapper>;
  using PartsPair = Portage::Parts<2, Wonton::Entity_kind::CELL,
                                      Wonton::Jali_Mesh_Wrapper,
                                      Wonton::Jali_State_Wrapper>;
  /**
   * @brief compute cell centroid.
   *
   * @param cell the given cell.
   * @param the mesh
   * @return centroid its centroid.
   */
  Wonton::Point<2> get_centroid(int cell, Wonton::Jali_Mesh_Wrapper const& mesh) {
    Wonton::Point<2> centroid;
    mesh.cell_centroid(cell, &centroid);
    return std::move(centroid);
  };

  /**
   * @brief compute cell density analytically.
   *
   * @param c    the given cell.
   * @param mesh a wrapper to the supporting mesh[source|target].
   * @return the computed cell density.
   */
  double compute_density(int c, Wonton::Jali_Mesh_Wrapper const& mesh) {
    double const x_gap = 0.4;
    double const t_min = 30.;
    double const t_max = 100.;

    auto centroid = get_centroid(c, mesh);
    return (centroid[0] < x_gap ? t_min : t_max);
  };

  /**
   * @brief Create a partition based on a threshold value.
   *
   * @param mesh     the current mesh to split
   * @param nb_cells its number of cells
   * @param thresh   x-axis threshold value
   * @param part     source or target mesh parts
   */
  void create_partition(Wonton::Jali_Mesh_Wrapper const& mesh, std::vector<int>* part) {
    assert(part != nullptr);
    assert(nb_parts == 2);

    auto const CELL = Wonton::Entity_kind::CELL;
    auto const ALL  = Wonton::Entity_type::ALL;

    int const nb_cells = mesh.num_entities(CELL, ALL);
    int const min_heap_size = static_cast<int>(nb_cells / 2);

    for (int i = 0; i < 2; ++i) {
      part[i].clear();
      part[i].reserve(min_heap_size);
    }

    for (int i = 0; i < nb_cells; ++i) {
      auto centroid = get_centroid(i, mesh);
      auto const& x = centroid[0];
      auto const& y = centroid[1];
      // inside block
      if (x > 0.2 and x < 0.6 and y > 0.2 and y < 0.6) {
        part[0].emplace_back(i);
      } else {
        part[1].emplace_back(i);
      }
    }
  }

public:
  /**
   * @brief Setup each test-case.
   *
   * It creates both source and target meshes,
   * then computes and assigns a density field on source mesh,
   * then creates parts couples for both source and target meshes.
   */
  PartDriverTest()
    : source_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5)),
      target_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 10, 10)),
      source_state(Jali::State::create(source_mesh)),
      target_state(Jali::State::create(target_mesh)),
      source_mesh_wrapper(*source_mesh),
      target_mesh_wrapper(*target_mesh),
      source_state_wrapper(*source_state),
      target_state_wrapper(*target_state)
  {
    // get rid of long namespaces
    auto const CELL = Wonton::Entity_kind::CELL;
    auto const ALL  = Wonton::Entity_type::ALL;

    // compute and add density field to the source mesh
    int const nb_source_cells = source_mesh_wrapper.num_entities(CELL, ALL);
    double source_density[nb_source_cells];
    for (int c = 0; c < nb_source_cells; c++) {
      source_density[c] = compute_density(c, source_mesh_wrapper);
    }

    source_state_wrapper.mesh_add_data(CELL, "density", source_density);
    target_state_wrapper.mesh_add_data<double>(CELL, "density", 0.);

    // create parts by picking entities within (0.2,0.2) and (0.6,0.6)
    create_partition(source_mesh_wrapper, source_cells);
    create_partition(target_mesh_wrapper, target_cells);

    // set parts
    for (int i = 0; i < 2; ++i) {
      parts.emplace_back(source_mesh_wrapper, target_mesh_wrapper,
                         source_state_wrapper,target_state_wrapper,
                         source_cells[i], target_cells[i], nullptr);
    }
  }

  // useful constants
  static constexpr double upper_bound = std::numeric_limits<double>::max();
  static constexpr double lower_bound = -upper_bound;
  static constexpr double epsilon = 1.E-10;

  // Source and target meshes and states
  std::shared_ptr<Jali::Mesh>  source_mesh;
  std::shared_ptr<Jali::Mesh>  target_mesh;
  std::shared_ptr<Jali::State> source_state;
  std::shared_ptr<Jali::State> target_state;

  // Wrappers for interfacing with the underlying mesh data structures
  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper;
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper;
  Wonton::Jali_State_Wrapper source_state_wrapper;
  Wonton::Jali_State_Wrapper target_state_wrapper;

  // source and target parts couples
  std::vector<PartsPair> parts;
  std::vector<int> source_cells[2];
  std::vector<int> target_cells[2];
};


TEST_F(PartDriverTest, Basic) {

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  // remap without redistribution, nor mismatch fixup
  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto source_weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);

  // check mismatch and compute cell volumes before
  for (int i = 0; i < 2; ++i) {
    // test for mismatch and compute volumes
    parts[i].test_mismatch(source_weights);
    assert(not partition.has_mismatch());

    // interpolate density part-by-part while fixing mismatched values
    remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
      "density", "density", source_weights, lower_bound, upper_bound,
      Portage::DEFAULT_LIMITER, Portage::DEFAULT_PARTIAL_FIXUP_TYPE,
      Portage::DEFAULT_EMPTY_FIXUP_TYPE, Portage::DEFAULT_CONSERVATION_TOL,
      Portage::DEFAULT_MAX_FIXUP_ITER, &(parts[i])
    );
  }

  // Finally check that we got the right target density values
  double* remapped;
  target_state_wrapper.mesh_get_data(Wonton::Entity_kind::CELL, "density", &remapped);

  for (int i = 0; i < 2; ++i) {
    for (auto&& c : target_cells[i]) {
      auto const& obtained = remapped[c];
      auto const  expected = compute_density(c, target_mesh_wrapper);
      #ifdef DEBUG_PART_BY_PART
        std::printf("target[%02d]: remapped: %7.3f, expected: %7.3f\n",
                    c, obtained, expected);
      #endif
      ASSERT_NEAR(obtained, expected, epsilon);
    }
  }
}
