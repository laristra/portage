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
#include "portage/interpolate/interpolate_2nd_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "portage/driver/coredriver.h"

/**
 * @brief Basic sanity tests for part-by-part interpolation.
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
 *   |   s1    :         |
 *   |    _____:__       |
 *   |   |     :  |      |
 *   |   |     :  |      |
 *   |   | s0  :  |      |
 *   |   |     :  |      |
 *   |   |     :  |      |
 *   |    -----:--       |
 *   | r=f(x,y): r=g(x,y)|
 *    ---------:---------
 *  0,0                 1,0
 */
class PartOrderTwoTest : public testing::Test {

protected:
  // useful shortcuts
  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                                          Wonton::Jali_Mesh_Wrapper,
                                          Wonton::Jali_State_Wrapper>;
  using PartPair = Portage::PartPair<2, Wonton::Entity_kind::CELL,
                                        Wonton::Jali_Mesh_Wrapper,
                                        Wonton::Jali_State_Wrapper>;

  /**
   * @brief Create a partition based on a threshold value.
   *
   * @tparam is_source: for the source mesh?
   * @param x_min: min centroid x-value threshold for first part.
   * @param x_max: max centroid x-value threshold for first part.
   * @param y_min: min centroid y-value threshold for first part.
   * @param y_max: max centroid y-value threshold for first part.
   */
  template<bool is_source>
  void create_partition(double x_min, double x_max, double y_min, double y_max) {

    auto size = is_source ? source_mesh_wrapper.num_owned_cells()
                          : target_mesh_wrapper.num_owned_cells();
    auto mesh = is_source ? source_mesh : target_mesh;
    auto part = is_source ? source_cells : target_cells;

    for (int i = 0; i < 2; ++i) {
      part[i].clear();
      part[i].reserve(size);
    }

    for (int i = 0; i < size; ++i) {
      auto centroid = mesh->cell_centroid(i);
      auto const& x = centroid[0];
      auto const& y = centroid[1];
      // check if inside block
      auto k = x > x_min and x < x_max and y > y_min and y < y_max ? 0 : 1;
      part[k].emplace_back(i);
    }

    for (int i = 0; i < 2; ++i)
      part[i].shrink_to_fit();
  }

public:
  /**
   * @brief Setup each test-case.
   *
   * It creates both source and target meshes,
   * then computes and assigns a density field on source mesh,
   * then creates parts couples for both source and target meshes.
   */
  PartOrderTwoTest()
    : source_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5)),
      target_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 10, 10)),
      source_state(Jali::State::create(source_mesh)),
      target_state(Jali::State::create(target_mesh)),
      source_mesh_wrapper(*source_mesh),
      target_mesh_wrapper(*target_mesh),
      source_state_wrapper(*source_state),
      target_state_wrapper(*target_state)
  {
    // save meshes sizes
    nb_source_cells = source_mesh_wrapper.num_owned_cells();
    nb_target_cells = target_mesh_wrapper.num_owned_cells();

    // add density field to both meshes
    source_state_wrapper.mesh_add_data<double>(CELL, "density", 0.);
    target_state_wrapper.mesh_add_data<double>(CELL, "density", 0.);
    double* original = nullptr;
    source_state_wrapper.mesh_get_data(CELL, "density", &original);
    std::printf("original[0]: %.3f\n", original[0]);

    // create parts by picking entities within (0.0,0.0) and (0.6,1.0)
    create_partition<true>(0.0, 0.6, 0.0, 1.0);
    create_partition<false>(0.0, 0.6, 0.0, 1.0);

    // set parts
    for (int i = 0; i < 2; ++i) {
      parts.emplace_back(source_mesh_wrapper, source_state_wrapper,
                         target_mesh_wrapper, target_state_wrapper,
                         source_cells[i], target_cells[i], nullptr);
    }
  }

protected:
  // useful constants and aliases
  static constexpr double upper_bound = std::numeric_limits<double>::max();
  static constexpr double lower_bound = std::numeric_limits<double>::min();
  static constexpr double epsilon = 1.E-10;
  static constexpr auto CELL = Wonton::Entity_kind::CELL;

  int nb_source_cells = 0;
  int nb_target_cells = 0;

  // source and target meshes and states
  std::shared_ptr<Jali::Mesh>  source_mesh;
  std::shared_ptr<Jali::Mesh>  target_mesh;
  std::shared_ptr<Jali::State> source_state;
  std::shared_ptr<Jali::State> target_state;

  // wrappers for interfacing with the underlying mesh data structures
  Wonton::Jali_Mesh_Wrapper  source_mesh_wrapper;
  Wonton::Jali_Mesh_Wrapper  target_mesh_wrapper;
  Wonton::Jali_State_Wrapper source_state_wrapper;
  Wonton::Jali_State_Wrapper target_state_wrapper;

  // source and target parts couples
  std::vector<PartPair> parts;
  std::vector<int> source_cells[2];
  std::vector<int> target_cells[2];
};

// sanity check 1: verify that part-by-part interpolation
// is strictly conservative for a piecewise constant field
// in absence of mismatch between source and target parts.
TEST_F(PartOrderTwoTest, PiecewiseLinearField) {

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  double* original = nullptr;
  double* remapped = nullptr;

  double const coef = 5.;
  double const x_max = 0.6;
  return;
  // assign a linear field on source mesh
  source_state_wrapper.mesh_get_data(CELL, "density", &original);
  for (int c = 0; c < nb_source_cells; c++) {
    auto centroid = source_mesh->cell_centroid(c);
    auto const& x = centroid[0];
    original[c] = (x < x_max ? coef * x : coef * (x - x_max));
  }

  // process remap
  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);
  bool has_mismatch = remapper.check_mesh_mismatch(weights);

  for (int i = 0; i < 2; ++i) {
    // test for mismatch and compute volumes
    parts[i].check_mismatch(weights);
    assert(not parts[i].has_mismatch());

    // interpolate density for current part
    remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
      "density", "density", weights, lower_bound, upper_bound,
      Portage::NOLIMITER, Portage::BND_NOLIMITER,
      Portage::DEFAULT_PARTIAL_FIXUP_TYPE, Portage::DEFAULT_EMPTY_FIXUP_TYPE,
      Portage::DEFAULT_CONSERVATION_TOL,
      Portage::DEFAULT_MAX_FIXUP_ITER, &(parts[i])
    );
  }


  // compare remapped values with analytically computed ones
  target_state_wrapper.mesh_get_data(CELL, "density", &remapped);

  for (auto const& part : target_cells) {
    for (auto&& c : part) {
      auto obtained = remapped[c];
      auto centroid = target_mesh->cell_centroid(c);
      auto const& x = centroid[0];
      auto expected = (x < x_max ? coef * x : coef * (x - x_max));
      #if DEBUG_PART_BY_PART
        std::printf("target[%02d]: remapped: %7.3f, expected: %7.3f\n",
                    c, obtained, expected);
      #endif
      ASSERT_NEAR(obtained, expected, epsilon);
    }
  }
}

/*
// sanity check 2: verify that both part-by-part and mesh-mesh
// interpolation schemes are equivalent for general fields
// in absence of mismatch between source and target parts.
TEST_F(PartOrderTwoTest, MeshMeshRemapComparison) {

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  double* original = nullptr;
  double* remapped = nullptr;
  double remapped_parts[nb_target_cells];

  // assign a gaussian density field on source mesh
  source_state_wrapper.mesh_get_data(CELL, "density", &original);
  for (int c = 0; c < nb_source_cells; c++) {
    auto centroid = source_mesh->cell_centroid(c);
    auto const& x = centroid[0];
    auto const& y = centroid[1];
    original[c] = std::exp(-10.*(x*x + y*y));
  }

  // process remap
  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);
  bool has_mismatch = remapper.check_mesh_mismatch(weights);

  for (int i = 0; i < 2; ++i) {
    // test for mismatch and compute volumes
    parts[i].check_mismatch(weights);
    assert(not parts[i].has_mismatch());

    // interpolate density for current part
    remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
      "density", "density", weights, lower_bound, upper_bound,
      Portage::DEFAULT_LIMITER, Portage::DEFAULT_BND_LIMITER, 
      Portage::DEFAULT_PARTIAL_FIXUP_TYPE, Portage::DEFAULT_EMPTY_FIXUP_TYPE, 
      Portage::DEFAULT_CONSERVATION_TOL,
      Portage::DEFAULT_MAX_FIXUP_ITER, &(parts[i])
    );
  }

  // store the part-by-part remapped values
  target_state_wrapper.mesh_get_data(CELL, "density", &remapped);
  std::copy(remapped, remapped + nb_target_cells, remapped_parts);
  std::fill(remapped, remapped + nb_target_cells, 0.);

  // interpolate density on whole source mesh
  remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
    "density", "density", weights, lower_bound, upper_bound
  );

  // now compare remapped value for each cell
  for (int c=0; c < nb_target_cells; ++c) {
    auto const& value_mesh_remap = remapped[c];
    auto const& value_part_remap = remapped_parts[c];
    #if DEBUG_PART_BY_PART
      std::printf("target[%02d]: value_mesh_remap: %7.3f, value_part_remap: %7.3f\n",
                  c, value_mesh_remap, value_part_remap);
    #endif
    ASSERT_NEAR(value_mesh_remap, value_part_remap, epsilon);
  }
}
 */