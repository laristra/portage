/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <memory>
#include "gtest/gtest.h"

#include "wonton/support/wonton.h"

#ifdef WONTON_ENABLE_MPI
  #include "mpi.h"
#endif

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/driver/coredriver.h"

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

/**
 * @brief Base fixture class for any order part-by-part remap tests.
 *
 * Here, source and target meshes are partitioned into two parts,
 * and each part is remapped independently. Remapped values are then
 * compared to the exact values given by an analytical function.
 * Here source and target parts are perfectly aligned (no mismatch),
 * but the target mesh is twice finer than the source mesh.
 */
class PartBaseTest : public testing::Test {

protected:
  // useful shortcuts
  using Remapper = Portage::CoreDriver<
    2, Wonton::Entity_kind::CELL,
    Wonton::Jali_Mesh_Wrapper,
    Wonton::Jali_State_Wrapper
  >;

  using PartPair = Portage::PartPair<2, Wonton::Jali_Mesh_Wrapper,
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
    assert(nb_parts == 2);

    auto size = is_source ? source_mesh_wrapper.num_owned_cells()
                          : target_mesh_wrapper.num_owned_cells();
    auto mesh = is_source ? source_mesh : target_mesh;
    auto part = is_source ? source_cells : target_cells;

    for (int i = 0; i < nb_parts; ++i) {
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

    for (int i = 0; i < nb_parts; ++i)
      part[i].shrink_to_fit();
  }

  /**
   * @brief Populate the part pairs partition.
   *
   */
  void populate_parts(double x_min, double x_max, double y_min, double y_max) {
    assert(std::abs(x_max - x_min) > epsilon);
    assert(std::abs(y_max - y_min) > epsilon);

    // create parts by picking entities within a certain range
    create_partition<true>(x_min, x_max, y_min, y_max);
    create_partition<false>(x_min, x_max, y_min, y_max);

    // set parts
    for (int i = 0; i < nb_parts; ++i) {
      parts.emplace_back(source_mesh_wrapper, source_state_wrapper,
                         target_mesh_wrapper, target_state_wrapper,
                         source_cells[i], target_cells[i], nullptr);
    }
  }

  /**
   * @brief Setup each test-case.
   *
   * It creates both source and target meshes,
   * then computes and assigns a density field on source mesh,
   * then creates parts couples for both source and target meshes.
   */
  PartBaseTest()
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
  }

protected:
  // useful constants and aliases
  static constexpr double const upper_bound = std::numeric_limits<double>::max();
  static constexpr double const lower_bound = std::numeric_limits<double>::min();
  static constexpr double const epsilon = 1.E-10;
  static constexpr int const nb_parts = 2;
  static constexpr auto const CELL = Wonton::Entity_kind::CELL;

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
  std::vector<int> source_cells[nb_parts];
  std::vector<int> target_cells[nb_parts];
};


/**
 * @brief Fixture class for first-order remap tests.
 *
 * Here we consider a piecewise constant field on a 2D
 * cartesian grid with:
 *
 *  f(x,y) = |r_min if (x < r)
 *           |r_max otherwise
 *           given [r_min,r_max]=[30,100]
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
 *   | r=30    : r=100   |
 *    ---------:---------
 *  0,0       r         1,0
 *  Notice that we do not have any part mismatch in this case.
 *  Hence we expect a strictly conservative remap using a first-order
 *  part-by-part scheme.
 */
class PartOrderOneTest : public PartBaseTest {

public:
  /**
   * @brief Setup each test-case.
   *
   * It creates both source and target meshes, then computes
   * and assigns a density field on source mesh, then creates
   * parts by picking entities within (0.2,0.2) and (0.6,0.6).
   */
  PartOrderOneTest() : PartBaseTest() {
    populate_parts(0.2, 0.6, 0.2, 0.6);
  }
};

/**
 * @brief Fixture class for second-order remap tests.
 *
 * Here we consider a piecewise linear field on a 2D
 * cartesian grid with:
 *
 *  f(x,y) = |c * x, if (x < r).
 *           |c * (x - r) otherwise.
 *           given c > 0 and 0 < r < 1.
 *
 *            f(x,y)
 *            /   :     /
 *          /     :    /
 *        /       :   /
 *      /    [0]  :  / [1]
 *    /           : /
 *    ------------:------
 *   0            r     1
 * Notice that we would obtain smoothed field values near 'r'
 * with a global remap since the field gradient would be smoothed
 * in that region. However, we should obtain a perfect remap with
 * a clear field discontinuity using a second-order part-by-part remap,
 * since the computed gradient near 'r' is distinct for each part.
 */
class PartOrderTwoTest : public PartBaseTest {

public:
  /**
   * @brief Setup each test-case.
   *
   * It creates both source and target meshes, then computes
   * and assigns a density field on source mesh, then creates
   * parts by picking entities within (0.0,0.0) and (0.6,1.0).
   */
  PartOrderTwoTest() : PartBaseTest() {
    populate_parts(0.0, 0.6, 0.0, 1.0);
  }
};


/**
 * verify that first-order part-by-part remap
 * is strictly conservative for a piecewise constant field
 * in absence of mismatch between source and target parts.
 */
TEST_F(PartOrderOneTest, PiecewiseConstantField) {

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  double* original = nullptr;
  double* remapped = nullptr;

  // assign a piecewise constant field on source mesh
  source_state_wrapper.mesh_get_data(CELL, "density", &original);
  for (int c = 0; c < nb_source_cells; c++) {
    auto centroid = source_mesh->cell_centroid(c);
    original[c] = (centroid[0] < 0.40 ? 30. : 100.);
  }

  // process remap
  auto candidates = remapper.search<Portage::SearchKDTree>();
  auto weights = remapper.intersect_meshes<Portage::IntersectR2D>(candidates);

  for (int i = 0; i < 2; ++i) {

    // interpolate density for current part
    remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
      "density", "density", weights, &(parts[i])
    );
      
    // test for mismatch and compute volumes
    parts[i].check_mismatch(weights);
    assert(not parts[i].has_mismatch());    
  }

  // compare remapped values with analytically computed ones
  target_state_wrapper.mesh_get_data(CELL, "density", &remapped);

  for (auto&& current_part : target_cells) {
    for (auto&& c : current_part) {
      auto obtained = remapped[c];
      auto centroid = target_mesh->cell_centroid(c);
      auto expected = (centroid[0] < 0.40 ? 30. : 100.);
      #if DEBUG_PART_BY_PART
        std::printf("target[%02d]: remapped: %7.3f, expected: %7.3f\n",
                    c, obtained, expected);
      #endif
      ASSERT_NEAR(obtained, expected, epsilon);
    }
  }
}


/**
 * verify that both first-order part-by-part and global remap
 * are equivalent for first-order interpolation of general fields.
 * in absence of mismatch between source and target parts.
 */
TEST_F(PartOrderOneTest, GlobalRemapComparison) {

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

  for (int i = 0; i < 2; ++i) {
    // interpolate density for current part
    remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
      "density", "density", weights, &(parts[i])
    );
    
    // test for mismatch and compute volumes
    parts[i].check_mismatch(weights);
    assert(not parts[i].has_mismatch());

    // this block won't run since no mismatch, but make sure syntactically correct
    if (parts[i].has_mismatch())
      parts[i].fix_mismatch("density", "density", lower_bound, upper_bound);
  }

  // store the part-by-part remapped values
  target_state_wrapper.mesh_get_data(CELL, "density", &remapped);
  std::copy(remapped, remapped + nb_target_cells, remapped_parts);
  std::fill(remapped, remapped + nb_target_cells, 0.);

  // interpolate density on whole source mesh
  remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
    "density", "density", weights
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


/**
 * verify that second-order part-by-part remap
 * is strictly conservative for a piecewise linear field
 * in absence of mismatch between source and target parts.
 * Notice that no gradient limiter is used here.
 */
TEST_F(PartOrderTwoTest, PiecewiseLinearField) {

  Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                    target_mesh_wrapper, target_state_wrapper);

  double* original = nullptr;
  double* remapped = nullptr;

  double const coef = 2.;
  double const x_max = 0.6;

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

  for (int i = 0; i < 2; ++i) {
    // test for mismatch and compute volumes
    parts[i].check_mismatch(weights);
    assert(not parts[i].has_mismatch());

    auto const& source_part = parts[i].source();
    auto gradients = remapper.compute_source_gradient("density",
                                                      Portage::NOLIMITER,
                                                      Portage::BND_NOLIMITER,
                                                      0, &source_part);

    // interpolate density for current part
    remapper.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
      "density", "density", weights, &(parts[i]), &gradients
    );


    // this block won't run since no mismatch, but make sure syntactically correct
    if (parts[i].has_mismatch())
      parts[i].fix_mismatch("density", "density", lower_bound, upper_bound,
      Portage::DEFAULT_NUMERIC_TOLERANCES<2>.relative_conservation_eps,
      Portage::DEFAULT_NUMERIC_TOLERANCES<2>.max_num_fixup_iter,
      Portage::DEFAULT_PARTIAL_FIXUP_TYPE, Portage::DEFAULT_EMPTY_FIXUP_TYPE);
  }

  // compare remapped values with analytically computed ones
  target_state_wrapper.mesh_get_data(CELL, "density", &remapped);

  for (auto&& current_part : target_cells) {
    for (auto&& c : current_part) {
      auto obtained = remapped[c];
      auto centroid = target_mesh->cell_centroid(c);
      auto const& x = centroid[0];
      auto expected = (x < x_max ? coef * x : coef * (x - x_max));
      #if DEBUG_PART_BY_PART
        auto const& y = centroid[1];
        std::printf("target[%02d]: x: %7.3f, y: %7.3f"
                    ", remapped: %7.3f, expected: %7.3f\n",
                    c, x, y, obtained, expected);
      #endif
      ASSERT_NEAR(obtained, expected, epsilon);
    }
  }
}
