/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <memory>
#include "gtest/gtest.h"

#ifdef WONTON_ENABLE_MPI
#include "mpi.h"
#endif

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "wonton/support/wonton.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/driver/coredriver.h"

/**
 * @brief Parts mismatch fixup tests for part-by-part interpolation.
 *
 * Here, the purpose is to check that parts mismatch are correctly
 * handled when interpolating each pair of source-target part independently.
 * More precisely, we aim to ensure that values within partially filled and
 * empty cells are correctly fixed according to the user-specified fixup
 * scheme which may be:
 * - constant and leave empty,
 * - constant and extrapolate,
 * - locally conservative and leave empty,
 * - locally conservative and extrapolate,
 * - shifted conservative and leave empty,
 * - shifted conservative and extrapolate.
 *
 * Again, source and target meshes are respectively split into two parts,
 * and each source-target part pair is remapped independently.
 * Remapped values are then compared to the exact expected values.
 *
 *  0,1     source      1,1
 *    ___________________      source mesh: 4x4 cartesian grid.
 *   |    |    :    |    |     Here, the domain is halved into two equal parts.
 *   |____|____:____|____|     The prescribed density is 100 for the first part
 *   |    |    :    |    |     and 1 for the second part.
 *   |____0____:____1____|
 *   |    |    :    |    |
 *   |____|____:____|____|
 *   |    |    :    |    |
 *   |____|____:____|____|
 *
 *    ___________________       target mesh: 5x5 cartesian grid.
 *   |   |   | . |   :   |      Here parts interface is shifted such that the
 *   |___|___|_._|___:___|      first part is 4/5 of the domain and the last
 *   |   |   | . |   :   |      part is only 1/5 of it.
 *   |___|___0_._|___:_1_|      Hence, we have a:
 *   |   |   | . |   :   |      - partially filled cells layer in [0.5-0.6]
 *   |___|___|_._|___:___|      - empty cells layer in [0.6-0.8]
 *   |   |   | . |   :   |
 *   |___|___|_._|___:___|
 *
 *  0,0     target      1,0
 */
class PartMismatchTest : public testing::Test {

protected:
  // useful shortcuts
  using Remapper = Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                                          Wonton::Jali_Mesh_Wrapper,
                                          Wonton::Jali_State_Wrapper>;
  using PartPair = Portage::PartPair<2, Wonton::Jali_Mesh_Wrapper,
                                        Wonton::Jali_State_Wrapper>;
  using Partial = Portage::Partial_fixup_type;
  using Empty   = Portage::Empty_fixup_type;

  /**
   * @brief Compute source cell density.
   *
   * @param c: the source cell local index
   * @return the exact density of the given cell.
   */
  double compute_source_density(int c) {
    double const rho_min = 1.;
    double const rho_max = 100.;
    auto centroid = source_mesh->cell_centroid(c);
    return (centroid[0] < 0.5 ? rho_max : rho_min);
  };

  /**
   * @brief Create a partition based on a threshold value.
   *
   * @tparam is_source source or target mesh?
   * @param thresh     x-axis threshold value
   */
  template<bool is_source>
  void create_partition(double thresh) {

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
      auto k = centroid[0] < thresh ? 0 : 1;
      part[k].emplace_back(i);
    }

    for (int i = 0; i < nb_parts; ++i)
      part[i].shrink_to_fit();
  }


  /**
   * @brief Compute the expected remapped density on target mesh.
   *
   * It depends on Partial/empty mismatch fixup schemes being used:
   * - locally-conservative/leave-empty: [100.0|100.0| 50.0|  0.0|1.0]
   * - locally-conservative/extrapolate: [100.0|100.0| 50.0| 50.0|1.0]
   * - constant/leave-empty:             [100.0|100.0|100.0|  0.0|1.0]
   * - constant/extrapolate:             [100.0|100.0|100.0|100.0|1.0]
   * - shifted-conservative/leave-empty: [ 83.3| 83.3| 83.3|  0.0|2.5]
   * - shifted-conservative/extrapolate: [ 62.5| 62.5| 62.5| 62.5|2.5]
   *
   * @param cell          the target cell local index
   * @param partial_fixup the Partial mismatch fixup scheme to use
   * @param empty_fixup   the empty cell fixup scheme to use
   * @return the expected remapped density value
   */
  double get_expected_remapped_density(int const cell,
                                       Partial partial_fixup,
                                       Empty   empty_fixup) {
    // get cell position
    auto centroid = target_mesh->cell_centroid(cell);
    auto const& x = centroid[0];

    switch (partial_fixup) {
      case Partial::LOCALLY_CONSERVATIVE: {
        if (empty_fixup == Empty::LEAVE_EMPTY) {
          if (x < 0.4) {
            return 100.;
          } else if (x < 0.6) {
            return 50.;
          } else if (x < 0.8) {
            return 0.;
          } else {
            return 1.;
          }
        } else /* EXTRAPOLATE */ {
          if (x < 0.4) {
            return 100.;
          } else if (x < 0.8) {
            return 50.;
          } else {
            return 1.;
          }
        }
      } case Partial::CONSTANT: {
        if (empty_fixup == Empty::LEAVE_EMPTY) {
          if (x < 0.6) {
            return 100.;
          } else if (x < 0.8) {
            return 0.;
          } else {
            return 1.;
          }
        } else /* EXTRAPOLATE */ {
          if (x < 0.8) {
            return 100.;
          } else {
            return 1.;
          }
        }
      } case Partial::GLOBALLY_CONSERVATIVE: {

        /* Correct target cell density for shifted-conservative fixup scheme. */
        auto compute_shifted_density = [](double mass_source,
                                          double unit_mass_target,
                                          double unit_volume_target,
                                          int nb_target_cells) {
          double mass_delta = (unit_mass_target * nb_target_cells) - mass_source;
          double discrepancy = mass_delta / nb_target_cells;
          return (unit_mass_target - discrepancy) / unit_volume_target;
        };

        if (empty_fixup == Empty::LEAVE_EMPTY) {
          if (x < 0.6) {
            return compute_shifted_density(50, 20, 0.2, 3);  // 83.35
          } else if(x < 0.8) {
            return 0.;
          } else {
            return compute_shifted_density(0.5, 0.2, 0.2, 1);  // 2.5
          }
        } else /* EXTRAPOLATE */ {
          if (x < 0.8) {
            return compute_shifted_density(50, 20, 0.2, 4);  // 62.5
          } else {
            return compute_shifted_density(0.5, 0.2, 0.2, 1);  // 2.5
          }
        } // end extrapolate
      } default: throw std::runtime_error("Invalid Partial fixup type");
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
  PartMismatchTest()
    : source_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 4, 4)),
      target_mesh(Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 5, 5)),
      source_state(Jali::State::create(source_mesh)),
      target_state(Jali::State::create(target_mesh)),
      source_mesh_wrapper(*source_mesh),
      target_mesh_wrapper(*target_mesh),
      source_state_wrapper(*source_state),
      target_state_wrapper(*target_state)
  {
    // get rid of long namespaces
    auto const CELL = Wonton::Entity_kind::CELL;

    // compute and add density field to the source mesh
    int const nb_cells = source_mesh_wrapper.num_owned_cells();
    double source_density[nb_cells];
    for (int c = 0; c < nb_cells; c++) {
      source_density[c] = compute_source_density(c);
    }

    source_state_wrapper.mesh_add_data(CELL, "density", source_density);
    target_state_wrapper.mesh_add_data<double>(CELL, "density", 0.);

    // create source and target mesh parts
    create_partition<true>(0.5);
    create_partition<false>(0.8);

    // set parts
    for (int i = 0; i < nb_parts; ++i) {
      parts.emplace_back(source_mesh_wrapper, source_state_wrapper,
                         target_mesh_wrapper, target_state_wrapper,
                         source_cells[i], target_cells[i], nullptr);
    }
  }

  /**
   * @brief Do the test.
   *
   * @param partial_fixup the partially filled cell fixup scheme to use.
   * @param empty_fixup   the empty cells fixup scheme to use.
   */
  void unitTest(Partial partial_fixup, Empty empty_fixup) {

    using Wonton::Entity_kind;
    using Wonton::Entity_type;

    // Perform remapping without redistribution
    Remapper remapper(source_mesh_wrapper, source_state_wrapper,
                      target_mesh_wrapper, target_state_wrapper);

    auto candidates = remapper.search<Portage::SearchKDTree>();
    auto weights = remapper.intersect_meshes<Portage::IntersectRnD>(candidates);

    for (int i = 0; i < nb_parts; ++i) {
      // interpolate density part-by-part while fixing mismatched values
      remapper.interpolate_mesh_var<double, Portage::Interpolate_1stOrder>(
        "density", "density", weights, &(parts[i])
      );
      
      // test for mismatch and compute volumes
      parts[i].check_mismatch(weights);

      if (parts[i].has_mismatch())
        parts[i].fix_mismatch("density", "density", lower_bound, upper_bound, 
        Portage::DEFAULT_NUMERIC_TOLERANCES<2>.relative_conservation_eps,
        Portage::DEFAULT_NUMERIC_TOLERANCES<2>.max_num_fixup_iter,
        partial_fixup, empty_fixup);
        
      }

    // Finally check that we got the right target density values
    double* remapped;
    target_state_wrapper.mesh_get_data(Entity_kind::CELL, "density", &remapped);

    for (int i = 0; i < nb_parts; ++i) {
      for (auto&& c : target_cells[i]) {
        auto obtained = remapped[c];
        auto expected = get_expected_remapped_density(c, partial_fixup, empty_fixup);
#if DEBUG_PART_BY_PART
        auto centroid = target_mesh->cell_centroid(c);
        std::printf("target[%02d]: (x=%.1f, y=%.1f), "
                    "remapped: %7.3f, expected: %7.3f\n",
                    c, centroid[0], centroid[1], obtained, expected);
#endif
        ASSERT_NEAR(obtained, expected, epsilon);
      }
    }
  }

protected:
  // useful constants
  static constexpr int const nb_parts = 2;
  static constexpr double upper_bound = std::numeric_limits<double>::max();
  static constexpr double lower_bound = std::numeric_limits<double>::min();
  static constexpr double epsilon = 1.E-10;

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


TEST_F(PartMismatchTest, LocallyConservative_LeaveEmpty) {
  unitTest(Partial::LOCALLY_CONSERVATIVE, Empty::LEAVE_EMPTY);
}

TEST_F(PartMismatchTest, LocallyConservative_Extrapolate) {
  unitTest(Partial::LOCALLY_CONSERVATIVE, Empty::EXTRAPOLATE);
}

TEST_F(PartMismatchTest, Constant_LeaveEmpty) {
  unitTest(Partial::CONSTANT, Empty::LEAVE_EMPTY);
}

TEST_F(PartMismatchTest, Constant_Extrapolate) {
  unitTest(Partial::CONSTANT, Empty::EXTRAPOLATE);
}

TEST_F(PartMismatchTest, ShiftedConservative_LeaveEmpty) {
  unitTest(Partial::GLOBALLY_CONSERVATIVE, Empty::LEAVE_EMPTY);
}

TEST_F(PartMismatchTest, ShiftedConservative_Extrapolate) {
  unitTest(Partial::GLOBALLY_CONSERVATIVE, Empty::EXTRAPOLATE);
}
