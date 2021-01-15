/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <string>
#include <vector>

#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/direct_product/direct_product_mesh.h"
#include "wonton/mesh/direct_product/direct_product_mesh_wrapper.h"
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh_wrapper.h"

#include "portage/intersect/intersect_boxes.h"
#include "portage/search/search_direct_product.h"

// ============================================================================

// Computes the integral of x^order dx from xbar - dx/2 to xbar + dx/2
double cartesian_moment_1d(
    int const order, double const xbar, const double dx) {
  double dr_2 = 0.5 * dx;
  double moment = 0.0;
  moment += std::pow(xbar + dr_2, order+1);
  moment -= std::pow(xbar - dr_2, order+1);
  moment /= (order + 1);
  return moment;
}

// ----------------------------------------------------------------------------

// Computes the moment specified by the exponents array, assuming that this
// cell is in Cartesian space (ignoring corrections from the coordinate system)
template<int D>
double cartesian_moment(std::array<int,D> const & exponents,
    Wonton::Point<D> const xbar, Wonton::Point<D> const dx) {
  double moment = 1.0;
  for (int d = 0; d < D; ++d) {
    moment *= cartesian_moment_1d(exponents[d], xbar[d], dx[d]);
  }
  return moment;
}

// ----------------------------------------------------------------------------

// Computes the moment specified by the exponents array, using the correct
// volume element for the coordinate system.
template<int D, typename CoordSys>
double analytic_moment(std::array<int,D> const & exponents,
    Wonton::Point<D> const xbar, Wonton::Point<D> const dx) {
  double moment = cartesian_moment_1d(
      exponents[0] + CoordSys::moment_shift, xbar[0], dx[0]);
  for (int d = 1; d < D; ++d) {
    moment *= cartesian_moment_1d(exponents[d], xbar[d], dx[d]);
  }
  moment *= CoordSys::moment_coefficient * CoordSys::inv_geo_fac;
  return moment;
}

// ----------------------------------------------------------------------------

// Computes the moment specified by the exponents array, using the correct
// volume element for the coordinate system.
// -- Not easy to compute for 3D spherical coordinates, so skip it
template<>
double analytic_moment<3,Wonton::Spherical3DCoordinates>(
    std::array<int,3> const & exponents,
    Wonton::Point<3> const xbar, Wonton::Point<3> const dx) {
  throw std::runtime_error("The analytic_moment method is not implemented for 3D spherical coordinates");
}

// ============================================================================

template<int D, typename CoordSys>
void intersect_test() {
  using SrcMesh_t = Wonton::Direct_Product_Mesh<D,CoordSys>;
  using TgtMesh_t = Wonton::Direct_Product_Mesh<D,CoordSys>;
  using SrcWrapper_t = Wonton::Direct_Product_Mesh_Wrapper<D,CoordSys>;
  using TgtWrapper_t = Wonton::Direct_Product_Mesh_Wrapper<D,CoordSys>;

  // Build meshes
  //     0.0   0.2   0.4   0.6   0.8   1.0
  // src: |     |                 |     |
  // tgt: |           |     |           |
  std::vector<double> src_axis_points_1d = {0.0, 0.2, 0.8, 1.0};
  std::array<std::vector<double>,D> src_axis_points;
  for (int d = 0; d < D; ++d) {
    src_axis_points[d] = src_axis_points_1d;
  }
  SrcMesh_t src_mesh(src_axis_points);
  std::vector<double> tgt_axis_points_1d = {0.0, 0.4, 0.6, 1.0};
  std::array<std::vector<double>,D> tgt_axis_points;
  for (int d = 0; d < D; ++d) {
    tgt_axis_points[d] = tgt_axis_points_1d;
  }
  TgtMesh_t tgt_mesh(tgt_axis_points);

  // Build mesh wrappers
  SrcWrapper_t src_wrapper(src_mesh);
  TgtWrapper_t tgt_wrapper(tgt_mesh);

  // Build the search object
  Portage::SearchDirectProduct<D,SrcWrapper_t,TgtWrapper_t>
    search(src_wrapper, tgt_wrapper);
  auto ekind = Wonton::Entity_kind::CELL;
  auto etype = Wonton::Entity_type::PARALLEL_OWNED;
  auto begin = tgt_wrapper.begin(ekind, etype);
  auto end   = tgt_wrapper.end(ekind, etype);

  // Perform search
  std::vector<std::vector<int>> candidates;
  candidates.resize(tgt_wrapper.num_owned_cells());
  for (auto iter = begin; iter != end; ++iter) {
    auto id = *iter;
    candidates[id] = search(id);
  }

  // Build the intersector object
  Portage::IntersectBoxes<D,SrcWrapper_t,TgtWrapper_t,CoordSys>
    intersect(src_wrapper, tgt_wrapper);

  // Perform intersect
  std::vector<std::vector<Portage::Weights_t>> weights(
      tgt_wrapper.num_owned_cells());
  for (auto iter = begin; iter != end; ++iter) {
    auto id = *iter;
    weights[id] = intersect(id, candidates[id]);
  }

  // How many moments are we going to check?
  constexpr int max_order = 1;
  constexpr int max_order_with_shift = max_order + CoordSys::moment_shift;
  constexpr int num_moments = Wonton::count_moments<D>(max_order);
  constexpr int num_moments_with_shift = Wonton::count_moments<D>(max_order_with_shift);

  // Verify weights
  for (auto t_iter = begin; t_iter != end; ++t_iter) {
    auto t_id = *t_iter;
    auto & tgt_cell_weights = weights[t_id];
    for (auto s_iter = tgt_cell_weights.begin();
              s_iter != tgt_cell_weights.end();
              s_iter++) {
      auto s_id = (*s_iter).entityID;
      auto & current_weights = (*s_iter).weights;

      // Verify that the search method provided enough moments (0th and 1st)
      ASSERT_GE(current_weights.size(), unsigned(num_moments));

      // Intersection for target cell t_id and source cell s_id.
      // Weights is list of weights in moment order (volume, M_x, ...).
      Wonton::Point<D> t_lo, t_hi;
      tgt_wrapper.cell_get_bounds(t_id, &t_lo, &t_hi);
      Wonton::Point<D> s_lo, s_hi;
      src_wrapper.cell_get_bounds(s_id, &s_lo, &s_hi);
      Wonton::Point<D> i_lo, i_hi;
      Wonton::Point<D> xbar, dx;
      for (int d = 0; d < D; ++d) {
        i_lo[d] = std::max(t_lo[d], s_lo[d]);
        i_hi[d] = std::min(t_hi[d], s_hi[d]);
        xbar[d] = 0.5 * (i_lo[d] + i_hi[d]);
        dx[d] = i_hi[d] - i_lo[d];
      }

      // Compute moments analytically
      std::vector<double> analytic_moments(num_moments);
      int const num_analytic_moments = analytic_moments.size();
      for (int idx = 0; idx < num_analytic_moments; ++idx) {
        int order;
        std::array<int,D> exponents{};
        std::tie(order,exponents) = Wonton::index_to_moment<D>(idx);
        analytic_moments[idx] = analytic_moment<D,CoordSys>(exponents, xbar, dx);
      }

      // Compute Cartesian moments
      std::vector<double> cartesian_moments(num_moments_with_shift);
      int const num_cartesian_moments = cartesian_moments.size();
      for (int idx = 0; idx < num_cartesian_moments; ++idx) {
        int order;
        std::array<int,D> exponents{};
        std::tie(order,exponents) = Wonton::index_to_moment<D>(idx);
        cartesian_moments[idx] = cartesian_moment<D>(exponents, xbar, dx);
      }

      // Shift moments to appropriate coordinate system
      std::vector<double> shifted_moments(cartesian_moments);
      CoordSys::template shift_moments_list<D>(shifted_moments);

      // Use axis-aligned-box optimizations to convert moments to appropriate
      // coordinate system
      // -- Unpack into volume and first moments
      auto volume = cartesian_moments[0];
      Wonton::Point<D> first_moments;
      for (int d = 0; d < D; ++d) {
        first_moments[d] = cartesian_moments[1+d];
      }
      // -- Modify according to the coordinate system
      CoordSys::template modify_volume<D>(volume, i_lo, i_hi);
      CoordSys::template modify_first_moments<D>(first_moments, i_lo, i_hi);
      // -- Repack into box_moments array
      std::vector<double> box_moments(1+D,0);
      box_moments[0] = volume;
      for (int d = 0; d < D; ++d) {
        box_moments[1+d] = first_moments[d];
      }

      // Check that results are equal from all methods
      // -- intersect (weights)
      // -- analytic result computed in coordinate system (analytic_moments)
      // -- analytic result with moment-shift method (shifted_moments)
      // -- analytic result with axis-aligned boxes (box_moments)
      for (int idx = 0; idx < num_moments; ++idx) {
        std::string name;
        switch (idx) {
          case 0 : name = "volume"; break;
          case 1 : name = "first moment x"; break;
          case 2 : name = "first moment y"; break;
          case 3 : name = "first moment z"; break;
          default: break;
        }
        constexpr double TOL = 5e-15;
        ASSERT_NEAR(current_weights[idx], analytic_moments[idx], TOL * current_weights[idx])
            << " With index of: " << idx << " (" << name << ")\n"
            << "          xbar: " << xbar << "\n"
            << "            dx: " << dx;
        ASSERT_NEAR(current_weights[idx], box_moments[idx], TOL * current_weights[idx])
            << " With index of: " << idx << " (" << name << ")\n"
            << "          xbar: " << xbar << "\n"
            << "            dx: " << dx;
        ASSERT_NEAR(current_weights[idx], shifted_moments[idx], TOL * current_weights[idx])
            << " With index of: " << idx << " (" << name << ")\n"
            << "          xbar: " << xbar << "\n"
            << "            dx: " << dx;
      }
    }
  }

}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest1DCartesian) {
  int constexpr D = 1;
  using CoordSys = Wonton::CartesianCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest2DCartesian) {
  int constexpr D = 2;
  using CoordSys = Wonton::CartesianCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest3DCartesian) {
  int constexpr D = 3;
  using CoordSys = Wonton::CartesianCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest1DCylindrical) {
  int constexpr D = 1;
  using CoordSys = Wonton::CylindricalRadialCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest2DCylindricalA) {
  int constexpr D = 2;
  using CoordSys = Wonton::CylindricalAxisymmetricCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest2DCylindricalP) {
  int constexpr D = 2;
  using CoordSys = Wonton::CylindricalPolarCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest3DCylindrical) {
  int constexpr D = 3;
  using CoordSys = Wonton::Cylindrical3DCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest1DSpherical) {
  int constexpr D = 1;
  using CoordSys = Wonton::SphericalRadialCoordinates;
  intersect_test<D,CoordSys>();
}

// ============================================================================

