/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <vector>

#include "gtest/gtest.h"

#include "portage/intersect/intersect_boxes.h"
#include "portage/search/search_direct_product.h"

#include "wonton/mesh/direct_product/direct_product_mesh.h"
#include "wonton/mesh/direct_product/direct_product_mesh_wrapper.h"
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh_wrapper.h"

// ============================================================================

int binomial(int n, int k) {
  int l = std::min(k, n-k); 
  int binomial; 
  switch(l) { 
  case 0:
    binomial = 1; 
    break;
  case 1:
    binomial = n;
    break;
  default: 
    std::vector<int> coefs(l+1,0);
    coefs[0] = 1; 
    for (int i = 1; i < n; ++i) {
      for (int j = std::min(l,i); j >= 1; --j) {
        coefs[j] += coefs[j-1]; 
      } 
    } 
    binomial = coefs[l] + coefs[l-1];
    break; 
  } 
  return binomial; 
}

// ----------------------------------------------------------------------------

double cartesian_moment_1d(
    int const order, double const xbar, const double dx) {
  double dr_2 = 0.5 * dx;
  double moment = 0.0;
  for (int j = 0; j <= order/2; ++j) {
    int jj = 2*j;
    moment += (binomial(order,jj) / (jj+1)) *
      std::pow(xbar,order-jj) * std::pow(dr_2,jj);
  }
  moment *= dx;
  return moment;
}

// ----------------------------------------------------------------------------

template<int D>
double cartesian_moment(std::array<int,D> const & exponents,
    Wonton::Point<D> const xbar, Wonton::Point<D> const dx) {
  double moment = 1.0;
  for (int d = 0; d < D; ++d) {
    moment *= cartesian_moment_1d(exponents[d], xbar[d], dx[d]);
  }
  return moment;
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest1D) {
  int constexpr D = 1;
  // TODO: Template on this then run in all coordinate systems?
  using CoordSys = Wonton::CartesianCoordinates;

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
  std::vector<std::vector<Portage::Weights_t>> weights;
  for (auto iter = begin; iter != end; ++iter) {
    auto id = *iter;
    weights[id] = intersect(id, candidates[id]);
  }

  // Verify weights
  for (auto t_iter = begin; t_iter != end; ++t_iter) {
    auto t_id = *t_iter;
    auto & tgt_cell_weights = weights[t_id];
    for (auto s_iter = tgt_cell_weights.begin();
        s_iter != tgt_cell_weights.end(); ++s_iter) {
      auto s_id = (*s_iter).entityID;
      auto & weights = (*s_iter).weights;
      // Verify that the search method provided enough moments (0th and 1st)
      ASSERT_GE(weights.size(), Wonton::count_moments<D>(1));
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
      // Compute Cartesian moments
      std::vector<double> cartesian_moments(
          Wonton::count_moments<D>(1 + CoordSys::moment_shift));
      for (int idx = 0; idx < cartesian_moments.size(); ++idx) {
        int order;
        std::array<int,D> exponents;
        std::tie(order,exponents) = Wonton::index_to_moment<D>(idx);
        cartesian_moments[idx] = cartesian_moment<D>(exponents, xbar, dx);
      }
      // Shift moments to appropriate coordinate system
      std::vector<double> shifted_moments(cartesian_moments);
      CoordSys::shift_moments_list<D>(shifted_moments);
      // Use axis-aligned-box optimizations to convert moments to appropriate
      // coordinate system
      auto volume = cartesian_moments[0];
      CoordSys::modify_volume<D>(volume, i_lo, i_hi);
      Wonton::Point<D> first_moments;
      for (int d = 0; d < D; ++d) {
        first_moments[d] = cartesian_moments[1+d];
      }
      CoordSys::modify_first_moments<D>(first_moments, i_lo, i_hi);
      std::vector<double> box_moments(1+D);
      box_moments[0] = volume;
      for (int d = 0; d < D; ++d) {
        box_moments[1+d] = first_moments[d];
      }
      // Check that results are equal from all three methods
      // -- intersect (weights)
      // -- analytic result with moment-shift method (shifted_moments)
      // -- analytic result with axis-aligned boxes (box_moments)
      for (int idx = 0; idx < Wonton::count_moments<D>(1); ++idx) {
        ASSERT_DOUBLE_EQ(weights[idx], shifted_moments[idx]);
        ASSERT_DOUBLE_EQ(weights[idx], box_moments[idx]);
      }
    }
  }

}

// ============================================================================

