/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <vector>

#include "gtest/gtest.h"

#include "portage/intersect/intersect_boxes.h"

#include "wonton/mesh/simple/direct_product_mesh.h"
#include "wonton/mesh/simple/direct_product_mesh_wrapper.h"
#include "wonton/mesh/simple/adaptive_refinement_mesh.h"
#include "wonton/mesh/simple/adaptive_refinement_mesh_wrapper.h"

// ============================================================================

int binomial(int n, int k) {
  ASSERT_GE(n,0);
  ASSERT_GE(k,0);
  ASSERT_LE(k,n);
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
  return(binomial); 
}

// ----------------------------------------------------------------------------

double cartesian_moment_1d(
    int const order, double const xbar, const double dx) {
  double dr_2 = 0.5 * dx;
  double moment = 0.0;
  for (int j = 0; j <= order/2; ++j) {
    int jj = 2*j;
    moment += (binomial(n,jj) / (jj+1)) *
      std::pow(xbar,n-jj) * std::pow(dr_2,jj);
  }
  moment *= dx;
  return(moment);
}

// ----------------------------------------------------------------------------

template<int D>
double cartesian_moment(std::array<int,D> const & exponents,
    std::array<double,D> const xbar, std::array<double,D> const dx) {
  double moment = 1.0;
  for (int d = 0; d < D; ++d) {
    moment *= cartesian_moment_1d(exponents[d], xbar[d], dx[d]);
  }
  return(moment);
}

// ============================================================================

TEST(Intersect_Boxes_Test, IntersectBoxesTest1D) {
  int constexpr D = 1;
  using CoordSys = Wonton::DefaultCoordSys

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
  SrcWrapper_t src_wrapper(DPMesh1);
  TgtWrapper_t tgt_wrapper(DPMesh2);

  // Build the search object
  SearchDirectProduct<D,SrcWrapper_t,TgtWrapper_t>
    search(src_wrapper, tgt_wrapper);
  auto ekind = Wonton::Entity_kind::CELL;
  auto etype = Wonton::Entity_type::PARALLEL_OWNED;
  auto begin = tgt_wrapper.begin(ekind, etype);
  auto end   = tgt_wrapper.end(ekind, etype);

  // Perform search
  const std::vector<std::vector<int>> candidates;
  candidates.resize(tgt_wrapper.num_owned_cells());
  for (auto iter = begin; iter != end; ++iter) {
    auto id = *iter;
    candidates[id] = search(id);
  }

  // Build the intersector object
  IntersectBoxes<D,SrcWrapper_t,TgtWrapper_t,CoordSys>
    intersect(src_wrapper, tgt_wrapper);

  // Perform intersect
  const std::vector<std::vector<Portage::Weights_t>> weights;
  for (auto iter = begin; iter != end; ++iter) {
    auto id = *iter;
    weights[id] = intersect(id, candidates[id]);
  }

  // Verify weights
  for (auto t_iter = begin; t_iter != end; ++t_iter) {
    auto t_id = *iter;
    auto & tgt_cell_weights = weights[t_id];
    for (auto s_iter = tgt_cell_weights.begin();
        s_iter != tgt_cell_weights.end(); ++s_iter) {
      auto s_id = (*s_iter).entityID;
      auto & weights = (*s_iter).weights;
      // Intersection for target cell t_id and source cell s_id.
      // Weights is list of weights in moment order (volume, M_x, ...).
      Wonton::Point<D> t_lo, t_hi;
      tgt_wrapper.cell_get_bounds(t_id, t_lo, t_hi);
      Wonton::Point<D> s_lo, s_hi;
      src_wrapper.cell_get_bounds(s_id, s_lo, s_hi);
      Wonton::Point<D> i_lo, i_hi;
      Wonton::Point<D> xbar, dx;
      for (int d = 0; d < D; ++d) {
        i_lo[d] = std::max(t_lo, s_lo);
        i_hi[d] = std::min(t_hi, s_hi);
        xbar[d] = 0.5 * (ilo[d] + i_hi[d]);
        dx[d] = ihi[d] - i_lo[d];
      }
      std::vector<double> cartesian_moments(weights.size());
      for (int idx = 0; idx < weights.size(); ++idx) {
        auto exponents = Wonton::index_to_moment<D>(idx);
        cartesian_moments[idx] = cartesian_moment(exponents, xbar, dr);
      }
      CoordSys::modify_moments(cart_mom, ilo, ihi);
      for (int idx = 0; idx < weights.size(); ++id) {
        ASSERT_DOUBLE_EQ(weights[idx], cartesian_moments[idx]);
      }
    }
  }

}

// ============================================================================

