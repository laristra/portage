/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

// portage includes
#include "portage/search/search_direct_product.h"

// wonton includes
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh_wrapper.h"
#include "wonton/mesh/direct_product/direct_product_mesh.h"
#include "wonton/mesh/direct_product/direct_product_mesh_wrapper.h"
#include "wonton/wonton/support/CellID.h"
#include "wonton/wonton/support/IntPoint.h"
#include "wonton/support/Point.h"

// ============================================================================

namespace search_direct_product_test {
  double refinement_function(Wonton::Point<3> coords) {
    // A simple refinement function: Assuming 3D coordinates and a grid that
    // extends from 0 to 1 along all axes, refine the full mesh into octants
    // then refine the (++-) octant one additional level.
    if ((coords[0] > 0.5) && (coords[1] > 0.5) && (coords[2] < 0.5))
      return(2.0);
    else
      return(1.0);
  }
}  // namespace search_direct_product_test

// ============================================================================

TEST(search_direct_product, DPtoDP) {
  /*
   * Overlaps: -------------------------------
   *      tgt(0,0):               src(0,0):
   *        src(0,0)                tgt(0,0)
   *        src(0,1)              src(0,1):
   *        src(1,0)                tgt(0,0)
   *        src(1,1)                tgt(0,1)
   *      tgt(0,1)                src(0,2):
   *        src(0,1)                tgt(0,1)
   *        src(0,2)              src(1,0):
   *        src(1,1)                tgt(0,0)
   *        src(1,2)                tgt(1,0)
   *      tgt(1,0)                src(1,1):
   *        src(1,0)                tgt(0,0)
   *        src(1,1)                tgt(0,1)
   *        src(2,0)                tgt(1,0)
   *        src(2,1)                tgt(1,1)
   *      tgt(1,1)                src(1,2):
   *        src(1,1)                tgt(0,1)
   *        src(1,2)                tgt(1,1)
   *        src(2,1)              src(2,0):
   *        src(2,2)                tgt(1,0)
   *                              src(2,1):
   *                                tgt(1,0)
   *                                tgt(1,1)
   *                              src(2,2):
   *                                tgt(1,1)
   */
  // dimensionality
  const int D = 2;
  // Create meshes
  const std::vector<double> x_tgt = {0.0, 0.5, 1.0};
  const std::vector<double> y_tgt = {0.0, 0.5, 1.0};
  const std::vector<double> edges_tgt[D] = {x_tgt, y_tgt};
  Wonton::Direct_Product_Mesh<D> tgt(edges_tgt);
  const std::vector<double> x_src = {0.00, 0.25, 0.75, 1.00};
  const std::vector<double> y_src = {0.00, 0.25, 0.75, 1.00};
  const std::vector<double> edges_src[D] = {x_src, y_src};
  Wonton::Direct_Product_Mesh<D> src(edges_src);

  // Create wrappers
  const Wonton::Direct_Product_Mesh_Wrapper<D> tgt_wrapper(tgt);
  const Wonton::Direct_Product_Mesh_Wrapper<D> src_wrapper(src);

  // Declare search
  Portage::SearchDirectProduct<D, Wonton::Direct_Product_Mesh_Wrapper<D>,
    Wonton::Direct_Product_Mesh_Wrapper<D>> search(src_wrapper, tgt_wrapper);

  // Verify overlaps
  for (int j = 0; j < y_tgt.size() - 1; ++j) {
    for (int i = 0; i < x_tgt.size() - 1; ++i) {
      Wonton::IntPoint<D> indices = {i,j};
      Wonton::CellID id = tgt_wrapper.indices_to_cellid(indices);
      const std::vector<Wonton::CellID> candidates = search(id);
      int n = 0;
      for (int j2 = j; j2 < j+2; ++j2) {
        for (int i2 = i; i2 < i+2; ++i2) {
          Wonton::IntPoint<D> candidate = {i2,j2};
          Wonton::CellID c_id = src_wrapper.indices_to_cellid(candidate);
          ASSERT_EQ(candidates[n], c_id);
          n++;
        }
      }
    }
  }

}  // TEST(search_direct_product, DPtoDP)

// ============================================================================

TEST(search_direct_product, DPtoAR) {
  /*
   * Overlaps: -------------------------------
   *      tgt(0):                 src(0,0,0):
   *        src(0,0,0)              tgt(0)
   *      tgt(1):                 src(0,0,1):
   *        src(1,0,0)              tgt(11)
   *      tgt(2):                 src(0,1,0):
   *        src(0,1,0)              tgt(2)
   *      tgt(3):                 src(0,1,1):
   *        src(1,1,0)              tgt(13)
   *      tgt(4):                 src(1,0,0):
   *        src(1,1,0)              tgt(1)
   *      tgt(5):                 src(1,0,1):
   *        src(1,1,0)              tgt(12)
   *      tgt(6):                 src(1,1,0):
   *        src(1,1,0)              tgt(3)
   *      tgt(7):                   tgt(4)
   *        src(1,1,0)              tgt(5)
   *      tgt(8):                   tgt(6)
   *        src(1,1,0)              tgt(7)
   *      tgt(9):                   tgt(8)
   *        src(1,1,0)              tgt(9)
   *      tgt(10):                  tgt(10)
   *        src(1,1,0)            src(1,1,1):
   *      tgt(11):                  tgt(14)
   *        src(0,0,1)
   *      tgt(12):
   *        src(1,0,1)
   *      tgt(13):
   *        src(0,1,1)
   *      tgt(14):
   *        src(1,1,1)
   */
  // dimensionality
  const int D = 3;
  // Create meshes
  const Wonton::Point<D> plo = {0,0,0};
  const Wonton::Point<D> phi = {1,1,1};
  Wonton::Adaptive_Refinement_Mesh<D> tgt(
      &search_direct_product_test::refinement_function, plo, phi);
  const std::vector<double> x_src = {0.00, 0.5, 1.00};
  const std::vector<double> y_src = {0.00, 0.5, 1.00};
  const std::vector<double> z_src = {0.00, 0.5, 1.00};
  const std::vector<double> edges_src[D] = {x_src, y_src, z_src};
  Wonton::Direct_Product_Mesh<D> src(edges_src);

  // Create wrappers
  const Wonton::Adaptive_Refinement_Mesh_Wrapper<D> tgt_wrapper(tgt);
  const Wonton::Direct_Product_Mesh_Wrapper<D> src_wrapper(src);

  // Declare search
  Portage::SearchDirectProduct<D, Wonton::Direct_Product_Mesh_Wrapper<D>,
    Wonton::Adaptive_Refinement_Mesh_Wrapper<D>> search(
        src_wrapper, tgt_wrapper);

  // Verify overlaps
  // -- Given the ID of a target cell, find the matching source cells
  /*
   * 0 --> 000  \
   * 1 --> 100   >  x = id%2, y = ((id-x)/2)%2, z = (((id-x)/2)-y)/2    CHECK THIS
   * 2 --> 010  /
   * 3-10 --> 110   -- always x=1, y=1, z=0
   * 11 --> 001   \
   * 12 --> 101    \    I think this is the same as 0-2, but first subtract 7 CHECK THIS
   * 13 --> 011    /
   * 14 --> 111   /
   */
  for (int n = 0; n < tgt_wrapper.total_num_cells(); ++n) {
    Wonton::CellID id = (Wonton::CellID) n;
    const std::vector<Wonton::CellID> candidates = search(id);
  }
  /*for (int k = 0; k < z_tgt.size() - 1; ++k) {
    for (int j = 0; j < y_tgt.size() - 1; ++j) {
      for (int i = 0; i < x_tgt.size() - 1; ++i) {
        Wonton::IntPoint<D> indices = {i,j,k};
        Wonton::CellID id = src_wrapper.indices_to_cellid(indices);
        const std::vector<Wonton::CellID> candidates = search(id);
        int n = 0;
        for (int j2 = j; j2 < j+2; ++j2) {
          for (int i2 = i; i2 < i+2; ++i2) {
            Wonton::IntPoint<D> candidate = {i2,j2};
            Wonton::CellID c_id = src_wrapper.indices_to_cellid(candidate);
            ASSERT_EQ(candidates[n], c_id);
            n++;
          }
        }
      }
    }
  }*/

}  // TEST(search_direct_product, DPtoAR)

