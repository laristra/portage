/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

// portage includes
#include "portage/search/search_direct_product.h"

// wonton includes
#include "wonton/mesh/direct_product/direct_product_mesh.h"
#include "wonton/mesh/direct_product/direct_product_mesh_wrapper.h"
#include "wonton/wonton/support/CellID.h"
#include "wonton/wonton/support/IntPoint.h"
#include "wonton/support/Point.h"

TEST(search_direct_product, case1) {
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
  Wonton::Direct_Product_Mesh tgt(x_tgt,y_tgt);
  const std::vector<double> x_src = {0.00, 0.25, 0.75, 1.00};
  const std::vector<double> y_src = {0.00, 0.25, 0.75, 1.00};
  Wonton::Direct_Product_Mesh src(x_src,y_src);

  // Create wrappers
  const Wonton::Direct_Product_Mesh_Wrapper tgt_wrapper(tgt);
  const Wonton::Direct_Product_Mesh_Wrapper src_wrapper(src);

  // Declare search
  Portage::SearchDirectProduct<D, Wonton::Direct_Product_Mesh_Wrapper,
    Wonton::Direct_Product_Mesh_Wrapper> search(src_wrapper, tgt_wrapper);

  // Verify overlaps
  for (int i = 0; i < x_tgt.size() - 1; ++i) {
    for (int j = 0; j < y_tgt.size() - 1; ++j) {
      Wonton::IntPoint<D> indices = {i,j};
      Wonton::CellID id = tgt_wrapper.indices_to_cellid<D>(indices);
      const std::vector<Wonton::CellID> candidates = search(id);
      int n = 0;
      for (int i2 = i; i2 < i+2; ++i2) {
        for (int j2 = j; j2 < j+2; ++j2) {
          Wonton::IntPoint<D> candidate = {i2,j2};
          Wonton::CellID c_id = src_wrapper.indices_to_cellid<D>(candidate);
          ASSERT_EQ(candidates[n], c_id);
          n++;
        }
      }
    }
  }

}  // TEST(search_direct_product, case1)

