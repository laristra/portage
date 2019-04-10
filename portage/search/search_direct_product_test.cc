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
#include "wonton/support/Point.h"


TEST(search_direct_product, Test1D) {
  /*
   * Overlaps: -------------------------------
   *      tgt(0):                 src(0):
   *        src(0)                  tgt(0)
   *        src(1)                src(1):
   *      tgt(1):                   tgt(0)
   *        src(1)                  tgt(1)
   *        src(2)                src(2):
   *                                tgt(1)
  */

  // dimensionality
  const int D = 1;

  // Create meshes
  const std::vector<double> x_tgt = {0.0, 0.5, 1.0};
  const std::array<std::vector<double>,D> edges_tgt = {x_tgt};
  Wonton::Direct_Product_Mesh<D> tgt(edges_tgt);
  const std::vector<double> x_src = {0.00, 0.25, 0.75, 1.00};
  const std::array<std::vector<double>,D> edges_src = {x_src};
  Wonton::Direct_Product_Mesh<D> src(edges_src);

  // Create wrappers
  const Wonton::Direct_Product_Mesh_Wrapper<D> tgt_wrapper(tgt);
  const Wonton::Direct_Product_Mesh_Wrapper<D> src_wrapper(src);

  // Declare search
  Portage::SearchDirectProduct<D, Wonton::Direct_Product_Mesh_Wrapper<D>,
    Wonton::Direct_Product_Mesh_Wrapper<D>> search(src_wrapper, tgt_wrapper);

  // Verify overlaps
  for (int i = 0; i < tgt_wrapper.axis_num_cells(0); ++i) {
    std::array<int,D> indices = {i};
    Wonton::CellID id = tgt_wrapper.indices_to_cellid(indices);
    const std::vector<Wonton::CellID> candidates = search(id);
    ASSERT_EQ(candidates.size(), 2);
    int n = 0;
    for (int i2 = i; i2 < i+2; ++i2) {
      ASSERT_TRUE(n < candidates.size());
      std::array<int,D> candidate = {i2};
      Wonton::CellID c_id = src_wrapper.indices_to_cellid(candidate);
      ASSERT_EQ(candidates[n], c_id);
      n++;
    }
  }

}  // TEST(search_direct_product, Test1D)

// ============================================================================

TEST(search_direct_product, Test2D) {
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
  const std::array<std::vector<double>,D> edges_tgt = {x_tgt, y_tgt};
  Wonton::Direct_Product_Mesh<D> tgt(edges_tgt);
  const std::vector<double> x_src = {0.00, 0.25, 0.75, 1.00};
  const std::vector<double> y_src = {0.00, 0.25, 0.75, 1.00};
  const std::array<std::vector<double>,D> edges_src = {x_src, y_src};
  Wonton::Direct_Product_Mesh<D> src(edges_src);

  // Create wrappers
  const Wonton::Direct_Product_Mesh_Wrapper<D> tgt_wrapper(tgt);
  const Wonton::Direct_Product_Mesh_Wrapper<D> src_wrapper(src);

  // Declare search
  Portage::SearchDirectProduct<D, Wonton::Direct_Product_Mesh_Wrapper<D>,
    Wonton::Direct_Product_Mesh_Wrapper<D>> search(src_wrapper, tgt_wrapper);

  // Verify overlaps
  for (int j = 0; j < tgt_wrapper.axis_num_cells(1); ++j) {
    for (int i = 0; i < tgt_wrapper.axis_num_cells(0); ++i) {
      std::array<int,D> indices = {i,j};
      Wonton::CellID id = tgt_wrapper.indices_to_cellid(indices);
      const std::vector<Wonton::CellID> candidates = search(id);
      ASSERT_EQ(candidates.size(), 4);
      int n = 0;
      for (int j2 = j; j2 < j+2; ++j2) {
        for (int i2 = i; i2 < i+2; ++i2) {
          ASSERT_TRUE(n < candidates.size());
          std::array<int,D> candidate = {i2,j2};
          Wonton::CellID c_id = src_wrapper.indices_to_cellid(candidate);
          ASSERT_EQ(candidates[n], c_id);
          n++;
        }
      }
    }
  }

}  // TEST(search_direct_product, Test2D)

// ============================================================================

TEST(search_direct_product, Test3D) {
  /*
   * Overlaps: -------------------------------
   *      tgt(0,0,0):             src(0,0,0):
   *        src(0,0,0)              tgt(0,0,0)
   *        src(0,0,1)            src(0,0,1):
   *        src(0,1,0)              tgt(0,0,0)
   *        src(0,1,1)              tgt(0,0,1)
   *        src(1,0,0)            src(0,0,2):
   *        src(1,0,1)              tgt(0,0,1)
   *        src(1,1,0)            src(0,1,0):
   *        src(1,1,1)              tgt(0,0,0)
   *      tgt(0,0,1):               tgt(0,1,0)
   *        src(0,0,1)            src(0,1,1):
   *        src(0,0,2)              tgt(0,0,0)
   *        src(0,1,1)              tgt(0,0,1)
   *        src(0,1,2)              tgt(0,1,0)
   *        src(1,0,1)              tgt(0,1,1)
   *        src(1,0,2)            src(0,1,2):
   *        src(1,1,1)              tgt(0,0,1)
   *        src(1,1,2)              tgt(0,1,1)
   *      tgt(0,1,0):             src(0,2,0):
   *        src(0,1,0)              tgt(0,1,0)
   *        src(0,1,1)            src(0,2,1):
   *        src(0,2,0)              tgt(0,1,0)
   *        src(0,2,1)              tgt(0,1,1)
   *        src(1,1,0)            src(0,2,2):
   *        src(1,1,1)              tgt(0,1,1)
   *        src(1,2,0)            src(1,0,0):
   *        src(1,2,1)              tgt(0,0,0)
   *      tgt(0,1,1):               tgt(1,0,0)
   *        src(0,1,1)            src(1,0,1):
   *        src(0,1,2)              tgt(0,0,0)
   *        src(0,2,1)              tgt(0,0,1)
   *        src(0,2,2)              tgt(1,0,0)
   *        src(1,1,1)              tgt(1,0,1)
   *        src(1,1,2)            src(1,0,2):
   *        src(1,2,1)              tgt(0,0,1)
   *        src(1,2,2)              tgt(1,0,1)
   *      tgt(1,0,0):             src(1,1,0):
   *        src(1,0,0)              tgt(0,0,0)
   *        src(1,0,1)              tgt(0,1,0)
   *        src(1,1,0)              tgt(1,0,0)
   *        src(1,1,1)              tgt(1,1,0)
   *        src(2,0,0)            src(1,1,1):
   *        src(2,0,1)              tgt(0,0,0)
   *        src(2,1,0)              tgt(0,0,1)
   *        src(2,1,1)              tgt(0,1,0)
   *      tgt(1,0,1):               tgt(0,1,1)
   *        src(1,0,1)              tgt(1,0,0)
   *        src(1,0,2)              tgt(1,0,1)
   *        src(1,1,1)              tgt(1,1,0)
   *        src(1,1,2)              tgt(1,1,1)
   *        src(2,0,1)            src(1,1,2):
   *        src(2,0,2)              tgt(0,0,1)
   *        src(2,1,1)              tgt(0,1,1)
   *        src(2,1,2)              tgt(1,0,1)
   *      tgt(1,1,0):               tgt(1,1,1)
   *        src(1,1,0)            src(1,2,0):
   *        src(1,1,1)              tgt(0,1,0)
   *        src(1,2,0)              tgt(1,1,0)
   *        src(1,2,1)            src(1,2,1):
   *        src(2,1,0)              tgt(0,1,0)
   *        src(2,1,1)              tgt(0,1,1)
   *        src(2,2,0)              tgt(1,1,0)
   *        src(2,2,1)              tgt(1,1,1)
   *      tgt(1,1,1):             src(1,2,2):
   *        src(1,1,1)              tgt(0,1,1)
   *        src(1,1,2)              tgt(1,1,1)
   *        src(1,2,1)            src(2,0,0):
   *        src(1,2,2)              tgt(1,0,0)
   *        src(2,1,1)            src(2,0,1):
   *        src(2,1,2)              tgt(1,0,0)
   *        src(2,2,1)              tgt(1,0,1)
   *        src(2,2,2)            src(2,0,2):
   *                                tgt(1,0,1)
   *                              src(2,1,0):
   *                                tgt(1,0,0)
   *                                tgt(1,1,0)
   *                              src(2,1,1):
   *                                tgt(1,0,0)
   *                                tgt(1,0,1)
   *                                tgt(1,1,0)
   *                                tgt(1,1,1)
   *                              src(2,1,2):
   *                                tgt(1,0,1)
   *                                tgt(1,1,1)
   *                              src(2,2,0):
   *                                tgt(1,1,0)
   *                              src(2,2,1):
   *                                tgt(1,1,0)
   *                                tgt(1,1,1)
   *                              src(2,2,2):
   *                                tgt(1,1,1)
   */

  // dimensionality
  const int D = 3;

  // Create meshes
  const std::vector<double> x_tgt = {0.0, 0.5, 1.0};
  const std::vector<double> y_tgt = {0.0, 0.5, 1.0};
  const std::vector<double> z_tgt = {0.0, 0.5, 1.0};
  const std::array<std::vector<double>,D> edges_tgt = {x_tgt, y_tgt, z_tgt};
  Wonton::Direct_Product_Mesh<D> tgt(edges_tgt);
  const std::vector<double> x_src = {0.00, 0.25, 0.75, 1.00};
  const std::vector<double> y_src = {0.00, 0.25, 0.75, 1.00};
  const std::vector<double> z_src = {0.00, 0.25, 0.75, 1.00};
  const std::array<std::vector<double>,D> edges_src = {x_src, y_src, z_src};
  Wonton::Direct_Product_Mesh<D> src(edges_src);

  // Create wrappers
  const Wonton::Direct_Product_Mesh_Wrapper<D> tgt_wrapper(tgt);
  const Wonton::Direct_Product_Mesh_Wrapper<D> src_wrapper(src);

  // Declare search
  Portage::SearchDirectProduct<D, Wonton::Direct_Product_Mesh_Wrapper<D>,
    Wonton::Direct_Product_Mesh_Wrapper<D>> search(src_wrapper, tgt_wrapper);

  // Verify overlaps
  for (int k = 0; k < tgt_wrapper.axis_num_cells(2); ++k) {
    for (int j = 0; j < tgt_wrapper.axis_num_cells(1); ++j) {
      for (int i = 0; i < tgt_wrapper.axis_num_cells(0); ++i) {
        std::array<int,D> indices = {i,j,k};
        Wonton::CellID id = tgt_wrapper.indices_to_cellid(indices);
        const std::vector<Wonton::CellID> candidates = search(id);
        ASSERT_EQ(candidates.size(), 8);
        int n = 0;
        for (int k2 = k; k2 < k+2; ++k2) {
          for (int j2 = j; j2 < j+2; ++j2) {
            for (int i2 = i; i2 < i+2; ++i2) {
              ASSERT_TRUE(n < candidates.size());
              std::array<int,D> candidate = {i2,j2,k2};
              Wonton::CellID c_id = src_wrapper.indices_to_cellid(candidate);
              ASSERT_EQ(candidates[n], c_id);
              n++;
            }
          }
        }
      }
    }
  }

}  // TEST(search_direct_product, Test3D)

