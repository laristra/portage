/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "search_kdtree.h"

#include <algorithm>

#include "gtest/gtest.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"

TEST(search_kdtree2, cell) {
  Portage::Simple_Mesh sm{0, 0, 1, 1, 3, 3};
  Portage::Simple_Mesh tm{0, 0, 1, 1, 2, 2};
  const Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(sm);
  const Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(tm);

  Portage::SearchKDTree<2, Portage::Entity_kind::CELL,
                        Wonton::Simple_Mesh_Wrapper,
                        Wonton::Simple_Mesh_Wrapper>
      search(source_mesh_wrapper, target_mesh_wrapper);

  for (int tc = 0; tc < 4; ++tc) {
    std::vector<int> candidates = search(tc);

    // there should be four candidate source cells, in a square
    // compute scbase = index of lower left source cell
    ASSERT_EQ(4, candidates.size());
    const int tx = tc % 2;
    const int ty = tc / 2;
    const int scbase = tx + ty * 3;
    // candidates might not be in order, so sort them
    std::sort(candidates.begin(), candidates.end());
    ASSERT_EQ(scbase, candidates[0]);
    ASSERT_EQ(scbase + 1, candidates[1]);
    ASSERT_EQ(scbase + 3, candidates[2]);
    ASSERT_EQ(scbase + 4, candidates[3]);
  }

}  // TEST(search_kdtree2, cell)

TEST(search_kdtree2, node) {
  Portage::Simple_Mesh sm{0, 0, 1, 1, 3, 3};
  Portage::Simple_Mesh tm{0, 0, 1, 1, 2, 2};
  const Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(sm);
  const Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(tm);

  Portage::SearchKDTree<2, Portage::Entity_kind::NODE,
                        Wonton::Simple_Mesh_Wrapper,
                        Wonton::Simple_Mesh_Wrapper>
      search(source_mesh_wrapper, target_mesh_wrapper);

  for (int tc = 0; tc < 9; ++tc) {
    std::vector<int> candidates = search(tc);

    // there should be four candidate source nodes, in a square
    // compute snbase = index of lower left source node
    ASSERT_EQ(4, candidates.size());
    const int tx = tc % 3;
    const int ty = tc / 3;
    const int snbase = tx + ty * 4;
    // candidates might not be in order, so sort them
    std::sort(candidates.begin(), candidates.end());
    ASSERT_EQ(snbase, candidates[0]);
    ASSERT_EQ(snbase + 1, candidates[1]);
    ASSERT_EQ(snbase + 4, candidates[2]);
    ASSERT_EQ(snbase + 5, candidates[3]);
  }

}  // TEST(search_kdtree2, node)
