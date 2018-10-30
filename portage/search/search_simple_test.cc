/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

// portage includes
#include "portage/search/search_simple.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

TEST(search_simple, case1) {
  Wonton::Simple_Mesh sm{0, 0, 1, 1, 3, 3};
  Wonton::Simple_Mesh tm{0, 0, 1, 1, 2, 2};
  const Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(sm);
  const Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(tm);

  Portage::SearchSimple<Wonton::Simple_Mesh_Wrapper,
                        Wonton::Simple_Mesh_Wrapper>
      search(source_mesh_wrapper, target_mesh_wrapper);

  for (int tc = 0; tc < 4; ++tc) {
    std::vector<int> candidates;
    search(tc, &candidates);

    // there should be four candidate source cells, in a square
    // compute scbase = index of lower left source cell
    ASSERT_EQ(4, candidates.size());
    const int tx = tc % 2;
    const int ty = tc / 2;
    const int scbase = tx + ty * 3;
    ASSERT_EQ(scbase, candidates[0]);
    ASSERT_EQ(scbase + 1, candidates[1]);
    ASSERT_EQ(scbase + 3, candidates[2]);
    ASSERT_EQ(scbase + 4, candidates[3]);
  }

}  // TEST(search_simple, case1)

class MeshWrapperDual {
 public:
  MeshWrapperDual(const Wonton::Simple_Mesh_Wrapper &w) : w_(w) {}
  int num_owned_cells() const { return w_.num_owned_nodes(); }
  int num_ghost_cells() const { return w_.num_ghost_nodes(); }
  void cell_get_coordinates(int const cellid,
                            std::vector<Portage::Point<2>> *pplist) const {
    w_.dual_cell_get_coordinates(cellid, pplist);
  }

 private:
  const Wonton::Simple_Mesh_Wrapper &w_;
};

TEST(search_simple, dual) {
  Wonton::Simple_Mesh sm{0, 0, 1, 1, 3, 3};
  Wonton::Simple_Mesh tm{0, 0, 1, 1, 2, 2};
  const Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(sm);
  const Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(tm);

  const MeshWrapperDual s2(source_mesh_wrapper);
  const MeshWrapperDual t2(target_mesh_wrapper);

  Portage::SearchSimple<MeshWrapperDual, MeshWrapperDual> search(s2, t2);

  for (int tc = 0; tc < 9; ++tc) {
    std::vector<int> candidates;
    search(tc, &candidates);

    // there should be four candidate source nodes, in a square
    // compute snbase = index of lower left source node
    ASSERT_EQ(4, candidates.size());
    const int tx = tc % 3;
    const int ty = tc / 3;
    const int snbase = tx + ty * 4;
    ASSERT_EQ(snbase, candidates[0]);
    ASSERT_EQ(snbase + 1, candidates[1]);
    ASSERT_EQ(snbase + 4, candidates[2]);
    ASSERT_EQ(snbase + 5, candidates[3]);
  }

}  // TEST(search_simple, dual)
