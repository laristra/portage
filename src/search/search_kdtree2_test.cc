/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "search_kdtree2.h"

#include <algorithm>

#include "gtest/gtest.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

TEST(search_kdtree2, case1)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    // overlay a 2x2 target mesh on a 3x3 source mesh
    // each target mesh cell gives four candidate source cells
    const std::unique_ptr<Jali::Mesh> smesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
    const std::unique_ptr<Jali::Mesh> tmesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
    const Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    const Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    Portage::SearchKDTree2<
        Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>
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
        // candidates might not be in order, so sort them
        std::sort(candidates.begin(), candidates.end());
        ASSERT_EQ(scbase,     candidates[0]);
        ASSERT_EQ(scbase + 1, candidates[1]);
        ASSERT_EQ(scbase + 3, candidates[2]);
        ASSERT_EQ(scbase + 4, candidates[3]);
    }

} // TEST(search_kdtree2, case1)

class MeshWrapperDual {
  public:
    MeshWrapperDual(const Portage::Jali_Mesh_Wrapper &w) : w_(w) {}
    int num_owned_cells() const { return w_.num_owned_nodes(); }
    int num_ghost_cells() const { return w_.num_ghost_nodes(); }
    void cell_get_coordinates(int const cellid,
            std::vector<std::pair<double, double>> *xylist) const {
        w_.dual_cell_get_coordinates(cellid, xylist);
    }
  private:
    const Portage::Jali_Mesh_Wrapper &w_;
};


TEST(search_kdtree2, dual)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    // overlay a 2x2 target mesh on a 3x3 source mesh
    // each target mesh node gives four candidate source nodes
    const std::unique_ptr<Jali::Mesh> smesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3, NULL, 
                                 true, true, true, true);
    const std::unique_ptr<Jali::Mesh> tmesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2, NULL, 
                                 true, true, true, true);
    const Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    const Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    const MeshWrapperDual s2(source_mesh_wrapper);
    const MeshWrapperDual t2(target_mesh_wrapper);

    Portage::SearchKDTree2<MeshWrapperDual, MeshWrapperDual>
        search(s2, t2);

    for (int tc = 0; tc < 9; ++tc) {
        std::vector<int> candidates;
        search(tc, &candidates);

        // there should be four candidate source nodes, in a square
        // compute snbase = index of lower left source node
        ASSERT_EQ(4, candidates.size());
        const int tx = tc % 3;
        const int ty = tc / 3;
        const int snbase = tx + ty * 4;
        // candidates might not be in order, so sort them
        std::sort(candidates.begin(), candidates.end());
        ASSERT_EQ(snbase,     candidates[0]);
        ASSERT_EQ(snbase + 1, candidates[1]);
        ASSERT_EQ(snbase + 4, candidates[2]);
        ASSERT_EQ(snbase + 5, candidates[3]);
    }

} // TEST(search_kdtree2, dual)

