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

#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"

TEST(search_kdtree3, case1)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    // overlay a 2x2x2 target mesh on a 3x3x3 source mesh
    // each target mesh cell gives eight candidate source cells
    const std::shared_ptr<Jali::Mesh> smesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
    const std::shared_ptr<Jali::Mesh> tmesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
    const Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    const Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    Portage::SearchKDTree<3,
        Portage::Jali_Mesh_Wrapper, Portage::Jali_Mesh_Wrapper>
        search(source_mesh_wrapper, target_mesh_wrapper);

    for (int tc = 0; tc < 8; ++tc) {
        std::vector<int> candidates;
        search(tc, &candidates);

        // there should be eight candidate source cells, in a cube
        // compute scbase = index of lower left source cell
        ASSERT_EQ(8, candidates.size());
        const int tx = tc % 2;
        const int ty = (tc / 2) % 2;
        const int tz = tc / 4;
        const int scbase = tx + ty * 3 + tz * 9;
        // candidates might not be in order, so sort them
        std::sort(candidates.begin(), candidates.end());
        ASSERT_EQ(scbase,      candidates[0]);
        ASSERT_EQ(scbase + 1,  candidates[1]);
        ASSERT_EQ(scbase + 3,  candidates[2]);
        ASSERT_EQ(scbase + 4,  candidates[3]);
        ASSERT_EQ(scbase + 9,  candidates[4]);
        ASSERT_EQ(scbase + 10, candidates[5]);
        ASSERT_EQ(scbase + 12, candidates[6]);
        ASSERT_EQ(scbase + 13, candidates[7]);
    }

} // TEST(search_kdtree3, case1)


