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
    std::unique_ptr<Jali::Mesh> smesh = mf(0.0,0.0,1.0,1.0,3,3);
    std::unique_ptr<Jali::Mesh> tmesh = mf(0.0,0.0,1.0,1.0,2,2);
    Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    Portage::SearchKDTree2<Portage::Jali_Mesh_Wrapper,Portage::Jali_Mesh_Wrapper> search(source_mesh_wrapper, target_mesh_wrapper);

    for (int tc = 0; tc < 4; ++tc) {
        std::vector<int> candidates;
        search.search(tc, &candidates);

        ASSERT_EQ(4, candidates.size());
        int tx = tc % 2; int ty = tc / 2;
        int scbase = tx + ty * 3;
        // candidates might not be in order, so sort them
        std::sort(candidates.begin(), candidates.end());
        ASSERT_EQ(scbase,     candidates[0]);
        ASSERT_EQ(scbase + 1, candidates[1]);
        ASSERT_EQ(scbase + 3, candidates[2]);
        ASSERT_EQ(scbase + 4, candidates[3]);
    }

} // TEST(search_kdtree2, ctor)

class MeshWrapperDual {
public:
    MeshWrapperDual(Portage::Jali_Mesh_Wrapper &w) : w_(w) {}
    int num_owned_cells() const { return w_.num_owned_nodes(); }
    int num_ghost_cells() const { return w_.num_ghost_nodes(); }
    void cell_get_coordinates(int const cellid,
            std::vector<std::pair<double,double> > *xylist) const {
        w_.dual_cell_get_coordinates(cellid, xylist);
    }
private:
    Portage::Jali_Mesh_Wrapper &w_;
};


TEST(search_kdtree2, dual)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    std::unique_ptr<Jali::Mesh> smesh = mf(0.0,0.0,1.0,1.0,3,3, NULL, true, true, true, true);
    std::unique_ptr<Jali::Mesh> tmesh = mf(0.0,0.0,1.0,1.0,2,2, NULL, true, true, true, true);
    Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    MeshWrapperDual s2(source_mesh_wrapper);
    MeshWrapperDual t2(target_mesh_wrapper);

    Portage::SearchKDTree2<MeshWrapperDual, MeshWrapperDual> search(s2, t2);

    for (int tc = 0; tc < 9; ++tc) {
        std::vector<int> candidates;
        search.search(tc, &candidates);

        ASSERT_EQ(4, candidates.size());
        int tx = tc % 3; int ty = tc / 3;
        int scbase = tx + ty * 4;
        // candidates might not be in order, so sort them
        std::sort(candidates.begin(), candidates.end());
        ASSERT_EQ(scbase,     candidates[0]);
        ASSERT_EQ(scbase + 1, candidates[1]);
        ASSERT_EQ(scbase + 4, candidates[2]);
        ASSERT_EQ(scbase + 5, candidates[3]);
    }

}

