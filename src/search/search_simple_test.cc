/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "gtest/gtest.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "search_simple.h"

TEST(search_simple, case1)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *smesh = mf(0.0,0.0,1.0,1.0,3,3);
    Jali::Mesh *tmesh = mf(0.0,0.0,1.0,1.0,2,2);
    Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    Portage::SearchSimple<Jali_Mesh_Wrapper,Jali_Mesh_Wrapper> search(source_mesh_wrapper, target_mesh_wrapper);

    for (int tc = 0; tc < 4; ++tc) {
        std::vector<int> candidates;
        search.search(tc, &candidates);

        ASSERT_EQ(candidates.size(), 4);
        int tx = tc % 2; int ty = tc / 2;
        int scbase = tx + ty * 3;
        ASSERT_EQ(candidates[0], scbase);
        ASSERT_EQ(candidates[1], scbase + 1);
        ASSERT_EQ(candidates[2], scbase + 3);
        ASSERT_EQ(candidates[3], scbase + 4);
    }

} // TEST(search_simple, ctor)

class MeshWrapper2 {
public:
    MeshWrapper2(Jali_Mesh_Wrapper &w) : w_(w) {}
    int num_owned_cells() const { return w_.num_owned_cells(); }
    int num_ghost_cells() const { return w_.num_ghost_cells(); }
    void cell_get_coordinates(int const cellid,
            std::vector<std::pair<double,double> > *xylist) const {
        w_.cell_get_coordinates(cellid, xylist);
    }
private:
    Jali_Mesh_Wrapper &w_;
};

class MeshWrapperDual {
public:
    MeshWrapperDual(Jali_Mesh_Wrapper &w) : w_(w) {}
    int num_owned_cells() const { return w_.num_owned_nodes(); }
    int num_ghost_cells() const { return w_.num_ghost_nodes(); }
    void cell_get_coordinates(int const cellid,
            std::vector<std::pair<double,double> > *xylist) const {
        w_.dual_cell_get_coordinates(cellid, xylist);
    }
private:
    Jali_Mesh_Wrapper &w_;
};

TEST(search_simple, wrapper2)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *smesh = mf(0.0,0.0,1.0,1.0,3,3);
    Jali::Mesh *tmesh = mf(0.0,0.0,1.0,1.0,2,2);
    Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    MeshWrapper2 s2(source_mesh_wrapper);
    MeshWrapper2 t2(target_mesh_wrapper);

    Portage::SearchSimple<MeshWrapper2, MeshWrapper2> search(s2, t2);

    for (int tc = 0; tc < 4; ++tc) {
        std::vector<int> candidates;
        search.search(tc, &candidates);

        ASSERT_EQ(candidates.size(), 4);
        int tx = tc % 2; int ty = tc / 2;
        int scbase = tx + ty * 3;
        ASSERT_EQ(candidates[0], scbase);
        ASSERT_EQ(candidates[1], scbase + 1);
        ASSERT_EQ(candidates[2], scbase + 3);
        ASSERT_EQ(candidates[3], scbase + 4);
    }

}

TEST(search_simple, dual)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *smesh = mf(0.0,0.0,1.0,1.0,3,3, NULL, true, true, true, true);
    Jali::Mesh *tmesh = mf(0.0,0.0,1.0,1.0,2,2, NULL, true, true, true, true);
    Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    MeshWrapperDual s2(source_mesh_wrapper);
    MeshWrapperDual t2(target_mesh_wrapper);

    Portage::SearchSimple<MeshWrapperDual, MeshWrapperDual> search(s2, t2);

    for (int tc = 0; tc < 9; ++tc) {
        std::vector<int> candidates;
        search.search(tc, &candidates);

        ASSERT_EQ(candidates.size(), 4);
        int tx = tc % 3; int ty = tc / 3;
        int scbase = tx + ty * 4;
        ASSERT_EQ(candidates[0], scbase);
        ASSERT_EQ(candidates[1], scbase + 1);
        ASSERT_EQ(candidates[2], scbase + 4);
        ASSERT_EQ(candidates[3], scbase + 5);
    }

}

/*-------------------------------------------------------------------------~--*
 * Formatting options for Emacs and vim.
 *
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:nil
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *-------------------------------------------------------------------------~--*/
