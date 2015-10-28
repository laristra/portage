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

    auto search = Portage::SearchSimple<Jali_Mesh_Wrapper,Jali_Mesh_Wrapper>(source_mesh_wrapper, target_mesh_wrapper);

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
