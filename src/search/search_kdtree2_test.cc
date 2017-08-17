/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



#include "search_kdtree.h"

#include <algorithm>

#include "gtest/gtest.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"

TEST(search_kdtree2, case1)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    // overlay a 2x2 target mesh on a 3x3 source mesh
    // each target mesh cell gives four candidate source cells
    const std::shared_ptr<Jali::Mesh> smesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
    const std::shared_ptr<Jali::Mesh> tmesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
    const Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    const Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    Portage::SearchKDTree<2,
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
            std::vector<Portage::Point<2>> *pplist) const {
        w_.dual_cell_get_coordinates(cellid, pplist);
    }
  private:
    const Portage::Jali_Mesh_Wrapper &w_;
};


TEST(search_kdtree2, dual)
{
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    mf.included_entities({Jali::Entity_kind::EDGE,
                          Jali::Entity_kind::FACE,
                          Jali::Entity_kind::WEDGE,
                          Jali::Entity_kind::CORNER});
    // overlay a 2x2 target mesh on a 3x3 source mesh
    // each target mesh node gives four candidate source nodes
    const std::shared_ptr<Jali::Mesh> smesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
    const std::shared_ptr<Jali::Mesh> tmesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
    const Portage::Jali_Mesh_Wrapper source_mesh_wrapper(*smesh);
    const Portage::Jali_Mesh_Wrapper target_mesh_wrapper(*tmesh);

    const MeshWrapperDual s2(source_mesh_wrapper);
    const MeshWrapperDual t2(target_mesh_wrapper);

    Portage::SearchKDTree<2, MeshWrapperDual, MeshWrapperDual>
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

