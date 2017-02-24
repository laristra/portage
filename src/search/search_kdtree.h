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



#ifndef SEARCH_KDTREE_H
#define SEARCH_KDTREE_H

#include <vector>

#include "portage/support/Point.h"
#include "BoundBox.h"
#include "kdtree.h"


namespace Portage {

/*!
  @class SearchKDTree "search_kdtree.h"
  @brief A search algorithm utilizing a k-d tree.
  @tparam D The dimension of the problem space.
  @tparam SourceMeshType The mesh type of the input mesh.
  @tparam TargetMeshType The mesh type of the output mesh.
 */
template <int D, typename SourceMeshType, typename TargetMeshType>
class SearchKDTree {
  public:

    //! Default constructor (disabled)
    SearchKDTree() = delete;
    
    /*!
      @brief Builds the k-d tree for searching for intersection
      candidates.
      @param[in] source_mesh Pointer to a mesh wrapper for getting the
      source mesh info
      @param[in] target_mesh Pointer to a mesh wrapper for getting the
      target mesh info 
      
      Constructor for k-d tree for finding cells from a source
      mesh that overlap the target mesh.
    */
    SearchKDTree(const SourceMeshType & source_mesh, 
                 const TargetMeshType & target_mesh)
            : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {

        const int numCells = sourceMesh_.num_owned_cells();
        std::vector<gk::IsotheticBBox<D>> bboxes;
        bboxes.reserve(numCells);

        // find bounding boxes for all cells
        for (int c = 0; c < numCells; ++c) {
            std::vector<Point<D>> cell_coord;
            sourceMesh_.cell_get_coordinates(c, &cell_coord);
            gk::IsotheticBBox<D> bb;
            for (const auto& cc : cell_coord) {
                bb.add(cc);
            }
            bboxes.emplace_back(bb);
        }

        // create the k-d tree
        tree_ = gk::KDTreeCreate(bboxes);

    } // SearchKDTree::SearchKDTree

    //! Copy constructor (disabled)
    SearchKDTree(const SearchKDTree &) = delete;

    //! Assignment operator (disabled)
    SearchKDTree & operator = (const SearchKDTree &) = delete;

    //! Destructor
    ~SearchKDTree() { if (tree_) delete tree_; }

    /*!
      @brief Find the source mesh cells potentially overlapping a given
      target cell.
      @param[in] cellId The index of the cell in the target mesh for
      which we wish to find the candidate overlapping cells in the
      source mesh.
      @param[in,out] candidates Pointer to a vector of potential candidate
      cells in the source mesh.
    */
    void operator() (const int cellId, std::vector<int> *candidates) const;

  private:

    // Aggregate data members
    const SourceMeshType & sourceMesh_;
    const TargetMeshType & targetMesh_;
    gk::KDTree<D>* tree_;

}; // class SearchKDTree



template<int D, typename SourceMeshType, typename TargetMeshType>
void SearchKDTree<D, SourceMeshType, TargetMeshType>::
operator() (const int cellId, std::vector<int> *candidates)
const {
    // find bounding box for target cell
    std::vector<Point<D>> cell_coord;
    targetMesh_.cell_get_coordinates(cellId, &cell_coord);
    gk::IsotheticBBox<D> bb;
    for (const auto& cc : cell_coord) {
        bb.add(cc);
    }

    // now see which sourceMesh cells have bounding boxes overlapping
    // with target cell, using the kdtree
    std::vector<long> lcandidates;
    gk::Intersect(bb, tree_, lcandidates);
    candidates->assign(lcandidates.begin(), lcandidates.end());

} // SearchKDTree::operator()


} // namespace Portage

#endif // SEARCH_KDTREE_H

