/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
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

