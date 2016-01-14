/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SEARCH_KDTREE2_H
#define SEARCH_KDTREE2_H

/*!
    \class SearchKDTree2 search_kdtree2.h
    \brief SearchKDTree2 provides...
 */

#include <vector>

#include "Point.h"
#include "BoundBox.h"
#include "kdtree.h"



namespace Portage {

template <typename SourceMeshType, typename TargetMeshType>
class SearchKDTree2 {
  public:

    //! Default constructor (disabled)
    SearchKDTree2() = delete;
    
    //! Constructor with Meshes
    /*!
      \brief Builds the kd-tree for finding intersection
      
      \param source_mesh   pointer to wrapper for getting the source mesh info
      \param target_mesh   pointer to wrapper for getting the target mesh info 
      
      Constructor for kd-tree for finding cells from a source
      mesh that overlap the target mesh

    */

    SearchKDTree2(const SourceMeshType & source_mesh, 
                 const TargetMeshType & target_mesh)
            : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {

        int numCells = sourceMesh_.num_owned_cells() + 
                sourceMesh_.num_ghost_cells();
        std::vector<gk::IsotheticBBox<2>> bboxes;
        bboxes.reserve(numCells);
        
        // find bounding boxes for all cells
        for (int c = 0; c < numCells; ++c) {
            std::vector<std::pair<double,double>> cell_coord;
            sourceMesh_.cell_get_coordinates(c, &cell_coord);
            gk::IsotheticBBox<2> bb;
            for (int n = 0; n < cell_coord.size(); ++n) {
                bb.add(gk::Point<2>(cell_coord[n].first,
                                    cell_coord[n].second));
            }
            bboxes.emplace_back(bb);
        }

        // create the kd-tree
        tree_ = gk::KDTreeCreate(bboxes);

    } // SearchKDTree2::SearchKDTree2

    //! Copy constructor (disabled)
    SearchKDTree2(const SearchKDTree2 &) = delete;

    //! Assignment operator (disabled)
    SearchKDTree2 & operator = (const SearchKDTree2 &) = delete;

    //! Destructor
    ~SearchKDTree2() { if (tree_) delete tree_; }

    /*!
      \brief returns source mesh cells potentially overlapping given target cell

      \param cellId the index of the target cell that I pass in...
      \param candidates pointer to vector of potential candidate cells in sourceMesh
    */
    void search(const int cellId, std::vector<int> *candidates) const;

  private:

    // Aggregate data members
    const SourceMeshType & sourceMesh_;
    const TargetMeshType & targetMesh_;
    gk::KDTree<2>* tree_;

}; // class SearchKDTree2



template<typename SourceMeshType, typename TargetMeshType>
void SearchKDTree2<SourceMeshType,TargetMeshType>::
search(const int cellId, std::vector<int> *candidates)
const {
    // find bounding box for target cell
    std::vector<std::pair<double,double>> cell_coord;
    targetMesh_.cell_get_coordinates(cellId, &cell_coord);
    gk::IsotheticBBox<2> bb;
    for (int n = 0; n < cell_coord.size(); ++n) {
        bb.add(gk::Point<2>(cell_coord[n].first,
                            cell_coord[n].second));
    }

    // now see which sourceMesh cells have bounding boxes overlapping
    // with target cell, using the kdtree
    std::vector<long> lcandidates;
    gk::Intersect(bb, tree_, lcandidates);
    candidates->assign(lcandidates.begin(), lcandidates.end());

} // SearchKDTree2::search




} // namespace Portage

#endif // SEARCH_KDTREE2_H

