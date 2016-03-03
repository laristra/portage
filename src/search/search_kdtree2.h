/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SEARCH_KDTREE2_H
#define SEARCH_KDTREE2_H

#include <vector>

#include "portage/support/Point.h"
#include "BoundBox.h"
#include "kdtree.h"


namespace Portage {

    /*!
      @class SearchKDTree2 "search_kdtree2.h"
      @brief A search algorithm utilizing a k-d tree in 2d.
      @tparam SourceMeshType The mesh type of the input mesh.
      @tparam TargetMeshType The mesh type of the output mesh.

      This search is only valid for 2d meshes.
      */
    template <typename SourceMeshType, typename TargetMeshType>
        class SearchKDTree2 {
        public:

           //! Default constructor (disabled)
           SearchKDTree2() = delete;

           /*!
             @brief Builds the k-d tree for searching for intersection candidates.
             @param[in] source_mesh Pointer to a mesh wrapper for getting the source
             mesh info.
             @param[in] target_mesh Pointer to a mesh wrapper for getting the target
             mesh info.

             Constructor for k-d tree for finding cells from a source
             mesh that overlap the target mesh.
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

                   // create the k-d tree
                   tree_ = gk::KDTreeCreate(bboxes);

               } // SearchKDTree2::SearchKDTree2

           //! Copy constructor (disabled)
           SearchKDTree2(const SearchKDTree2 &) = delete;

           //! Assignment operator (disabled)
           SearchKDTree2 & operator = (const SearchKDTree2 &) = delete;

           //! Destructor
           ~SearchKDTree2() { if (tree_) delete tree_; }

           /*!
             @brief Find the source mesh cells potentially overlapping a given target cell
             @param[in] cellId The index of the cell in the target mesh for which we wish
             to find the candidate overlapping cells in the source mesh.
             @param[in,out] candidates Pointer to a vector of potential candidate cells in
             the source mesh.
             */
           void operator() (const int cellId, std::vector<int> *candidates) const;

       private:

           // Aggregate data members
           const SourceMeshType & sourceMesh_;
           const TargetMeshType & targetMesh_;
           gk::KDTree<2>* tree_;

        }; // class SearchKDTree2



    template<typename SourceMeshType, typename TargetMeshType>
        void SearchKDTree2<SourceMeshType,TargetMeshType>::
        operator() (const int cellId, std::vector<int> *candidates)
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

        } // SearchKDTree2::operator()




} // namespace Portage

#endif // SEARCH_KDTREE2_H

