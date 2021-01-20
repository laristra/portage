/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_SEARCH_SEARCH_KDTREE_H_
#define PORTAGE_SEARCH_SEARCH_KDTREE_H_

#include <vector>
#include <memory>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

// portage includes
#include "portage/support/portage.h"
#include "portage/search/BoundBox.h"
#include "portage/search/kdtree.h"

namespace Portage {

/*!
  @class SearchKDTree "search_kdtree.h"
  @brief A k-d tree search class that allows us to search for control
  volumes of entities from one mesh (source) that potentially overlap
  the control volume of an entity from the second mesh (target)
  @tparam D The dimension of the problem space.
  @tparam on_what  The kind of entity we are doing a search on (NODE, CELL)
  @tparam SourceMeshType The mesh type of the source mesh.
  @tparam TargetMeshType The mesh type of the target mesh.
*/
template <int D, Entity_kind on_what,
          typename SourceMeshType, typename TargetMeshType>
class SearchKDTree {
 public:

  //! Default constructor (disabled)
  SearchKDTree() = delete;

  /*!
    @brief Builds the k-d tree for searching for intersection
    candidates.
    @param[in] source_mesh Mesh in which we search for candidates
    @param[in] target_mesh Mesh containing entity for which we search

    Constructor for k-d tree for finding entities from a source mesh
    whose control volumes overlap the control volume of an entity of
    the target mesh.
  */
  SearchKDTree(const SourceMeshType & source_mesh,
               const TargetMeshType & target_mesh)
      : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {}

  /*!
    @brief Find the source mesh entities whose control volumes
    potentially overlap control volumes of a given target entity

    @param[in] entityId The index of the cell in the target mesh for
    which we wish to find the candidate overlapping cells in the
    source mesh.
    @param[in,out] candidates Pointer to a vector of potential candidate
    cells in the source mesh.
  */
  std::vector<int> operator() (const int entityId) const {
    std::vector<int> candidates;
    throw std::runtime_error("Search not implemented for generic entity kind");
  }

 private:
  const SourceMeshType & sourceMesh_;
  const TargetMeshType & targetMesh_;
  std::shared_ptr<Portage::KDTree<D>> tree_;
};  // class SearchKDTree




//////////////////////////////////////////////////////////////////////////////
/*!
  @brief A k-d tree search class (specialization) that allows us to
  search for cells from one mesh (source) that potentially overlap a
  cell from the second mesh (target)

  @tparam D The dimension of the problem space.
  @tparam SourceMeshType The mesh type of the source mesh.
  @tparam TargetMeshType The mesh type of the target mesh.
*/
template <int D, typename SourceMeshType, typename TargetMeshType>
class SearchKDTree<D, Entity_kind::CELL, SourceMeshType, TargetMeshType> {
 public:

  //! Default constructor (disabled)
  SearchKDTree() = delete;

  /*!
    @brief Builds the k-d tree for searching for intersection
    candidates.
    @param[in] source_mesh Mesh in which we search for candidates
    @param[in] target_mesh Mesh containing entity for which we search

    Constructor for k-d tree for finding cells from a source
    mesh that overlap the target mesh.
  */
  SearchKDTree(const SourceMeshType & source_mesh,
               const TargetMeshType & target_mesh)
      : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {

    const int numCells = sourceMesh_.num_owned_cells() + sourceMesh_.num_ghost_cells();
    std::vector<Portage::IsotheticBBox<D>> bboxes;
    bboxes.reserve(numCells);

    // find bounding boxes for all cells
    for (int c = 0; c < numCells; ++c) {
      std::vector<Wonton::Point<D>> cell_coord;
      sourceMesh_.cell_get_coordinates(c, &cell_coord);
      Portage::IsotheticBBox<D> bb;
      for (const auto& cc : cell_coord) {
        bb.add(cc);
      }
      bboxes.emplace_back(bb);
    }

    // create the k-d tree
    tree_ = std::make_shared<Portage::KDTree<D>>(*Portage::KDTreeCreate(bboxes));

  }  // SearchKDTree::SearchKDTree

  /*!  @brief Find the source mesh entities whose control volumes
    potentially overlap control volumes of a given target entity
    @param[in] cellId The index of the cell in the target mesh for
    which we wish to find the candidate overlapping cells in the
    source mesh.
    @param[in,out] candidates Pointer to a vector of potential candidate
    cells in the source mesh.
  */
  std::vector<int> operator() (const int cellId) const {
    // find bounding box for target cell
    std::vector<Wonton::Point<D>> cell_coord;
    targetMesh_.cell_get_coordinates(cellId, &cell_coord);
    Portage::IsotheticBBox<D> bb;
    for (const auto& cc : cell_coord)
      bb.add(cc);

    // now see which sourceMesh cells have bounding boxes overlapping
    // with target cell, using the kdtree - since Portage::Intersect does
    // not take a shared_ptr, we have have to dereference and take
    // address of tree_
    std::vector<int> lcandidates;
    Portage::Intersect(bb, &(*tree_), lcandidates);

    std::vector<int> candidates(lcandidates.begin(), lcandidates.end());
    return candidates;
  }  // SearchKDTree::operator()

 private:
  const SourceMeshType & sourceMesh_;
  const TargetMeshType & targetMesh_;
  std::shared_ptr<Portage::KDTree<D>> tree_;
};  // class SearchKDTree (CELL specialization)




//////////////////////////////////////////////////////////////////////////////
/*!
  @brief A k-d tree search class (specialization) that allows us to
  search for nodes from one mesh (source) whose control volumes
  potentially overlap the control volumes of a node from the second
  mesh (target)

  @tparam D The dimension of the problem space.
  @tparam SourceMeshType The mesh type of the source mesh.
  @tparam TargetMeshType The mesh type of the target mesh.
*/
template <int D, typename SourceMeshType, typename TargetMeshType>
class SearchKDTree<D, Entity_kind::NODE, SourceMeshType, TargetMeshType> {
 public:

  //! Default constructor (disabled)
  SearchKDTree() = delete;

  /*!
    @brief Builds the k-d tree for searching for intersection
    candidates.
    @param[in] source_mesh Mesh in which we search for candidates
    @param[in] target_mesh Mesh containing entity for which we search

    Constructor for k-d tree for finding nodes from a source
    mesh whose control volumes overlap the control volume of a node from
    the target mesh.
  */
  SearchKDTree(const SourceMeshType & source_mesh,
               const TargetMeshType & target_mesh)
      : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {

    const int numNodes = sourceMesh_.num_owned_nodes() + sourceMesh_.num_ghost_nodes();
    std::vector<Portage::IsotheticBBox<D>> bboxes;
    bboxes.reserve(numNodes);

    // find bounding boxes for all cells
    for (int n = 0; n < numNodes; ++n) {
      std::vector<Wonton::Point<D>> dual_cell_coord;
      sourceMesh_.dual_cell_get_coordinates(n, &dual_cell_coord);
      Portage::IsotheticBBox<D> bb;
      for (const auto& cc : dual_cell_coord)
        bb.add(cc);
      bboxes.emplace_back(bb);
    }

    // create the k-d tree
    tree_ = std::make_shared<Portage::KDTree<D>>(*Portage::KDTreeCreate(bboxes));

  }  // SearchKDTree::SearchKDTree

  //! Destructor
  //  ~SearchKDTree() { if (tree_) delete tree_; }

  /*!
    @brief Find the source mesh entities whose control volumes
    potentially overlap control volumes of a given target entity

    @param[in] nodeId The index of the node in the target mesh for
    which we wish to find the candidate "overlapping" nodes in the
    source mesh.
    @param[in,out] candidates Pointer to a vector of potential candidate
    nodes in the source mesh.
  */
  std::vector<int> operator() (const int nodeId) const {
    // find bounding box for dual cell of target node
    std::vector<Wonton::Point<D>> dual_cell_coord;
    targetMesh_.dual_cell_get_coordinates(nodeId, &dual_cell_coord);
    Portage::IsotheticBBox<D> bb;
    for (const auto& cc : dual_cell_coord)
      bb.add(cc);

    // now see which sourceMesh dual cells have bounding boxes
    // overlapping with dual cell of targetMesh, using the kdtree -
    // since Portage::Intersect does not take a shared_ptr, we have have to
    // dereference and take address of tree_
    std::vector<int> lcandidates;
    Portage::Intersect(bb, &(*tree_), lcandidates);

    std::vector<int> candidates(lcandidates.begin(), lcandidates.end());
    return candidates;
  }  // SearchKDTree::operator()

 private:
  const SourceMeshType & sourceMesh_;
  const TargetMeshType & targetMesh_;
  std::shared_ptr<Portage::KDTree<D>> tree_;
};  // class SearchKDTree (NODE specialization)

}  // namespace Portage

#endif  // PORTAGE_SEARCH_SEARCH_KDTREE_H_
