/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_SEARCH_SEARCH_SWEPT_FACE_H_
#define PORTAGE_SEARCH_SEARCH_SWEPT_FACE_H_

#include <vector>
#include <algorithm>

// portage includes
#include "portage/support/portage.h"

namespace Portage {

/*!
  @class SearchSweptFace "search_swept_face.h"
  @brief The search algorithm that is used for the swept face remap.
  @tparam D  The dimension of the problem space.
  @tparam on_what  The kind of entity we are doing a search on (NODE, CELL)    
  @tparam SourceMeshType  The mesh type of the input mesh.
  @tparam TargetMeshType  The mesh type of the output mesh.
*/
template <int D, Entity_kind on_what,
          typename SourceMeshType, typename TargetMeshType>
class SearchSweptFace {
public:
  //! Default constructor (disabled)
  SearchSweptFace() = delete;

  // Constructor with Meshes
  /*!
    @brief Builds the search structure for the swept face algorithm.
    @param[in] source_mesh_wrapper Pointer to a mesh wrapper for
    getting the source mesh info.
    @param[in] target_mesh_wrapper Pointer to a mesh wrapper for
    getting the target mesh info.

    Constructor for the search structure that finds mesh entities that
    are neighbors through the faces. Topology and indexing of the source 
    and target meshes are assumed to be identical.
  */
  SearchSweptFace(const SourceMeshType & source_mesh,
                  const TargetMeshType & target_mesh)
    : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {
#ifdef DEBUG
    int numSourceCells = sourceMesh_.num_owned_cells() +
      sourceMesh_.num_ghost_cells();
    int numTargetCells = targetMesh_.num_owned_cells() +
      targetMesh_.num_ghost_cells();
    assert(numSourceCells == numTargetCells);
#endif    
  }  // SearchSweptFace::SearchSweptFace

  //! Copy constructor (disabled)
  SearchSweptFace(const SearchSweptFace&) = delete;

  //! Assignment operator (disabled)
  SearchSweptFace & operator = (const SearchSweptFace &) = delete;

  //! Destructor
  ~SearchSweptFace() = default;

  /*!
    @brief Find the source mesh entities that are adjacent to the given target
    mesh entity through the faces. Topology and indexing of the source and target
    meshes are assumed to be identical.
    @param[in] entityId The index of the entity in the target mesh for
    which we wish to find the candidate overlapping entities in the
    source mesh.
    @return Vector of potential candidate entities in the source mesh.
  */
  std::vector<int> operator() (const int entityId) const {
    std::vector<int> candidates;
    std::cerr << "Swept face search is not implemented for a generic entity kind" << std::endl;
    return candidates;
  }

private:
  const SourceMeshType &sourceMesh_;
  const TargetMeshType &targetMesh_;
};  // class SearchSweptFace

/*!
  @brief Specialization of the search algorithm that is used for the swept face
  remap of fields associatred with nodal control volumes.
  @tparam D  The dimension of the problem space.
  @tparam SourceMeshType  The mesh type of the input mesh.
  @tparam TargetMeshType  The mesh type of the output mesh.
*/
template <int D, typename SourceMeshType, typename TargetMeshType>
class SearchSweptFace<D, Entity_kind::NODE, SourceMeshType, TargetMeshType> {
public:
  //! Default constructor (disabled)
  SearchSweptFace() = delete;

  // Constructor with Meshes
  /*!
    @brief Builds the search structure for the swept face algorithm.
    @param[in] source_mesh_wrapper Pointer to a mesh wrapper for
    getting the source mesh info.
    @param[in] target_mesh_wrapper Pointer to a mesh wrapper for
    getting the target mesh info.

    Constructor for the search structure that finds source mesh nodes
    for which their control volumes are adjacent to the control volume
    of the target node. Topology and indexing of the source and
    target meshes are assumed to be identical.
  */
  SearchSweptFace(const SourceMeshType & source_mesh,
                  const TargetMeshType & target_mesh)
    : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {
#ifdef DEBUG
    int numSourceCells = sourceMesh_.num_owned_cells() +
      sourceMesh_.num_ghost_cells();
    int numTargetCells = targetMesh_.num_owned_cells() +
      targetMesh_.num_ghost_cells();
    assert(numSourceCells == numTargetCells);
#endif    
  }  // SearchSweptFace::SearchSweptFace

  //! Assignment operator (disabled)
  SearchSweptFace & operator = (const SearchSweptFace &) = delete;

  //! Destructor
  ~SearchSweptFace() = default;

  /*!
    @brief Find the source mesh nodes for which their control volumes
    are adjacent to the control volume of the given target node.
    Topology and indexing of the source and target meshes are
    assumed to be identical.
    @param[in] entityId The index of the entity in the target mesh for
    which we wish to find the candidate overlapping entities in the
    source mesh.
    @return Vector of potential candidate entities in the source mesh.
  */
  std::vector<int> operator() (const int entityId) const {
    std::vector<int> candidates;
    std::cerr << "Swept face search is not implemented for nodal control volumes" << std::endl;
    return candidates;
  }

private:
  const SourceMeshType &sourceMesh_;
  const TargetMeshType &targetMesh_;
};  // class SearchSweptFace

/*!
  @brief Specialization of the search algorithm that is used for the swept face
  remap of fields associated with cells.
  @tparam D  The dimension of the problem space.
  @tparam SourceMeshType  The mesh type of the input mesh.
  @tparam TargetMeshType  The mesh type of the output mesh.
*/
template <int D, typename SourceMeshType, typename TargetMeshType>
class SearchSweptFace<D, Entity_kind::CELL, SourceMeshType, TargetMeshType> {
public:
  //! Default constructor (disabled)
  SearchSweptFace() = delete;

  // Constructor with Meshes
  /*!
    @brief Builds the search structure for the swept face algorithm.
    @param[in] source_mesh_wrapper Pointer to a mesh wrapper for
    getting the source mesh info.
    @param[in] target_mesh_wrapper Pointer to a mesh wrapper for
    getting the target mesh info.

    Constructor for the search structure that finds cells that are
    neighbors through the faces. Topology and indexing of the source 
    and target meshes are assumed to be identical.
  */
  SearchSweptFace(const SourceMeshType & source_mesh,
                  const TargetMeshType & target_mesh)
    : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {
#ifdef DEBUG
    int numSourceCells = sourceMesh_.num_owned_cells() +
      sourceMesh_.num_ghost_cells();
    int numTargetCells = targetMesh_.num_owned_cells() +
      targetMesh_.num_ghost_cells();
    assert(numSourceCells == numTargetCells);
#endif    
  }  // SearchSweptFace::SearchSweptFace

  //! Assignment operator (disabled)
  SearchSweptFace & operator = (const SearchSweptFace &) = delete;

  //! Destructor
  ~SearchSweptFace() = default;

  /*!
    @brief Find the source mesh cells that are adjacent to the given target cell
    through the faces. Topology and indexing of the source and target meshes are
    assumed to be identical.
    @param[in] entityId The index of the entity in the target mesh for
    which we wish to find the candidate overlapping entities in the
    source mesh.
    @return Vector of potential candidate entities in the source mesh.
  */
  std::vector<int> operator() (const int entityId) const {
    std::vector<int> candidates;
    sourceMesh_.cell_get_face_adj_cells(entityId, Entity_type::ALL, &candidates);
#ifdef DEBUG
    std::vector<int> alt_candidates;
    targetMesh_.cell_get_face_adj_cells(entityId, Entity_type::ALL, &alt_candidates);
    assert(candidates.size() == alt_candidates.size());
    for (int icc = 0; icc < candidates.size(); icc++)
      assert(candidates[icc] == alt_candidates[icc]);
#endif  

    return candidates;
  }

private:
  const SourceMeshType &sourceMesh_;
  const TargetMeshType &targetMesh_;
};  // class SearchSweptFace

}  // namespace Portage

#endif  // PORTAGE_SEARCH_SEARCH_SWEPT_FACE_H_
