/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef FLAT_MESH_WRAPPER_H_
#define FLAT_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <map>

#include "portage/wrappers/mesh/AuxMeshTopology.h"
#include "portage/support/portage.h"
#include "portage/support/Point.h"

namespace Portage {

/*!
  \class Flat_Mesh_Wrapper flat_mesh_wrapper.h
  \brief Flat_Mesh_Wrapper implements mesh methods

  Flat_Mesh_Wrapper stores mesh coordinates in a flat vector.
  It supports arbitrary polygons in 2D and arbitrary polyhedra
  in 3D.
*/

template <class T=double>
class Flat_Mesh_Wrapper : public AuxMeshTopology<Flat_Mesh_Wrapper<>> {
 public:

  //! Constructor
  Flat_Mesh_Wrapper() {};

  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Flat_Mesh_Wrapper & operator=(Flat_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Flat_Mesh_Wrapper() {};

  template<class Mesh_Wrapper>
  void initialize(Mesh_Wrapper& input)
  {
    dim_ = input.space_dimension();
    int numCells = input.num_owned_cells() + input.num_ghost_cells();
    numOwnedCells_ = input.num_owned_cells();
    int numNodes = input.num_owned_nodes() + input.num_ghost_nodes();
    numOwnedNodes_ = input.num_owned_nodes();
    int numFaces = -1;
    if (dim_ == 3)
    {
      numFaces = input.num_owned_faces() + input.num_ghost_faces();
      numOwnedFaces_ = input.num_owned_faces();
    }
    nodeCoords_.clear();
    cellNodeCounts_.clear();
    cellToNodeList_.clear();
    cellFaceCounts_.clear();
    cellToFaceList_.clear();
    faceNodeCounts_.clear();
    faceToNodeList_.clear();
    cellGlobalIds_.clear();
    nodeCoords_.reserve(numNodes*dim_);
    cellNodeCounts_.reserve(numCells);
    // reserve (dim_+1) nodes per cell (lower bound)
    cellToNodeList_.reserve(numCells*(dim_+1));
    if (dim_ == 3)
    {
      cellFaceCounts_.reserve(numCells);
      // reserve (dim_+1) faces per cell (lower bound)
      cellToFaceList_.reserve(numCells*(dim_+1));
      faceNodeCounts_.reserve(numFaces);
      // reserve dim_ nodes per face (lower bound)
      faceToNodeList_.reserve(numFaces*dim_);
    }
    cellGlobalIds_.reserve(numCells);

    for (unsigned int c=0; c<numCells; c++)
    {
      cellGlobalIds_.push_back(input.get_global_id(c, Entity_kind::CELL));

      std::vector<int> cellNodes;
      input.cell_get_nodes(c, &cellNodes);
      int cellNumNodes = cellNodes.size();
      cellNodeCounts_.push_back(cellNumNodes);
      cellToNodeList_.insert(cellToNodeList_.end(),
                            cellNodes.begin(), cellNodes.end());


      std::vector<int> cellNeighbors;
      input.cell_get_node_adj_cells(c, ALL, &(cellNeighbors));
      cellNeighborCounts_.push_back(cellNeighbors.size());
      for (unsigned int j=0; j<cellNeighbors.size(); j++)
        neighbors_.push_back(
            input.get_global_id(cellNeighbors[j], Entity_kind::CELL));
    }

    if (dim_ == 3)
    {
      for (unsigned int c=0; c<numCells; ++c)
      {
        std::vector<int> cellFaces, cfDirs;
        input.cell_get_faces_and_dirs(c, &cellFaces, &cfDirs);
        int cellNumFaces = cellFaces.size();
        cellFaceCounts_.push_back(cellNumFaces);
        cellToFaceList_.insert(cellToFaceList_.end(),
                              cellFaces.begin(), cellFaces.end());
        for (unsigned int j=0; j<cellNumFaces; ++j)
          cellToFaceDirs_.push_back(cfDirs[j] >= 0);
      }

      for (unsigned int f=0; f<numFaces; ++f)
      {
        std::vector<int> faceNodes;
        input.face_get_nodes(f, &faceNodes);
        int faceNumNodes = faceNodes.size();
        faceNodeCounts_.push_back(faceNumNodes);
        faceToNodeList_.insert(faceToNodeList_.end(),
                              faceNodes.begin(), faceNodes.end());
      }
    } // if dim_ == 3

    for (unsigned int n=0; n<numNodes; ++n) {
      // ugly hack, since dim_ is not known at compile time
      if (dim_ == 3)
      {
        Portage::Point<3> nodeCoord;
        input.node_get_coordinates(n, &nodeCoord);
        for (unsigned int j=0; j<3; ++j)
          nodeCoords_.push_back(nodeCoord[j]);
      }
      else if (dim_ == 2)
      {
        Portage::Point<2> nodeCoord;
        input.node_get_coordinates(n, &nodeCoord);
        for (unsigned int j=0; j<2; ++j)
          nodeCoords_.push_back(nodeCoord[j]);
      }
    }

    finish_init();
  }

  //! Finish mesh initialization, after initialize or MPI distribute
  void finish_init()
  {
    // Create all index maps
    make_index_maps();

    // Redo AuxMeshTopology info
    build_aux_entities();
  }

  //! Global to local index conversion
  int global_to_local(int globalId) const {
    std::map<int, int>::const_iterator it = globalCellMap_.find(globalId);
    if (it != globalCellMap_.end()) return (it->second);
    std::cout << "Global id " << globalId << " not found" << std::endl;
    return -1;
  }

  //! Create maps for index space conversions
  void make_index_maps() {

    // Global to local map
    globalCellMap_.clear();
    for (unsigned int i=0; i<cellGlobalIds_.size(); i++)
      if (globalCellMap_.find(cellGlobalIds_[i]) == globalCellMap_.end())
        globalCellMap_[cellGlobalIds_[i]] = i;

    // Node, face and neighbor offsets
    compute_offsets(cellNodeCounts_, &cellNodeOffsets_);
    if (dim_ == 3)
    {
      compute_offsets(faceNodeCounts_, &faceNodeOffsets_);
      compute_offsets(cellFaceCounts_, &cellFaceOffsets_);
    }
    compute_offsets(cellNeighborCounts_, &cellNeighborOffsets_);
  }

  //! Compute offsets from counts
  int compute_offsets(const std::vector<int>& counts,
                      std::vector<int>* offsets) {
    offsets->resize(counts.size());
    (*offsets)[0] = 0;
    std::partial_sum(counts.begin(),
                     counts.end()-1,
                     offsets->begin()+1);
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return numOwnedCells_;
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return (cellNeighborCounts_.size() - numOwnedCells_);
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return numOwnedNodes_;
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return (nodeCoords_.size()/dim_ - numOwnedNodes_);
  }

  //! Number of owned faces in the mesh
  int num_owned_faces() const {
    if (dim_ == 2)
      return num_owned_nodes();
    else
      return numOwnedFaces_;
  }

  //! Number of ghost faces in the mesh
  int num_ghost_faces() const {
    if (dim_ == 2)
      return num_ghost_nodes();
    else
      return (faceNodeCounts_.size() - numOwnedFaces_);
  }

  //! Coords of a node
  template <long D>
  void node_get_coordinates(int const nodeid, Point<D>* pp) const {
    assert(D == dim_);
    for (unsigned int j=0; j<dim_; j++)
      (*pp)[j] = nodeCoords_[nodeid*dim_+j];
  }

  //! Get the type of the cell - PARALLEL_OWNED or PARALLEL_GHOST
  Portage::Entity_type cell_get_type(int const cellid) const {
    return (cellid < numOwnedCells_ ? PARALLEL_OWNED : PARALLEL_GHOST);
  }

  //! Get the element type of a cell - TRI, QUAD, POLYGON, TET, HEX,
  //! PRISM OR POLYHEDRON
  Portage::Element_type cell_get_element_type(int const cellid) const {
    // FIXME
    return (dim_ == 2 ? POLYGON : POLYHEDRON);
  }

  //! Get cell faces and the directions in which they are used
  void cell_get_faces_and_dirs(int const cellid, std::vector<int> *cfaces,
                               std::vector<int> *cfdirs) const {
    if (dim_ == 2)
    {
      int count = cellNodeCounts_[cellid];
      int offset = cellNodeOffsets_[cellid];
      cfaces->resize(count);
      std::iota(cfaces->begin(), cfaces->end(), offset);
      cfdirs->assign(count, 1);
    }
    else
    {
      int count = cellFaceCounts_[cellid];
      int offset = cellFaceOffsets_[cellid];
      cfaces->assign(&cellToFaceList_[offset],
                     &cellToFaceList_[offset+count]);
      cfdirs->clear();
      for (unsigned int j=0; j<count; ++j)
        cfdirs->push_back(cellToFaceDirs_[offset + j] ? 1 : -1);
    }
  }

  //! Get list of nodes for a cell
  void cell_get_nodes(int const cellid, std::vector<int> *nodes) const {
    int offset = cellNodeOffsets_[cellid];
    int count = cellNodeCounts_[cellid];
    nodes->assign(&cellToNodeList_[offset],
                  &cellToNodeList_[offset+count]);
  }

  //! Get nodes of a face
  void face_get_nodes(int const faceid, std::vector<int> *fnodes) const {
    if (dim_ == 2)
    {
      auto itr = std::upper_bound(cellNodeOffsets_.begin(),
                                  cellNodeOffsets_.end(), faceid);
      int cellid = itr - cellNodeOffsets_.begin() - 1;
      int offset = cellNodeOffsets_[cellid];
      int count = cellNodeCounts_[cellid];
      int fid = faceid - offset;
      fnodes->resize(2);
      (*fnodes)[0] = cellToNodeList_[faceid];
      (*fnodes)[1] = cellToNodeList_[offset + ((fid + 1) % count)];
    }
    else
    {
      int offset = faceNodeOffsets_[faceid];
      int count = faceNodeCounts_[faceid];
      fnodes->assign(&faceToNodeList_[offset],
                     &faceToNodeList_[offset+count]);
    }
  }

  //! Coords of nodes of a cell
  template<long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Portage::Point<D>> *pplist) const {

    assert(D == dim_);
    std::vector<int> nodes;
    cell_get_nodes(cellid, &nodes);
    int cellNumNodes = nodes.size();
    pplist->resize(cellNumNodes);
    for (unsigned int i=0; i<cellNumNodes; ++i)
      node_get_coordinates(nodes[i], &((*pplist)[i]));
  }

  //! Compute the centroid of the cell
  template<long D>
  void cell_centroid(int const cellid,
                     Point<D> *centroid) const {
    assert(D == dim_);
    std::vector<Portage::Point<D>> cellCoord;
    cell_get_coordinates(cellid, &cellCoord);

    for (unsigned int i=0; i<D; i++) (*centroid)[i] = 0.0;
    for (unsigned int i=0; i<cellCoord.size(); i++)
      for (unsigned int j=0; j<D; j++)
        (*centroid)[j] = (*centroid)[j] + cellCoord[i][j];
    for (unsigned int i=0; i<D; i++) (*centroid)[i] /= cellCoord.size();

  }

  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const cellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {

    int index = cellNeighborOffsets_[cellid];
    adjcells->resize(cellNeighborCounts_[cellid]);
    for (unsigned int i=0; i<adjcells->size(); i++)
      (*adjcells)[i] = global_to_local(neighbors_[index+i]);
  }

  //! get coordinates
  std::vector<T>& get_coords() { return nodeCoords_; }

  //! get node counts
  std::vector<int>& get_node_counts() { return cellNodeCounts_; }

  //! get global cell ids
  std::vector<int>& get_global_cell_ids() { return cellGlobalIds_; }

  //! set the number of owned cells
  void set_num_owned_cells(int numOwnedCells) { numOwnedCells_ = numOwnedCells; }

  //! set the number of owned nodes
  void set_num_owned_nodes(int numOwnedNodes) { numOwnedNodes_ = numOwnedNodes; }

  //! get neighbor counts
  std::vector<int>& get_neighbor_counts() { return cellNeighborCounts_; }

  //! get neighbors
  std::vector<int>& get_neighbors() { return neighbors_; }

  //! get spatial dimension
  int space_dimension() const { return dim_; }

private:
  friend class MPI_Bounding_Boxes;
  std::vector<T> nodeCoords_;
  std::vector<int> cellToNodeList_;
  std::vector<int> cellNodeCounts_;
  std::vector<int> cellNodeOffsets_;
  // cell-to-face only used in 3D
  std::vector<int> cellToFaceList_;
  std::vector<bool> cellToFaceDirs_;
  std::vector<int> cellFaceCounts_;
  std::vector<int> cellFaceOffsets_;
  // face-to-node only used in 3D
  std::vector<int> faceToNodeList_;
  std::vector<int> faceNodeCounts_;
  std::vector<int> faceNodeOffsets_;
  std::vector<int> neighbors_; // node-connected neighbors of each cell
  std::vector<int> cellNeighborCounts_;
  std::vector<int> cellNeighborOffsets_;
  std::vector<int> cellGlobalIds_;
  // map of global to local cell IDs
  std::map<int, int> globalCellMap_;
  int dim_;
  int numOwnedCells_;
  int numOwnedFaces_; // only used in 3D
  int numOwnedNodes_;

}; // class Flat_Mesh_Wrapper


} // end namespace Portage

#endif // FLAT_MESH_WRAPPER_H_
