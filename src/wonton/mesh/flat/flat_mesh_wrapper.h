/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#ifndef FLAT_MESH_WRAPPER_H_
#define FLAT_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <numeric>
#include <array>
#include <limits>
#include <map>
#include <set>
#include <vector>

#include "portage/wonton/mesh/AuxMeshTopology.h"
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
    nodeGlobalIds_.clear();
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
    nodeGlobalIds_.reserve(numNodes);

    for (unsigned int c=0; c<numCells; c++)
    {
      cellGlobalIds_.push_back(input.get_global_id(c, Entity_kind::CELL));

      std::vector<int> cellNodes;
      input.cell_get_nodes(c, &cellNodes);
      int cellNumNodes = cellNodes.size();
      cellNodeCounts_.push_back(cellNumNodes);
      cellToNodeList_.insert(cellToNodeList_.end(),
                             cellNodes.begin(), cellNodes.end());
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
      nodeGlobalIds_.push_back(input.get_global_id(n, Entity_kind::NODE));
    }

    // ugly hack, since dim_ is not known at compile time
    if (dim_ == 3) {
      for (unsigned int n=0; n<numNodes; ++n) {
        Portage::Point<3> nodeCoord;
        input.node_get_coordinates(n, &nodeCoord);
        for (unsigned int j=0; j<3; ++j)
          nodeCoords_.push_back(nodeCoord[j]);
      }
    }
    else if (dim_ == 2) {
      for (unsigned int n=0; n<numNodes; ++n) {
        Portage::Point<2> nodeCoord;
        input.node_get_coordinates(n, &nodeCoord);
        for (unsigned int j=0; j<2; ++j)
          nodeCoords_.push_back(nodeCoord[j]);
      }
    }

    // TODO:  Most of what's done in this call isn't needed
    // before calling distribute - can we simplify it?
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

  //! Create maps for index space conversions
  void make_index_maps() {

    // Global to local maps for cells and nodes
    std::map<int, int> globalCellMap;
    std::vector<int> cellUniqueRep;
    cellUniqueRep.reserve(cellGlobalIds_.size());
    for (unsigned int i=0; i<cellGlobalIds_.size(); ++i) {
      auto itr = globalCellMap.find(cellGlobalIds_[i]);
      if (itr == globalCellMap.end()) {
        globalCellMap[cellGlobalIds_[i]] = i;
        cellUniqueRep[i] = i;
      }
      else {
        cellUniqueRep[i] = itr->second;
      }
    }

    std::map<int, int> globalNodeMap;
    std::vector<int> nodeUniqueRep;
    nodeUniqueRep.reserve(nodeGlobalIds_.size());
    for (unsigned int i=0; i<nodeGlobalIds_.size(); ++i) {
      auto itr = globalNodeMap.find(nodeGlobalIds_[i]);
      if (itr == globalNodeMap.end()) {
        globalNodeMap[nodeGlobalIds_[i]] = i;
        nodeUniqueRep[i] = i;
      }
      else {
        nodeUniqueRep[i] = itr->second;
      }
    }

    // Compute node and face offsets
    if (dim_ == 2)
    {
      compute_offsets(cellNodeCounts_, &cellNodeOffsets_);
    }
    if (dim_ == 3)
    {
      compute_offsets(faceNodeCounts_, &faceNodeOffsets_);
      compute_offsets(cellFaceCounts_, &cellFaceOffsets_);
    }

    // Remove duplicate nodes
    if (dim_ == 2)
    {
      for (unsigned int i=0; i<cellToNodeList_.size(); ++i)
        cellToNodeList_[i] = nodeUniqueRep[cellToNodeList_[i]];
    }
    if (dim_ == 3)
    {
      for (unsigned int i=0; i<faceToNodeList_.size(); ++i)
        faceToNodeList_[i] = nodeUniqueRep[faceToNodeList_[i]];
    }

    // Compute cell-to-node adjacency lists (3D only)
    if (dim_ == 3)
    {
      cellNodeCounts_.clear();
      cellToNodeList_.clear();
      cellNodeCounts_.reserve(cellFaceCounts_.size());
      cellToNodeList_.reserve(cellFaceCounts_.size() * 4);
      for (unsigned int c=0; c<cellFaceCounts_.size(); ++c) {
        std::vector<int> cellfaces, dummy;
        cell_get_faces_and_dirs(c, &cellfaces, &dummy);
        std::set<int> cellnodes;

        for (auto const& f : cellfaces) {
          std::vector<int> facenodes;
          face_get_nodes(f, &facenodes);
          cellnodes.insert(facenodes.begin(), facenodes.end());
        }

        cellNodeCounts_.emplace_back(cellnodes.size());
        cellToNodeList_.insert(cellToNodeList_.end(),
                               cellnodes.begin(), cellnodes.end());
      }
      compute_offsets(cellNodeCounts_, &cellNodeOffsets_);
    }

    // Compute node-to-cell adjacency lists
    int numNodes = nodeCoords_.size() / dim_;
    std::vector<std::set<int>> nodeToCellTmp(numNodes);
    for (unsigned int c=0; c<cellNodeCounts_.size(); ++c) {
      int offset = cellNodeOffsets_[c];
      int count = cellNodeCounts_[c];
      for (unsigned int i=0; i<count; ++i) {
        int n = cellToNodeList_[offset+i];
        nodeToCellTmp[n].insert(cellUniqueRep[c]);
      }
    }
    nodeCellCounts_.clear();
    nodeToCellList_.clear();
    nodeCellCounts_.reserve(numNodes);
    nodeToCellList_.reserve(cellToNodeList_.size());
    for (unsigned int n=0; n<numNodes; ++n) {
      const std::set<int>& nodes = nodeToCellTmp[n];
      nodeCellCounts_.emplace_back(nodes.size());
      nodeToCellList_.insert(nodeToCellList_.end(), nodes.begin(), nodes.end());
    }
    compute_offsets(nodeCellCounts_, &nodeCellOffsets_);

    // Compute cell-to-face and face-to-node (2D only)
    // This isn't quite accurate:  some faces at rank boundaries
    // may end up being "owned" by two ranks.  But for what we are
    // doing, I don't think this will make a difference.
    if (dim_ == 2) {
      cellToFaceList_.clear();
      cellToFaceDirs_.clear();
      faceToNodeList_.clear();
      cellToFaceList_.reserve(cellNodeCounts_.size());
      cellToFaceDirs_.reserve(cellNodeCounts_.size());
      faceToNodeList_.reserve(cellToNodeList_.size()); // slight underestimate
      std::map<std::pair<int, int>, int> nodeToFaceTmp;
      int facecount = 0;
      for (unsigned int c=0; c<cellNodeCounts_.size(); ++c) {
        int offset = cellNodeOffsets_[c];
        int count = cellNodeCounts_[c];
        for (unsigned int i=0; i<count; ++i) {
          int n0 = cellToNodeList_[offset+i];
          int n1 = cellToNodeList_[offset+((i+1)%count)];
          // put nodes in canonical order
          int p0 = std::min(n0, n1);
          int p1 = std::max(n0, n1);
          auto npair = std::make_pair(p0, p1);
          auto it = nodeToFaceTmp.find(npair);
          int face;
          if (it == nodeToFaceTmp.end()) {
            // new face
            face = facecount;
            nodeToFaceTmp.insert(std::make_pair(npair, face));
            faceToNodeList_.emplace_back(p0);
            faceToNodeList_.emplace_back(p1);
            ++facecount;
          }
          else {
            // existing face
            face = it->second;
          }
          cellToFaceList_.emplace_back(face);
          cellToFaceDirs_.push_back(n0 == p0);
        }  // for i

        if (c == numOwnedCells_ - 1) numOwnedFaces_ = facecount;

      }  // for c
    }  // if dim == 2

    // Delete adjacency lists for inactive cells
    // (Note that this leaves unused spaces in the list; at some later
    // time we could compress the list, but skip that for now)
    for (unsigned int c=numOwnedCells_; c<cellNodeCounts_.size(); ++c) {
      if (cellUniqueRep[c] != c) {
        cellNodeCounts_[c] = 0;
      }
    }

  } // make_index_maps

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
    return (cellNodeCounts_.size() - numOwnedCells_);
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
    return numOwnedFaces_;
  }

  //! Number of ghost faces in the mesh
  int num_ghost_faces() const {
    if (dim_ == 2) {
      return (faceToNodeList_.size() / 2 - numOwnedFaces_);
    }
    else {
      return (faceNodeCounts_.size() - numOwnedFaces_);
    }
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

  //! Get the type of the node - PARALLEL_OWNED or PARALLEL_GHOST
  Portage::Entity_type node_get_type(int const nodeid) const {
    return (nodeid < numOwnedNodes_ ? PARALLEL_OWNED : PARALLEL_GHOST);
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
    int offset, count;
    if (dim_ == 2)
    {
      offset = cellNodeOffsets_[cellid];
      count = cellNodeCounts_[cellid];
    }
    else
    {
      offset = cellFaceOffsets_[cellid];
      count = cellFaceCounts_[cellid];
    }
    cfaces->assign(&cellToFaceList_[offset],
                   &cellToFaceList_[offset+count]);
    cfdirs->clear();
    for (unsigned int j=0; j<count; ++j)
      cfdirs->push_back(cellToFaceDirs_[offset + j] ? 1 : -1);
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
    int offset, count;
    if (dim_ == 2)
    {
      offset = 2 * faceid;
      count = 2;
    }
    else
    {
      offset = faceNodeOffsets_[faceid];
      count = faceNodeCounts_[faceid];
    }
    fnodes->assign(&faceToNodeList_[offset],
                   &faceToNodeList_[offset+count]);
  }

  //! Get list of cells for a node
  void node_get_cells(int const nodeid, Entity_type const ptype,
                      std::vector<int> *cells) const {
    int offset = nodeCellOffsets_[nodeid];
    int count = nodeCellCounts_[nodeid];
    if (ptype == Entity_type::ALL)
      cells->assign(&nodeToCellList_[offset],
                    &nodeToCellList_[offset+count]);
    else {
      cells->clear();
      cells->reserve(count);
      for (int i = 0; i < count; i++) {
        int c = nodeToCellList_[offset+i];
        if (cell_get_type(c) == ptype)
          cells->push_back(c);
      }
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

  //! get spatial dimension
  int space_dimension() const { return dim_; }

private:
  friend class MPI_Bounding_Boxes;
  std::vector<T> nodeCoords_;
  std::vector<int> cellToNodeList_;
  std::vector<int> cellNodeCounts_;
  std::vector<int> cellNodeOffsets_;
  std::vector<int> cellToFaceList_;
  std::vector<bool> cellToFaceDirs_;
  // unused in 2D (identical to cellNodeCounts_)
  std::vector<int> cellFaceCounts_;
  // unused in 2D (identical to cellNodeOffsets_)
  std::vector<int> cellFaceOffsets_;
  std::vector<int> faceToNodeList_;
  std::vector<int> faceNodeCounts_; // unused in 2D (always 2)
  std::vector<int> faceNodeOffsets_; // unused in 2D (can be computed)
  std::vector<int> nodeToCellList_;
  std::vector<int> nodeCellCounts_;
  std::vector<int> nodeCellOffsets_;
  std::vector<int> cellGlobalIds_;
  std::vector<int> nodeGlobalIds_;
  int dim_;
  int numOwnedCells_;
  int numOwnedFaces_;
  int numOwnedNodes_;

}; // class Flat_Mesh_Wrapper

//! Get radius of minimum-enclosing-sphere of a cell centered at the centroid
template<size_t D>
void cell_radius(Portage::Flat_Mesh_Wrapper<double> &wrapper,
                 int const cellid, double *radius) {
  Point<D> centroid;
  wrapper.cell_centroid<D>(cellid, &centroid);
  std::vector<int> nodes;
  wrapper.cell_get_nodes(cellid, &nodes);
  Point<D> arm, node;
  *radius = 0.;
  for (int i=0; i<nodes.size(); i++) {
    wrapper.node_get_coordinates(nodes[i], &node);
    for (int j=0; j<D; j++) arm[j] = node[j]-centroid[j];
    double distance = 0.0;
    for (int j=0; j<D; j++) distance += arm[j]*arm[j];
    distance = sqrt(distance);
    if (distance > *radius) *radius = distance;
  }
}


//! Get radius of minimum-enclosing-sphere of all nodes connected by a cell
template<size_t D>
void node_radius(Portage::Flat_Mesh_Wrapper<double> &wrapper,
                 int const nodeid, double *radius) {
  std::vector<int> nodes;
  wrapper.node_get_cell_adj_nodes(nodeid, ALL, &nodes);
  Point<D> center;
  wrapper.node_get_coordinates(nodeid, &center);
  *radius = 0.;
  Point<D> arm, node;
  for (int i=0; i<nodes.size(); i++) {
    wrapper.node_get_coordinates(nodes[i], &node);
    for (int j=0; j<D; j++) arm[j] = node[j]-center[j];
    double distance = 0.0;
    for (int j=0; j<D; j++) distance += arm[j]*arm[j];
    distance = sqrt(distance);
    if (distance > *radius) *radius = distance;
  }
}

} // end namespace Portage

#endif // FLAT_MESH_WRAPPER_H_
