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
  In 2D, it supports arbitrary polygons.  In 3D, it currently
  only handles the cases of axis-aligned hexahedra and general
  tetrahedra when decomposing into tets.
*/

template <class T=double>
class Flat_Mesh_Wrapper : public AuxMeshTopology<Flat_Mesh_Wrapper<>> {
 public:

  //! Constructor
  Flat_Mesh_Wrapper() {};

  template<class Mesh_Wrapper>
  void initialize(int nodes_per_cell, Mesh_Wrapper& input)
  {
    nodesPerCell_ = nodes_per_cell;
    dim_ = input.space_dimension();
    int numCells = input.num_owned_cells() + input.num_ghost_cells();
    numOwnedCells_ = input.num_owned_cells();
    coords_.clear();
    nodeCounts_.clear();
    // reserve (dim_+1) nodes per cell (lower bound)
    coords_.reserve(numCells*(dim_+1)*dim_);
    nodeCounts_.reserve(numCells);
      
    for (unsigned int c=0; c<numCells; c++)
    {
      globalCellIds_.push_back(input.get_global_id(c, Entity_kind::CELL));
      int cellNumNodes;
      // ugly hack, since dim_ is not known at compile time
      if (dim_ == 3)
      {
        std::vector<Portage::Point<3>> cellCoord;
        input.cell_get_coordinates(c, &cellCoord);
        cellNumNodes = cellCoord.size();
        for (unsigned int i=0; i<cellNumNodes; i++)
          for (unsigned int j=0; j<dim_; j++)
            coords_.push_back(cellCoord[i][j]);
      }
      else if (dim_ == 2)
      {
        std::vector<Portage::Point<2>> cellCoord;
        input.cell_get_coordinates(c, &cellCoord);
        cellNumNodes = cellCoord.size();
        for (unsigned int i=0; i<cellNumNodes; i++)
          for (unsigned int j=0; j<dim_; j++)
            coords_.push_back(cellCoord[i][j]);
      }
      nodeCounts_.push_back(cellNumNodes);

      std::vector<int> cellNeighbors;
      input.cell_get_node_adj_cells(c, ALL, &(cellNeighbors));
      neighborCounts_.push_back(cellNeighbors.size());
      for (unsigned int j=0; j<cellNeighbors.size(); j++)
        neighbors_.push_back(input.get_global_id(cellNeighbors[j], Entity_kind::CELL));
    }
    numOwnedNodes_ = std::accumulate(
        &nodeCounts_[0], &nodeCounts_[numOwnedCells_], 0);

    make_index_maps();
  }

  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Flat_Mesh_Wrapper & operator=(Flat_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Flat_Mesh_Wrapper() {};

  //! Cell area or volume
  double cell_volume(int cellID) const {

    // TODO:  Generalize from regular grid to arbitrary polyhedron
    std::vector<T> extrema(2*dim_); 
    for (unsigned int i=0; i<2*dim_; i+=2) extrema[i] = std::numeric_limits<T>::max();
    for (unsigned int i=1; i<2*dim_; i+=2) extrema[i] = -std::numeric_limits<T>::max();
    int count = nodeCounts_[cellID];
    int offset = nodeOffsets_[cellID];
    for (unsigned int i=0; i<count; i++)
    {
      for (unsigned int j=0; j<dim_; j++)
      {
        T v = coords_[offset*dim_+i*dim_+j];
        if (v < extrema[2*j]) extrema[2*j] = v;
        if (v > extrema[2*j+1]) extrema[2*j+1] = v;
      }
    }
    double volume = 1.0f;
    for (unsigned int i=0; i<dim_; i++) volume *= (extrema[2*i+1] - extrema[2*i]);

    return volume;
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
    for (unsigned int i=0; i<globalCellIds_.size(); i++)
      if (globalCellMap_.find(globalCellIds_[i]) == globalCellMap_.end())
        globalCellMap_[globalCellIds_[i]] = i; 

    // Node offsets
    nodeOffsets_.clear();
    nodeOffsets_.resize(nodeCounts_.size());
    nodeOffsets_[0] = 0;
    std::partial_sum(nodeCounts_.begin(), nodeCounts_.end()-1,
                     nodeOffsets_.begin()+1);

    // Neighbor offsets
    neighborOffsets_.clear();
    neighborOffsets_.resize(neighborCounts_.size());
    neighborOffsets_[0] = 0;
    std::partial_sum(neighborCounts_.begin(), neighborCounts_.end()-1,
                     neighborOffsets_.begin()+1);
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return num_entities(Entity_kind::CELL, Entity_type::PARALLEL_OWNED);
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return num_entities(Entity_kind::CELL, Entity_type::PARALLEL_GHOST);
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return num_entities(Entity_kind::NODE, Entity_type::PARALLEL_OWNED);
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return num_entities(Entity_kind::NODE, Entity_type::PARALLEL_GHOST);
  }

  //! Number of owned faces in the mesh
  int num_owned_faces() const {
    return num_entities(Entity_kind::FACE, Entity_type::PARALLEL_OWNED);
  }

  //! Number of ghost faces in the mesh
  int num_ghost_faces() const {
    return num_entities(Entity_kind::FACE, Entity_type::PARALLEL_GHOST);
  }

  //! coords of nodes of a cell
  template<long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Portage::Point<D>> *pplist) const {

    int count = nodeCounts_[cellid];
    int offset = nodeOffsets_[cellid];
    pplist->resize(count);
    for (unsigned int i=0; i<count; i++)
      for (unsigned int j=0; j<dim_; j++)
        (*pplist)[i][j] = coords_[offset*dim_+i*dim_+j];
  }

  //! Compute the centroid of the cell
  template <long D>
  void cell_centroid(int const cellid,
                     Point<D> *centroid) const {
 
    std::vector<Portage::Point<D> > cellCoord;
    cell_get_coordinates(cellid, &cellCoord);
 
    for (unsigned int i=0; i<D; i++) (*centroid)[i] = 0.0;
    for (unsigned int i=0; i<cellCoord.size(); i++)
      for (unsigned int j=0; j<D; j++)
        (*centroid)[j] = (*centroid)[j] + cellCoord[i][j];
    for (unsigned int i=0; i<D; i++) (*centroid)[i] /= cellCoord.size();

  }

  //! Get the simplest possible decomposition of a 3D cell into tets.
  //! This currently only handles the cases of planar hexahedra or tetrahedra
  //! cells
  void decompose_cell_into_tets(const int cellID,
      std::vector<std::array<Portage::Point<3>, 4>> *tcoords,
      const bool planar_hex) const {

    int count = nodeCounts_[cellID];
    int offset = nodeOffsets_[cellID];
    if (/*planar_hex &&*/ (count == 8) && (dim_ == 3))
    {
      std::vector<Portage::Point<3>> vertices(count);
      std::array<T, 6> extrema;
      extrema[0] = extrema[2] = extrema[4] = std::numeric_limits<T>::max();
      extrema[1] = extrema[3] = extrema[5] = -std::numeric_limits<T>::max();

      for (unsigned int i=0; i<count; i++)
      {
        for (unsigned int j=0; j<dim_; j++) 
        {
          vertices[i][j] = coords_[offset*dim_+i*dim_+j];
          if (vertices[i][j] < extrema[2*j]) extrema[2*j] = vertices[i][j];
          if (vertices[i][j] > extrema[2*j+1]) extrema[2*j+1] = vertices[i][j];
        }
      }
      Portage::Point<3> center;
      for (unsigned int i=0; i<dim_; i++) center[i] = extrema[2*i] + (extrema[2*i+1] - extrema[2*i])/2.0;
      vertices.push_back(center);

      int indexes[12][4] = { {0,2,8,1}, {0,2,3,8}, {1,6,8,5}, {1,6,2,8}, {4,3,8,0}, {4,3,7,8}, 
                             {3,6,8,2}, {3,6,7,8}, {4,1,8,5}, {4,1,0,8}, {5,7,8,4}, {5,7,6,8} }; 
      for (unsigned int t=0; t<12; t++)
      {
        std::array<Portage::Point<3>, 4> tmp;
        for (int i=0; i<4; i++)
          tmp[i] = vertices[indexes[t][i]];
        tcoords->push_back(tmp);
      }
    }

    if ((count == 4) && (dim_ == 3))
    {
      std::array<Portage::Point<3>, 4> tmp;
      for (int i=0; i<4; i++)
        for (int j=0; j<3; j++)
          tmp[i][j] = coords_[offset*dim_+i*dim_+j];
      tcoords->push_back(tmp);
    }
  }

  //! Number of items of given entity
  int num_entities(Entity_kind const entity, Entity_type const etype=Entity_type::ALL) const {
    if (entity == Entity_kind::CELL)
    {
      if (etype == Entity_type::PARALLEL_OWNED)  return numOwnedCells_;
      else if (etype == Entity_type::PARALLEL_GHOST) return (neighborCounts_.size() - numOwnedCells_);
      else if (etype == Entity_type::ALL) return (neighborCounts_.size());
    }
    else if (entity == Entity_kind::NODE) 
    {
      if (etype == Entity_type::PARALLEL_OWNED)  return numOwnedNodes_;
      else if (etype == Entity_type::PARALLEL_GHOST)
        return (coords_.size()/dim_ - numOwnedNodes_);
      else if (etype == Entity_type::ALL) return coords_.size()/dim_;
    }
    else return 0;
  }

  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const cellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {
    
    int index = neighborOffsets_[cellid];
    adjcells->resize(neighborCounts_[cellid]);
    for (unsigned int i=0; i<adjcells->size(); i++)
      (*adjcells)[i] = global_to_local(neighbors_[index+i]);
  }

  //! get coordinates
  std::vector<T>& get_coords() { return coords_; }

  //! get node counts
  std::vector<int>& get_node_counts() { return nodeCounts_; }

  //! get global cell ids
  std::vector<int>& get_global_cell_ids() { return globalCellIds_; }

  //! set the number of owned cells
  void set_num_owned_cells(int numOwnedCells) { numOwnedCells_ = numOwnedCells; }

  //! set the number of owned nodes
  void set_num_owned_nodes(int numOwnedNodes) { numOwnedNodes_ = numOwnedNodes; }

  //! get neighbor counts
  std::vector<int>& get_neighbor_counts() { return neighborCounts_; }

  //! get neighbors
  std::vector<int>& get_neighbors() { return neighbors_; }

  //! get spatial dimension
  int space_dimension() const { return dim_; }

private:
  std::vector<T> coords_;
  std::vector<int> nodeCounts_;
  std::vector<int> nodeOffsets_;
  std::vector<int> neighbors_;
  std::vector<int> neighborCounts_;
  std::vector<int> neighborOffsets_;
  std::vector<int> globalCellIds_;
  std::map<int, int> globalCellMap_;
  int nodesPerCell_;
  int dim_;
  int numOwnedCells_;
  int numOwnedNodes_;

}; // class Flat_Mesh_Wrapper


} // end namespace Portage

#endif // FLAT_HEX_MESH_WRAPPER_H_
