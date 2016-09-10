/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef FLAT_MESH_WRAPPER_H_
#define FLAT_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <array>
#include <limits>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

namespace Portage {

/*!
  \class Flat_Mesh_Wrapper flat_hex_mesh_wrapper.h
  \brief Flat_Mesh_Wrapper implements mesh methods

  Flat_Mesh_Wrapper stores mesh coordinates in a flat vector.
  It currently assumes all cells are of the same type (i.e.,
  same number of nodes per cell). It also currently only
  handles the cases of hexahedra and tetrahedra when
  decomposing into tets.
*/

template <class T=double>
class Flat_Mesh_Wrapper {
 public:

  //! Constructor
  template<class Mesh_Wrapper>
  Flat_Mesh_Wrapper(int nodes_per_cell, Mesh_Wrapper& input) :
                    nodesPerCell_(nodes_per_cell), dim_(input.space_dimension())
  {
    int numCells = input.num_owned_cells() + input.num_ghost_cells();
    coords_.resize(numCells*nodesPerCell_*dim_);
      
    for (unsigned int c=0; c<numCells; c++)
    {
      if (c < input.num_owned_cells()) ownedCellIndexes_.push_back(c);
      globalCellIds_.push_back(input.get_global_id(c, Entity_kind::CELL));
      virtualCellIds_.push_back(c);
      std::vector<Portage::Point<3>> cellCoord;
      input.cell_get_coordinates(c, &cellCoord);
      for (unsigned int j=0; j<nodesPerCell_; j++)
      {
        coords_[c*nodesPerCell_*dim_+j*dim_+0] = cellCoord[j][0];
        coords_[c*nodesPerCell_*dim_+j*dim_+1] = cellCoord[j][1];
        coords_[c*nodesPerCell_*dim_+j*dim_+2] = cellCoord[j][2];
      }

      std::vector<int> cellNeighbors;
      input.cell_get_node_adj_cells(c, ALL, &(cellNeighbors));
      neighborCounts_.push_back(cellNeighbors.size());
      for (unsigned int j=0; j<cellNeighbors.size(); j++)
        neighbors_.push_back(input.get_global_id(cellNeighbors[j], Entity_kind::CELL));
    }
  }

  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Flat_Mesh_Wrapper & operator=(Flat_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Flat_Mesh_Wrapper() {};

  //! Cell area/volume
  double cell_volume(int cellID) const {

    cellID = virtual_to_local(cellID);

    std::vector<T> extrema(2*dim_); 
    for (unsigned int i=0; i<2*dim_; i+=2) extrema[i] = std::numeric_limits<T>::max();
    for (unsigned int i=1; i<2*dim_; i+=2) extrema[i] = -std::numeric_limits<T>::max();
    for (unsigned int i=0; i<nodesPerCell_; i++)
    {
      for (unsigned int j=0; j<dim_; j++)
      {
        T v = coords_[cellID*nodesPerCell_*dim_+i*dim_+j];
        if (v < extrema[2*j]) extrema[2*j] = v;
        if (v > extrema[2*j+1]) extrema[2*j+1] = v;
      }
    }
    double volume = 1.0f;
    for (unsigned int i=0; i<dim_; i++) volume *= (extrema[2*i+1] - extrema[2*i]);

    return volume;
  }

  //! Virtual to local
  int virtual_to_local(int virtualId) const {
    for (unsigned int i=0; i<virtualCellIds_.size(); i++)
      if (virtualCellIds_[i] == virtualId)
        return i;
    std::cout << "Virtual id " << virtualId << " not found" << std::endl;
    return 0;
  }

  //! Global to virtual
  int global_to_virtual(int globalId) const {
    for (unsigned int i=0; i<globalCellIds_.size(); i++)
      if (globalCellIds_[i] == globalId)
        return virtualCellIds_[i];
    std::cout << "Global id " << globalId << " not found" << std::endl;
    return 0;
  }

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return num_entities(Entity_kind::CELL, Entity_type::PARALLEL_OWNED);
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return 0;
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return num_entities(Entity_kind::CELL, Entity_type::PARALLEL_GHOST);
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return 0;
  }

  //! coords of nodes of a cell
  template<long D>
  void cell_get_coordinates(int const vcellid,
                            std::vector<Portage::Point<D>> *pplist) const {

    int cellid = virtual_to_local(vcellid);

    pplist->resize(nodesPerCell_);
    for (unsigned int i=0; i<nodesPerCell_; i++)
      for (unsigned int j=0; j<dim_; j++)
        (*pplist)[i][j] = coords_[cellid*nodesPerCell_*dim_+i*dim_+j];
  }

  //! Compute the centroid of the cell
  template <long D>
  void cell_centroid(int const cellid,
                     Point<D> *centroid) const {
    int dim = space_dimension();
    double boundingBoxes[2*dim];
    for (unsigned int i=0; i<2*dim; i+=2)
    {
      boundingBoxes[i+0] = std::numeric_limits<double>::max();
      boundingBoxes[i+1] = -std::numeric_limits<double>::max();
    }
    std::vector<Portage::Point<D>> cellCoord;
    cell_get_coordinates(cellid, &cellCoord);
    for (unsigned int j=0; j<cellCoord.size(); j++)
    {
      for (unsigned int k=0; k<dim; k++)
        if (cellCoord[j][k] < boundingBoxes[2*k])
          boundingBoxes[2*k] = cellCoord[j][k];
      for (unsigned int k=0; k<dim; k++)
        if (cellCoord[j][k] > boundingBoxes[2*k+1])
          boundingBoxes[2*k+1] = cellCoord[j][k];
    }

    for (int i = 0; i < dim; ++i)
      (*centroid)[i] = boundingBoxes[2*i] + (boundingBoxes[2*i+1] - boundingBoxes[2*i])/2.0f; 
  }

  //! Get the simplest possible decomposition of a 3D cell into tets.
  //! This currently only handles the cases of planar hexahedra or tetrahedra
  //! cells
  void decompose_cell_into_tets(const int vcellID,
      std::vector<std::array<Portage::Point<3>, 4>> *tcoords,
      const bool planar_hex) const {
    
    int cellID = virtual_to_local(vcellID);

    if (/*planar_hex &&*/ (nodesPerCell_ == 8) && (dim_ == 3))
    {
      std::vector<Portage::Point<3>> vertices(nodesPerCell_);
      std::array<T, 6> extrema;
      extrema[0] = extrema[2] = extrema[4] = std::numeric_limits<T>::max();
      extrema[1] = extrema[3] = extrema[5] = -std::numeric_limits<T>::max();

      for (unsigned int i=0; i<nodesPerCell_; i++)
      {
        for (unsigned int j=0; j<dim_; j++) 
        {
          vertices[i][j] = coords_[cellID*nodesPerCell_*dim_+i*dim_+j];
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

    if ((nodesPerCell_ == 4) && (dim_ == 3))
    {
      std::array<Portage::Point<3>, 4> tmp;
      for (int i=0; i<4; i++)
        for (int j=0; j<3; j++)
          tmp[i][j] = coords_[cellID*nodesPerCell_*dim_+i*dim_+j];
      tcoords->push_back(tmp);
    }
  }

  //! Number of items of given entity
  int num_entities(Entity_kind const entity, Entity_type const etype=Entity_type::ALL) const {
    if (entity == Entity_kind::CELL)
    {
      if (etype == Entity_type::PARALLEL_OWNED)  return ownedCellIndexes_.size();
      else if (etype == Entity_type::PARALLEL_GHOST) return (neighborCounts_.size() - ownedCellIndexes_.size());
      else if (etype == Entity_type::ALL) return (neighborCounts_.size());
    }
    else if (entity == Entity_kind::NODE) 
    {
      return coords_.size();
    }
    else return 0;
  }

  //! Get node connected neighbors of cell
  void cell_get_node_adj_cells(int const vcellid,
                               Entity_type const ptype,
                               std::vector<int> *adjcells) const {
    // TODO: Compute the scan of the neighborCounts_ vector once so that we don't redo this
    // loop everytime we access neighbor information
    
    int cellid = virtual_to_local(vcellid);

    int index = 0;
    for (unsigned int i=0; i<cellid; i++)
      index += neighborCounts_[i];
    adjcells->resize(neighborCounts_[cellid]);
    for (unsigned int i=0; i<adjcells->size(); i++)
      (*adjcells)[i] = global_to_virtual(neighbors_[index+i]);
  }

  //! get coordinates
  std::vector<T>& get_coords() { return coords_; }

  //! get global cell ids
  std::vector<int>& get_global_cell_ids() { return globalCellIds_; }

  //! get owned indexes
  std::vector<int>& get_owned_cell_indexes() { return ownedCellIndexes_; }

  //! get virtual cell ids
  std::vector<int>& get_virtual_cell_ids() { return virtualCellIds_; }

  //! get neighbor counts
  std::vector<int>& get_neighbor_counts() { return neighborCounts_; }

  //! get neighbors
  std::vector<int>& get_neighbors() { return neighbors_; }

  //! get nodes per cell
  int get_nodes_per_cell() { return nodesPerCell_; }

  //! get spatial dimension
  int space_dimension() const { return dim_; }

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity) const {
    int start_index = 0;
    return make_counting_iterator(start_index);
  }

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity, Entity_type const etype=Entity_type::ALL) const {
    int start_index = 0;
    return (make_counting_iterator(start_index) + num_entities(entity, etype));
  }


private:
  std::vector<T> coords_;
  std::vector<int> neighbors_;
  std::vector<int> neighborCounts_;
  std::vector<int> globalCellIds_;
  std::vector<int> virtualCellIds_;
  std::vector<int> ownedCellIndexes_;
  const int nodesPerCell_;
  const int dim_;

}; // class Flat_Mesh_Wrapper


} // end namespace Portage

#endif // FLAT_HEX_MESH_WRAPPER_H_
