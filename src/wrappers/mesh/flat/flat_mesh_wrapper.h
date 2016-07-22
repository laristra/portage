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

namespace {  // unnamned
//! helper function:  convert Jali point to Portage point

}  // namespace unnamed

template <class T=double>
class Flat_Mesh_Wrapper {
 public:

  //! Constructor
  template<class Mesh_Wrapper>
  Flat_Mesh_Wrapper(int nodes_per_cell, Mesh_Wrapper& input) :
                    nodesPerCell_(nodes_per_cell), dim_(input.space_dimension())
  {
    int numCells = input.num_owned_cells();
    coords_.resize(numCells*nodesPerCell_*dim_);
      
    for (unsigned int c=0; c<numCells; c++)
    {
      std::vector<Portage::Point<3>> cellCoord;
      input.cell_get_coordinates(c, &cellCoord);
      for (unsigned int j=0; j<nodesPerCell_; j++)
      {
        coords_[c*nodesPerCell_*dim_+j*dim_+0] = cellCoord[j][0];
        coords_[c*nodesPerCell_*dim_+j*dim_+1] = cellCoord[j][1];
        coords_[c*nodesPerCell_*dim_+j*dim_+2] = cellCoord[j][2];
      }
    }
  }

  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Flat_Mesh_Wrapper & operator=(Flat_Mesh_Wrapper const &) = delete;

  //! Empty destructor
  ~Flat_Mesh_Wrapper() {};

  //! Cell area/volume
  double cell_volume(int cellID) const {

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

  //! Number of owned cells in the mesh
  int num_owned_cells() const {
    return (coords_.size() / (nodesPerCell_*dim_));
  }

  //! Number of owned nodes in the mesh
  int num_owned_nodes() const {
    return (coords_.size() / dim_);
  }

  //! Number of ghost cells in the mesh
  int num_ghost_cells() const {
    return 0;
  }

  //! Number of ghost nodes in the mesh
  int num_ghost_nodes() const {
    return 0;
  }

  //! coords of nodes of a cell
  template<long D>
  void cell_get_coordinates(int const cellid,
                            std::vector<Portage::Point<D>> *pplist) const {
    pplist->resize(nodesPerCell_);
    for (unsigned int i=0; i<nodesPerCell_; i++)
      for (unsigned int j=0; j<dim_; j++)
        (*pplist)[i][j] = coords_[cellid*nodesPerCell_*dim_+i*dim_+j];
  }

  //! Get the simplest possible decomposition of a 3D cell into tets.
  //! This currently only handles the cases of hexahedra or tetrahedra cells
  void decompose_cell_into_tets(const int cellID,
      std::vector<std::array<Portage::Point<3>, 4>> *tcoords) const {
    
    if ((nodesPerCell_ == 8) && (dim_ == 3))
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

  //! get coordinates
  std::vector<T>& get_coords() { return coords_; }

  //! get nodes per cell
  int get_nodes_per_cell() { return nodesPerCell_; }

  //! get spatial dimension
  int space_dimension() { return dim_; }

private:
  std::vector<T> coords_;
  const int nodesPerCell_;
  const int dim_;

}; // class Flat_Mesh_Wrapper


} // end namespace Portage

#endif // FLAT_HEX_MESH_WRAPPER_H_
