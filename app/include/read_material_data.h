/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <stdlib.h>

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>

#include "portage/support/portage.h"

/*!
 @brief Reads material data for a mesh from a provided binary file.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.
 @tparam D            Dimensionality of problem

 @param[in] mesh Mesh wrapper.
 @param[in] mesh_data_fname Name of the file with material data.
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
*/
template <class Mesh_Wrapper, int D>
void read_material_data(const Mesh_Wrapper& mesh,
                        const std::string& mesh_data_fname,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Portage::Point<D>>& cell_mat_centroids) {

  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();

  std::ifstream os(mesh_data_fname.c_str(), std::ifstream::binary);
  if (!os.good()) {
    std::ostringstream os;
    os << std::endl << "Cannot open " << mesh_data_fname <<
      " for binary input" << std::endl;
    throw std::runtime_error(os.str());
  }

  int data_dim;
  os.read(reinterpret_cast<char *>(&data_dim), sizeof(int));
  assert(data_dim == D);
  int ncells;
  os.read(reinterpret_cast<char *>(&ncells), sizeof(int));
  assert(ncells = mesh.num_owned_cells());
  cell_num_mats.resize(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    os.read(reinterpret_cast<char *>(&cell_num_mats[icell]), sizeof(int));
    for (int im = 0; im < cell_num_mats[icell]; im++) {
      int imat;
      os.read(reinterpret_cast<char *>(&imat), sizeof(int));
      cell_mat_ids.push_back(imat);
    }
  }
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      cell_mat_volfracs.push_back(1.0);
      continue;
    }
    for (int im = 0; im < cell_num_mats[icell]; im++) {
      double vfrac;
      os.read(reinterpret_cast<char *>(&vfrac), sizeof(double));
      cell_mat_volfracs.push_back(vfrac);
    }
  }
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      Portage::Point<D> cur_cell_cen;
      mesh.cell_centroid(icell, &cur_cell_cen);
      cell_mat_centroids.push_back(Portage::Point<D>(cur_cell_cen));
      continue;
    }
    for (int im = 0; im < cell_num_mats[icell]; im++) {
      double cen_x, cen_y, cen_z;
      os.read(reinterpret_cast<char *>(&cen_x), sizeof(double));
      os.read(reinterpret_cast<char *>(&cen_y), sizeof(double));
      if (D == 2)
        cell_mat_centroids.push_back(Portage::Point<D>(cen_x, cen_y));
      else if (D == 3) {
        os.read(reinterpret_cast<char *>(&cen_z), sizeof(double));
        cell_mat_centroids.push_back(Portage::Point<D>(cen_x, cen_y, cen_z));
      }
    }
  }
  os.close();
}
