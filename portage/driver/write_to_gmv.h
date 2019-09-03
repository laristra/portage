/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef PORTAGE_WRITE_TO_GMV_H_
#define PORTAGE_WRITE_TO_GMV_H_

#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "portage/support/portage.h"

/* Do we need to write a variant of this which will work if we don't
 * have TANGRAM and multi-material state? */

#include "tangram/driver/CellMatPoly.h"

namespace Portage {

/*!
  @brief Method to write out material polygons from an interface
  reconstruction and any associated material fields

  @tparam D   Dimension of problem
  @tparam Mesh_Wrapper   A Mesh class
  @tparam State_Wrapper  A State Manager class
  @tparam InterfaceReconstructor  An Interface Reconstruction class

  @param mesh   A mesh object  (for mesh topology/geometry queries)
  @param state  A state object (for field and material queries)
  @param ir     An interface reconstructor to create material polygons from volume fractions (and optionally material centroids)
  @fieldnames   List of fields to write out
  @filename     Name of file to write to
*/

template<int D, class Mesh_Wrapper, class State_Wrapper,
         class InterfaceReconstructor>
void write_to_gmv(Mesh_Wrapper const& mesh,
                  State_Wrapper const& state,
                  std::shared_ptr<InterfaceReconstructor> ir,
                  std::vector<std::string> &fieldnames,
                  std::string filename) {

  std::cerr << "Writing material and field data to GMV file\n";

  int nmats = state.num_materials();
  int nc = mesh.num_entities(Entity_kind::CELL, Entity_type::PARALLEL_OWNED);
  int nallc = mesh.num_entities(Entity_kind::CELL, Entity_type::ALL);

  std::vector<int> cell_num_mats(nc, 0);
  std::vector<std::vector<int>> cell_mat_ids(nc);

  std::vector<int> icell_owned2all(nc, -1);
  for (int ic = 0, iowned = 0; ic < nallc; ic++)
    if (mesh.cell_get_type(ic) == Wonton::PARALLEL_OWNED) {
      state.cell_get_mats(iowned, &(cell_mat_ids[iowned]));
      cell_num_mats[iowned] = cell_mat_ids[iowned].size();
      icell_owned2all[iowned] = ic;
      iowned++;
    }

  std::ofstream fout(filename);
  fout << std::scientific;
  fout.precision(17);

  fout << "gmvinput ascii" << std::endl;
  fout << "codename Portage" << std::endl;
  fout << "simdate 01/01/01" << std::endl;

  // Make a list of unique points in the mesh and cell mat poly
  // structures combined

  int np = mesh.num_entities(Entity_kind::NODE, Entity_type::PARALLEL_OWNED);
  int nallp = mesh.num_entities(Entity_kind::NODE, Entity_type::ALL);
  int nmatpnts = np;

  std::vector<int> inode_owned2all(np, -1);
  std::vector<int> inode_all2owned(nallp, -1);
  for (int ip = 0, iowned = 0; ip < nallp; ip++)
    if (mesh.node_get_type(ip) == Wonton::PARALLEL_OWNED) {
      inode_owned2all[iowned] = ip;
      inode_all2owned[ip] = iowned;
      iowned++;
    }

  for (int c = 0; c < nc; c++) {
    if (cell_num_mats[c] > 1) {  // Mixed cell
      Tangram::CellMatPoly<D> const& cellmatpoly = 
        ir->cell_matpoly_data(icell_owned2all[c]);

      int ncp = cellmatpoly.num_matvertices();
      for (int i = 0; i < ncp; i++)
        if (cellmatpoly.matvertex_parent_kind(i) != Tangram::Entity_kind::NODE)
          nmatpnts++;  // point does not correspond to mesh node - add to list
    }
  }

  std::vector<Point<D>> points(nmatpnts);

  for (int i = 0; i < np; i++)
    mesh.node_get_coordinates(inode_owned2all[i], &(points[i]));

  nmatpnts = np;  // reset nmatpnts so it can be used as a counter
  for (int c = 0; c < nc; c++) {
    if (cell_num_mats[c] > 1) {  // Mixed cell
      Tangram::CellMatPoly<D> const& cellmatpoly = 
        ir->cell_matpoly_data(icell_owned2all[c]);

      int ncp = cellmatpoly.num_matvertices();
      for (int i = 0; i < ncp; i++) {
        if (cellmatpoly.matvertex_parent_kind(i) != Tangram::Entity_kind::NODE) {
          points[nmatpnts] = Point<D>(cellmatpoly.matvertex_point(i));
          nmatpnts++;
        }
      }
    }
  }

  fout << "nodev " << nmatpnts << std::endl;
  for (int ip = 0; ip < nmatpnts; ip++) {
    fout << points[ip];
    if (D == 2)  // GMV requires us to output 3 coordinates
      fout << " " << 0.000000;
    fout << std::endl;
  }

  // Count the number of polygons we will write out - assume one
  // polygon per material per cell

  int npoly = 0;
  for (int c = 0; c < nc; c++)
    npoly += cell_num_mats[c];
  fout << "cells " << npoly << std::endl;

  // Now write out material polygons one material at a time (this
  // makes it easier to also write out corresponding field data

  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    state.mat_get_cells(m, &matcells);

    for (int c : matcells) {
      if (cell_num_mats[c] > 1) {  // multi-material cell

        Tangram::CellMatPoly<D> const& cellmatpoly = 
          ir->cell_matpoly_data(icell_owned2all[c]);

        // Write out the polyhedra corresponding to this material
        int nmp = cellmatpoly.num_matpolys();
        for (int i = 0; i < nmp; i++) {

          if (cellmatpoly.matpoly_matid(i) != m) continue;

          if (D == 1 || D == 2) {
            std::vector<int> mverts = cellmatpoly.matpoly_vertices(i);
            fout << "general 1 " << mverts.size() << " ";
            for (auto n : mverts) {
              if (cellmatpoly.matvertex_parent_kind(n) == Tangram::Entity_kind::NODE)
                fout << inode_all2owned[cellmatpoly.matvertex_parent_id(n)] + 1 << " ";
              else {
                Point<D> pnt = cellmatpoly.matvertex_point(n);
                for (int j = np; j < nmatpnts; j++) {
                  if (points[j] == pnt) {
                    fout << j+1 << " ";
                    break;
                  }
                }
              }
            }
            fout << std::endl;

          } else if (D == 3) {

            std::vector<int> const& mfaces = cellmatpoly.matpoly_faces(i);
            fout << "general " << mfaces.size() << std::endl;
            for (auto f : mfaces) {
              std::vector<int> const& mfverts = cellmatpoly.matface_vertices(f);
              fout << mfverts.size() << " ";
            }
            fout << std::endl;
            for (auto f : mfaces) {
              std::vector<int> const& mfverts = cellmatpoly.matface_vertices(f);
              int nfv = mfverts.size();
              int dir, mp0, mp1;
              cellmatpoly.matface_matpolys(f, &mp0, &mp1);
              if (mp0 == i) {  // Natural order of vertices
                for (int j = 0; j < nfv; j++) {
                  int n = mfverts[j];
                  if (cellmatpoly.matvertex_parent_kind(n) == Tangram::Entity_kind::NODE)
                    fout << inode_all2owned[cellmatpoly.matvertex_parent_id(n)] + 1 << " ";
                  else {
                    Point<D> pnt = cellmatpoly.matvertex_point(n);
                    for (int j = np; j < nmatpnts; j++) {
                      if (points[j] == pnt) {
                        fout << j+1 << " ";
                        break;
                      }
                    }
                  }
                }
                fout << std::endl;
              } else {  // Reverse order of vertices
                for (int j = 0; j < nfv; j++) {
                  int n = mfverts[nfv-j-1];
                  if (cellmatpoly.matvertex_parent_kind(n) == Tangram::Entity_kind::NODE)
                    fout << inode_all2owned[cellmatpoly.matvertex_parent_id(n)] + 1 << " ";
                  else {
                    Point<D> pnt = cellmatpoly.matvertex_point(n);
                    for (int j = np; j < nmatpnts; j++) {
                      if (points[j] == pnt) {
                        fout << j+1 << " ";
                        break;
                      }
                    }
                  }
                }
                fout << std::endl;
              }
            }  // for (auto f : mfaces)
          }  // else if (D == 3)
        }  // for (i = 0; i < nmp; i++)

      } else {  // single material cell

        if (cell_mat_ids[c][0] == m) {
          if (D == 1 || D == 2) {
            // Write out the cell
            std::vector<int> cverts;
            mesh.cell_get_nodes(icell_owned2all[c], &cverts);

            // write cell out as polygons (even if its a quad or tri)
            fout << "general 1 " << cverts.size() << " ";
            for (auto n : cverts)
              fout << n+1 << " ";
            fout << std::endl;
          } else if (D == 3) {
            int nf;
            std::vector<int> cfaces;
            std::vector<int> cfdirs;
            mesh.cell_get_faces_and_dirs(icell_owned2all[c], &cfaces, &cfdirs);
            fout << "general " << cfaces.size() << std::endl;
            for (auto f : cfaces) {
              std::vector<int> fverts;
              mesh.face_get_nodes(f, &fverts);
              fout << fverts.size() << " ";
            }
            fout << std::endl;

            int j = 0;
            for (auto f : cfaces) {
              std::vector<int> fverts;
              mesh.face_get_nodes(f, &fverts);
              int nfverts = fverts.size();
              if (cfdirs[j] == 1) {
                for (int k = 0; k < nfverts; k++)
                  fout << fverts[k]+1 << " ";
                fout << std::endl;
              } else {
                for (int k = 0; k < nfverts; k++)
                  fout << fverts[nfverts-k-1]+1 << " ";
                fout << std::endl;
              }
            }
          }
        }  // if (cell_mat_ids[c][0] == m)

      }  // if (cell_num_mats[c] > 1) ... else ...
    }  // for (c : matcells)
  }  // for (m = 0; m < nmats; m++)


  // Write out material ID of polygons

  fout << "material " << std::endl;
  fout << nmats << " 0" << std::endl;
  for (int m = 0; m < nmats; m++)
    fout << "mat" << m+1 << std::endl;

  for (int m = 0; m < nmats; m++) {
    std::vector<int> matcells;
    state.mat_get_cells(m, &matcells);

    for (int c : matcells) {
      for (int cm : cell_mat_ids[c]) {
        if (cm == m) {
          fout << m+1 << " ";
          break;
        }
      }
    }
  }
  fout << std::endl;

  // Write out any requested fields
  if (fieldnames.size()) {
    fout << "variable" << std::endl;
    for (auto & fieldname : fieldnames) {
      if (state.field_type(Portage::Entity_kind::CELL, fieldname) ==
          Portage::Field_type::UNKNOWN_TYPE_FIELD)
        continue;

      fout << fieldname << " 0 " << std::endl;
      for (int m = 0; m < nmats; m++) {
        std::vector<int> matcells;
        state.mat_get_cells(m, &matcells);
        int nmatcells = matcells.size();

        double const *matvec;
        state.mat_get_celldata(fieldname, m, &matvec);

        for (int i = 0; i < nmatcells; i++)
          fout << matvec[i] << std::endl;
      }
    }
    fout << "endvars" << std::endl;
  }

  fout << "endgmv" << std::endl;
}


}  // namespace Portage

#endif  // PORTAGE_WRITE_TO_GMV_H_
