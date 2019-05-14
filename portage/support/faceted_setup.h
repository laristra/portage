/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef FACETED_SETUP_H
#define FACETED_SETUP_H

#include <vector>
#include <algorithm>

#include "wonton/mesh/AuxMeshTopology.h"
#include "wonton/support/Point.h"
#include "portage/support/weight.h"
#include "portage/support/portage.h"

namespace Portage {
  namespace Meshfree {
    namespace Weight {

      template<int DIM, class Mesh_Wrapper> void faceted_setup_cell
        (Mesh_Wrapper &mesh, 
         vector<std::vector<std::vector<double>>> &smoothing_lengths,
         vector<Wonton::Point<DIM>> &extents,
         double smoothing_factor)
      {
        int ncells=mesh.num_owned_cells();
        smoothing_lengths.resize(ncells);
        extents.resize(ncells);
        
        // get face normal and distance information
        for (int i=0; i<ncells; i++) {
          std::vector<int> faces, fdirs;
          mesh.cell_get_faces_and_dirs(i, &faces, &fdirs);

          int nfaces = faces.size();
          smoothing_lengths[i].resize(nfaces);

          Wonton::Point<DIM> fcent, ccent, normal, distance;
          std::vector<int> fnodes;
          std::vector<Wonton::Point<DIM>> fncoord;
          for (int j=0; j<nfaces; j++) {
            // get face data
            mesh.face_centroid(faces[j], &fcent);
            mesh.cell_centroid(i, &ccent);
            mesh.face_get_nodes(faces[j], &fnodes);
            fncoord.resize(fnodes.size()+1);
            for (int k=0; k<fnodes.size(); k++) {
              mesh.node_get_coordinates(fnodes[k], &(fncoord[k]));
            }
            fncoord[fnodes.size()] = fncoord[0]; 

            if (DIM==2) {
              // Get 2d face normal 
              std::vector<double> delta(2);
              for (int k=0;k<2;k++) delta[k] = fncoord[1][k] - fncoord[0][k];
              normal[0] =  delta[1];
              normal[1] = -delta[0];
              if (fdirs[j]<0) normal*=-1.0;
              double norm = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
              normal /= norm;
              smoothing_lengths[i][j].resize(3);
              smoothing_lengths[i][j][0] = normal[0];
              smoothing_lengths[i][j][1] = normal[1];
            } else if (DIM==3) {
              // Get 3d face normal using average cross product, in case face is not flat
              // or has nodes close together.
              std::vector<std::vector<double>> bivec(2,std::vector<double>(3,0.));
              for (int k=0;k<3;k++) {
                bivec[0][k] = fncoord[0][k] - fcent[k];
                normal[k] = 0.;
              }
              for (int k=1; k<=fnodes.size(); k++) {
                for (int m=0;m<3;m++) bivec[1][m] = fncoord[k][m] - fcent[m];
                std::vector<double> cross(3,0.);
                cross[0] =  bivec[1][2]*bivec[0][1] - bivec[1][1]*bivec[0][2];
                cross[1] = -bivec[1][2]*bivec[0][0] + bivec[1][0]*bivec[0][2];
                cross[2] =  bivec[1][1]*bivec[0][0] - bivec[1][0]*bivec[0][1];
                double norm=sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
                for (int m=0;m<3;m++) {
                  cross[m]/=norm;
                }
                normal += cross;
                bivec[0] = bivec[1];
              }
              if (fdirs[j]<0) normal*=-1.0;
              double norm=0.;
              for (int k=0;k<3;k++) norm+=normal[k]*normal[k];
              norm=sqrt(norm);
              normal /= norm;
              smoothing_lengths[i][j].resize(4);
              for (int k=0;k<3;k++) smoothing_lengths[i][j][k] = normal[k];
            }

            distance = Wonton::Point<DIM>(fcent-ccent);
            double smoothing = 0.0;
            for (int k=0; k<DIM; k++) smoothing += 2*distance[k]*normal[k];
            if (smoothing<0.0) {  // when nodes are ordered backwards
              normal = -normal;
              smoothing = -smoothing;
            }
            smoothing_lengths[i][j][DIM] = smoothing_factor*smoothing;
          }
        }

        // get extent information
        for (int i=0; i<ncells; i++) {
          // calculate the min and max of coordinates
          std::vector<int> node_indices;
          mesh.cell_get_nodes(i, &node_indices);
          Wonton::Point<DIM> cmin, cmax, first_coords; 
          mesh.node_get_coordinates(node_indices[0], &first_coords);
          for (int k=0; k<DIM; k++) {cmin[k] = first_coords[k]; cmax[k] = first_coords[k];}
          for (int j=1; j<node_indices.size(); j++) {
            Wonton::Point<DIM> node_coords;
            mesh.node_get_coordinates(node_indices[j], &node_coords);
            for (int k=0; k<DIM; k++) {
              cmin[k] = std::min(node_coords[k], cmin[k]);
              cmax[k] = std::max(node_coords[k], cmax[k]);
            }
          }

          // subtract to get extents
          for (int k=0; k<DIM; k++) extents[i][k] = cmax[k] - cmin[k];

          // multiply by smoothing_factor
          for (int k=0; k<DIM; k++) extents[i][k] *= 2.*smoothing_factor;
        }
      }

    } // namespace Weight
  } // namespace Meshfree
} // namespace Portage
#endif
