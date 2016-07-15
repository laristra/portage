/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef MPI_BOUNDING_BOXES_H_
#define MPI_BOUNDING_BOXES_H_

#define DEBUG_MPI

#include "portage/support/portage.h"

/*!
  @file mpi_bounding_boxes.h
  @brief Distributes source data using MPI based on bounding boxes
 */

namespace Portage {

/*!
  @class MPI_Bounding_Boxes
  @brief Distributes source data using MPI based on bounding boxes

         Currently assumes coordinates and all fields are doubles.
*/
class MPI_Bounding_Boxes {
 public:

  /*!
    @brief Constructor of MPI_Bounding_Boxes
   */
  MPI_Bounding_Boxes() {}

  
  /*!
    @brief Destructor of MPI_Bounding_Boxes
   */
  ~MPI_Bounding_Boxes() {};


  /*!
    @brief Compute bounding boxes for all partitions, and send source mesh and state
           information to all target partitions with an overlapping bounding box using MPI
    @param[in] source_mesh_flat  Input mesh (must be flat representation)
    @param[in] source_state_flat Input state (must be flat representation)
    @param[in] target_mesh       Target mesh
    @param[in] target_state      Target state (not actually used for now)
   */
  template <class Source_Mesh, class Source_State, class Target_Mesh, class Target_State>
  void distribute(Source_Mesh &source_mesh_flat, Source_State &source_state_flat, 
                  Target_Mesh &target_mesh, Target_State &target_state)
  {
    // Get the MPI communicator size and rank information
    int commSize, commRank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    // Compute the bounding box for the target mesh on this rank, and put it in an
    // array that will later store the target bounding box for each rank 
    int targetNumCells = target_mesh.num_owned_cells();
    int targetDim = target_mesh.space_dimension();
    double targetBoundingBoxes[2*targetDim*commSize];
    for (unsigned int i=0; i<2*targetDim; i+=2)
    {
      targetBoundingBoxes[2*targetDim*commRank+i+0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2*targetDim*commRank+i+1] = -std::numeric_limits<double>::max();
    }
    for (unsigned int c=0; c<targetNumCells; c++)
    {
      std::vector<Portage::Point<3>> cellCoord;
      target_mesh.cell_get_coordinates(c, &cellCoord);
      for (unsigned int j=0; j<cellCoord.size(); j++)
      {
        for (unsigned int k=0; k<targetDim; k++)
          if (cellCoord[j][k] < targetBoundingBoxes[2*targetDim*commRank+2*k])
            targetBoundingBoxes[2*targetDim*commRank+2*k] = cellCoord[j][k];
        for (unsigned int k=0; k<targetDim; k++)
          if (cellCoord[j][k] > targetBoundingBoxes[2*targetDim*commRank+2*k+1])
            targetBoundingBoxes[2*targetDim*commRank+2*k+1] = cellCoord[j][k];
      }
    }

    // Get the source mesh properties and cordinates
    int sourceNumCells = source_mesh_flat.num_owned_cells();
    int sourceDim = source_mesh_flat.space_dimension();
    int sourceNodesPerCell = source_mesh_flat.get_nodes_per_cell();
    int sourceCellStride = sourceNodesPerCell*sourceDim;
    Portage::vector<double>& sourceCoords = source_mesh_flat.get_coords();

    // Compute the bounding box for the source mesh on this rank
    double sourceBoundingBox[2*sourceDim];
    for (unsigned int i=0; i<2*sourceDim; i+=2)
    {
      sourceBoundingBox[i+0] = std::numeric_limits<double>::max();
      sourceBoundingBox[i+1] = -std::numeric_limits<double>::max();
    }
    for (unsigned int c=0; c<sourceNumCells; c++)
    {
      for (unsigned int j=0; j<sourceNodesPerCell; j++)
      {
        for (unsigned int k=0; k<sourceDim; k++)
          if (sourceCoords[c*sourceCellStride+j*sourceDim+k] < sourceBoundingBox[2*k])
            sourceBoundingBox[2*k] = sourceCoords[c*sourceCellStride+j*sourceDim+k];
        for (unsigned int k=0; k<sourceDim; k++)
          if (sourceCoords[c*sourceCellStride+j*sourceDim+k] > sourceBoundingBox[2*k+1])
            sourceBoundingBox[2*k+1] = sourceCoords[c*sourceCellStride+j*sourceDim+k];
      }
    }

    // Broadcast the target bounding boxes so that each rank knows the bounding boxes for all ranks
    for (unsigned int i=0; i<commSize; i++)
      MPI_Bcast(&(targetBoundingBoxes[0])+i*2*targetDim, 2*targetDim, MPI_DOUBLE, i, MPI_COMM_WORLD);

#ifdef DEBUG_MPI
    std::cout << "Target boxes: ";
    for (unsigned int i=0; i<2*targetDim*commSize; i++) std::cout << targetBoundingBoxes[i] << " ";
    std::cout << std::endl;

    std::cout << "Source box : " << commRank << " ";
    for (unsigned int i=0; i<2*sourceDim; i++) std::cout << sourceBoundingBox[i] << " ";
    std::cout << std::endl;
#endif

    // Offset the source bounding boxes by a fudge factor so that we don't send source data to a rank
    // with a target bounding box that is incident to but not overlapping the source bounding box
    const double boxOffset = 2.0*std::numeric_limits<double>::epsilon();
    double min_x2 = sourceBoundingBox[0]+boxOffset;  double max_x2 = sourceBoundingBox[1]-boxOffset;
    double min_y2 = sourceBoundingBox[2]+boxOffset;  double max_y2 = sourceBoundingBox[3]-boxOffset;
    double min_z2 = sourceBoundingBox[4]+boxOffset;  double max_z2 = sourceBoundingBox[5]-boxOffset;

    // For each target rank with a bounding box that overlaps the bounding box for this rank's partition
    // of the source mesh, we will send it all our source cells; otherwise, we will send it nothing
    std::vector<int> sendCounts(commSize, 0);
    std::vector<int> recvCounts(commSize); 
    for (unsigned int i=0; i<commSize; i++)
    {
      double min_x1 = targetBoundingBoxes[2*targetDim*i+0]+boxOffset;  
      double max_x1 = targetBoundingBoxes[2*targetDim*i+1]-boxOffset;
      double min_y1 = targetBoundingBoxes[2*targetDim*i+2]+boxOffset;  
      double max_y1 = targetBoundingBoxes[2*targetDim*i+3]-boxOffset;
      double min_z1 = targetBoundingBoxes[2*targetDim*i+4]+boxOffset;  
      double max_z1 = targetBoundingBoxes[2*targetDim*i+5]-boxOffset;

      sendCounts[i] = ( ((min_x1 <= min_x2 && min_x2 <= max_x1) || (min_x2 <= min_x1 && min_x1 <= max_x2)) &&
                      ((min_y1 <= min_y2 && min_y2 <= max_y1) || (min_y2 <= min_y1 && min_y1 <= max_y2)) &&
                      ((min_z1 <= min_z2 && min_z2 <= max_z1) || (min_z2 <= min_z1 && min_z1 <= max_z2)) ) ? sourceNumCells : 0;
    }

    // Each rank will tell each other rank how many cells it is going to send it
    for (unsigned int i=0; i<commSize; i++)
        MPI_Scatter(&(sendCounts[0]), 1, MPI_INT, &(recvCounts[0])+i, 1, MPI_INT, i, MPI_COMM_WORLD);

    // Compute the total number of source cells this rank will receive from all ranks, 
    // and allocate a vector to hold them
    int totalRecvSize = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvSize += recvCounts[i]; 
    std::vector<double> newData(sourceCellStride*totalRecvSize);

#ifdef DEBUG_MPI
    std::cout << "Received " << commRank << " " << recvCounts[0] << " " << recvCounts[1] << " " << totalRecvSize
              << " " << (totalRecvSize - recvCounts[commRank]) << std::endl;
#endif

    // Copy source cells that will stay on this rank into the proper place in the new vector
    int localOffset = 0;
    for (unsigned int i=0; i<commRank; i++) localOffset += recvCounts[i]; 
    if (recvCounts[commRank] > 0)
      std::copy(sourceCoords.begin(), sourceCoords.end(), newData.begin() + sourceCellStride*localOffset);

    // Each rank will send and receive the appropriate source cells
    MPI_Status stat;
    int writeOffset = 0;

    for (unsigned int i=0; i<commSize; i++)
    {
      if (commRank == i)
      {
        for (unsigned int j=0; j<commSize; j++)
        {
          if ((j != i) && (sendCounts[j] > 0))
          {
            //std::cout << commRank << " sending " << sourceCoords[10] << " " << sendCounts[j] << " to " << j << std::endl;
            MPI_Send(&(sourceCoords[0]), sourceCellStride*sendCounts[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
          }
        }
      }
      else
      {
        if (recvCounts[i] > 0)
        {
          MPI_Recv(&(newData[0])+writeOffset, sourceCellStride*recvCounts[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
          //std::cout << commRank << " receiving " << recvCounts[i] << " from " << i << " at " << writeOffset << std::endl;
        }
      }
      writeOffset += sourceCellStride*recvCounts[i];
      MPI_Barrier(MPI_COMM_WORLD);
    }

    // We will now use the received source data as our new source mesh on this partition
    sourceCoords = newData;

    // Each rank will send and receive the appropriate source state
    for (unsigned int s=0; s<source_state_flat.get_num_vectors(); s++)
    {
      Portage::vector<double>& sourceState = source_state_flat.get_vector(s);

      std::vector<double> newField(totalRecvSize);
      std::copy(sourceState.begin(), sourceState.end(), newField.begin() + localOffset);
      writeOffset = 0;
      for (unsigned int i=0; i<commSize; i++)
      {
        if (commRank == i)
        {
          for (unsigned int j=0; j<commSize; j++)
          {
            if ((j != i) && (sendCounts[j] > 0))
            {
              //std::cout << commRank << " sending " << sourceState[0] << " " << sendCounts[j] << " to " << j << std::endl;
              MPI_Send(&(sourceState[0]), sendCounts[j], MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            }
          }
        }
        else
        {
          if (recvCounts[i] > 0)
          {
            MPI_Recv(&(newField[0])+writeOffset, recvCounts[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stat);
            //std::cout << commRank << " receiving " << newField[writeOffset] << " " << recvCounts[i] << " from " << i << std::endl;
          }
        }
        writeOffset += recvCounts[i];
        MPI_Barrier(MPI_COMM_WORLD);
      }

      // We will now use the received source state as our new source state on this partition
      sourceState = newField;

#ifdef DEBUG_MPI
      if (commRank == 1)
      {
        std::cout << "Sizes: " << commRank << " " << totalRecvSize << " " << localOffset << " " << targetNumCells
                  << " " << sourceNumCells << std::endl;
        for (unsigned int i=0; i<sourceCoords.size(); i++)
          std::cout << sourceCoords[i] << " " << newData[i] << " ";
        std::cout << std::endl;
        for (unsigned int i=0; i<sourceState.size(); i++)
          std::cout << sourceState[i] << " ";
        std::cout << std::endl;
      }
#endif
    }
  }

  
}; // MPI_Bounding_Boxes

} // namespace Portage

#endif // MPI_Bounding_Boxes_H_
