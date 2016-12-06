/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef MPI_BOUNDING_BOXES_H_
#define MPI_BOUNDING_BOXES_H_

//#define DEBUG_MPI

#include <cassert>
#include <algorithm>
#include <memory>

#include "portage/support/portage.h"

#include "mpi.h"

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

    int dim = source_mesh_flat.space_dimension();
    assert(dim == target_mesh.space_dimension());

    // Compute the bounding box for the target mesh on this rank, and put it in an
    // array that will later store the target bounding box for each rank 
    int targetNumCells = target_mesh.num_owned_cells();
    double targetBoundingBoxes[2*dim*commSize];
    for (unsigned int i=0; i<2*dim; i+=2)
    {
      targetBoundingBoxes[2*dim*commRank+i+0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2*dim*commRank+i+1] = -std::numeric_limits<double>::max();
    }
    for (unsigned int c=0; c<targetNumCells; c++)
    {
      // ugly hack, since dim is not known at compile time
      if (dim == 3)
      {
        std::vector<Portage::Point<3>> cellCoord;
        target_mesh.cell_get_coordinates(c, &cellCoord);
        for (unsigned int j=0; j<cellCoord.size(); j++)
        {
          for (unsigned int k=0; k<dim; k++)
            if (cellCoord[j][k] < targetBoundingBoxes[2*dim*commRank+2*k])
              targetBoundingBoxes[2*dim*commRank+2*k] = cellCoord[j][k];
          for (unsigned int k=0; k<dim; k++)
            if (cellCoord[j][k] > targetBoundingBoxes[2*dim*commRank+2*k+1])
              targetBoundingBoxes[2*dim*commRank+2*k+1] = cellCoord[j][k];
        }
      }
      else if (dim == 2)
      {
        std::vector<Portage::Point<2>> cellCoord;
        target_mesh.cell_get_coordinates(c, &cellCoord);
        for (unsigned int j=0; j<cellCoord.size(); j++)
        {
          for (unsigned int k=0; k<dim; k++)
            if (cellCoord[j][k] < targetBoundingBoxes[2*dim*commRank+2*k])
              targetBoundingBoxes[2*dim*commRank+2*k] = cellCoord[j][k];
          for (unsigned int k=0; k<dim; k++)
            if (cellCoord[j][k] > targetBoundingBoxes[2*dim*commRank+2*k+1])
              targetBoundingBoxes[2*dim*commRank+2*k+1] = cellCoord[j][k];
        }
      }
    }

    // Get the source mesh properties and cordinates
    int sourceNumOwnedCells = source_mesh_flat.num_owned_cells();
    int sourceNumCells = source_mesh_flat.num_owned_cells() + source_mesh_flat.num_ghost_cells();
    int sourceNumOwnedNodes = source_mesh_flat.num_owned_nodes();
    int sourceNumNodes = source_mesh_flat.num_owned_nodes() + source_mesh_flat.num_ghost_nodes();
    std::vector<double>& sourceCoords = source_mesh_flat.get_coords();
    std::vector<int>& sourceNodeCounts = source_mesh_flat.get_node_counts();
    std::vector<int>& sourceGlobalCellIds = source_mesh_flat.get_global_cell_ids();
    std::vector<int>& sourceNeighborCounts = source_mesh_flat.get_neighbor_counts();
    std::vector<int>& sourceNeighbors = source_mesh_flat.get_neighbors();

    // Compute the bounding box for the source mesh on this rank
    double sourceBoundingBox[2*dim];
    for (unsigned int i=0; i<2*dim; i+=2)
    {
      sourceBoundingBox[i+0] = std::numeric_limits<double>::max();
      sourceBoundingBox[i+1] = -std::numeric_limits<double>::max();
    }
    // Since we're doing a combined bounding box for all cells,
    // don't need to distinguish which nodes belong to which cell
    for (unsigned int n=0; n<sourceNumOwnedNodes; n++)
    {
      for (unsigned int k=0; k<dim; k++)
        if (sourceCoords[n*dim+k] < sourceBoundingBox[2*k])
          sourceBoundingBox[2*k] = sourceCoords[n*dim+k];
      for (unsigned int k=0; k<dim; k++)
        if (sourceCoords[n*dim+k] > sourceBoundingBox[2*k+1])
          sourceBoundingBox[2*k+1] = sourceCoords[n*dim+k];
    }

    // Broadcast the target bounding boxes so that each rank knows the bounding boxes for all ranks
    for (unsigned int i=0; i<commSize; i++)
      MPI_Bcast(&(targetBoundingBoxes[0])+i*2*dim, 2*dim, MPI_DOUBLE, i, MPI_COMM_WORLD);

#ifdef DEBUG_MPI
    std::cout << "Target boxes: ";
    for (unsigned int i=0; i<2*dim*commSize; i++) std::cout << targetBoundingBoxes[i] << " ";
    std::cout << std::endl;

    std::cout << "Source box : " << commRank << " ";
    for (unsigned int i=0; i<2*dim; i++) std::cout << sourceBoundingBox[i] << " ";
    std::cout << std::endl;
#endif

    // Offset the source bounding boxes by a fudge factor so that we don't send source data to a rank
    // with a target bounding box that is incident to but not overlapping the source bounding box
    const double boxOffset = 2.0*std::numeric_limits<double>::epsilon();
    double min2[dim], max2[dim];
    for (unsigned int k=0; k<dim; k++)
    {
      min2[k] = sourceBoundingBox[2*k]+boxOffset;
      max2[k] = sourceBoundingBox[2*k+1]-boxOffset;
    }

    // For each target rank with a bounding box that overlaps the bounding box for this rank's partition
    // of the source mesh, we will send it all our source cells; otherwise, we will send it nothing
    std::vector<int> sendCellCounts(commSize, 0);
    std::vector<int> recvCellCounts(commSize); 
    for (unsigned int i=0; i<commSize; i++)
    {
      double min1[dim], max1[dim];
      bool sendThis = true;
      for (unsigned int k=0; k<dim; k++)
      {
        min1[k] = targetBoundingBoxes[2*dim*i+2*k]+boxOffset;  
        max1[k] = targetBoundingBoxes[2*dim*i+2*k+1]-boxOffset;
        sendThis = sendThis &&
            ((min1[k] <= min2[k] && min2[k] <= max1[k]) ||
             (min2[k] <= min1[k] && min1[k] <= max2[k]));
      }

      sendCellCounts[i] = sendThis ? sourceNumCells : 0;
    }

    // Each rank will tell each other rank how many cells it is going to send it
    MPI_Alltoall(&(sendCellCounts[0]), 1, MPI_INT, &(recvCellCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Send how many cells are locally owned
    std::vector<int> sendOwnedCellCounts(commSize, 0);
    std::vector<int> recvOwnedCellCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendOwnedCellCounts[i] = (sendCellCounts[i] > 0) ? sourceNumOwnedCells : 0;
    MPI_Alltoall(&(sendOwnedCellCounts[0]), 1, MPI_INT, &(recvOwnedCellCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Compute the total number of source cells this rank will receive from all ranks
    int totalRecvCells = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvCells += recvCellCounts[i]; 
    int totalRecvOwnedCells = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvOwnedCells += recvOwnedCellCounts[i];

#ifdef DEBUG_MPI
    std::cout << "Received " << commRank << " " << recvCellCounts[0] << " " << recvCellCounts[1] << " " << totalRecvCells
              << " " << (totalRecvCells - recvCellCounts[commRank]) << std::endl;
#endif

    // Tell each other rank how many nodes it is going to send it
    std::vector<int> sendNodeCounts(commSize, 0);
    std::vector<int> recvNodeCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendNodeCounts[i] = (sendCellCounts[i] > 0) ? sourceNumNodes : 0;
    MPI_Alltoall(&(sendNodeCounts[0]), 1, MPI_INT, &(recvNodeCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Send how many nodes are locally owned
    std::vector<int> sendOwnedNodeCounts(commSize, 0);
    std::vector<int> recvOwnedNodeCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendOwnedNodeCounts[i] = (sendCellCounts[i] > 0) ? sourceNumOwnedNodes : 0;
    MPI_Alltoall(&(sendOwnedNodeCounts[0]), 1, MPI_INT, &(recvOwnedNodeCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Compute the total number of source nodes this rank will receive from all ranks
    int totalRecvNodes = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvNodes += recvNodeCounts[i]; 
    int totalRecvOwnedNodes = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvOwnedNodes += recvOwnedNodeCounts[i];

    // Compute the total number of neighbor indexes and owned neighbor indexes on this rank
    int sourceNumNeighbors = 0;
    for (unsigned int i=0; i<sourceNeighborCounts.size(); i++)
      sourceNumNeighbors += sourceNeighborCounts[i];
    int sourceNumOwnedNeighbors = 0;
    for (unsigned int i=0; i<sourceNumOwnedCells; i++)
      sourceNumOwnedNeighbors += sourceNeighborCounts[i];

    // Each rank will tell each other rank how many neighbor indexes it is going to send it
    std::vector<int> sendNeighborsCounts(commSize, 0);
    std::vector<int> recvNeighborsCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendNeighborsCounts[i] = (sendCellCounts[i] > 0) ? sourceNumNeighbors : 0;
    MPI_Alltoall(&(sendNeighborsCounts[0]), 1, MPI_INT, &(recvNeighborsCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Each rank will tell each other rank how many owned neighbor indexes it is going to send it
    std::vector<int> sendOwnedNeighborsCounts(commSize, 0);
    std::vector<int> recvOwnedNeighborsCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendOwnedNeighborsCounts[i] = (sendOwnedCellCounts[i] > 0) ? sourceNumOwnedNeighbors : 0;
    MPI_Alltoall(&(sendOwnedNeighborsCounts[0]), 1, MPI_INT, &(recvOwnedNeighborsCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Compute the total number of neighbor indexes this rank will receive from all ranks
    int totalRecvNeighbors = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvNeighbors += recvNeighborsCounts[i];
    int totalRecvOwnedNeighbors = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvOwnedNeighbors += recvOwnedNeighborsCounts[i];

    // Data structures to hold mesh data received from other ranks
    std::vector<double> newCoords(dim*totalRecvNodes);
    std::vector<int> newNodeCounts(totalRecvCells);
    std::vector<int> newGlobalCellIds(totalRecvCells);
    std::vector<int> newNeighborCounts(totalRecvCells);
    std::vector<int> newNeighbors(totalRecvNeighbors);

    // Perform two rounds of sends: the first for owned cells, and the second for ghost cells
    for (unsigned int r=0; r<2; r++)
    {
      // Compute all the offsets depending on whether this is the owned round or the ghost round
      std::vector<int> curCellRecvCounts(commSize);
      std::vector<int> curCellSendCounts(commSize);
      std::vector<int> curNodeRecvCounts(commSize);
      std::vector<int> curNodeSendCounts(commSize);
      std::vector<int> curNeighborRecvCounts(commSize);
      std::vector<int> curNeighborSendCounts(commSize);

      int sourceCellStart = 0;
      int sourceCellEnd = sourceNumOwnedCells;
      int localCellStart = 0;
      int localCellEnd = totalRecvOwnedCells;
      int sourceNodeStart = 0;
      int sourceNodeEnd = sourceNumOwnedNodes;
      int localNodeStart = 0;
      int localNodeEnd = totalRecvOwnedNodes;
      int sourceNeighborStart = 0;
      int sourceNeighborEnd = sourceNumOwnedNeighbors;
      int localNeighborStart = 0;
      int localNeighborEnd = totalRecvOwnedNeighbors;

      if (r == 1)
      {
        sourceCellStart = sourceNumOwnedCells;
        sourceCellEnd = sourceNumCells;
        localCellStart = totalRecvOwnedCells;
        localCellEnd = totalRecvCells;
        sourceNodeStart = sourceNumOwnedNodes;
        sourceNodeEnd = sourceNumNodes;
        localNodeStart = totalRecvOwnedNodes;
        localNodeEnd = totalRecvNodes;
        sourceNeighborStart = sourceNumOwnedNeighbors;
        sourceNeighborEnd = sourceNumNeighbors;
        localNeighborStart = totalRecvOwnedNeighbors;
        localNeighborEnd = totalRecvNeighbors;
      }

      for (unsigned int i=0; i<commSize; i++)
      {
        if (r == 0)
        {
          curCellRecvCounts[i] = recvOwnedCellCounts[i];
          curCellSendCounts[i] = sendOwnedCellCounts[i];
          curNodeRecvCounts[i] = recvOwnedNodeCounts[i];
          curNodeSendCounts[i] = sendOwnedNodeCounts[i];
          curNeighborRecvCounts[i] = recvOwnedNeighborsCounts[i];
          curNeighborSendCounts[i] = sendOwnedNeighborsCounts[i];
        }
        else 
        {
          curCellRecvCounts[i] = recvCellCounts[i] - recvOwnedCellCounts[i];
          curCellSendCounts[i] = sendCellCounts[i] - sendOwnedCellCounts[i];
          curNodeRecvCounts[i] = recvNodeCounts[i] - recvOwnedNodeCounts[i];
          curNodeSendCounts[i] = sendNodeCounts[i] - sendOwnedNodeCounts[i];
          curNeighborRecvCounts[i] = recvNeighborsCounts[i] - recvOwnedNeighborsCounts[i];
          curNeighborSendCounts[i] = sendNeighborsCounts[i] - sendOwnedNeighborsCounts[i];
        }
      }

      // SEND NUMBER OF NODES FOR EACH CELL 

      moveData(commRank, commSize, MPI_INT, 1,
               sourceCellStart, sourceCellEnd,
               localCellStart,
               curCellSendCounts, curCellRecvCounts,
               sourceNodeCounts, newNodeCounts);

      // SEND NODE COORDINATES

      moveData(commRank, commSize, MPI_DOUBLE, dim,
               sourceNodeStart, sourceNodeEnd,
               localNodeStart,
               curNodeSendCounts, curNodeRecvCounts,
               sourceCoords, newCoords);

      // SEND GLOBAL CELL IDS

      moveData(commRank, commSize, MPI_INT, 1,
               sourceCellStart, sourceCellEnd,
               localCellStart,
               curCellSendCounts, curCellRecvCounts,
               sourceGlobalCellIds, newGlobalCellIds);

      // SEND NUMBER OF NEIGHBORS FOR EACH CELL 
 
      moveData(commRank, commSize, MPI_INT, 1,
               sourceCellStart, sourceCellEnd,
               localCellStart,
               curCellSendCounts, curCellRecvCounts,
               sourceNeighborCounts, newNeighborCounts);

      // SEND LIST OF NEIGHBOR GLOBAL IDS FOR EACH CELL

      moveData(commRank, commSize, MPI_INT, 1,
               sourceNeighborStart, sourceNeighborEnd,
               localNeighborStart,
               curNeighborSendCounts, curNeighborRecvCounts,
               sourceNeighbors, newNeighbors);

#ifdef DEBUG_MPI
      if (commRank == 1)
      {
        std::cout << "newNeighbors: ";
        for (unsigned int i=0; i<newNeighbors.size(); i++) std::cout << newNeighbors[i] << " ";
        std::cout << std::endl;
      }
#endif
  
    } // for r

    
    // SEND FIELD VALUES

    // Send and receive each field to be remapped (might be more efficient to consolidate sends)
    for (int s=0; s<source_state_flat.get_num_vectors(); s++)
    {
      std::vector<double> newField(totalRecvCells);

      // Perform two rounds of sends: the first for owned cells, and the second for ghost cells
      for (unsigned int r=0; r<2; r++)
      {
        // Set up the offsets depending on whether this is the owned round or the ghost round
        std::vector<int> curCellRecvCounts(commSize);
        std::vector<int> curCellSendCounts(commSize);

        int sourceCellStart = 0;
        int sourceCellEnd = sourceNumOwnedCells;
        int localCellStart = 0;
        int localCellEnd = totalRecvOwnedCells;

        if (r == 1)
        {
          sourceCellStart = sourceNumOwnedCells;
          sourceCellEnd = sourceNumCells;
          localCellStart = totalRecvOwnedCells;
          localCellEnd = totalRecvCells;
        }

        for (unsigned int i=0; i<commSize; i++)
        {
          if (r == 0)
          {
            curCellRecvCounts[i] = recvOwnedCellCounts[i];
            curCellSendCounts[i] = sendOwnedCellCounts[i];
          }
          else
          {
            curCellRecvCounts[i] = recvCellCounts[i] - recvOwnedCellCounts[i];
            curCellSendCounts[i] = sendCellCounts[i] - sendOwnedCellCounts[i];
          }
        }

        std::shared_ptr<std::vector<double>> sourceField = source_state_flat.get_vector(s);
        int sourceFieldStride = source_state_flat.get_field_dim(s);

        moveData(commRank, commSize, MPI_DOUBLE, sourceFieldStride,
                 sourceCellStart, sourceCellEnd,
                 localCellStart,
                 curCellSendCounts, curCellRecvCounts,
                 *sourceField, newField);
    
#ifdef DEBUG_MPI
        if ((r == 1) && (commRank == 1))
        {
          std::cout << "Test: field " << newField.size() << std::endl;
          for (unsigned int f=0; f<newField.size(); f++)
            std::cout << newField[f] << " ";
          std::cout << std::endl << std::endl;
        }
#endif
      } // for r

      // We will now use the received source state as our new source state on this partition
      std::shared_ptr<std::vector<double>> sourceField = source_state_flat.get_vector(s);
      sourceField->resize(newField.size());
      std::copy(newField.begin(), newField.end(), sourceField->begin());
    } // for s

    // We will now use the received source mesh data as our new source mesh on this partition
    sourceCoords = newCoords;
    sourceNodeCounts = newNodeCounts;
    sourceGlobalCellIds = newGlobalCellIds;
    sourceNeighborCounts = newNeighborCounts;
    sourceNeighbors = newNeighbors;
    source_mesh_flat.set_num_owned_cells(totalRecvOwnedCells);
    source_mesh_flat.set_num_owned_nodes(totalRecvOwnedNodes);

    // Create all index maps
    source_mesh_flat.make_index_maps();

#ifdef DEBUG_MPI
    if (commRank == 1)
    {
      std::cout << "Sizes: " << commRank << " " << totalRecvCells << " " << targetNumCells
                << " " << sourceNumCells << " " << sourceCoords.size() << std::endl;
      
      for (unsigned int i=0; i<sourceGlobalCellIds.size(); i++)
        std::cout << sourceGlobalCellIds[i] << " ";
      std::cout << std::endl;
      for (unsigned int i=0; i<sourceCoords.size(); i++)
      {
        std::cout << sourceCoords[i] << " " ; 
        if (i % 24 == 23) std::cout << std::endl;
      }
      std::cout << std::endl << std::endl;
    }
#endif

  } // distribute

 private:
  /*!
    @brief Move values for a single data field to all ranks as needed
    @tparam[in] T                C++ type of data to be moved
    @param[in] commRank          MPI rank of this PE
    @param[in] commSize          Total number of MPI ranks
    @param[in] mpiType           MPI type of data (MPI_???) to be moved
    @param[in] stride            Stride of data field
    @param[in] sourceStart       Start location in source (send) data
    @param[in] sourceEnd         End location in source (send) data
    @param[in] localStart        Start location in local (recv) data
    @param[in] curSendCounts     Array of send sizes from me to all PEs
    @param[in] curRecvCounts     Array of recv sizes to me from all PEs
    @param[in] sourceData        Array of (old) source data
    @param[in] newData           Array of new source data
   */
  template<typename T>
  void moveData(int commRank, int commSize,
                MPI_Datatype mpiType, int stride,
                int sourceStart, int sourceEnd,
                int localStart,
                const std::vector<int>& curSendCounts,
                const std::vector<int>& curRecvCounts,
                const std::vector<T>& sourceData,
                std::vector<T>& newData)
  {
    // Each rank will do a non-blocking receive from each rank from
    // which it will receive data values
    int writeOffset = localStart;
    int myOffset;
    std::vector<MPI_Request> requests;
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (curRecvCounts[i] > 0))
      {
        MPI_Request request;
        MPI_Irecv((void *)&(newData[stride*writeOffset]),
                  stride*curRecvCounts[i], mpiType, i,
                  MPI_ANY_TAG, MPI_COMM_WORLD, &request);
        requests.push_back(request);
      }
      else if (i == commRank)
      {
        myOffset = writeOffset;
      }
      writeOffset += curRecvCounts[i];
    }

    // Copy data values that will stay on this rank into the
    // proper place in the new vector
    if (curRecvCounts[commRank] > 0)
      std::copy(sourceData.begin() + stride*sourceStart,
                sourceData.begin() + stride*sourceEnd, 
                newData.begin() + stride*myOffset);

    // Each rank will send its data values to appropriate ranks
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (curSendCounts[i] > 0))
      {
        MPI_Send((void *)&(sourceData[stride*sourceStart]),
                 stride*curSendCounts[i], mpiType, i, 0, MPI_COMM_WORLD);
      }
    }

    // Wait for all receives to complete
    if (requests.size() > 0)
    {
      std::vector<MPI_Status> statuses(requests.size());
      MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
    }
  } // moveData
  
}; // MPI_Bounding_Boxes

} // namespace Portage

#endif // MPI_Bounding_Boxes_H_
