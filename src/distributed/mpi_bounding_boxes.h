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
    std::vector<int> sendCounts(commSize, 0);
    std::vector<int> recvCounts(commSize); 
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

      sendCounts[i] = sendThis ? sourceNumCells : 0;
    }

    // Each rank will tell each other rank how many cells it is going to send it
    MPI_Alltoall(&(sendCounts[0]), 1, MPI_INT, &(recvCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Send how many cells are locally owned
    std::vector<int> sendOwnedCounts(commSize, 0);
    std::vector<int> recvOwnedCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendOwnedCounts[i] = (sendCounts[i] > 0) ? sourceNumOwnedCells : 0;
    MPI_Alltoall(&(sendOwnedCounts[0]), 1, MPI_INT, &(recvOwnedCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Compute the total number of source cells this rank will receive from all ranks
    int totalRecvSize = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvSize += recvCounts[i]; 
    int totalOwnedRecvSize = 0;
    for (unsigned int i=0; i<commSize; i++) totalOwnedRecvSize += recvOwnedCounts[i];

#ifdef DEBUG_MPI
    std::cout << "Received " << commRank << " " << recvCounts[0] << " " << recvCounts[1] << " " << totalRecvSize
              << " " << (totalRecvSize - recvCounts[commRank]) << std::endl;
#endif

    // Tell each other rank how many nodes it is going to send it
    std::vector<int> sendNodes(commSize, 0);
    std::vector<int> recvNodes(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendNodes[i] = (sendCounts[i] > 0) ? sourceNumNodes : 0;
    MPI_Alltoall(&(sendNodes[0]), 1, MPI_INT, &(recvNodes[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Send how many nodes are locally owned
    std::vector<int> sendOwnedNodes(commSize, 0);
    std::vector<int> recvOwnedNodes(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendOwnedNodes[i] = (sendCounts[i] > 0) ? sourceNumOwnedNodes : 0;
    MPI_Alltoall(&(sendOwnedNodes[0]), 1, MPI_INT, &(recvOwnedNodes[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Compute the total number of source nodes this rank will receive from all ranks
    int totalRecvNodes = 0;
    for (unsigned int i=0; i<commSize; i++) totalRecvNodes += recvNodes[i]; 
    int totalOwnedRecvNodes = 0;
    for (unsigned int i=0; i<commSize; i++) totalOwnedRecvNodes += recvOwnedNodes[i];

    // Compute the total number of neighbor indexes and owned neighbor indexes on this rank
    int neighborSize = 0;
    for (unsigned int i=0; i<sourceNeighborCounts.size(); i++)
      neighborSize += sourceNeighborCounts[i];
    int ownedNeighborSize = 0;
    for (unsigned int i=0; i<sourceNumOwnedCells; i++)
      ownedNeighborSize += sourceNeighborCounts[i];

    // Each rank will tell each other rank how many neighbor indexes it is going to send it
    std::vector<int> sendNeighborsCounts(commSize, 0);
    std::vector<int> recvNeighborsCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendNeighborsCounts[i] = (sendCounts[i] > 0) ? neighborSize : 0;
    MPI_Alltoall(&(sendNeighborsCounts[0]), 1, MPI_INT, &(recvNeighborsCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);
    int neighborsRecvSize = 0;
    for (unsigned int i=0; i<commSize; i++) neighborsRecvSize += recvNeighborsCounts[i];

    // Each rank will tell each other rank how many owned neighbor indexes it is going to send it
    std::vector<int> sendOwnedNeighborsCounts(commSize, 0);
    std::vector<int> recvOwnedNeighborsCounts(commSize);
    for (unsigned int i=0; i<commSize; i++)
      sendOwnedNeighborsCounts[i] = (sendOwnedCounts[i] > 0) ? ownedNeighborSize : 0;
    MPI_Alltoall(&(sendOwnedNeighborsCounts[0]), 1, MPI_INT, &(recvOwnedNeighborsCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);
    int ownedNeighborsRecvSize = 0;
    for (unsigned int i=0; i<commSize; i++) ownedNeighborsRecvSize += recvOwnedNeighborsCounts[i];

    // Data structures to hold mesh data received from other ranks
    std::vector<double> newCoords(dim*totalRecvNodes);
    std::vector<int> newNodeCounts(totalRecvSize);
    std::vector<int> newGlobalCellIds(totalRecvSize);
    std::vector<int> newNeighborCounts(totalRecvSize);
    std::vector<int> newNeighbors(neighborsRecvSize);

    // Perform two rounds of sends: the first for owned cells, and the second for ghost cells
    for (unsigned int r=0; r<2; r++)
    {
      // Compute all the offsets depending on whether this is the owned round or the ghost round
      std::vector<int> curRecvCounts(commSize);
      std::vector<int> curSendCounts(commSize);
      std::vector<int> curNodeRecvCounts(commSize);
      std::vector<int> curNodeSendCounts(commSize);
      std::vector<int> curNeighborRecvCounts(commSize);
      std::vector<int> curNeighborSendCounts(commSize);

      int sourceStart = 0;
      int sourceEnd = sourceNumOwnedCells;
      int localStart = 0;
      int localEnd = totalOwnedRecvSize;
      int sourceNodeStart = 0;
      int sourceNodeEnd = sourceNumOwnedNodes;
      int localNodeStart = 0;
      int localNodeEnd = totalOwnedRecvNodes;
      int sourceNeighborStart = 0;
      int sourceNeighborEnd = ownedNeighborSize;
      int localNeighborStart = 0;
      int localNeighborEnd = ownedNeighborsRecvSize;

      if (r == 1)
      {
        sourceStart = sourceNumOwnedCells;
        sourceEnd = sourceNumCells;
        localStart = totalOwnedRecvSize;
        localEnd = totalRecvSize;
        sourceNodeStart = sourceNumOwnedNodes;
        sourceNodeEnd = sourceNumNodes;
        localNodeStart = totalOwnedRecvNodes;
        localNodeEnd = totalRecvNodes;
        sourceNeighborStart = ownedNeighborSize;
        sourceNeighborEnd = neighborSize;
        localNeighborStart = ownedNeighborsRecvSize;
        localNeighborEnd = neighborsRecvSize;
      }

      for (unsigned int i=0; i<commSize; i++)
      {
        if (r == 0)
        {
          curRecvCounts[i] = recvOwnedCounts[i];
          curSendCounts[i] = sendOwnedCounts[i];
          curNodeRecvCounts[i] = recvOwnedNodes[i];
          curNodeSendCounts[i] = sendOwnedNodes[i];
          curNeighborRecvCounts[i] = recvOwnedNeighborsCounts[i];
          curNeighborSendCounts[i] = sendOwnedNeighborsCounts[i];
        }
        else 
        {
          curRecvCounts[i] = recvCounts[i] - recvOwnedCounts[i];
          curSendCounts[i] = sendCounts[i] - sendOwnedCounts[i];
          curNodeRecvCounts[i] = recvNodes[i] - recvOwnedNodes[i];
          curNodeSendCounts[i] = sendNodes[i] - sendOwnedNodes[i];
          curNeighborRecvCounts[i] = recvNeighborsCounts[i] - recvOwnedNeighborsCounts[i];
          curNeighborSendCounts[i] = sendNeighborsCounts[i] - sendOwnedNeighborsCounts[i];
        }
      }

      int localOffset = localStart;
      for (unsigned int i=0; i<commRank; i++) localOffset += curRecvCounts[i]; 

      int localNodeOffset = localNodeStart;
      for (unsigned int i=0; i<commRank; i++) localNodeOffset += curNodeRecvCounts[i]; 

      int localNeighborsOffset = localNeighborStart;
      for (unsigned int i=0; i<commRank; i++) localNeighborsOffset += curNeighborRecvCounts[i];
      
      // SEND NUMBER OF NODES FOR EACH CELL 

      if (curRecvCounts[commRank] > 0)
        std::copy(sourceNodeCounts.begin()+sourceStart, sourceNodeCounts.begin()+sourceEnd, 
                  newNodeCounts.begin() + localOffset);

      // Each rank will do a non-blocking receive from each rank from which it will receive node counts
      MPI_Status stat;
      int writeOffset = localStart;
      std::vector<MPI_Request> requests;
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curRecvCounts[i] > 0))
        {
          MPI_Request request;
          MPI_Irecv(&(newNodeCounts[0])+writeOffset, curRecvCounts[i], MPI_INT, i,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &request);
          requests.push_back(request);
        }
        writeOffset += curRecvCounts[i];
      }

      // Each rank will send its node counts to appropriate ranks
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curSendCounts[i] > 0))
        {
          MPI_Send(&(sourceNodeCounts[sourceStart]), curSendCounts[i], MPI_INT, i, 0, MPI_COMM_WORLD);
        }
      }

      // Wait for all receives to complete
      if (requests.size() > 0)
      {
        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
      }


      // SEND NODE COORDINATES

      // Copy source nodes that will stay on this rank into the proper place in the new vector
      if (curRecvCounts[commRank] > 0)
        std::copy(sourceCoords.begin()+dim*sourceNodeStart, sourceCoords.begin()+dim*sourceNodeEnd, 
                  newCoords.begin() + dim*localNodeOffset);

      // Each rank will do a non-blocking receive from each rank from which it will receive source node coordinates
      writeOffset = dim * localNodeStart;
      requests.clear();
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curRecvCounts[i] > 0))
        {
          MPI_Request request;
          MPI_Irecv(&(newCoords[0])+writeOffset, dim*curNodeRecvCounts[i], MPI_DOUBLE, i, 
                    MPI_ANY_TAG, MPI_COMM_WORLD, &request);
          requests.push_back(request);
        }
        writeOffset += dim * curNodeRecvCounts[i];
      }

      // Each rank will send its source node coordinates to appropriate ranks
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curSendCounts[i] > 0))
        {
          MPI_Send(&(sourceCoords[dim*sourceNodeStart]), dim*curNodeSendCounts[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
      }
    
      // Wait for all receives to complete
      if (requests.size() > 0)
      {
        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));	
      }
    

      // SEND GLOBAL CELL IDS

      // Copy source cells that will stay on this rank into the proper place in the new vector
      if (curRecvCounts[commRank] > 0)
        std::copy(sourceGlobalCellIds.begin()+sourceStart, sourceGlobalCellIds.begin()+sourceEnd, 
                  newGlobalCellIds.begin() + localOffset);

      // Each rank will do a non-blocking receive from each rank from which it will receive global cell ids
      writeOffset = localStart;
      requests.clear();
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curRecvCounts[i] > 0))
        {
          MPI_Request request;
          MPI_Irecv(&(newGlobalCellIds[0])+writeOffset, curRecvCounts[i], MPI_INT, i,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &request);
          requests.push_back(request);
        }
        writeOffset += curRecvCounts[i];
      }

      // Each rank will send its global cell ids to appropriate ranks
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curSendCounts[i] > 0))
        {
          MPI_Send(&(sourceGlobalCellIds[sourceStart]), curSendCounts[i], MPI_INT, i, 0, MPI_COMM_WORLD);
        }
      }

      // Wait for all receives to complete
      if (requests.size() > 0)
      {
        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
      }


      // SEND NUMBER OF NEIGHBORS FOR EACH CELL 
 
      // Copy source cells that will stay on this rank into the proper place in the new vector
      if (curRecvCounts[commRank] > 0)
        std::copy(sourceNeighborCounts.begin()+sourceStart, sourceNeighborCounts.begin()+sourceEnd, 
                  newNeighborCounts.begin() + localOffset);

      // Each rank will do a non-blocking receive from each rank from which it will receive neighbor counts
      writeOffset = localStart;
      requests.clear();
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curRecvCounts[i] > 0))
        {
          MPI_Request request;
          MPI_Irecv(&(newNeighborCounts[0])+writeOffset, curRecvCounts[i], MPI_INT, i,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &request);
          requests.push_back(request);
        }
        writeOffset += curRecvCounts[i];
      }

      // Each rank will send its neighbor counts to appropriate ranks
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curSendCounts[i] > 0))
        {
          MPI_Send(&(sourceNeighborCounts[sourceStart]), curSendCounts[i], MPI_INT, i, 0, MPI_COMM_WORLD);
        }
      }

      // Wait for all receives to complete
      if (requests.size() > 0)
      {
        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
      }


      // SEND LIST OF NEIGHBOR GLOBAL IDS FOR EACH CELL

      // Copy neighbor ids that will stay on this rank into the proper place in the new vector
      if (curNeighborRecvCounts[commRank] > 0)
        std::copy(sourceNeighbors.begin()+sourceNeighborStart, sourceNeighbors.begin()+sourceNeighborEnd, 
                  newNeighbors.begin() + localNeighborsOffset);

      // Each rank will do a non-blocking receive from each rank from which it will receive neighbor ids
      writeOffset = localNeighborStart;  
      requests.clear();
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curNeighborRecvCounts[i] > 0))
        {
          MPI_Request request;
          MPI_Irecv(&(newNeighbors[0])+writeOffset, curNeighborRecvCounts[i], MPI_INT, i,
                    MPI_ANY_TAG, MPI_COMM_WORLD, &request);
          requests.push_back(request);
        }
        writeOffset += curNeighborRecvCounts[i];
      }

      // Each rank will send its neighbor global ids to appropriate ranks
      for (unsigned int i=0; i<commSize; i++)
      {
        if ((i != commRank) && (curNeighborSendCounts[i] > 0))
        {
          MPI_Send(&(sourceNeighbors[sourceNeighborStart]), curNeighborSendCounts[i], MPI_INT, i, 0, MPI_COMM_WORLD);
        }
      }

      // Wait for all receives to complete
      if (requests.size() > 0)
      {
        std::vector<MPI_Status> statuses(requests.size());
        MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
      }

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
      std::vector<double> newField(totalRecvSize);

      // Perform two rounds of sends: the first for owned cells, and the second for ghost cells
      for (unsigned int r=0; r<2; r++)
      {
        // Set up the offsets depending on whether this is the owned round or the ghost round
        std::vector<int> curRecvCounts(commSize);
        std::vector<int> curSendCounts(commSize);

        int sourceStart = 0;
        int sourceEnd = sourceNumOwnedCells;
        int localStart = 0;
        int localEnd = totalOwnedRecvSize;

        if (r == 1)
        {
          sourceStart = sourceNumOwnedCells;
          sourceEnd = sourceNumCells;
          localStart = totalOwnedRecvSize;
          localEnd = totalRecvSize;
        }

        for (unsigned int i=0; i<commSize; i++)
        {
          if (r == 0)
          {
            curRecvCounts[i] = recvOwnedCounts[i];
            curSendCounts[i] = sendOwnedCounts[i];
          }
          else
          {
            curRecvCounts[i] = recvCounts[i] - recvOwnedCounts[i];
            curSendCounts[i] = sendCounts[i] - sendOwnedCounts[i];
          }
        }

        int localOffset = localStart;
        for (unsigned int i=0; i<commRank; i++) localOffset += curRecvCounts[i];

        std::shared_ptr<std::vector<double>> sourceState = source_state_flat.get_vector(s);
        int sourceFieldStride = source_state_flat.get_field_dim(s);

        // Copy field values that will stay on this rank into the proper place in the new vector
        if (curRecvCounts[commRank] > 0)
          std::copy(sourceState->begin()+sourceStart, sourceState->begin()+sourceEnd, 
                    newField.begin() + sourceFieldStride*localOffset);
        int writeOffset = localStart;

        // Each rank will do a non-blocking receive from each rank from which it will receive source state
        std::vector<MPI_Request> requests; 
        for (unsigned int i=0; i<commSize; i++)
        {
          if ((i != commRank) && (curRecvCounts[i] > 0))
          {
            MPI_Request request;
            MPI_Irecv(&(newField[0])+writeOffset, sourceFieldStride*curRecvCounts[i], MPI_DOUBLE, i,
                      MPI_ANY_TAG, MPI_COMM_WORLD, &request);
            requests.push_back(request);
          }
          writeOffset += sourceFieldStride*curRecvCounts[i];
        }

        // Each rank will send its source fields to appropriate ranks
        for (unsigned int i=0; i<commSize; i++)
        {
          if ((i != commRank) && (curSendCounts[i] > 0))
          {
            MPI_Send(&((*sourceState)[sourceStart]), sourceFieldStride*curSendCounts[i], MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
          }
        }

        // Wait for all receives to complete
        if (requests.size() > 0)
        {
          std::vector<MPI_Status> statuses(requests.size());
          MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
        }
    
#ifdef DEBUG_MPI
        if ((r == 1) && (commRank == 1))
        {
          std::cout << "Test: field " << newField.size() << std::endl;
          for (unsigned int f=0; f<newField.size(); f++)
            std::cout << newField[f] << " ";
          std::cout << std::endl << std::endl;
        }
#endif
      }

      // We will now use the received source state as our new source state on this partition
      std::shared_ptr<std::vector<double>> sourceState = source_state_flat.get_vector(s);
      sourceState->resize(newField.size());
      std::copy(newField.begin(), newField.end(), sourceState->begin());
    }

    // We will now use the received source mesh data as our new source mesh on this partition
    sourceCoords = newCoords;
    sourceNodeCounts = newNodeCounts;
    sourceGlobalCellIds = newGlobalCellIds;
    sourceNeighborCounts = newNeighborCounts;
    sourceNeighbors = newNeighbors;
    source_mesh_flat.set_num_owned_cells(totalOwnedRecvSize);
    source_mesh_flat.set_num_owned_nodes(totalOwnedRecvNodes);

    // Create all index maps
    source_mesh_flat.make_index_maps();

#ifdef DEBUG_MPI
    if (commRank == 1)
    {
      std::cout << "Sizes: " << commRank << " " << totalRecvSize << " " << targetNumCells
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

  }
  
}; // MPI_Bounding_Boxes

} // namespace Portage

#endif // MPI_Bounding_Boxes_H_
