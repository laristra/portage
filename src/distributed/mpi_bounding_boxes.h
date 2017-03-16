/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



#ifndef MPI_BOUNDING_BOXES_H_
#define MPI_BOUNDING_BOXES_H_

//#define DEBUG_MPI

#include <cassert>
#include <algorithm>
#include <numeric>
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
    @brief Helper structure containg comms info for a given entity type
   */
  struct comm_info_t {
    //!< Number of total/owned entities in source field
    int sourceNum = 0, sourceNumOwned = 0;
    //!< Number of total/owned entities in new field
    int newNum = 0, newNumOwned = 0;
    //! Array of total/owned send sizes from me to all PEs
    std::vector<int> sendCounts, sendOwnedCounts;
    //! Array of total/owned recv sizes to me from all PEs
    std::vector<int> recvCounts, recvOwnedCounts;
  };


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
    int targetNumOwnedCells = target_mesh.num_owned_cells();
    std::vector<double> targetBoundingBoxes(2*dim*commSize);
    for (unsigned int i=0; i<2*dim; i+=2)
    {
      targetBoundingBoxes[2*dim*commRank+i+0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2*dim*commRank+i+1] = -std::numeric_limits<double>::max();
    }
    for (unsigned int c=0; c<targetNumOwnedCells; ++c)
    {
      std::vector<int> nodes;
      target_mesh.cell_get_nodes(c, &nodes);
      int cellNumNodes = nodes.size();
      for (unsigned int j=0; j<cellNumNodes; ++j)
      {
        int n = nodes[j];
        // ugly hack, since dim is not known at compile time
        if (dim == 3)
        {
          Portage::Point<3> nodeCoord;
          target_mesh.node_get_coordinates(n, &nodeCoord);
          for (unsigned int k=0; k<dim; ++k)
          {
            if (nodeCoord[k] < targetBoundingBoxes[2*dim*commRank+2*k])
              targetBoundingBoxes[2*dim*commRank+2*k] = nodeCoord[k];
            if (nodeCoord[k] > targetBoundingBoxes[2*dim*commRank+2*k+1])
              targetBoundingBoxes[2*dim*commRank+2*k+1] = nodeCoord[k];
          }
        }
        else if (dim == 2)
        {
          Portage::Point<2> nodeCoord;
          target_mesh.node_get_coordinates(n, &nodeCoord);
          for (unsigned int k=0; k<dim; ++k)
          {
            if (nodeCoord[k] < targetBoundingBoxes[2*dim*commRank+2*k])
              targetBoundingBoxes[2*dim*commRank+2*k] = nodeCoord[k];
            if (nodeCoord[k] > targetBoundingBoxes[2*dim*commRank+2*k+1])
              targetBoundingBoxes[2*dim*commRank+2*k+1] = nodeCoord[k];
          }
        }
      } // for j
    } // for c

    // Get the source mesh properties and cordinates
    int sourceNumOwnedCells = source_mesh_flat.num_owned_cells();
    int sourceNumCells = source_mesh_flat.num_owned_cells() + source_mesh_flat.num_ghost_cells();
    int sourceNumOwnedNodes = source_mesh_flat.num_owned_nodes();
    int sourceNumNodes = source_mesh_flat.num_owned_nodes() + source_mesh_flat.num_ghost_nodes();
    // only used in 3D:
    int sourceNumOwnedFaces = -1;
    int sourceNumFaces = -1;
    if (dim == 3)
    {
      sourceNumOwnedFaces = source_mesh_flat.num_owned_faces();
      sourceNumFaces = source_mesh_flat.num_owned_faces() + source_mesh_flat.num_ghost_faces();
    }
    std::vector<double>& sourceCoords = source_mesh_flat.get_coords();
    std::vector<int>& sourceNodeCounts = source_mesh_flat.get_node_counts();
    std::vector<int>& sourceCellGlobalIds = source_mesh_flat.cellGlobalIds_;
    std::vector<int>& sourceNodeGlobalIds = source_mesh_flat.nodeGlobalIds_;

    // Compute the bounding box for the source mesh on this rank
    std::vector<double> sourceBoundingBox(2*dim);
    for (unsigned int i=0; i<2*dim; i+=2)
    {
      sourceBoundingBox[i+0] = std::numeric_limits<double>::max();
      sourceBoundingBox[i+1] = -std::numeric_limits<double>::max();
    }
    for (unsigned int c=0; c<sourceNumOwnedCells; ++c)
    {
      std::vector<int> nodes;
      source_mesh_flat.cell_get_nodes(c, &nodes);
      int cellNumNodes = nodes.size();
      for (unsigned int j=0; j<cellNumNodes; ++j)
      {
        int n = nodes[j];
        for (unsigned int k=0; k<dim; ++k)
        {
          if (sourceCoords[n*dim+k] < sourceBoundingBox[2*k])
            sourceBoundingBox[2*k] = sourceCoords[n*dim+k];
          if (sourceCoords[n*dim+k] > sourceBoundingBox[2*k+1])
            sourceBoundingBox[2*k+1] = sourceCoords[n*dim+k];
        }
      } // for j
    } // for c

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
    for (unsigned int k=0; k<dim; ++k)
    {
      min2[k] = sourceBoundingBox[2*k]+boxOffset;
      max2[k] = sourceBoundingBox[2*k+1]-boxOffset;
    }

    // For each target rank with a bounding box that overlaps the bounding box for this rank's partition
    // of the source mesh, we will send it all our source cells; otherwise, we will send it nothing
    std::vector<bool> sendFlags(commSize);
    for (unsigned int i=0; i<commSize; ++i)
    {
      double min1[dim], max1[dim];
      bool sendThis = true;
      for (unsigned int k=0; k<dim; ++k)
      {
        min1[k] = targetBoundingBoxes[2*dim*i+2*k]+boxOffset;
        max1[k] = targetBoundingBoxes[2*dim*i+2*k+1]-boxOffset;
        sendThis = sendThis &&
            ((min1[k] <= min2[k] && min2[k] <= max1[k]) ||
             (min2[k] <= min1[k] && min1[k] <= max2[k]));
      }

      sendFlags[i] = sendThis;
    }

    comm_info_t cellInfo, nodeInfo, cellToNodeInfo;
    // only used in 3D:
    comm_info_t faceInfo, cellToFaceInfo, faceToNodeInfo;

    setInfo(&cellInfo, commSize, sendFlags,
            sourceNumCells, sourceNumOwnedCells);

#ifdef DEBUG_MPI
    std::cout << "Received " << commRank << " " << cellInfo.recvCounts[0] << " " << cellInfo.recvCounts[1] << " " << cellInfo.newNum
              << " " << (cellInfo.newNum - cellInfo.recvCounts[commRank]) << std::endl;
#endif

    setInfo(&nodeInfo, commSize, sendFlags,
            sourceNumNodes, sourceNumOwnedNodes);

    // TODO:  in 3D, don't communicate this - compute it from
    // cell-to-face and face-to-node
    int sizeCellToNodeList = source_mesh_flat.cellToNodeList_.size();
    int sizeOwnedCellToNodeList = (
        sourceNumCells == sourceNumOwnedCells ? sizeCellToNodeList :
        source_mesh_flat.cellNodeOffsets_[sourceNumOwnedCells]);
    setInfo(&cellToNodeInfo, commSize, sendFlags,
            sizeCellToNodeList, sizeOwnedCellToNodeList);

    if (dim == 3)
    {
      setInfo(&faceInfo, commSize, sendFlags,
              sourceNumFaces, sourceNumOwnedFaces);

      int sizeCellToFaceList = source_mesh_flat.cellToFaceList_.size();
      int sizeOwnedCellToFaceList = (
          sourceNumCells == sourceNumOwnedCells ? sizeCellToFaceList :
          source_mesh_flat.cellFaceOffsets_[sourceNumOwnedCells]);

      setInfo(&cellToFaceInfo, commSize, sendFlags,
              sizeCellToFaceList, sizeOwnedCellToFaceList);

      int sizeFaceToNodeList = source_mesh_flat.faceToNodeList_.size();
      int sizeOwnedFaceToNodeList = (
          sourceNumFaces == sourceNumOwnedFaces ? sizeFaceToNodeList :
          source_mesh_flat.faceNodeOffsets_[sourceNumOwnedFaces]);

      setInfo(&faceToNodeInfo, commSize, sendFlags,
              sizeFaceToNodeList, sizeOwnedFaceToNodeList);
    }

    // Data structures to hold mesh data received from other ranks
    std::vector<double> newCoords(dim*nodeInfo.newNum);
    std::vector<int> newCellNodeCounts(cellInfo.newNum);
    std::vector<int> newCellToNodeList(cellToNodeInfo.newNum);
    std::vector<int> newCellFaceCounts;
    std::vector<int> newCellToFaceList;
    std::vector<bool> newCellToFaceDirs;
    std::vector<int> newFaceNodeCounts;
    std::vector<int> newFaceToNodeList;
    if (dim == 3)
    {
      newCellFaceCounts.resize(cellInfo.newNum);
      newCellToFaceList.resize(cellToFaceInfo.newNum);
      newCellToFaceDirs.resize(cellToFaceInfo.newNum);
      newFaceNodeCounts.resize(faceInfo.newNum);
      newFaceToNodeList.resize(faceToNodeInfo.newNum);
    }
    std::vector<int> newCellGlobalIds(cellInfo.newNum);
    std::vector<int> newNodeGlobalIds(nodeInfo.newNum);

    // SEND NUMBER OF NODES FOR EACH CELL

    moveField(cellInfo, commRank, commSize, MPI_INT, 1,
              sourceNodeCounts, &newCellNodeCounts);

    // SEND CELL-TO-NODE MAP

    moveField(cellToNodeInfo, commRank, commSize, MPI_INT, 1,
              source_mesh_flat.cellToNodeList_, &newCellToNodeList);

    // SEND NODE COORDINATES

    moveField(nodeInfo, commRank, commSize, MPI_DOUBLE, dim,
              sourceCoords, &newCoords);

    if (dim == 3)
    {
      // SEND NUMBER OF FACES FOR EACH CELL

      moveField(cellInfo, commRank, commSize, MPI_INT, 1,
                source_mesh_flat.cellFaceCounts_, &newCellFaceCounts);

      // SEND CELL-TO-FACE MAP
      // For this array only, pack up face IDs + dirs and send together
      int size = source_mesh_flat.cellToFaceList_.size();
      for (unsigned int j=0; j<size; ++j)
      {
        int f = source_mesh_flat.cellToFaceList_[j];
        int dir = static_cast<int>(source_mesh_flat.cellToFaceDirs_[j]);
        source_mesh_flat.cellToFaceList_[j] = (f << 1) | dir;
      }
      moveField(cellToFaceInfo, commRank, commSize, MPI_INT, 1,
                source_mesh_flat.cellToFaceList_, &newCellToFaceList);
      // Unpack face IDs and dirs
      for (unsigned int j=0; j<newCellToFaceList.size(); ++j)
      {
        int fd = newCellToFaceList[j];
        newCellToFaceList[j] = fd >> 1;
        newCellToFaceDirs[j] = fd & 1;
      }

      // SEND NUMBER OF NODES FOR EACH FACE

      moveField(faceInfo, commRank, commSize, MPI_INT, 1,
                source_mesh_flat.faceNodeCounts_, &newFaceNodeCounts);

      // SEND FACE-TO-NODE MAP

      moveField(faceToNodeInfo, commRank, commSize, MPI_INT, 1,
                source_mesh_flat.faceToNodeList_, &newFaceToNodeList);
    }

    // SEND GLOBAL CELL IDS

    moveField(cellInfo, commRank, commSize, MPI_INT, 1,
              sourceCellGlobalIds, &newCellGlobalIds);

    // SEND GLOBAL NODE IDS

    moveField(nodeInfo, commRank, commSize, MPI_INT, 1,
              sourceNodeGlobalIds, &newNodeGlobalIds);

    // SEND FIELD VALUES

    // Send and receive each field to be remapped
    for (int s=0; s<source_state_flat.get_num_vectors(); s++)
    {
      std::shared_ptr<std::vector<double>> sourceField = source_state_flat.get_vector(s);
      int sourceFieldStride = source_state_flat.get_field_stride(s);

      // Currently only cell and node fields are supported
      comm_info_t& info = (source_state_flat.get_entity(s) == NODE ?
                           nodeInfo : cellInfo);
      std::vector<double> newField(info.newNum);

      moveField(info, commRank, commSize,
                MPI_DOUBLE, sourceFieldStride,
                *sourceField, &newField);

#ifdef DEBUG_MPI
      if ((r == 1) && (commRank == 1))
      {
        std::cout << "Test: field " << newField.size() << std::endl;
        for (unsigned int f=0; f<newField.size(); f++)
          std::cout << newField[f] << " ";
        std::cout << std::endl << std::endl;
      }
#endif

      // We will now use the received source state as our new source state on this partition
      sourceField->resize(newField.size());
      std::copy(newField.begin(), newField.end(), sourceField->begin());

    } // for s

    // We will now use the received source mesh data as our new source mesh on this partition
    // TODO:  replace with swap
    sourceCoords = newCoords;
    sourceNodeCounts = newCellNodeCounts;
    sourceCellGlobalIds = newCellGlobalIds;
    std::swap(source_mesh_flat.nodeGlobalIds_, newNodeGlobalIds);
    source_mesh_flat.set_num_owned_cells(cellInfo.newNumOwned);
    source_mesh_flat.set_num_owned_nodes(nodeInfo.newNumOwned);

    fixListIndices(cellToNodeInfo, nodeInfo, commSize, &newCellToNodeList);
    std::swap(source_mesh_flat.cellToNodeList_, newCellToNodeList);

    if (dim == 3)
    {
      std::swap(source_mesh_flat.cellFaceCounts_, newCellFaceCounts);
      std::swap(source_mesh_flat.faceNodeCounts_, newFaceNodeCounts);

      source_mesh_flat.numOwnedFaces_ = faceInfo.newNumOwned;

      fixListIndices(cellToFaceInfo, faceInfo, commSize, &newCellToFaceList);
      std::swap(source_mesh_flat.cellToFaceList_, newCellToFaceList);
      std::swap(source_mesh_flat.cellToFaceDirs_, newCellToFaceDirs);

      fixListIndices(faceToNodeInfo, nodeInfo, commSize, &newFaceToNodeList);
      std::swap(source_mesh_flat.faceToNodeList_, newFaceToNodeList);
    }

    // Finish initialization using redistributed data
    source_mesh_flat.finish_init();

#ifdef DEBUG_MPI
    if (commRank == 1)
    {
      std::cout << "Sizes: " << commRank << " " << cellInfo.newNum << " " << targetNumOwnedCells
                << " " << cellInfo.sourceNum << " " << sourceCoords.size() << std::endl;

      for (unsigned int i=0; i<sourceCellGlobalIds.size(); i++)
        std::cout << sourceCellGlobalIds[i] << " ";
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
    @brief Compute fields needed to do comms for a given entity type
    @param[in] info              Info data structure to be filled
    @param[in] commSize          Total number of MPI ranks
    @param[in] sendFlags         Array of flags:  do I send to PE n?
    @param[in] sourceNum         Number of entities (total) on this rank
    @param[in] sourceNumOwned    Number of owned entities on this rank
   */
  void setInfo(comm_info_t* info,
               const int commSize,
               const std::vector<bool>& sendFlags,
               const int sourceNum,
               const int sourceNumOwned)
  {
    // Set my counts of all entities and owned entities
    info->sourceNum = sourceNum;
    info->sourceNumOwned = sourceNumOwned;

    // Each rank will tell each other rank how many indexes it is going to send it
    info->sendCounts.resize(commSize);
    info->recvCounts.resize(commSize);
    for (unsigned int i=0; i<commSize; i++)
      info->sendCounts[i] = sendFlags[i] ? info->sourceNum : 0;
    MPI_Alltoall(&(info->sendCounts[0]), 1, MPI_INT,
                 &(info->recvCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Each rank will tell each other rank how many owned indexes it is going to send it
    info->sendOwnedCounts.resize(commSize);
    info->recvOwnedCounts.resize(commSize);
    for (unsigned int i=0; i<commSize; i++)
      info->sendOwnedCounts[i] = sendFlags[i] ? info->sourceNumOwned : 0;
    MPI_Alltoall(&(info->sendOwnedCounts[0]), 1, MPI_INT,
                 &(info->recvOwnedCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);

    // Compute the total number of indexes this rank will receive from all ranks
    for (unsigned int i=0; i<commSize; i++)
      info->newNum += info->recvCounts[i];
    for (unsigned int i=0; i<commSize; i++)
      info->newNumOwned += info->recvOwnedCounts[i];
  } // setInfo


  /*!
    @brief Move values for a single data field to all ranks as needed
    @tparam[in] T                C++ type of data to be moved
    @param[in] info              Info struct for entity type of field
    @param[in] commRank          MPI rank of this PE
    @param[in] commSize          Total number of MPI ranks
    @param[in] mpiType           MPI type of data (MPI_???) to be moved
    @param[in] stride            Stride of data field
    @param[in] sourceData        Array of (old) source data
    @param[in] newData           Array of new source data
   */
  template<typename T>
  void moveField(const comm_info_t& info,
                 const int commRank, const int commSize,
                 const MPI_Datatype mpiType, const int stride,
                 const std::vector<T>& sourceData,
                 std::vector<T>* newData)
  {
    // Perform two rounds of sends: the first for owned cells, and the second for ghost cells
    std::vector<int> recvGhostCounts(commSize);
    std::vector<int> sendGhostCounts(commSize);

    for (unsigned int i=0; i<commSize; i++)
    {
      recvGhostCounts[i] = info.recvCounts[i] - info.recvOwnedCounts[i];
      sendGhostCounts[i] = info.sendCounts[i] - info.sendOwnedCounts[i];
    }

    moveData(commRank, commSize, mpiType, stride,
             0, info.sourceNumOwned,
             0,
             info.sendOwnedCounts, info.recvOwnedCounts,
             sourceData, newData);
    moveData(commRank, commSize, mpiType, stride,
             info.sourceNumOwned, info.sourceNum,
             info.newNumOwned,
             sendGhostCounts, recvGhostCounts,
             sourceData, newData);

#ifdef DEBUG_MPI
    if (commRank == 1)
    {
      std::cout << "Test: field " << newData.size() << std::endl;
      for (unsigned int f=0; f<newData.size(); f++)
        std::cout << newData[f] << " ";
      std::cout << std::endl << std::endl;
    }
#endif
  } // moveField


  /*!
    @brief Move values for a single range of data to all ranks as needed
    @tparam[in] T                C++ type of data to be moved
    @param[in] commRank          MPI rank of this PE
    @param[in] commSize          Total number of MPI ranks
    @param[in] mpiType           MPI type of data (MPI_???) to be moved
    @param[in] stride            Stride of data field
    @param[in] sourceStart       Start location in source (send) data
    @param[in] sourceEnd         End location in source (send) data
    @param[in] newStart          Start location in new (recv) data
    @param[in] curSendCounts     Array of send sizes from me to all PEs
    @param[in] curRecvCounts     Array of recv sizes to me from all PEs
    @param[in] sourceData        Array of (old) source data
    @param[in] newData           Array of new source data
   */
  template<typename T>
  void moveData(const int commRank, const int commSize,
                const MPI_Datatype mpiType, const int stride,
                const int sourceStart, const int sourceEnd,
                const int newStart,
                const std::vector<int>& curSendCounts,
                const std::vector<int>& curRecvCounts,
                const std::vector<T>& sourceData,
                std::vector<T>* newData)
  {
    // Each rank will do a non-blocking receive from each rank from
    // which it will receive data values
    int writeOffset = newStart;
    int myOffset;
    std::vector<MPI_Request> requests;
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (curRecvCounts[i] > 0))
      {
        MPI_Request request;
        MPI_Irecv((void *)&((*newData)[stride*writeOffset]),
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
                newData->begin() + stride*myOffset);

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

  //! Correct a map to account for concatenated lists
  void fixListIndices(const comm_info_t& mapInfo,
              const comm_info_t& rangeInfo,
              const int commSize,
              std::vector<int>* newMap)
  {
    // compute corrections for each rank
    std::vector<int> ownedOffsets(commSize), ghostOffsets(commSize);
    ownedOffsets[0] = 0;
    std::partial_sum(rangeInfo.recvOwnedCounts.begin(),
                     rangeInfo.recvOwnedCounts.end()-1,
                     ownedOffsets.begin()+1);
    std::vector<int> recvGhostCounts(commSize);
    for (unsigned int i=0; i<commSize; ++i)
      recvGhostCounts[i] =
          rangeInfo.recvCounts[i] - rangeInfo.recvOwnedCounts[i];
    ghostOffsets[0] = 0;
    std::partial_sum(recvGhostCounts.begin(), recvGhostCounts.end()-1,
                     ghostOffsets.begin()+1);
    for (unsigned int i=0; i<commSize; ++i)
      ghostOffsets[i] += rangeInfo.newNumOwned - rangeInfo.recvOwnedCounts[i];

    // correct owned entities, one rank at a time
    int base = 0;
    for (unsigned int i=0; i<commSize; ++i) {
      int ownedCount = mapInfo.recvOwnedCounts[i];
      int ownedNodeCount = rangeInfo.recvOwnedCounts[i];
      for (int j=0; j<ownedCount; ++j) {
        int n = (*newMap)[base + j];
        (*newMap)[base + j] +=
            (n < ownedNodeCount ? ownedOffsets[i] : ghostOffsets[i]);
      }
      base += ownedCount;
    } // for i

    // correct ghost entities, one rank at a time
    for (unsigned int i=0; i<commSize; ++i) {
      int ownedCount = mapInfo.recvOwnedCounts[i];
      int ghostCount = mapInfo.recvCounts[i] - ownedCount;
      int ownedNodeCount = rangeInfo.recvOwnedCounts[i];
      for (int j=0; j<ghostCount; ++j) {
        int n = (*newMap)[base + j];
        (*newMap)[base + j] +=
            (n < ownedNodeCount ? ownedOffsets[i] : ghostOffsets[i]);
      }
      base += ghostCount;
    } // for i

  } // fixListIndices

}; // MPI_Bounding_Boxes

} // namespace Portage

#endif // MPI_Bounding_Boxes_H_
