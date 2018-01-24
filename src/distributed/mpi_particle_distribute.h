/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef MPI_PARTICLE_DISTRIBUTE_H_
#define MPI_PARTICLE_DISTRIBUTE_H_

//#define DEBUG_MPI

#include <cassert>
#include <algorithm>
#include <numeric>
#include <memory>

#include "portage/accumulate/accumulate.h"
#include "portage/support/portage.h"

#include "mpi.h"

/*!
  @file mpi_bounding_boxes.h
  @brief Distributes source data using MPI based on bounding boxes
 */

namespace Portage {

/*!
  @class MPI_Particle_Distribute
  @brief Distributes source data using MPI based on bounding boxes

         Currently assumes coordinates and all fields are doubles.
*/
template<size_t dim> 
class MPI_Particle_Distribute {
 public:

  /*!
    @brief Constructor of MPI_Particle_Distribute
   */
  MPI_Particle_Distribute() {}


  /*!
    @brief Destructor of MPI_Particle_Distribute
   */
  ~MPI_Particle_Distribute() {};


  /*!
    @brief Helper structure containg comms info for a given entity type
   */
  struct comm_info_t {
    //!< Number of entities in source swarm to be sent to each rank
    std::vector<int> sourceNum;
    //!< Number of total entities in new field
    int newNum = 0;
    //! Array of send sizes from me to all PEs
    std::vector<int> sendCounts;
    //! Array of total recv sizes to me from all PEs
    std::vector<int> recvCounts;
  };


  /*!
    @brief Compute bounding boxes for all partitions, and send source mesh and state
           information to all target partitions with an overlapping bounding box using MPI
    @param[in] source_mesh_flat  Input mesh (must be flat representation)
    @param[in] source_state_flat Input state (must be flat representation)
    @param[in] target_mesh       Target mesh
    @param[in] target_state      Target state (not actually used for now)
   */
  template <class SourceSwarm, class SourceState, class TargetSwarm, class TargetState>
  void distribute(SourceSwarm &source_swarm, SourceState &source_state,
                  TargetSwarm &target_swarm, TargetState &target_state,
                  std::vector<std::vector<std::vector<double>>>& smoothing_lengths,
                  Meshfree::WeightCenter center=Meshfree::WeightCenter::Gather)
  {
    // Get the MPI communicator size and rank information
    int commSize, commRank;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);

    assert(dim == source_swarm.space_dimension());
    assert(dim == target_swarm.space_dimension());
    assert(center == Meshfree::WeightCenter::Gather);

    /************************************************************************** 
    * Step 1: Compute bounding box for target swarm based on weight center    *
    *         for the current rank, and put it in a vector that will be       *
    *         used later to store the target bounding box for each rank       *
    **************************************************************************/

    size_t targetNumOwnedPts = target_swarm.num_particles(PARALLEL_OWNED);
    
    std::vector<double> targetBoundingBoxes(2*dim*commSize);

    for (size_t i=0; i<2*dim; i+=2)
    {
      targetBoundingBoxes[2*dim*commRank+i+0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2*dim*commRank+i+1] = -std::numeric_limits<double>::max();
    }

    for (size_t c=0; c<targetNumOwnedPts; ++c)
    {
      Point<dim> coord = target_swarm.get_particle_coordinates(c);
      Point<dim> ext = Point<dim>(smoothing_lengths[c][0]);
     
      for (size_t k=0; k < dim; ++k)
      {
        double val0 = coord[k]-ext[k];
        double val1 = coord[k]+ext[k];

        if (val0 < targetBoundingBoxes[2*dim*commRank+2*k])
           targetBoundingBoxes[2*dim*commRank+2*k] = val0;

        if (val1 > targetBoundingBoxes[2*dim*commRank+2*k+1])
           targetBoundingBoxes[2*dim*commRank+2*k+1] = val1;
      }      
     }// for c

    /************************************************************************** 
    * Step 2: Broadcast the target bounding boxes so that each                *
    *         rank knows the bounding boxes for all ranks                     *
    **************************************************************************/

    for (size_t i=0; i<commSize; i++)
      MPI_Bcast(&(targetBoundingBoxes[0])+2*dim*i, 2*dim, MPI_DOUBLE, i, MPI_COMM_WORLD);

#ifdef DEBUG_MPI
    std::cout << "Target boxes: ";
    for (size_t i=0; i<2*dim*commSize; i++) std::cout << targetBoundingBoxes[i] << " ";
    std::cout << std::endl;
#endif
 
    /************************************************************************** 
    * Step 3: Collect the source particles on the current rank that           *
    *         lie in the target bounding box for each rank in the             *
    *         communicator                                                    *
    **************************************************************************/

    std::vector<bool> sendFlags(commSize, false);
    std::vector<std::vector<int>> sourcePtsToSend(commSize); 
    std::vector<int> sourcePtsToSendSize(commSize); 


    size_t sourceNumPts = source_swarm.num_particles(PARALLEL_OWNED);

    for (size_t i=0; i<commSize; ++i)
    {
      if (i == commRank)
        continue;
      else
      {
        double min[dim], max[dim];

        for (size_t k=0; k<dim; ++k)
        {
          min[k] = targetBoundingBoxes[2*dim*i+2*k];
          max[k] = targetBoundingBoxes[2*dim*i+2*k+1];
        }

        for (size_t c = 0; c < sourceNumPts; ++c)
        {
          Point<dim> coord = source_swarm.get_particle_coordinates(c);

          bool thisPt = true; 
          for (size_t k=0; k<dim; ++k)
          {
            thisPt = thisPt && (coord[k] <= max[k] && coord[k] >= min[k]);
          } 

          if (thisPt)
            sourcePtsToSend[i].emplace_back(c);
        }
      }

      sourcePtsToSendSize[i] = sourcePtsToSend[i].size();
      sendFlags[i] = (sourcePtsToSend[i].size() > 0);
    }

    /************************************************************************** 
    * Step 4: Set up communication info and collect coordinate data for       *
    *         source particles that need to be sent to other ranks            *
    **************************************************************************/
    
    comm_info_t src_info;
    setInfo(&src_info, commSize, sendFlags, sourcePtsToSendSize);
 
    std::vector<std::vector<double>> sourceSendCoords(commSize);
    for (size_t i = 0; i < commSize; ++i)
    {
      if ((!sendFlags[i]) || (i==commRank))
        continue;
      else 
      {
        for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
        {
          Point<dim> coord = source_swarm.get_particle_coordinates(sourcePtsToSend[i][j]);
          for (size_t d = 0 ; d < dim; ++d)
            sourceSendCoords[i].insert(sourceSendCoords[i].end(), coord[d]);
        }
      }
    }

    //move this coordinate data
    std::vector<double> sourceRecvCoords(src_info.newNum*dim);
    moveField<double>(&src_info, commRank, commSize, MPI_DOUBLE, dim, sourceSendCoords, &sourceRecvCoords);

    //DEBUG
    for (size_t i =0; i < commSize; ++i)
    {
      for (size_t j = 0 ; j < sourceSendCoords[i].size(); ++j)
       std::cout<<"Rank-"<<commRank<<"::sourceSendCoords["<<i<<"]["<<j<<"] = "<<sourceSendCoords[i][j]<<std::endl;
      //for (size_t j = 0 ; j < sourceRecvCoords[i].size(); ++j)
      // std::cout<<"Rank-"<<commRank<<"::sourceRecvCoords["<<i<<"]["<<j<<"] = "<<sourceRecvCoords[i][j]<<std::endl;
    }
     for (size_t j = 0 ; j < sourceRecvCoords.size(); ++j)
       std::cout<<"Rank-"<<commRank<<"::sourceRecvCoords["<<j<<"] = "<<sourceRecvCoords[j]<<std::endl;

    //update local source field data with the received data
  /*  int writeOffset = 0;
    for (size_t i = 0; i < commSize; ++i)
    {
      source_swarm.extend_particle_list({sourceRecvCoords[writeOffset], sourceRecvCoords[writeOffset+src_info[i].recvCounts[i]*dim-1]});
      writeOffset += src_info[i].recvCounts[i]*dim; 
    }
*/

      source_swarm.extend_particle_list(sourceRecvCoords);

    /************************************************************************** 
    * Step 5: Collect integer field data from source swarm to be sent         * 
    *         to other ranks                                                  *
    **************************************************************************/
  /*
    std::vector<std::string> int_field_names = source_state->field_names_int();

    for (size_t nvars = 0; nvars < int_field_names.size(); ++nvars)
    {
      // Get field data from source state
      std::shared_ptr<std::vector<int>> srcdata; 
      source_state->get_field(int_field_names[nvars], srcdata);

      // Collect field data for source particles that need to be sent to other ranks
      std::vector<std::vector<int>> sourceSendData(commSize);
      for (size_t i = 0; i < commSize; ++i)
      {
        if ((!sendFlags[i]) || (i==commRank))
          continue;
        else 
        {
          for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
          {
            int sid = sourcePtsToSend[i][j];
            sourceSendData[i].emplace_back((*srcdata)[sid]);
          }
        }
      }

      //move this field data
      std::vector<std::vector<int>> sourceRecvData(commSize);
      moveField<int>(&src_info, commRank, commSize, MPI_INT, 1, sourceSendData, sourceRecvData);

      //update local source field data with the received data

      for (size_t i = 0; i < commSize; ++i)
      {
        source_state->extend_field(int_field_names[nvars], sourceRecvData[i]);
      }
    }

    */
    /************************************************************************** 
    * Step 6: Collect double field data from source swarm to be sent          * 
    *         to other ranks                                                  *
    **************************************************************************/
/*
    std::vector<std::string> dbl_field_names = source_state->field_names_double();

    for (size_t nvars = 0; nvars < dbl_field_names.size(); ++nvars)
    {
      // Get field data from source state
      std::shared_ptr<std::vector<double>> srcdata; 
      source_state->get_field(dbl_field_names[nvars], srcdata);

      // Collect field data for source particles that need to be sent to other ranks
      std::vector<std::vector<double>> sourceSendData(commSize);
      for (size_t i = 0; i < commSize; ++i)
      {
        if ((!sendFlags[i]) || (i==commRank))
          continue;
        else 
        {
          for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
          {
            int sid = sourcePtsToSend[i][j];
            sourceSendData[i].emplace_back((*srcdata)[sid]);
          }
        }
      }

      //move this field data
      std::vector<std::vector<double>> sourceRecvData(commSize);
      moveField<double>(&src_info, commRank, commSize, MPI_DOUBLE, 1, sourceSendData, sourceRecvData);

      //update local source field data with the received data

      for (size_t i = 0; i < commSize; ++i)
      {
        source_state->extend_field(dbl_field_names[nvars], sourceRecvData[i]);
      }
    }
   */

  } // distribute

 private:

  /*!
    @brief Compute fields needed to do comms 
    @param[in] info              Info data structure to be filled
    @param[in] commSize          Total number of MPI ranks
    @param[in] sendFlags         Array of flags:  do I send to PE n?
    @param[in] sourceNum         List of number of particles to be sent to each rank
   */
  void setInfo(comm_info_t* info,
               const int commSize,
               const std::vector<bool>& sendFlags,
               const std::vector<int>& sourceNum)
  {
    // Set my counts of all entities and owned entities
    info->sourceNum.resize(commSize);

    // Each rank will tell each other rank how many indexes it is going to send it
    info->sendCounts.resize(commSize);
    info->recvCounts.resize(commSize);
    
    for (unsigned int i=0; i<commSize; i++)
    {
      info->sourceNum[i] = sourceNum[i];
      info->sendCounts[i] = sendFlags[i] ? info->sourceNum[i] : 0;
    }

    MPI_Alltoall(&(info->sendCounts[0]), 1, MPI_INT,
                 &(info->recvCounts[0]), 1, MPI_INT, MPI_COMM_WORLD);
    
    // Compute the total number of indexes this rank will receive from all ranks
    for (unsigned int i=0; i<commSize; i++)
      info->newNum += info->recvCounts[i];
  } // setInfo

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
  void moveField(comm_info_t* info,
                const int commRank, const int commSize,
                const MPI_Datatype mpiType, const size_t nvals,
                const std::vector<std::vector<T>>& sourceData,
                std::vector<T>* newData)
  {
    // Each rank will do a non-blocking receive from each rank from
    // which it will receive data values
    int writeOffset = 0;
    std::vector<MPI_Request> requests;
   // std::vector<T> buffer(info->newNum*nvals);
    size_t rcount = 0; 
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (info->recvCounts[i] > 0))
      {
        MPI_Request request;
        //std::vector<T> buffer(info->recvCounts[i]*nvals);
        MPI_Irecv((void *)&((*newData)[writeOffset]),
                  info->recvCounts[i]*nvals, mpiType, i,
                  MPI_ANY_TAG, MPI_COMM_WORLD, &request);
        requests.push_back(request);

        //DEBUG
        //for (size_t j = 0; j < buffer.size(); ++j)
        // std::cout<<"buffer["<<j<<"] = "<<buffer[j]<<std::endl;

        //newData[i].insert(newData[i].end(), buffer.begin(), buffer.end());
       // rcount += buffer.size();
      }
     writeOffset += info->recvCounts[i]*nvals;
    }
    //assert(rcount == info->newNum*nvals);
    assert(writeOffset == info->newNum*nvals);

   // for (size_t i = 0; i < commSize; ++i)
   //{
   //     newData[i].insert(newData[i].end(), buffer.begin(), buffer.end());
  // }


    // Each rank will send its data values to appropriate ranks
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (info->sendCounts[i] > 0))
      {
        MPI_Send((void *)&(sourceData[i][0]),
                 info->sendCounts[i]*nvals, mpiType, i, 0, MPI_COMM_WORLD);
      }
    }

    // Wait for all receives to complete
    if (requests.size() > 0)
    {
      std::vector<MPI_Status> statuses(requests.size());
      MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
    }

    std::cout<<"placeholder"<<std::endl;
  } // moveField


}; // MPI_Particle_Distribute
} // namespace Portage

#endif // MPI_Particle_Distribute
