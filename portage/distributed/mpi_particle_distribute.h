/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef MPI_PARTICLE_DISTRIBUTE_H_
#define MPI_PARTICLE_DISTRIBUTE_H_

#ifdef PORTAGE_ENABLE_MPI

#include <cassert>
#include <algorithm>
#include <numeric>
#include <memory>
#include <vector>

#include "portage/accumulate/accumulate.h"
#include "portage/support/portage.h"
#include "portage/support/weight.h"
#include "wonton/support/Point.h"

#include "mpi.h"

/*!
  @file mpi_particle_distribute.h
  @brief Distributes source particles and data using MPI.
 */
using Wonton::Point;

namespace Portage {

/*!
  @class MPI_Particle_Distribute
  @brief Distributes source particles and their data using MPI.
*/
template<size_t dim>
class MPI_Particle_Distribute {
 public:

  /*!
    @brief Constructor of MPI_Particle_Distribute
   */
  MPI_Particle_Distribute(Wonton::MPIExecutor_type const *mpiexecutor) {
    assert(mpiexecutor);
    comm_ = mpiexecutor->mpicomm;
  }


  /*!
    @brief Destructor of MPI_Particle_Distribute
   */
  ~MPI_Particle_Distribute() {};


  /*!
    @brief Helper structure containg comms info
   */
  struct comm_info_t {
    //!< Number of particles in source swarm to be sent to each rank
    std::vector<int> sourceNum;
    //!< Number of total particles in new field
    int newNum = 0;
    //! Array of send sizes from me to all PEs
    std::vector<int> sendCounts;
    //! Array of total recv sizes to me from all PEs
    std::vector<int> recvCounts;
  };


  /*!
    @brief Compute bounding boxes for target swarm on all partitions, and
           send source particles and their data to all target partitions
           within a bounding box using MP. Currently, only Gather scheme
           is supported.
    @param[in] source_swarm       Input swarm
    @param[in] source_state       Input swarm state
    @param[in] target_swarm       Target swarm
    @param[in] target_state       Target swarm state
    @param[in] smoothing_lengths  Extents on source(Scatter-form). Not used for Gather-form.
    @param[in] center             Weight center
   */
  template <class SourceSwarm, class SourceState, class TargetSwarm, class TargetState>
  void distribute(SourceSwarm &source_swarm, SourceState &source_state,
                  TargetSwarm &target_swarm, TargetState &target_state,
                  vector<std::vector<std::vector<double>>> &smoothing_lengths,
                  vector<Point<dim>> &extents,
                  vector<Meshfree::Weight::Kernel>& kernel_types,
                  vector<Meshfree::Weight::Geometry>& geom_types,
                  Meshfree::WeightCenter center=Meshfree::WeightCenter::Gather)
  {
    // Get the MPI communicator size and rank information
    int commSize, commRank;
    MPI_Comm_size(comm_, &commSize);
    MPI_Comm_rank(comm_, &commRank);

    assert(dim == source_swarm.space_dimension());
    assert(dim == target_swarm.space_dimension());

    if (center == Meshfree::WeightCenter::Gather) {
      assert(smoothing_lengths.size() == target_swarm.num_particles(Entity_type::PARALLEL_OWNED));
      assert(kernel_types.size() == target_swarm.num_particles(Entity_type::PARALLEL_OWNED));
      assert(geom_types.size() == target_swarm.num_particles(Entity_type::PARALLEL_OWNED));
    } else if (center == Meshfree::WeightCenter::Scatter) {
      assert(smoothing_lengths.size() == source_swarm.num_particles(Entity_type::PARALLEL_OWNED));
      assert(extents.size() == source_swarm.num_particles(Entity_type::PARALLEL_OWNED));
      assert(kernel_types.size() == source_swarm.num_particles(Entity_type::PARALLEL_OWNED));
      assert(geom_types.size() == source_swarm.num_particles(Entity_type::PARALLEL_OWNED));
    }

    /**************************************************************************
    * Step 1: Compute bounding box for target swarm based on weight center    *
    *         for the current rank, and put it in a vector that will be       *
    *         used later to store the target bounding box for each rank       *
    **************************************************************************/
    size_t targetNumOwnedPts = target_swarm.num_particles(Entity_type::PARALLEL_OWNED);

    std::vector<double> targetBoundingBoxes(2*dim*commSize);

    for (size_t i=0; i<2*dim; i+=2)
    {
      targetBoundingBoxes[2*dim*commRank+i+0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2*dim*commRank+i+1] = -std::numeric_limits<double>::max();
    }

    for (size_t c=0; c<targetNumOwnedPts; ++c)
      {
        Point<dim> coord = target_swarm.get_particle_coordinates(c);
        Point<dim> ext;
        if (center == Meshfree::WeightCenter::Scatter)
          {
            ext = extents[c];
          }
        for (size_t k=0; k < dim; ++k)
          {
            double val0, val1;
            if (center == Meshfree::WeightCenter::Gather)
              {
                val0 = coord[k];
                val1 = coord[k];
              }
            else if (center == Meshfree::WeightCenter::Scatter)
              {
                val0 = coord[k]-ext[k];
                val1 = coord[k]+ext[k];
              }
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
      MPI_Bcast(&(targetBoundingBoxes[0])+2*dim*i, 2*dim, MPI_DOUBLE, i, comm_);

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

    size_t sourceNumPts = source_swarm.num_particles(Entity_type::PARALLEL_OWNED);

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
          Point<dim> ext;
          if (center == Meshfree::WeightCenter::Scatter)
          {
            ext = extents[c];
          }
          bool thisPt = true;
          for (size_t k=0; k<dim; ++k)
          {
            double val0, val1;

            if (center == Meshfree::WeightCenter::Gather)
            {
              val0 = coord[k];
              val1 = coord[k];
            }
            else if (center == Meshfree::WeightCenter::Scatter)
            {
              val0 = coord[k]-ext[k];
              val1 = coord[k]+ext[k];
            }

            //check if the coordinates or the bnds of the current
            //source point either inside or intersecting with the
            //bounding box of the target swarm.
            thisPt = thisPt && (val0 <= max[k] && val1 >= min[k]);
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
            sourceSendCoords[i].push_back(coord[d]);
        }
      }
    }

    //move this coordinate data
    std::vector<double> sourceRecvCoords(src_info.newNum*dim);
    moveField<double>(&src_info, commRank, commSize, MPI_DOUBLE, dim, sourceSendCoords, &sourceRecvCoords);

    // update local source particle list with received new particles
    std::vector<Point<dim>> RecvCoords;
    for (size_t i = 0; i < src_info.newNum; ++i)
    {
      Point<dim> coord;
      for (size_t d = 0 ; d < dim; ++d)
        coord[d] = sourceRecvCoords[dim*i+d];
      RecvCoords.emplace_back(coord);
    }
    source_swarm.extend_particle_list(RecvCoords);

    /**************************************************************************
    * Step 5: Set up communication info and collect smoothing length data     *
    *         for source particles that need to be sent to other ranks for    *
    *         the Scatter scheme                                              *
    ***************************************************************************/
    if (center == Meshfree::WeightCenter::Scatter)
    {
      //collect smoothing length sizes and get max
      size_t smlen_dim = smoothing_lengths[0][0].size();
      size_t nsources = source_swarm.num_particles(Entity_type::PARALLEL_OWNED);
      std::vector<int> smoothing_length_sizes(nsources);
      size_t max_slsize=0;
      for (size_t i=0; i<nsources; i++) {
	smoothing_length_sizes[i] = smoothing_lengths[i].size();
	max_slsize = std::max<size_t>(smoothing_length_sizes[i], max_slsize);
	for (size_t j=0; j<smoothing_lengths[i].size(); j++)
	  assert(smoothing_lengths[i][j].size() == smlen_dim);	  
      }

      //-----------------------------------------------
      //communicate smoothing length sizes
      std::vector<std::vector<int>> sourceSendSmoothSizes(commSize);
      for (size_t i = 0; i < commSize; ++i)
      {
        if ((!sendFlags[i]) || (i==commRank))
          continue;
        else
        {
          for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
          {
            int smlensize = smoothing_length_sizes[sourcePtsToSend[i][j]];
            sourceSendSmoothSizes[i].push_back(smlensize);
          }
        }
      }
      //move this smoothing length size data
      std::vector<int> sourceRecvSmoothSizes(src_info.newNum);
      moveField<int>(&src_info, commRank, commSize, MPI_INT, 1,
                     sourceSendSmoothSizes,&sourceRecvSmoothSizes);

      // update local source particle list with received new particles
      for (size_t i = 0; i < src_info.newNum; ++i)
      {
	int smlensize = sourceRecvSmoothSizes[i];
        smoothing_length_sizes.push_back(smlensize);
      }

      // create cumulative received sizes
      std::vector<int> recvSLSizeCum(src_info.newNum, 0);
      for (size_t i = 1; i < src_info.newNum; ++i) recvSLSizeCum[i] += recvSLSizeCum[i-1];

      //-----------------------------------------------
      //communicate smoothing_lengths
      //this is inefficient if doing facets and there is large variation in number of facets 
      std::vector<std::vector<double>> sourceSendSmoothLengths(commSize);
      for (size_t i = 0; i < commSize; ++i)
      {
        if ((!sendFlags[i]) || (i==commRank))
          continue;
        else
        {
          for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
          {
            std::vector<std::vector<double>> smlen(max_slsize,std::vector<double>(smlen_dim,0.));
	    std::copy(smoothing_lengths[sourcePtsToSend[i][j]].begin(),
		      smoothing_lengths[sourcePtsToSend[i][j]].end(),
		      smlen.begin());
	    for (size_t k = 0 ; k < max_slsize; ++k)
              for (size_t d = 0; d < smlen_dim; d++) 
		sourceSendSmoothLengths[i].push_back(smlen[k][d]);
          }
        }
      }
      //move this smoothing length data
      size_t sl_unit_size = max_slsize*smlen_dim;
      std::vector<double> sourceRecvSmoothLengths(src_info.newNum*max_slsize*smlen_dim);
      moveField<double>(&src_info, commRank, commSize, MPI_DOUBLE, max_slsize*smlen_dim,
                        sourceSendSmoothLengths,&sourceRecvSmoothLengths);

      // update local source particle list with received new particles
      for (size_t i = 0; i < src_info.newNum; ++i)
      {
	std::vector<std::vector<double>> smlen
          (smoothing_length_sizes[nsources+i],std::vector<double>(smlen_dim));
        for (size_t k=0; k < smoothing_length_sizes[nsources+i]; k++) 
          for (size_t d = 0 ; d < smlen_dim; ++d)
            smlen[k][d] = sourceRecvSmoothLengths[recvSLSizeCum[i]+k*smlen_dim+d];

        smoothing_lengths.push_back(smlen);
      }

      //-----------------------------------------------
      //communicate extents
      if (center == Meshfree::WeightCenter::Scatter)
        {
          std::vector<std::vector<double>> sourceSendExtents(commSize);
          for (size_t i = 0; i < commSize; ++i)
            {
              if ((!sendFlags[i]) || (i==commRank))
                continue;
              else
                {
                  for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
                    {
                      Point<dim> ext = extents[sourcePtsToSend[i][j]];
                      for (size_t d = 0 ; d < dim; ++d)
                        sourceSendExtents[i].push_back(ext[d]);
                    }
                }
            }
          //move this extents data
          std::vector<double> sourceRecvExtents(src_info.newNum*dim);
          moveField<double>(&src_info, commRank, commSize, MPI_DOUBLE, dim,
                            sourceSendExtents,&sourceRecvExtents);

          // update local source particle list with received new particles
          for (size_t i = 0; i < src_info.newNum; ++i)
            {
              Point<dim> ext;
              for (size_t d = 0 ; d < dim; ++d)
                ext[d] = sourceRecvExtents[dim*i+d];

              extents.push_back(ext);
            }
        }

      //-----------------------------------------------
      //communicate kernel types
      std::vector<std::vector<int>> sourceSendKernels(commSize);
      for (size_t i = 0; i < commSize; ++i)
      {
        if ((!sendFlags[i]) || (i==commRank))
          continue;
        else
        {
          for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
          {
              sourceSendKernels[i].push_back(kernel_types[sourcePtsToSend[i][j]]);
          }
        }
      }
      //move this kernel type data
      std::vector<int> sourceRecvKernels(src_info.newNum);
      moveField<int>(&src_info, commRank, commSize, MPI_INT, 1,
                     sourceSendKernels,&sourceRecvKernels);

      // update local source swarm with received kernel types
      for (size_t i = 0; i < src_info.newNum; ++i)
      {
        kernel_types.push_back(static_cast<Meshfree::Weight::Kernel>(sourceRecvKernels[i]));
      }

      //-----------------------------------------------
      //communicate geom_types
      std::vector<std::vector<int>> sourceSendGeoms(commSize);
      for (size_t i = 0; i < commSize; ++i)
      {
        if ((!sendFlags[i]) || (i==commRank))
          continue;
        else
        {
          for (size_t j = 0; j < sourcePtsToSendSize[i]; ++j)
          {
              sourceSendGeoms[i].push_back(geom_types[sourcePtsToSend[i][j]]);
          }
        }
      }
      //move this geom type data
      std::vector<int> sourceRecvGeoms(src_info.newNum);
      moveField<int>(&src_info, commRank, commSize, MPI_INT, 1,
                        sourceSendGeoms,&sourceRecvGeoms);

      // update local source swarm with received geom types
      for (size_t i = 0; i < src_info.newNum; ++i)
      {
        geom_types.push_back(static_cast<Meshfree::Weight::Geometry>(sourceRecvGeoms[i]));
      }
    }
    /**************************************************************************
    * Step 6: Collect integer field data from source swarm to be sent         *
    *         to other ranks                                                  *
    **************************************************************************/
    std::vector<std::string> int_field_names = source_state.field_names_int();

    for (size_t nvars = 0; nvars < int_field_names.size(); ++nvars)
    {
      // Get field data from source state
      std::shared_ptr<vector<int>> srcdata;
      source_state.get_field(int_field_names[nvars], srcdata);

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
      std::vector<int> sourceRecvData(src_info.newNum);
      moveField<int>(&src_info, commRank, commSize, MPI_INT, 1, sourceSendData, &sourceRecvData);

      //update local source field data with the received data
      {vector<int> recvtmp(sourceRecvData);
       source_state.extend_field(int_field_names[nvars], recvtmp);}
    }

    /**************************************************************************
    * Step 7: Collect double field data from source swarm to be sent          *
    *         to other ranks                                                  *
    **************************************************************************/
    std::vector<std::string> dbl_field_names = source_state.field_names_double();

    for (size_t nvars = 0; nvars < dbl_field_names.size(); ++nvars)
    {
      // Get field data from source state
      std::shared_ptr<vector<double>> srcdata;
      source_state.get_field(dbl_field_names[nvars], srcdata);

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
      std::vector<double> sourceRecvData(src_info.newNum);
      moveField<double>(&src_info, commRank, commSize, MPI_DOUBLE, 1, sourceSendData, &sourceRecvData);

      //update local source field data with the received data
      {vector<double> recvtmp(sourceRecvData);
       source_state.extend_field(dbl_field_names[nvars], recvtmp);}
    }


  } // distribute

 private:

  MPI_Comm comm_;
  
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
                 &(info->recvCounts[0]), 1, MPI_INT, comm_);

    // Compute the total number of indexes this rank will receive from all ranks
    for (unsigned int i=0; i<commSize; i++)
      info->newNum += info->recvCounts[i];
  } // setInfo

  /*!
    @brief Move values for a single range of data to all ranks as needed
    @tparam[in] T                C++ type of data to be moved
    @param[in] info              Structure with send/recv counts
    @param[in] commRank          MPI rank of this PE
    @param[in] commSize          Total number of MPI ranks
    @param[in] mpiType           MPI type of data (MPI_???) to be moved
    @param[in] nvals             Number of values per particle
    @param[in] sourceData        Source data
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
    size_t rcount = 0;
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (info->recvCounts[i] > 0))
      {
        MPI_Request request;
        MPI_Irecv((void *)&((*newData)[writeOffset]),
                  info->recvCounts[i]*nvals, mpiType, i,
                  MPI_ANY_TAG, comm_, &request);
        requests.push_back(request);
      }
     writeOffset += info->recvCounts[i]*nvals;
    }
    assert(writeOffset == info->newNum*nvals);

    // Each rank will send its data values to appropriate ranks
    for (unsigned int i=0; i<commSize; i++)
    {
      if ((i != commRank) && (info->sendCounts[i] > 0))
      {
        MPI_Send((void *)&(sourceData[i][0]),
                 info->sendCounts[i]*nvals, mpiType, i, 0, comm_);
      }
    }

    // Wait for all receives to complete
    if (requests.size() > 0)
    {
      std::vector<MPI_Status> statuses(requests.size());
      MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
    }

  } // moveField


}; // MPI_Particle_Distribute
} // namespace Portage

#endif  // PORTAGE_ENABLE_MPI

#endif // MPI_Particle_Distribute
