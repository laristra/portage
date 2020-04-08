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
template<int dim>
class MPI_Particle_Distribute {
public:

  /*!
    @brief Constructor of MPI_Particle_Distribute
   */
  explicit MPI_Particle_Distribute(Wonton::MPIExecutor_type const *mpiexecutor) {
    assert(mpiexecutor);
    comm_ = mpiexecutor->mpicomm;
  }


  /*!
    @brief Destructor of MPI_Particle_Distribute
   */
  ~MPI_Particle_Distribute() = default;;


  /*!
    @brief Helper structure containg comms info
   */
  struct comm_info_t {
    //!< Number of particles in source swarm to be sent to each rank
    std::vector<int> source_num;
    //!< Number of total particles in new field
    int new_num = 0;
    //! Array of send sizes from me to all PEs
    std::vector<int> send_counts;
    //! Array of total recv sizes to me from all PEs
    std::vector<int> recv_counts;
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
    @param[in] smoothing_lengths  Smoothing lengths, sized by target for Gather, source for Scatter
    @param[in] source_extents     Extents for source partices in Scatter mode
    @param[in] smoothing_lengths  Extents for target particles in Gather mode
    @param[in] kernel_types       Kernel types
    @param[in] geom_types         Geometry types
    @param[in] center             Weight center
   */
  template<class SourceSwarm, class SourceState, class TargetSwarm, class TargetState>
  void distribute(SourceSwarm& source_swarm, SourceState& source_state,
                  TargetSwarm& target_swarm, TargetState& target_state,
                  Portage::vector<std::vector<std::vector<double>>>& smoothing_lengths,
                  Portage::vector<Point<dim>>& source_extents,
                  Portage::vector<Point<dim>>& target_extents,
                  Portage::vector<Meshfree::Weight::Kernel>& kernel_types,
                  Portage::vector<Meshfree::Weight::Geometry>& geom_types,
                  Meshfree::WeightCenter center = Meshfree::WeightCenter::Gather) {
    // Get the MPI communicator size and rank information
    int nb_ranks, rank;
    MPI_Comm_size(comm_, &nb_ranks);
    MPI_Comm_rank(comm_, &rank);

    assert(dim == source_swarm.space_dimension());
    assert(dim == target_swarm.space_dimension());

    int nb_target_points = target_swarm.num_particles(Wonton::PARALLEL_OWNED);
    int nb_source_points = source_swarm.num_particles(Wonton::PARALLEL_OWNED);

#ifdef DEBUG
    if (center == Meshfree::WeightCenter::Gather) {
      assert(smoothing_lengths.size() == unsigned(nb_target_points));
      assert(target_extents.size()    == unsigned(nb_target_points));
      assert(kernel_types.size()      == unsigned(nb_target_points));
      assert(geom_types.size()        == unsigned(nb_target_points));
    } else if (center == Meshfree::WeightCenter::Scatter) {
      assert(smoothing_lengths.size() == unsigned(nb_source_points));
      assert(source_extents.size()    == unsigned(nb_source_points));
      assert(kernel_types.size()      == unsigned(nb_source_points));
      assert(geom_types.size()        == unsigned(nb_source_points));
    }
#endif

    /**************************************************************************
    * Step 1: Compute bounding box for target swarm based on weight center    *
    *         for the current rank, and put it in a vector that will be       *
    *         used later to store the target bounding box for each rank       *
    **************************************************************************/

    std::vector<double> targetBoundingBoxes(2 * dim * nb_ranks);

    for (int i = 0; i < 2 * dim; i += 2) {
      targetBoundingBoxes[2 * dim * rank + i + 0] = std::numeric_limits<double>::max();
      targetBoundingBoxes[2 * dim * rank + i + 1] = -std::numeric_limits<double>::max();
    }

    for (int c = 0; c < nb_target_points; ++c) {
      Point<dim> coord = target_swarm.get_particle_coordinates(c);
      Point<dim> ext;
      if (center == Meshfree::WeightCenter::Gather) {
        ext = target_extents[c];
      }
      for (int k = 0; k < dim; ++k) {
        double val0 = 0., val1 = 0.;
        if (center == Meshfree::WeightCenter::Gather) {
          val0 = coord[k] - ext[k];
          val1 = coord[k] + ext[k];
        } else if (center == Meshfree::WeightCenter::Scatter) {
          val0 = coord[k];
          val1 = coord[k];
        }
        if (val0 < targetBoundingBoxes[2 * dim * rank + 2 * k])
          targetBoundingBoxes[2 * dim * rank + 2 * k] = val0;

        if (val1 > targetBoundingBoxes[2 * dim * rank + 2 * k + 1])
          targetBoundingBoxes[2 * dim * rank + 2 * k + 1] = val1;
      }
    }// for c

    /**************************************************************************
    * Step 2: Broadcast the target bounding boxes so that each                *
    *         rank knows the bounding boxes for all ranks                     *
    **************************************************************************/
    for (int i = 0; i < nb_ranks; i++)
      MPI_Bcast(&(targetBoundingBoxes[0]) + 2 * dim * i, 2 * dim, MPI_DOUBLE, i, comm_);

#ifdef DEBUG_MPI
    std::cout << "Target boxes: ";
    for (int i=0; i<2*dim*nb_ranks; i++) std::cout << targetBoundingBoxes[i] << " ";
    std::cout << std::endl;
#endif

    /**************************************************************************
    * Step 3: Collect the source particles on the current rank that           *
    *         lie in the target bounding box for each rank in the             *
    *         communicator                                                    *
    **************************************************************************/
    std::vector<bool> sendFlags(nb_ranks, false);
    std::vector<std::vector<int>> sourcePtsToSend(nb_ranks);
    std::vector<int> sourcePtsToSendSize(nb_ranks);

    for (int i = 0; i < nb_ranks; ++i) {
      if (i == rank)
        continue;
      else {
        double min[dim], max[dim];

        for (int k = 0; k < dim; ++k) {
          min[k] = targetBoundingBoxes[2 * dim * i + 2 * k];
          max[k] = targetBoundingBoxes[2 * dim * i + 2 * k + 1];
        }

        for (int c = 0; c < nb_source_points; ++c) {
          Point<dim> coord = source_swarm.get_particle_coordinates(c);
          Point<dim> ext;
          if (center == Meshfree::WeightCenter::Scatter) {
            ext = source_extents[c];
          }
          bool thisPt = true;
          for (int k = 0; k < dim; ++k) {
            double val0 = 0., val1 = 0.;

            if (center == Meshfree::WeightCenter::Gather) {
              val0 = coord[k];
              val1 = coord[k];
            } else if (center == Meshfree::WeightCenter::Scatter) {
              val0 = coord[k] - ext[k];
              val1 = coord[k] + ext[k];
            }

            //check if the coordinates or the bnds of the current
            //source point either inside or intersecting with the
            //bounding box of the target swarm.
            thisPt = thisPt && (val0 <= max[k] && val1 >= min[k]);
          }

          if (thisPt) sourcePtsToSend[i].push_back(c);
        }
      }

      sourcePtsToSendSize[i] = sourcePtsToSend[i].size();
      sendFlags[i] = not sourcePtsToSend[i].empty();
    }

    /**************************************************************************
    * Step 4: Set up communication info and collect coordinate data for       *
    *         source particles that need to be sent to other ranks            *
    **************************************************************************/
    comm_info_t src_info;
    setInfo(&src_info, nb_ranks, sendFlags, sourcePtsToSendSize);

    std::vector<std::vector<double>> sourceSendCoords(nb_ranks);
    for (int i = 0; i < nb_ranks; ++i) {
      if ((!sendFlags[i]) || (i == rank))
        continue;
      else {
        for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
          Point<dim> coord = source_swarm.get_particle_coordinates(sourcePtsToSend[i][j]);
          for (int d = 0; d < dim; ++d)
            sourceSendCoords[i].push_back(coord[d]);
        }
      }
    }

    //move this coordinate data
    std::vector<double> sourceRecvCoords(src_info.new_num * dim);
    moveField<double>(&src_info, rank, nb_ranks, MPI_DOUBLE, dim, sourceSendCoords, &sourceRecvCoords);

    // update local source particle list with received new particles
    std::vector<Point<dim>> RecvCoords;
    for (int i = 0; i < src_info.new_num; ++i) {
      Point<dim> coord;
      for (int d = 0; d < dim; ++d)
        coord[d] = sourceRecvCoords[dim * i + d];
      RecvCoords.push_back(coord);
    }
    source_swarm.extend_particle_list(RecvCoords);

    /**************************************************************************
    * Step 5: Set up communication info and collect smoothing length data     *
    *         for source particles that need to be sent to other ranks for    *
    *         the Scatter scheme                                              *
    ***************************************************************************/
    if (center == Meshfree::WeightCenter::Scatter) {
      //collect smoothing length sizes and get max
      std::vector<std::vector<double>> h0 = smoothing_lengths[0];
      int smlen_dim = h0[0].size();
      std::vector<int> smoothing_length_sizes(nb_source_points);
      int max_slsize = 0;
      for (int i = 0; i < nb_source_points; i++) {
        std::vector<std::vector<double>> h = smoothing_lengths[i];
        smoothing_length_sizes[i] = h.size();
        max_slsize = std::max<int>(smoothing_length_sizes[i], max_slsize);
        for (int j = 0; j < smoothing_length_sizes[i]; j++)
          assert(h[j].size() == unsigned(smlen_dim));
      }

      //-----------------------------------------------
      //communicate smoothing length sizes
      std::vector<std::vector<int>> sourceSendSmoothSizes(nb_ranks);
      for (int i = 0; i < nb_ranks; ++i) {
        if ((!sendFlags[i]) || (i == rank))
          continue;
        else {
          for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
            int smlensize = smoothing_length_sizes[sourcePtsToSend[i][j]];
            sourceSendSmoothSizes[i].push_back(smlensize);
          }
        }
      }
      //move this smoothing length size data
      std::vector<int> sourceRecvSmoothSizes(src_info.new_num);
      moveField<int>(&src_info, rank, nb_ranks, MPI_INT, 1,
                     sourceSendSmoothSizes, &sourceRecvSmoothSizes);

      // update local source particle list with received new particles
      for (int i = 0; i < src_info.new_num; ++i) {
        int smlensize = sourceRecvSmoothSizes[i];
        smoothing_length_sizes.push_back(smlensize);
      }

      // create cumulative received sizes
      std::vector<int> recvSLSizeCum(src_info.new_num, 0);
      for (int i = 1; i < src_info.new_num; ++i) recvSLSizeCum[i] += recvSLSizeCum[i - 1];

      //-----------------------------------------------
      //communicate smoothing_lengths
      //this is inefficient if doing facets and there is large variation in number of facets
      std::vector<std::vector<double>> sourceSendSmoothLengths(nb_ranks);
      for (int i = 0; i < nb_ranks; ++i) {
        if ((!sendFlags[i]) || (i == rank))
          continue;
        else {
          for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
            std::vector<std::vector<double>> smlen(max_slsize, std::vector<double>(smlen_dim, 0.));
            std::vector<std::vector<double>> h = smoothing_lengths[sourcePtsToSend[i][j]];
            std::copy(h.begin(), h.end(), smlen.begin());
            for (int k = 0; k < max_slsize; ++k)
              for (int d = 0; d < smlen_dim; d++)
                sourceSendSmoothLengths[i].push_back(smlen[k][d]);
          }
        }
      }
      //move this smoothing length data
      std::vector<double> sourceRecvSmoothLengths(src_info.new_num * max_slsize * smlen_dim);
      moveField<double>(&src_info, rank, nb_ranks, MPI_DOUBLE, max_slsize * smlen_dim,
                        sourceSendSmoothLengths, &sourceRecvSmoothLengths);

      // update local source particle list with received new particles
      for (int i = 0; i < src_info.new_num; ++i) {
        std::vector<std::vector<double>> smlen
          (smoothing_length_sizes[nb_source_points + i], std::vector<double>(smlen_dim));
        for (int k = 0; k < smoothing_length_sizes[nb_source_points + i]; k++)
          for (int d = 0; d < smlen_dim; ++d)
            smlen[k][d] = sourceRecvSmoothLengths[recvSLSizeCum[i] + k * smlen_dim + d];

        smoothing_lengths.push_back(smlen);
      }

      //-----------------------------------------------
      //communicate extents
      if (center == Meshfree::WeightCenter::Scatter) {
        std::vector<std::vector<double>> sourceSendExtents(nb_ranks);
        for (int i = 0; i < nb_ranks; ++i) {
          if ((!sendFlags[i]) || (i == rank))
            continue;
          else {
            for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
              Point<dim> ext = source_extents[sourcePtsToSend[i][j]];
              for (int d = 0; d < dim; ++d)
                sourceSendExtents[i].push_back(ext[d]);
            }
          }
        }
        //move this extents data
        std::vector<double> sourceRecvExtents(src_info.new_num * dim);
        moveField<double>(&src_info, rank, nb_ranks, MPI_DOUBLE, dim,
                          sourceSendExtents, &sourceRecvExtents);

        // update local source particle list with received new particles
        for (int i = 0; i < src_info.new_num; ++i) {
          Point<dim> ext;
          for (int d = 0; d < dim; ++d)
            ext[d] = sourceRecvExtents[dim * i + d];

          source_extents.push_back(ext);
        }
      }

      //-----------------------------------------------
      //communicate kernel types
      std::vector<std::vector<int>> sourceSendKernels(nb_ranks);
      for (int i = 0; i < nb_ranks; ++i) {
        if ((!sendFlags[i]) || (i == rank))
          continue;
        else {
          for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
            sourceSendKernels[i].push_back(kernel_types[sourcePtsToSend[i][j]]);
          }
        }
      }
      //move this kernel type data
      std::vector<int> sourceRecvKernels(src_info.new_num);
      moveField<int>(&src_info, rank, nb_ranks, MPI_INT, 1,
                     sourceSendKernels, &sourceRecvKernels);

      // update local source swarm with received kernel types
      for (int i = 0; i < src_info.new_num; ++i) {
        kernel_types.push_back(static_cast<Meshfree::Weight::Kernel>(sourceRecvKernels[i]));
      }

      //-----------------------------------------------
      //communicate geom_types
      std::vector<std::vector<int>> sourceSendGeoms(nb_ranks);
      for (int i = 0; i < nb_ranks; ++i) {
        if ((!sendFlags[i]) || (i == rank))
          continue;
        else {
          for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
            sourceSendGeoms[i].push_back(geom_types[sourcePtsToSend[i][j]]);
          }
        }
      }
      //move this geom type data
      std::vector<int> sourceRecvGeoms(src_info.new_num);
      moveField<int>(&src_info, rank, nb_ranks, MPI_INT, 1,
                     sourceSendGeoms, &sourceRecvGeoms);

      // update local source swarm with received geom types
      for (int i = 0; i < src_info.new_num; ++i) {
        geom_types.push_back(static_cast<Meshfree::Weight::Geometry>(sourceRecvGeoms[i]));
      }
    }
    /**************************************************************************
    * Step 6: Collect integer field data from source swarm to be sent         *
    *         to other ranks                                                  *
    **************************************************************************/
    auto int_field_names = source_state.template get_field_names<int>();
    int const num_int_fields = int_field_names.size();
    for (int nvars = 0; nvars < num_int_fields; ++nvars) {
      // Get field data from source state
      auto& srcdata = source_state.get_field_int(int_field_names[nvars]);

      // Collect field data for source particles that need to be sent to other ranks
      std::vector<std::vector<int>> sourceSendData(nb_ranks);
      for (int i = 0; i < nb_ranks; ++i) {
        if ((!sendFlags[i]) || (i == rank))
          continue;
        else {
          for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
            int sid = sourcePtsToSend[i][j];
            sourceSendData[i].push_back(srcdata[sid]);
          }
        }
      }

      //move this field data
      std::vector<int> sourceRecvData(src_info.new_num);
      moveField<int>(&src_info, rank, nb_ranks, MPI_INT, 1, sourceSendData, &sourceRecvData);

      //update local source field data with the received data
      {
        vector<int> recvtmp(sourceRecvData);
        source_state.extend_field(int_field_names[nvars], recvtmp);
      }
    }

    /**************************************************************************
    * Step 7: Collect double field data from source swarm to be sent          *
    *         to other ranks                                                  *
    **************************************************************************/
    auto dbl_field_names = source_state.template get_field_names<double>();
    int const num_dbl_fields = dbl_field_names.size();
    for (int nvars = 0; nvars < num_dbl_fields; ++nvars) {
      // Get field data from source state
      auto& srcdata = source_state.get_field(dbl_field_names[nvars]);

      // Collect field data for source particles that need to be sent to other ranks
      std::vector<std::vector<double>> sourceSendData(nb_ranks);
      for (int i = 0; i < nb_ranks; ++i) {
        if ((!sendFlags[i]) || (i == rank))
          continue;
        else {
          for (int j = 0; j < sourcePtsToSendSize[i]; ++j) {
            int sid = sourcePtsToSend[i][j];
            sourceSendData[i].push_back(srcdata[sid]);
          }
        }
      }

      //move this field data
      std::vector<double> sourceRecvData(src_info.new_num);
      moveField<double>(&src_info, rank, nb_ranks, MPI_DOUBLE, 1, sourceSendData, &sourceRecvData);

      //update local source field data with the received data
      {
        vector<double> recvtmp(sourceRecvData);
        source_state.extend_field(dbl_field_names[nvars], recvtmp);
      }
    }


  } // distribute

private:

  MPI_Comm comm_ = MPI_COMM_NULL;

  /*!
    @brief Compute fields needed to do comms
    @param[in] info              Info data structure to be filled
    @param[in] nb_ranks          Total number of MPI ranks
    @param[in] sendFlags         Array of flags:  do I send to PE n?
    @param[in] sourceNum         List of number of particles to be sent to each rank
   */
  void setInfo(comm_info_t *info,
               const int nb_ranks,
               const std::vector<bool>& sendFlags,
               const std::vector<int>& sourceNum) {
    // Set my counts of all entities and owned entities
    info->source_num.resize(nb_ranks);

    // Each rank will tell each other rank how many indexes it is going to send it
    info->send_counts.resize(nb_ranks);
    info->recv_counts.resize(nb_ranks);

    for (int i = 0; i < nb_ranks; i++) {
      info->source_num[i] = sourceNum[i];
      info->send_counts[i] = sendFlags[i] ? info->source_num[i] : 0;
    }

    MPI_Alltoall(&(info->send_counts[0]), 1, MPI_INT,
                 &(info->recv_counts[0]), 1, MPI_INT, comm_);

    // Compute the total number of indexes this rank will receive from all ranks
    for (int i = 0; i < nb_ranks; i++)
      info->new_num += info->recv_counts[i];
  } // setInfo

  /*!
    @brief Move values for a single range of data to all ranks as needed
    @tparam[in] T                C++ type of data to be moved
    @param[in] info              Structure with send/recv counts
    @param[in] rank          MPI rank of this PE
    @param[in] nb_ranks          Total number of MPI ranks
    @param[in] datatype           MPI type of data (MPI_???) to be moved
    @param[in] nvals             Number of values per particle
    @param[in] sourceData        Source data
    @param[in] newData           Array of new source data
   */
  template<typename T>
  void moveField(comm_info_t *info,
                 int rank, int nb_ranks,
                 const MPI_Datatype datatype, const int nvals,
                 const std::vector<std::vector<T>>& sourceData,
                 std::vector<T> *newData) {
    // Each rank will do a non-blocking receive from each rank from
    // which it will receive data values
    int writeOffset = 0;
    std::vector<MPI_Request> requests;

    for (int i = 0; i < nb_ranks; i++) {
      if ((i != rank) && (info->recv_counts[i] > 0)) {
        MPI_Request request;
        MPI_Irecv((void *) &((*newData)[writeOffset]),
                  info->recv_counts[i] * nvals, datatype, i,
                  MPI_ANY_TAG, comm_, &request);
        requests.push_back(request);
      }
      writeOffset += info->recv_counts[i] * nvals;
    }
    assert(writeOffset == info->new_num * nvals);

    // Each rank will send its data values to appropriate ranks
    for (int i = 0; i < nb_ranks; i++) {
      if ((i != rank) && (info->send_counts[i] > 0)) {
        MPI_Send((void *) &(sourceData[i][0]),
                 info->send_counts[i] * nvals, datatype, i, 0, comm_);
      }
    }

    // Wait for all receives to complete
    if (not requests.empty()) {
      std::vector<MPI_Status> statuses(requests.size());
      MPI_Waitall(requests.size(), &(requests[0]), &(statuses[0]));
    }

  } // moveField


}; // MPI_Particle_Distribute
} // namespace Portage

#endif  // PORTAGE_ENABLE_MPI

#endif // MPI_Particle_Distribute
