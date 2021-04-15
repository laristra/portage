/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#pragma once

#include "portage/support/timer.h"
#include "params.h"

/**
 * @class Time manager for each step.
 *
 */
class Step {

public:
  /**
   * @brief Create an instance.
   *
   * @param params: input params.
   */
  explicit Step(Params const& params)
    : rank(params.mpi.rank),
      comm(params.mpi.comm)
  {
    tic = timer::now();
  }

  /**
   * @brief End step.
   *
   */
  void stop() {
    MPI_Barrier(comm);
    if (rank == 0) {
      std::printf(" done. \e[32m(%.3f s)\e[0m\n", timer::elapsed(tic));
    }
  }

  /**
   * @brief Start step.
   *
   * @param description: its description.
   * @param first: whether it is the initial step or not.
   */
  void start(std::string const& description) {
    // show elapsed time of previous step
    if (count > 0) { stop(); }

    if (rank == 0) {
      tic = timer::now();
      std::printf("%s ... ", description.data());
      std::fflush(stdout);
    }
    count++;
  }

  /**
   * @brief Reinitialize
   *
   */
  void reset() { count = 0; }

private:
  /** current time point */
  std::chrono::high_resolution_clock::time_point tic;
  /** current step count */
  int count = 0;
  /** current MPI rank */
  int rank = 0;
  /** MPI communicator */
  MPI_Comm comm = MPI_COMM_WORLD;
};