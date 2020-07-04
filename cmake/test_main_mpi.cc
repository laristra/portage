/*
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include <mpi.h>

#include <gtest/gtest.h>

//----------------------------------------------------------------------------//
// Main
//----------------------------------------------------------------------------//

int main(int argc, char ** argv) {

  // Initialize GTest
  ::testing::InitGoogleTest(&argc, argv);

  // Initialize MPI
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Make sure listeners are present only on rank 0 - otherwise the
  // output could be quite garbled
  ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
  if (rank != 0)
    delete listeners.Release(listeners.default_result_printer());
  
  int result(0);

  result = RUN_ALL_TESTS();

  // Shutdown MPI
  MPI_Finalize();

  return result;
}

