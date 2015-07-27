/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "gtest/gtest.h"
#include "mpi.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc,argv);

  MPI_Init(&argc,&argv);

  int status = RUN_ALL_TESTS();

  MPI_Finalize();

  return status;
}
