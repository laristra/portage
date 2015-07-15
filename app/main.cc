#include <cstdio>
#include <cstdlib>
#include <vector>

#include "mpi.h"

#include "driver.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  if (numpe > 1) {
      std::printf("error - only 1 mpi rank is allowed\n");
      std::exit(1);
  }

  std::printf("starting portageapp...\n");

  Driver d;
  d.run();

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


