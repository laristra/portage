#include <cstdio>
#include <cstdlib>
#include <vector>

#include "mpi.h"

#include "driver.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

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

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  // Create a 2d quad input mesh from (0,0) to (1,1) with 2x2 zones
  Jali::Mesh* inputMesh = mf.create(0.0, 0.0,
				    1.0, 1.0,
				    2, 2);
  // Create a 2d quad output mesh from (0,0) to (1,1) with 2x2 zones
  Jali::Mesh* targetMesh = mf.create(0.0, 0.0,
				     1.0, 1.0,
				     2, 2);

  // TODO: populate inputMesh with data using Rao's StateManager

  Driver d(inputMesh, targetMesh);
  d.run();

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


