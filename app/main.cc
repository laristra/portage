#include <cstdio>
#include <cstdlib>
#include <vector>

#include "mpi.h"

#include "driver.h"
#include "state_vector.h"
#include "state.h"

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
  Jali::Mesh* inputMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  // Create a 2d quad output mesh from (0,0) to (1,1) with 2x2 zones
  Jali::Mesh* targetMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  // TODO: populate inputMesh with data using Rao's StateManager
  
  Portage::State inputState(inputMesh);
  std::vector<double> inputData = {0.0,1.0,2.0,3.0};
  inputState.add("celldata",Jali::CELL,&(inputData[0]));

  Portage::State targetState(targetMesh);
  std::vector<double> targetData = {0.0,0.0,0.0,0.0};
  targetState.add("celldata",Jali::CELL,&(targetData[1]));

  Portage::Driver d(*inputMesh, inputState, *targetMesh, targetState);
  d.run();

  // When done, the "celldata" vector on the target mesh cells should
  // be identical to the "celldata" vector on the destination mesh

  Portage::State::const_iterator it = targetState.find("celldata",Jali::CELL);
  if (it == targetState.end()) {
    std::cerr << "Could not find vector with name celldata in targetState" << std::endl;
    std::exit(-1);
  }
  Portage::StateVector outvec = *it;
  
  std::cerr << "celldata vector on target mesh after remapping is:" << std::endl;
  std::cerr << "   " << outvec[0] << ", " << outvec[1] << ", " << outvec[2] << ", " << outvec[3] << std::endl;

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


