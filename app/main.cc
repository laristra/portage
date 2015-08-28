#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>

#include "mpi.h"

#include "portage/driver/driver.h"
#include "portage/state/state_vector.h"
#include "portage/state/state.h"

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

  // Create a 2d quad input mesh from (0,0) to (1,1) with 3x3 zones
  Jali::Mesh* inputMesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
  // Create a 2d quad output mesh from (0,0) to (1,1) with 4x4 zones
  Jali::Mesh* targetMesh = mf(0.0, 0.0, 1.0, 1.0, 4, 4);

  Portage::State inputState(inputMesh);
  std::vector<double> inputData = {0.0,1.0,2.0,1.0,2.0,3.0,2.0,3.0,4.0};
  Portage::StateVector & cellvecin = 
      inputState.add("celldata",Jali::CELL,&(inputData[0]));

  Portage::State targetState(targetMesh);
  std::vector<double> targetData(16,0.0);
  Portage::StateVector & cellvecout = 
      targetState.add("celldata",Jali::CELL,&(targetData[0]));

  Portage::Driver d(*inputMesh, inputState, *targetMesh, targetState);
  std::vector<std::string> remap_fields;
  remap_fields.push_back("celldata");
  d.set_remap_var_names(remap_fields);
  d.run();

  // When done, the "remapped_data" vector on the target mesh cells should
  // be identical to the "celldata" vector on the destination mesh

  std::cerr << "celldata vector on target mesh after remapping is:" << std::endl;
  std::cerr << cellvecout << std::endl;

  std::printf("finishing portageapp...\n");

  MPI_Finalize();
}


