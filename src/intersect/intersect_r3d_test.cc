#include "intersectClipper.h" // TODO: remove this include
#include "intersect_r3d.h"
#include "gtest/gtest.h"
#include "Mesh.hh"
#include "../driver/driver.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "MeshFactory.hh"

TEST(intersectR3D, simple) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(0,0,0, 2,2,2, 1,1,1, NULL, true, true, true, true);
  Jali::Mesh* tm = mf(1,1,1, 2,2,2, 1,1,1, NULL, true, true, true, true);
  Portage::Jali_Mesh_Wrapper s(*sm);
  Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  std::vector<std::vector<double> > moments = isect(0, 0);
  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j] << std::endl;
    }
  }

  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      ASSERT_TRUE(std::abs(moments[i][0] - 1/48.) < eps);
    }
  }
}
