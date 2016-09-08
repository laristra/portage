#include "intersect_r2d.h"
#include "gtest/gtest.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

/*! 
 * @brief Intersect two cells on two single cell meshes to compute moments.
 * Intersect two cells contained in the mesh: 
 * (0, 0) (2, 0) (2, 2) (0,2) with (1,1) (2,1) (2,2) (1,2)
 * Results should be an area of 1 and a centroid of 1.5, 1.5.
 */
TEST(intersectR2D, simple1) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  //Create mesh from 0, 0 to 2.4, 2 1x1
  std::shared_ptr<Jali::Mesh> sm = mf(0, 0, 2, 2, 1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(1,1,2, 2, 1, 1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);
  Portage::IntersectR2D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  std::vector<std::vector<double> > moments = isect(0, 0); 
  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j] << std::endl;
    }   
  }

  ASSERT_EQ(moments[0][0], 1);
  ASSERT_EQ(moments[0][1], 1.5);
  ASSERT_EQ(moments[0][2], 1.5);
}

