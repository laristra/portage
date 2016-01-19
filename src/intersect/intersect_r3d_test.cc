#include "intersect_r3d.h"
#include "gtest/gtest.h"
#include "../driver/driver.h"
#include "MeshFactory.hh"

TEST(intersectR3D, simple1) {
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 1.5) < eps);
}

TEST(intersectR3D, simple2) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(0,0,0, 2,2,2, 1,1,1, NULL, true, true, true, true);
  Jali::Mesh* tm = mf(0,0,0, 2,2,2, 1,1,1, NULL, true, true, true, true);
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 8)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 1) < eps);
}

TEST(intersectR3D, simple3) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(0,0,0, 3,3,3, 1,1,1, NULL, true, true, true, true);
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 1.5) < eps);
}

TEST(intersectR3D, simple4) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(0,0,0, 2,2,2, 1,1,1, NULL, true, true, true, true);
  Jali::Mesh* tm = mf(1,1,1, 3,3,3, 1,1,1, NULL, true, true, true, true);
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 1.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 1.5) < eps);
}

TEST(intersectR3D, simple5) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(0,0,0, 10,10,10, 1,1,1, NULL, true, true, true, true);
  Jali::Mesh* tm = mf(-5,-5,-5, 5,5,5, 1,1,1, NULL, true, true, true, true);
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 125)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 2.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 2.5) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 2.5) < eps);
}

TEST(intersectR3D, simple6) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(0,0,0, 10,10,10, 5,5,5, NULL, true, true, true, true);
  Jali::Mesh* tm = mf(0,0,0, 10,10,10, 2,2,2, NULL, true, true, true, true);
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 8)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 1) < eps);


  moments = isect(1, 0);
  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j] << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 8)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 3) < eps);

  moments = isect(2, 0);
  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j] << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 4)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - 1) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - 4.5) < eps);
}

TEST(intersectR3D, simple7) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(-2,-2,-2, 0,0,0, 1,1,1, NULL, true, true, true, true);
  Jali::Mesh* tm = mf(-1,-1,-1, 0,0,0, 1,1,1, NULL, true, true, true, true);
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 1)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - (-0.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - (-0.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - (-0.5)) < eps);
}

TEST(intersectR3D, simple8) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh* sm = mf(-4,-4,-4, 0,0,0, 1,1,1, NULL, true, true, true, true);
  Jali::Mesh* tm = mf(-3,-3,-3, 0,0,0, 1,1,1, NULL, true, true, true, true);
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

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 27)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - (-1.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - (-1.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - (-1.5)) < eps);
}
