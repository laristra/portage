/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#include "intersect_r3d.h"
#include "gtest/gtest.h"
#include "MeshFactory.hh"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"

TEST(intersectR3D, simple1) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0,0,0, 2,2,2, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(1,1,1, 2,2,2, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
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
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0,0,0, 2,2,2, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(0,0,0, 2,2,2, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
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
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0,0,0, 3,3,3, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(1,1,1, 2,2,2, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
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
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0,0,0, 2,2,2, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(1,1,1, 3,3,3, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
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
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0,0,0, 10,10,10, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(-5,-5,-5, 5,5,5, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
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
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(0,0,0, 10,10,10, 5,5,5);
  std::shared_ptr<Jali::Mesh> tm = mf(0,0,0, 10,10,10, 2,2,2);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
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
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(-2,-2,-2, 0,0,0, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(-1,-1,-1, 0,0,0, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
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
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(-4,-4,-4, 0,0,0, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(-3,-3,-3, 0,0,0, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
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

TEST(intersectR3D, simple9) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> sm = mf(-4,-3,-2, 0,1,2, 1,1,1);
  std::shared_ptr<Jali::Mesh> tm = mf(-3,-2,-1, 0,1,2, 1,1,1);
  const Portage::Jali_Mesh_Wrapper s(*sm);
  const Portage::Jali_Mesh_Wrapper t(*tm);

  const double eps = 1e-12;
  const Portage::IntersectR3D<Portage::Jali_Mesh_Wrapper> isect{s , t};
  const std::vector<std::vector<double> > moments = isect(0, 0);
  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j] << std::endl;
    }
  }

  ASSERT_TRUE(moments.size() == 1);
  ASSERT_TRUE(moments[0].size() == 4);

  ASSERT_TRUE(std::abs(moments[0][0] - 27)   < eps);
  ASSERT_TRUE(std::abs(moments[0][1]/moments[0][0] - (-1.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][2]/moments[0][0] - (-0.5)) < eps);
  ASSERT_TRUE(std::abs(moments[0][3]/moments[0][0] - (+0.5)) < eps);
}
