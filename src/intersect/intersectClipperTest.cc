#include "intersect.h"
#include "intersectClipper.h"
#include "gtest/gtest.h"
TEST(intersectClipper, simple){
  std::vector<std::pair<double, double> > polyA, polyB;
  polyA.emplace_back(.6,4);
  polyA.emplace_back(3,4);
  polyA.emplace_back(3,2);

  polyB.emplace_back(2,5);
  polyB.emplace_back(4.4,3);
  polyB.emplace_back(2, 3);
  
  IntersectClipper isect;
  std::vector<std::vector<double> > moments = isect(polyA, polyB); 
  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j] << std::endl;
    }   
  }
  ASSERT_EQ(moments[0][0], 1);
  ASSERT_EQ(moments[0][1], 2.5);
  ASSERT_EQ(moments[0][2], 3.5);
}

TEST(intersectClipper, convex){
  std::vector<std::pair<double, double> > polyA, polyB;
  polyA.emplace_back(2,5);
  polyA.emplace_back(10.5,5);
  polyA.emplace_back(10.5, 3.5);
  polyA.emplace_back(2,3.5);

  polyB.emplace_back(2.5,0);
  polyB.emplace_back(8.5,0);
  polyB.emplace_back(8.5,5);
  polyB.emplace_back(5.5,2.5);
  polyB.emplace_back(2.5,4);
  
  IntersectClipper isect;
  std::vector<std::vector<double> > moments = isect(polyA, polyB); 
  for(int i=0;i<moments.size();i++){
    for(int j=0;j<moments[i].size();j++){
      std::cout << "i, j, m " << i << ", " << j << ", " << moments[i][j] << std::endl;
    }   
  }
  double eps = 1e-6;
  EXPECT_NEAR(moments[0][0], 1.35, eps);
  EXPECT_NEAR(moments[0][1], 10.665, eps);
  EXPECT_NEAR(moments[0][2], 5.4, eps);
  EXPECT_NEAR(moments[1][0], .25, eps);
  EXPECT_NEAR(moments[1][1], .708333, eps);
  EXPECT_NEAR(moments[1][2], .916667, eps);
}

