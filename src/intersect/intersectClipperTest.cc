#include "intersectClipper.h"
#include "gtest/gtest.h"
#include "Mesh.hh"
#include "../driver/driver.h"

TEST(intersectClipper, simple){

  std::vector<JaliGeometry::Point> cellA;
  cellA.emplace_back(.6,4);
  cellA.emplace_back(3,4);
  cellA.emplace_back(3,2);

  std::vector<JaliGeometry::Point> cellB;
  cellB.emplace_back(2,5);
  cellB.emplace_back(4.4,3);
  cellB.emplace_back(2,3);

  IntersectClipper<std::vector<JaliGeometry::Point> > isect = IntersectClipper<std::vector<JaliGeometry::Point> >(Portage::pointToXY());
  std::vector<std::vector<double> > moments = isect(cellA, cellB); 
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
  std::vector<JaliGeometry::Point> cellA, cellB;
  cellA.emplace_back(2,5);
  cellA.emplace_back(10.5,5);
  cellA.emplace_back(10.5, 3.5);
  cellA.emplace_back(2,3.5);

  cellB.emplace_back(2.5,0);
  cellB.emplace_back(8.5,0);
  cellB.emplace_back(8.5,5);
  cellB.emplace_back(5.5,2.5);
  cellB.emplace_back(2.5,4);
  
  IntersectClipper<std::vector<JaliGeometry::Point> > isect = IntersectClipper<std::vector<JaliGeometry::Point> >(Portage::pointToXY());
  std::vector<std::vector<double> > moments = isect(cellA, cellB); 
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

