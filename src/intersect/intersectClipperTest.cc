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
