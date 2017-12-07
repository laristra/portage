/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

#include "portage/support/operator.h"
#include "portage/support/Point.h"

using Portage::Point;

using Portage::Meshfree::Basis::Unitary;

using Portage::Meshfree::Operator::VolumeIntegral;
using Portage::Meshfree::Operator::Interval;

using std::vector;

class OperatorTest : public ::testing::Test 
{
public:
  vector<Point<1>> interval_points_;

  virtual void SetUp(){
    interval_points_.resize(2);
    interval_points_[0][0] = 0.0;
    interval_points_[1][0] = 1.0;
  }

  virtual void TearDown(){
  }

  void check_operator()
  {
    vector<vector<double>> result;
    Portage::Meshfree::Operator::apply<1>(VolumeIntegral, Unitary, Interval, interval_points_, result);
    ASSERT_EQ(result[0][0], 1.0);
  }
};

TEST_F(OperatorTest, VolumeIntegralUnitaryInterval) {
  check_operator();
}
