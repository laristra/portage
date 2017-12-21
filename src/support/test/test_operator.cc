/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

#include <vector>

#include "portage/support/operator.h"
#include "portage/support/Point.h"

using Portage::Point;

using Portage::Meshfree::Basis::Unitary;
using Portage::Meshfree::Basis::Linear;
using Portage::Meshfree::Basis::Quadratic;

using Portage::Meshfree::Operator::Operator;
using Portage::Meshfree::Operator::VolumeIntegral;
using Portage::Meshfree::Operator::Interval;
using Portage::Meshfree::Operator::Quadrilateral;
using Portage::Meshfree::Operator::Triangle;
using Portage::Meshfree::Operator::Hexahedron;
using Portage::Meshfree::Operator::Wedge;
using Portage::Meshfree::Operator::Tetrahedron;

using std::vector;

vector<Point<1>> interval_points_=
  {{0},{1}};
vector<Point<2>> quadrilateral_points_=
  {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
vector<Point<2>> triangle_points_=
  {{0, 0}, {1, 0}, {0, 1}};
vector<Point<3>> hexahedron_points_=
  {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
vector<Point<3>> wedge_points_=
  {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
vector<Point<3>> tetrahedron_points_=
  {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

// Basic tests on unit interval, square, cube

TEST(VolumeIntegral, UnitaryIntervalBasic) {
vector<vector<double>> result;
Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Interval>>(interval_points_, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.0);
}

TEST(VolumeIntegral, LinearIntervalBasic) {
vector<vector<double>> result;
Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Interval>>(interval_points_, result);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.0);
  ASSERT_EQ(result[1][0], 0.5);
}

TEST(VolumeIntegral, QuadraticIntervalBasic) {
vector<vector<double>> result;
Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Interval>>(interval_points_, result);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.0);
  ASSERT_EQ(result[1][0], 0.5);
  ASSERT_EQ(result[2][0], 1.0/6.0);
}

TEST(VolumeIntegral, UnitaryQuadrilateralBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Quadrilateral>>(quadrilateral_points_, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.0);
}

TEST(VolumeIntegral, LinearQuadrilateralBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Quadrilateral>>(quadrilateral_points_, result);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.0);
  ASSERT_EQ(result[1][0], 0.5);
  ASSERT_EQ(result[2][0], 0.5);
}

TEST(VolumeIntegral, QuadraticQuadrilateralBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Quadrilateral>>(quadrilateral_points_, result);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.0);
  ASSERT_EQ(result[1][0], 0.5);
  ASSERT_EQ(result[2][0], 0.5);
  ASSERT_EQ(result[3][0], 1.0/6.0);
  ASSERT_EQ(result[4][0], 0.25);
  ASSERT_EQ(result[5][0], 1.0/6.0);
}

TEST(VolumeIntegral, UnitaryTriangleBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Triangle>>(triangle_points_, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 0.5);
}

TEST(VolumeIntegral, LinearTriangleBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Triangle>>(triangle_points_, result);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 0.5);
  ASSERT_EQ(result[1][0], 1./6.);
  ASSERT_EQ(result[2][0], 1./6.);
}

TEST(VolumeIntegral, QuadraticTriangleBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Triangle>>(triangle_points_, result);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 0.5);
  ASSERT_EQ(result[1][0], 1./6.);
  ASSERT_EQ(result[2][0], 1./6.);
  ASSERT_EQ(result[3][0], 1./24.);
  ASSERT_EQ(result[4][0], 1./24.);
  ASSERT_EQ(result[5][0], 1./24.);
}

TEST(VolumeIntegral, UnitaryHexahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Hexahedron>>(hexahedron_points_, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.);
}

TEST(VolumeIntegral, LinearHexahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Hexahedron>>(hexahedron_points_, result);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1.0);
  ASSERT_EQ(result[1][0], 0.5);
  ASSERT_EQ(result[2][0], 0.5);
  ASSERT_EQ(result[3][0], 0.5);
}

/*
TEST(VolumeIntegral, UnitaryWedgeBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Wedge>>(wedge_points_, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 0.5);
}

TEST(VolumeIntegral, LinearWedgeBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Wedge>>(wedge_points_, result);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 0.5);
  ASSERT_EQ(result[1][0], 1./6.);
  ASSERT_EQ(result[2][0], 1./6.);
  ASSERT_EQ(result[3][0], 0.25);
}

TEST(VolumeIntegral, UnitaryTetrahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Tetrahedron>>(tetrahedron_points_, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1./6.);
}

TEST(VolumeIntegral, LinearTetrahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Tetrahedron>>(tetrahedron_points_, result);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_EQ(result[0][0], 1./6.);
  ASSERT_EQ(result[1][0], 1./24.);
  ASSERT_EQ(result[2][0], 1./24.);
  ASSERT_EQ(result[3][0], 1./24.);
}
*/
