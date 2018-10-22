/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <cmath>

#include "gtest/gtest.h"

#include "portage/support/portage.h"
#include "portage/support/operator.h"
#include "portage/support/test_operator_data.h"

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

using Portage::Meshfree::Operator::size_info;

using std::vector;

// Basic tests on unit interval, square, cube, wedge, tet

TEST(VolumeIntegral, UnitaryIntervalBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactUnitaryInterval[i], 1.e-12);
}

TEST(VolumeIntegral, LinearIntervalBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactLinearInterval[i], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearInterval[1], 1.e-12);
}

TEST(VolumeIntegral, QuadraticIntervalBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactQuadraticInterval[i], 1.e-12);
}

TEST(VolumeIntegral, UnitaryQuadrilateralBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactUnitaryQuadrilateral[i], 1.e-12);
}

TEST(VolumeIntegral, LinearQuadrilateralBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactLinearQuadrilateral[i], 1.e-12);
}

TEST(VolumeIntegral, QuadraticQuadrilateralBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactQuadraticQuadrilateral[i], 1.e-12);
}

TEST(VolumeIntegral, UnitaryTriangleBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactUnitaryTriangle[i], 1.e-12);
}

TEST(VolumeIntegral, LinearTriangleBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactLinearTriangle[i], 1.e-12);
}

TEST(VolumeIntegral, QuadraticTriangleBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactQuadraticTriangle[i], 1.e-12);
}

TEST(VolumeIntegral, UnitaryHexahedronBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Hexahedron>>(hexahedron_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactUnitaryHexahedron[i], 1.e-12);
}

TEST(VolumeIntegral, LinearHexahedronBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Hexahedron>>(hexahedron_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactLinearHexahedron[i], 1.e-12);
}

TEST(VolumeIntegral, UnitaryWedgeBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Wedge>>(wedge_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactUnitaryWedge[i], 1.e-12);
}

TEST(VolumeIntegral, LinearWedgeBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Wedge>>(wedge_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactLinearWedge[i], 1.e-12);
}

TEST(VolumeIntegral, UnitaryTetrahedronBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactUnitaryTetrahedron[i], 1.e-12);
}

TEST(VolumeIntegral, LinearTetrahedronBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactLinearTetrahedron[i], 1.e-12);
}

TEST(VolumeIntegral, QuadraticTetrahedronBasic) {
  vector<vector<double>> result;
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<result.size(); i++) ASSERT_NEAR(result[i][0], exactQuadraticTetrahedron[i], 1.e-12);
}

// Test transfactor operation on unit interval, square, cube, wedge, tet

TEST(VolumeIntegral, UnitaryIntervalTF) {
  vector<vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Interval>>(points, result);
  auto tex = makeTranslatedExact<Unitary,1>(exactUnitaryInterval, shift1d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, LinearIntervalTF) {
  vector<vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Interval>>(points, result);
  auto tex = makeTranslatedExact<Linear,1>(exactLinearInterval, shift1d);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, QuadraticIntervalTF) {
  vector<vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Interval>>(points, result);
  auto tex = makeTranslatedExact<Quadratic,1>(exactQuadraticInterval, shift1d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, UnitaryQuadrilateralTF) {
  vector<vector<double>> result;
  auto points = shift_points<2>(quadrilateral_points_, shift2d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Quadrilateral>>(points, result);
  auto tex = makeTranslatedExact<Unitary,2>(exactUnitaryQuadrilateral, shift2d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, LinearQuadrilateralTF) {
  vector<vector<double>> result;
  auto points = shift_points<2>(quadrilateral_points_, shift2d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Quadrilateral>>(points, result);
  auto tex = makeTranslatedExact<Linear,2>(exactLinearQuadrilateral, shift2d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, QuadraticQuadrilateralTF) {
  vector<vector<double>> result;
  auto points = shift_points<2>(quadrilateral_points_, shift2d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Quadrilateral>>(points, result);
  auto tex = makeTranslatedExact<Quadratic,2>(exactQuadraticQuadrilateral, shift2d);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, UnitaryTriangleTF) {
  vector<vector<double>> result;
  auto points = shift_points<2>(triangle_points_, shift2d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Triangle>>(points, result);
  auto tex = makeTranslatedExact<Unitary,2>(exactUnitaryTriangle, shift2d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, LinearTriangleTF) {
  vector<vector<double>> result;
  auto points = shift_points<2>(triangle_points_, shift2d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Triangle>>(points, result);
  auto tex = makeTranslatedExact<Linear,2>(exactLinearTriangle, shift2d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, QuadraticTriangleTF) {
  vector<vector<double>> result;
  auto points = shift_points<2>(triangle_points_, shift2d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Triangle>>(points, result);
  auto tex = makeTranslatedExact<Quadratic,2>(exactQuadraticTriangle, shift2d);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, UnitaryHexahedronTF) {
  vector<vector<double>> result;
  auto points = shift_points<3>(hexahedron_points_, shift3d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Hexahedron>>(points, result);
  auto tex = makeTranslatedExact<Unitary,3>(exactUnitaryHexahedron, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, LinearHexahedronTF) {
  vector<vector<double>> result;
  auto points = shift_points<3>(hexahedron_points_, shift3d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Hexahedron>>(points, result);
  auto tex = makeTranslatedExact<Linear,3>(exactLinearHexahedron, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, UnitaryWedgeTF) {
  vector<vector<double>> result;
  auto points = shift_points<3>(wedge_points_, shift3d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Wedge>>(points, result);
  auto tex = makeTranslatedExact<Unitary,3>(exactUnitaryWedge, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, LinearWedgeTF) {
  vector<vector<double>> result;
  auto points = shift_points<3>(wedge_points_, shift3d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Wedge>>(points, result);
  auto tex = makeTranslatedExact<Linear,3>(exactLinearWedge, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, UnitaryTetrahedronTF) {
  vector<vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Tetrahedron>>(points, result);
  auto tex = makeTranslatedExact<Unitary,3>(exactUnitaryTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, LinearTetrahedronTF) {
  vector<vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Tetrahedron>>(points, result);
  auto tex = makeTranslatedExact<Linear,3>(exactLinearTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, QuadraticTetrahedronTF) {
  vector<vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Tetrahedron>>(points, result);
  auto tex = makeTranslatedExact<Quadratic,3>(exactQuadraticTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

// Test correctness of volume integrals under rotations and stretches

TEST(VolumeIntegral, UnitaryQuadrilateralDeform) {
  vector<vector<double>> result;
  auto points = deform_points<2>(quadrilateral_points_, matrix2);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Quadrilateral>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryQuadrilateral[0]*determinant2, 1.e-12);
}

TEST(VolumeIntegral, UnitaryTriangleDeform) {
  vector<vector<double>> result;
  auto points = deform_points<2>(triangle_points_, matrix2);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Triangle>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryTriangle[0]*determinant2, 1.e-12);
}

TEST(VolumeIntegral, UnitaryHexahedronDeform) {
  vector<vector<double>> result;
  auto points = deform_points<3>(hexahedron_points_, matrix3);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Hexahedron>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryHexahedron[0]*determinant3, 1.e-12);
}

TEST(VolumeIntegral, UnitaryWedgeDeform) {
  vector<vector<double>> result;
  auto points = deform_points<3>(wedge_points_, matrix3);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Wedge>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryWedge[0]*determinant3, 1.e-12);
}

TEST(VolumeIntegral, UnitaryTetrahedronDeform) {
  vector<vector<double>> result;
  auto points = deform_points<3>(tetrahedron_points_, matrix3);
  Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Tetrahedron>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryTetrahedron[0]*determinant3, 1.e-12);
}

// Test dynamic interface.

TEST(VolumeIntegral, QuadraticIntervalTFDynamic) {
  vector<vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  Portage::Meshfree::Operator::apply<1>(VolumeIntegral, Quadratic, Interval, points, result);
  auto tex = makeTranslatedExact<Quadratic,1>(exactQuadraticInterval, shift1d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}

TEST(VolumeIntegral, UnitaryQuadrilateralDeformDynamic) {
  vector<vector<double>> result;
  auto points = deform_points<2>(quadrilateral_points_, matrix2);
  Portage::Meshfree::Operator::apply<2>(VolumeIntegral, Unitary, Quadrilateral, points, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryQuadrilateral[0]*determinant2, 1.e-12);
}

TEST(VolumeIntegral, UnitaryHexahedronDeformDynamic) {
  vector<vector<double>> result;
  auto points = deform_points<3>(hexahedron_points_, matrix3);
  Portage::Meshfree::Operator::apply<3>(VolumeIntegral, Unitary, Hexahedron, points, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryHexahedron[0]*determinant3, 1.e-12);
}

TEST(VolumeIntegral, QuadraticTetrahedronTFDynamic) {
  vector<vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  Portage::Meshfree::Operator::apply<3>(VolumeIntegral, Quadratic, Tetrahedron, points, result);
  auto tex = makeTranslatedExact<Quadratic,3>(exactQuadraticTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i=0; i<tex.size(); i++) ASSERT_NEAR(result[i][0], tex[i], 1.e-7*fabs(tex[i]));
}
