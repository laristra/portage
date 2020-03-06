/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "portage/support/test_operator_data.h"

using namespace Portage::swarm::oper;
using namespace Portage::swarm::basis;
using namespace Portage::swarm::reference;

class OperatorTest : public testing::Test {
protected:
  OperatorTest() {
    interval_      = points<Interval>();
    quadrilateral_ = points<Quadrilateral>();
    triangle_      = points<Triangle>();
    hexahedron_    = points<Hexahedron>();
    wedge_         = points<Wedge>();
    tetrahedron_   = points<Tetrahedron>();
  }

  std::vector<Point<1>> interval_;
  std::vector<Point<2>> quadrilateral_;
  std::vector<Point<2>> triangle_;
  std::vector<Point<3>> hexahedron_;
  std::vector<Point<3>> wedge_;
  std::vector<Point<3>> tetrahedron_;
};

// Basic tests on unit interval, square, cube, wedge, tet
TEST_F(OperatorTest, UnitaryIntervalBasic) {
  auto result = get_result<OP<VolumeIntegral, Unitary, Interval>>(interval_, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], unitary_interval[i], 1.e-12);
}

TEST_F(OperatorTest, LinearIntervalBasic) {
  auto result = get_result<OP<VolumeIntegral, Linear, Interval>>(interval_, false);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], linear_interval[i], 1.e-12);
  ASSERT_NEAR(result[1][0], linear_interval[1], 1.e-12);
}

TEST_F(OperatorTest, QuadraticIntervalBasic) {
  auto result = get_result<OP<VolumeIntegral, Quadratic, Interval>>(interval_, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], quadratic_interval[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryQuadrilateralBasic) {
  auto result = get_result<OP<VolumeIntegral, Unitary, Quadrilateral>>(quadrilateral_, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], unitary_quadrilateral[i], 1.e-12);
}

TEST_F(OperatorTest, LinearQuadrilateralBasic) {
  auto result = get_result<OP<VolumeIntegral, Linear, Quadrilateral>>(quadrilateral_, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], linear_quadrilateral[i], 1.e-12);
}

TEST_F(OperatorTest, QuadraticQuadrilateralBasic) {
  auto result = get_result<OP<VolumeIntegral, Quadratic, Quadrilateral>>(quadrilateral_, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], quadratic_quadrilateral[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryTriangleBasic) {
  auto result = get_result<OP<VolumeIntegral, Unitary, Triangle>>(triangle_, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], unitary_triangle[i], 1.e-12);
}

TEST_F(OperatorTest, LinearTriangleBasic) {
  auto result = get_result<OP<VolumeIntegral, Linear, Triangle>>(triangle_, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], linear_triangle[i], 1.e-12);
}

TEST_F(OperatorTest, QuadraticTriangleBasic) {
  auto result = get_result<OP<VolumeIntegral, Quadratic, Triangle>>(triangle_, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], quadratic_triangle[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryHexahedronBasic) {
  auto result = get_result<OP<VolumeIntegral, Unitary, Hexahedron>>(hexahedron_, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], unitary_hexahedron[i], 1.e-12);
}

TEST_F(OperatorTest, LinearHexahedronBasic) {
  auto result = get_result<OP<VolumeIntegral, Linear, Hexahedron>>(hexahedron_, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], linear_hexahedron[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryWedgeBasic) {
  auto result = get_result<OP<VolumeIntegral, Unitary, Wedge>>(wedge_, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], unitary_wedge[i], 1.e-12);
}

TEST_F(OperatorTest, LinearWedgeBasic) {
  auto result = get_result<OP<VolumeIntegral, Linear, Wedge>>(wedge_, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], linear_wedge[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryTetrahedronBasic) {
  auto result = get_result<OP<VolumeIntegral, Unitary, Tetrahedron>>(tetrahedron_, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], unitary_tetrahedron[i], 1.e-12);
}

TEST_F(OperatorTest, LinearTetrahedronBasic) {
  auto result = get_result<OP<VolumeIntegral, Linear, Tetrahedron>>(tetrahedron_, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], linear_tetrahedron[i], 1.e-12);
}

TEST_F(OperatorTest, QuadraticTetrahedronBasic) {
  auto result = get_result<OP<VolumeIntegral, Quadratic, Tetrahedron>>(tetrahedron_, false);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], quadratic_tetrahedron[i], 1.e-12);
}

// Test transfactor operation on unit interval, square, cube, wedge, tet

TEST_F(OperatorTest, UnitaryIntervalTF) {
  auto points = shift_points<1>(interval_, shift1d);
  auto result = get_result<OP<VolumeIntegral, Unitary, Interval>>(points);
  auto tex = make_translated_exact<Unitary, 1>(unitary_interval, shift1d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearIntervalTF) {
  auto points = shift_points<1>(interval_, shift1d);
  auto result = get_result<OP<VolumeIntegral, Linear, Interval>>(points);
  auto tex = make_translated_exact<Linear, 1>(linear_interval, shift1d);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticIntervalTF) {
  auto points = shift_points<1>(interval_, shift1d);
  auto result = get_result<OP<VolumeIntegral, Quadratic, Interval>>(points);
  auto tex = make_translated_exact<Quadratic, 1>(quadratic_interval, shift1d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryQuadrilateralTF) {
  auto points = shift_points<2>(quadrilateral_, shift2d);
  auto result = get_result<OP<VolumeIntegral, Unitary, Quadrilateral>>(points);
  auto tex = make_translated_exact<Unitary, 2>(unitary_quadrilateral, shift2d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearQuadrilateralTF) {
  auto points = shift_points<2>(quadrilateral_, shift2d);
  auto result = get_result<OP<VolumeIntegral, Linear, Quadrilateral>>(points);
  auto tex = make_translated_exact<Linear, 2>(linear_quadrilateral, shift2d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticQuadrilateralTF) {
  auto points = shift_points<2>(quadrilateral_, shift2d);
  auto result = get_result<OP<VolumeIntegral, Quadratic, Quadrilateral>>(points);
  auto tex = make_translated_exact<Quadratic, 2>(quadratic_quadrilateral, shift2d);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryTriangleTF) {
  auto points = shift_points<2>(triangle_, shift2d);
  auto result = get_result<OP<VolumeIntegral, Unitary, Triangle>>(points);
  auto tex = make_translated_exact<Unitary, 2>(unitary_triangle, shift2d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearTriangleTF) {
  auto points = shift_points<2>(triangle_, shift2d);
  auto result = get_result<OP<VolumeIntegral, Linear, Triangle>>(points);
  auto tex = make_translated_exact<Linear, 2>(linear_triangle, shift2d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticTriangleTF) {
  auto points = shift_points<2>(triangle_, shift2d);
  auto result = get_result<OP<VolumeIntegral, Quadratic, Triangle>>(points);
  auto tex = make_translated_exact<Quadratic, 2>(quadratic_triangle, shift2d);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryHexahedronTF) {
  auto points = shift_points<3>(hexahedron_, shift3d);
  auto result = get_result<OP<VolumeIntegral, Unitary, Hexahedron>>(points);
  auto tex = make_translated_exact<Unitary, 3>(unitary_hexahedron, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearHexahedronTF) {
  auto points = shift_points<3>(hexahedron_, shift3d);
  auto result = get_result<OP<VolumeIntegral, Linear, Hexahedron>>(points);
  auto tex = make_translated_exact<Linear, 3>(linear_hexahedron, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryWedgeTF) {
  auto points = shift_points<3>(wedge_, shift3d);
  auto result = get_result<OP<VolumeIntegral, Unitary, Wedge>>(points);
  auto tex = make_translated_exact<Unitary, 3>(unitary_wedge, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearWedgeTF) {
  auto points = shift_points<3>(wedge_, shift3d);
  auto result = get_result<OP<VolumeIntegral, Linear, Wedge>>(points);
  auto tex = make_translated_exact<Linear, 3>(linear_wedge, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryTetrahedronTF) {
  auto points = shift_points<3>(tetrahedron_, shift3d);
  auto result = get_result<OP<VolumeIntegral, Unitary, Tetrahedron>>(points);
  auto tex = make_translated_exact<Unitary, 3>(unitary_tetrahedron, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearTetrahedronTF) {
  auto points = shift_points<3>(tetrahedron_, shift3d);
  auto result = get_result<OP<VolumeIntegral, Linear, Tetrahedron>>(points);
  auto tex = make_translated_exact<Linear, 3>(linear_tetrahedron, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticTetrahedronTF) {
  auto points = shift_points<3>(tetrahedron_, shift3d);
  auto result = get_result<OP<VolumeIntegral, Quadratic, Tetrahedron>>(points);
  auto tex = make_translated_exact<Quadratic, 3>(quadratic_tetrahedron, shift3d);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

// Test correctness of volume integrals under rotations and stretches

TEST_F(OperatorTest, UnitaryQuadrilateralDeform) {
  auto points = deform_points<2>(quadrilateral_, matrix2);
  auto result = get_result<OP<VolumeIntegral, Unitary, Quadrilateral>>(points);
  ASSERT_NEAR(result[0][0], unitary_quadrilateral[0] * determinant2, 1.e-12);
}

TEST_F(OperatorTest, UnitaryTriangleDeform) {
  auto points = deform_points<2>(triangle_, matrix2);
  auto result = get_result<OP<VolumeIntegral, Unitary, Triangle>>(points);
  ASSERT_NEAR(result[0][0], unitary_triangle[0] * determinant2, 1.e-12);
}

TEST_F(OperatorTest, UnitaryHexahedronDeform) {
  auto points = deform_points<3>(hexahedron_, matrix3);
  auto result = get_result<OP<VolumeIntegral, Unitary, Hexahedron>>(points);
  ASSERT_NEAR(result[0][0], unitary_hexahedron[0] * determinant3, 1.e-12);
}

TEST_F(OperatorTest, UnitaryWedgeDeform) {
  auto points = deform_points<3>(wedge_, matrix3);
  auto result = get_result<OP<VolumeIntegral, Unitary, Wedge>>(points);
  ASSERT_NEAR(result[0][0], unitary_wedge[0] * determinant3, 1.e-12);
}

TEST_F(OperatorTest, UnitaryTetrahedronDeform) {
  auto points = deform_points<3>(tetrahedron_, matrix3);
  auto result = get_result<OP<VolumeIntegral, Unitary, Tetrahedron>>(points);
  ASSERT_NEAR(result[0][0], unitary_tetrahedron[0] * determinant3, 1.e-12);
}

// Test dynamic interface.

TEST_F(OperatorTest, QuadraticIntervalTFDynamic) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<1>(interval_, shift1d);
  apply<1>(VolumeIntegral, Quadratic, Interval, points, result);
  auto tex = make_translated_exact<Quadratic, 1>(quadratic_interval, shift1d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryQuadrilateralDeformDynamic) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<2>(quadrilateral_, matrix2);
  apply<2>(VolumeIntegral, Unitary, Quadrilateral, points, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], unitary_quadrilateral[0] * determinant2, 1.e-12);
}

TEST_F(OperatorTest, UnitaryHexahedronDeformDynamic) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<3>(hexahedron_, matrix3);
  apply<3>(VolumeIntegral, Unitary, Hexahedron, points, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], unitary_hexahedron[0] * determinant3, 1.e-12);
}

TEST_F(OperatorTest, QuadraticTetrahedronTFDynamic) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(tetrahedron_, shift3d);
  apply<3>(VolumeIntegral, Quadratic, Tetrahedron, points, result);
  auto tex = make_translated_exact<Quadratic, 3>(quadratic_tetrahedron, shift3d);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}
