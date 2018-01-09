/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

#include <vector>
#include <cmath>

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
using std::array;

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

template<size_t dim>
vector<Point<dim>> shift_points(const vector<Point<dim>> &points, 
				const Point<dim> &shift) {
  auto result(points);
  for (int i=0; i<points.size(); i++) {
    for (int j=0; j<dim; j++) {
      result[i][j] = points[i][j] + shift[j];
    }
  }
  return result;
}
/*
array<array<double,2>,2> matrix2 = 
  {{2.8549186535902855, -4.346198279115005}, {0.14208725143260595, 0.0933338790596824}};
array<array<double,3>,3> matrix3 = 
  {{3.1611889909865836, -3.1727215693209625, -2.6421056009990864}, 
   {0.0636728533375156, 0.13338461548842906, -0.0839899523685015}, 
   {1.4212135017018008, 0.22338659728810717, 1.4321838606591486}};
*/
vector<vector<double>> matrix2 = 
  {{2.8549186535902855, -4.346198279115005}, {0.14208725143260595, 0.0933338790596824}};
vector<vector<double>> matrix3 = 
  {{3.1611889909865836, -3.1727215693209625, -2.6421056009990864}, 
   {0.0636728533375156, 0.13338461548842906, -0.0839899523685015}, 
   {1.4212135017018008, 0.22338659728810717, 1.4321838606591486}};
double determinant2 = 0.884;
double determinant3 = 1.79452;

template<size_t dim>
vector<Point<dim>> deform_points(const vector<Point<dim>> &points, 
				 const vector<vector<double>> &matrix) {
  auto result(points);
  for (int i=0; i<points.size(); i++) {
    for (int j=0; j<dim; j++) {
      result[i][j] = 0.;
    }
  }
  for (int i=0; i<points.size(); i++) {
    for (int j=0; j<dim; j++) {
      for (int k=0; k<dim; k++) {
	result[i][j] += matrix[j][k]*points[i][k];
      }
    }
  }
  return result;
}

Point<1> shift1d={1.e8};
Point<2> shift2d={1.e8,2.e8};
Point<3> shift3d={1.e8,2.e8,3.e8};

typename Portage::Meshfree::Basis::Traits<Unitary,1>::array_t exactUnitaryInterval={1.0};
typename Portage::Meshfree::Basis::Traits<Unitary,2>::array_t exactUnitaryQuadrilateral={1.0};
typename Portage::Meshfree::Basis::Traits<Unitary,2>::array_t exactUnitaryTriangle={0.5};
typename Portage::Meshfree::Basis::Traits<Unitary,3>::array_t exactUnitaryHexahedron={1.0};
typename Portage::Meshfree::Basis::Traits<Unitary,3>::array_t exactUnitaryWedge={0.5};
typename Portage::Meshfree::Basis::Traits<Unitary,3>::array_t exactUnitaryTetrahedron={1./6.};

typename Portage::Meshfree::Basis::Traits<Linear,1>::array_t exactLinearInterval=
  {1.0, 0.5};
typename Portage::Meshfree::Basis::Traits<Linear,2>::array_t exactLinearQuadrilateral=
  {1.0, 0.5, 0.5};
typename Portage::Meshfree::Basis::Traits<Linear,2>::array_t exactLinearTriangle=
  {0.5, 1./6., 1./6.};
typename Portage::Meshfree::Basis::Traits<Linear,3>::array_t exactLinearHexahedron=
  {1.0, 0.5, 0.5, 0.5};
typename Portage::Meshfree::Basis::Traits<Linear,3>::array_t exactLinearWedge=
  {0.5, 1./6., 1./6., 1./4.};
typename Portage::Meshfree::Basis::Traits<Linear,3>::array_t exactLinearTetrahedron=
  {1./6., 1./24., 1./24., 1./24.};

typename Portage::Meshfree::Basis::Traits<Quadratic,1>::array_t exactQuadraticInterval=
  {1.0, 0.5, 1./6.};
typename Portage::Meshfree::Basis::Traits<Quadratic,2>::array_t exactQuadraticQuadrilateral=
  {1.0, 0.5, 0.5, 1./6., 1./4., 1./6.};
typename Portage::Meshfree::Basis::Traits<Quadratic,2>::array_t exactQuadraticTriangle=
  {0.5, 1./6., 1./6., 1./24., 1./24., 1./24.};
typename Portage::Meshfree::Basis::Traits<Quadratic,3>::array_t exactQuadraticTetrahedron=
  {1./6., 1./24., 1./24., 1./24., 1./120., 1./120., 1./120., 1./120, 1./120., 1./120.};

template<Portage::Meshfree::Basis::Type type, size_t dim>
typename Portage::Meshfree::Basis::Traits<type,dim>::array_t 
makeTranslatedExact(
  typename Portage::Meshfree::Basis::Traits<type,dim>::array_t &values,
  Point<dim> &point) 
{
  auto tf = Portage::Meshfree::Basis::transfactor<dim>(type, point);
  typename Portage::Meshfree::Basis::Traits<type,dim>::array_t tex;
  for (int i=0; i<tf.size(); i++) tex[i] = 0.;
  for (int i=0; i<tf.size(); i++) 
    for (int j=0; j<tf.size(); j++) 
      tex[i] += tf[i][j]*values[j];
  return tex;
}

// Basic tests on unit interval, square, cube, wedge, tet

TEST(VolumeIntegral, UnitaryIntervalBasic) {
vector<vector<double>> result;
Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryInterval[0], 1.e-12);
}

TEST(VolumeIntegral, LinearIntervalBasic) {
vector<vector<double>> result;
Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactLinearInterval[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearInterval[1], 1.e-12);
}

TEST(VolumeIntegral, QuadraticIntervalBasic) {
vector<vector<double>> result;
Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactQuadraticInterval[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactQuadraticInterval[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactQuadraticInterval[2], 1.e-12);
}

TEST(VolumeIntegral, UnitaryQuadrilateralBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryQuadrilateral[0], 1.e-12);
}

TEST(VolumeIntegral, LinearQuadrilateralBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactLinearQuadrilateral[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearQuadrilateral[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactLinearQuadrilateral[2], 1.e-12);
}

TEST(VolumeIntegral, QuadraticQuadrilateralBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactQuadraticQuadrilateral[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactQuadraticQuadrilateral[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactQuadraticQuadrilateral[2], 1.e-12);
  ASSERT_NEAR(result[3][0], exactQuadraticQuadrilateral[3], 1.e-12);
  ASSERT_NEAR(result[4][0], exactQuadraticQuadrilateral[4], 1.e-12);
  ASSERT_NEAR(result[5][0], exactQuadraticQuadrilateral[5], 1.e-12);
}

TEST(VolumeIntegral, UnitaryTriangleBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryTriangle[0], 1.e-12);
}

TEST(VolumeIntegral, LinearTriangleBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactLinearTriangle[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearTriangle[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactLinearTriangle[2], 1.e-12);
}

TEST(VolumeIntegral, QuadraticTriangleBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactQuadraticTriangle[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactQuadraticTriangle[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactQuadraticTriangle[2], 1.e-12);
  ASSERT_NEAR(result[3][0], exactQuadraticTriangle[3], 1.e-12);
  ASSERT_NEAR(result[4][0], exactQuadraticTriangle[4], 1.e-12);
  ASSERT_NEAR(result[5][0], exactQuadraticTriangle[5], 1.e-12);
}

TEST(VolumeIntegral, UnitaryHexahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Hexahedron>>(hexahedron_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryHexahedron[0], 1.e-12);
}

TEST(VolumeIntegral, LinearHexahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Hexahedron>>(hexahedron_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactLinearHexahedron[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearHexahedron[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactLinearHexahedron[2], 1.e-12);
  ASSERT_NEAR(result[3][0], exactLinearHexahedron[3], 1.e-12);
}

TEST(VolumeIntegral, UnitaryWedgeBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Wedge>>(wedge_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryWedge[0], 1.e-12);
}

TEST(VolumeIntegral, LinearWedgeBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Wedge>>(wedge_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactLinearWedge[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearWedge[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactLinearWedge[2], 1.e-12);
  ASSERT_NEAR(result[3][0], exactLinearWedge[3], 1.e-12);
}

TEST(VolumeIntegral, UnitaryTetrahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Unitary, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryTetrahedron[0], 1.e-12);
}

TEST(VolumeIntegral, LinearTetrahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Linear, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactLinearTetrahedron[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearTetrahedron[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactLinearTetrahedron[2], 1.e-12);
  ASSERT_NEAR(result[3][0], exactLinearTetrahedron[3], 1.e-12);
}

TEST(VolumeIntegral, QuadraticTetrahedronBasic) {
vector<vector<double>> result;
 Portage::Meshfree::Operator::get_result<Operator<VolumeIntegral, Quadratic, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactQuadraticTetrahedron[0], 1.e-12);
  ASSERT_NEAR(result[1][0], exactQuadraticTetrahedron[1], 1.e-12);
  ASSERT_NEAR(result[2][0], exactQuadraticTetrahedron[2], 1.e-12);
  ASSERT_NEAR(result[3][0], exactQuadraticTetrahedron[3], 1.e-12);
  ASSERT_NEAR(result[4][0], exactQuadraticTetrahedron[4], 1.e-12);
  ASSERT_NEAR(result[5][0], exactQuadraticTetrahedron[5], 1.e-12);
  ASSERT_NEAR(result[6][0], exactQuadraticTetrahedron[6], 1.e-12);
  ASSERT_NEAR(result[7][0], exactQuadraticTetrahedron[7], 1.e-12);
  ASSERT_NEAR(result[8][0], exactQuadraticTetrahedron[8], 1.e-12);
  ASSERT_NEAR(result[9][0], exactQuadraticTetrahedron[9], 1.e-12);
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
