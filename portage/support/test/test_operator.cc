/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <cmath>
#include "gtest/gtest.h"

#include "portage/support/operator.h"
#include "wonton/support/Point.h"

using Wonton::Point;

// avoid compiler confusion on 'Operator'
template<Portage::Meshfree::Operator::Type type,
         Portage::Meshfree::Basis::Type basis_type,
         Portage::Meshfree::Operator::Domain domain_type>
using OP = Portage::Meshfree::Operator::Operator<type, basis_type, domain_type>;

// avoid very long namespaces
using namespace Portage::Meshfree;
using Basis::Unitary;
using Basis::Linear;
using Basis::Quadratic;
using Operator::VolumeIntegral;
using Operator::Interval;
using Operator::Quadrilateral;
using Operator::Triangle;
using Operator::Hexahedron;
using Operator::Tetrahedron;
using Operator::Wedge;
using Operator::dimension;
using Operator::get_result;
using Operator::apply;

/**
 * @brief Fixture class for operator tests.
 */
class OperatorTest : public testing::Test {
protected:
  // reference points for different domains
  std::vector<Point<1>> const interval_points_ =
    {{0}, {1}};
  std::vector<Point<2>> const quadrilateral_points_ =
    {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  std::vector<Point<2>> const triangle_points_ =
    {{0, 0}, {1, 0}, {0, 1}};
  std::vector<Point<3>> const hexahedron_points_ =
    {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
     {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  std::vector<Point<3>> const wedge_points_ =
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
  std::vector<Point<3>> const tetrahedron_points_ =
    {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

  // accessor for different reference points
  template<Operator::Domain domain>
  std::vector<Point<dimension(domain)>> reference_points() {
    switch (domain) {
      case Interval:      return interval_points_;
      case Quadrilateral: return quadrilateral_points_;
      case Triangle:      return triangle_points_;
      case Hexahedron:    return hexahedron_points_;
      case Wedge:         return wedge_points_;
      case Tetrahedron:   return tetrahedron_points_;
      default: throw std::runtime_error("invalid domain");
    }
  };

  // function to shift reference points
  template<int dim>
  std::vector<Point<dim>> shift_points(const std::vector<Point<dim>>& points,
                                       const Point<dim>& shift) {
    auto result(points);
    for (int i = 0; i < points.size(); i++) {
      for (int j = 0; j < dim; j++) {
        result[i][j] = points[i][j] + shift[j];
      }
    }
    return result;
  }

  // standard shift std::vectors
  Point<1> const shift1d = {1.e8};
  Point<2> const shift2d = {1.e8, 2.e8};
  Point<3> const shift3d = {1.e8, 2.e8, 3.e8};

  // standard deformation std::vectors
  std::vector<std::vector<double>> const matrix2 = {
    {2.8549186535902855,  -4.346198279115005},
    {0.14208725143260595, 0.0933338790596824}
  };

  std::vector<std::vector<double>> const matrix3 = {
    {3.1611889909865836, -3.1727215693209625, -2.6421056009990864},
    {0.0636728533375156, 0.13338461548842906, -0.0839899523685015},
    {1.4212135017018008, 0.22338659728810717, 1.4321838606591486}
  };

  double const determinant2 = 0.884;
  double const determinant3 = 1.79452;

  // function to deform reference points
  template<int dim>
  std::vector<Point<dim>> deform_points(const std::vector<Point<dim>>& points,
                                        const std::vector<std::vector<double>>& matrix) {
    auto result(points);
    for (int i = 0; i < points.size(); i++) {
      for (int j = 0; j < dim; j++) {
        result[i][j] = 0.;
        for (int k = 0; k < dim; k++) {
          result[i][j] += matrix[j][k] * points[i][k];
        }
      }
    }
    return result;
  }

  // exact results for integrations of bases over different domains
  typename Basis::Traits<Unitary, 1>::array_t exactUnitaryInterval = {1.0};
  typename Basis::Traits<Unitary, 2>::array_t exactUnitaryQuadrilateral = {1.0};
  typename Basis::Traits<Unitary, 2>::array_t exactUnitaryTriangle = {0.5};
  typename Basis::Traits<Unitary, 3>::array_t exactUnitaryHexahedron = {1.0};
  typename Basis::Traits<Unitary, 3>::array_t exactUnitaryWedge = {0.5};
  typename Basis::Traits<Unitary, 3>::array_t exactUnitaryTetrahedron = {1./6.};

  typename Basis::Traits<Linear, 1>::array_t exactLinearInterval = {1.0, 0.5};
  typename Basis::Traits<Linear, 2>::array_t exactLinearQuadrilateral = {1.0, 0.5, 0.5};
  typename Basis::Traits<Linear, 2>::array_t exactLinearTriangle = {0.5, 1./6., 1./6.};
  typename Basis::Traits<Linear, 3>::array_t exactLinearHexahedron = {1.0, 0.5, 0.5, 0.5};
  typename Basis::Traits<Linear, 3>::array_t exactLinearWedge = {0.5, 1./6., 1./6., 1./4.};
  typename Basis::Traits<Linear, 3>::array_t exactLinearTetrahedron = {1./6., 1./24., 1./24., 1./24.};

  typename Basis::Traits<Quadratic, 1>::array_t exactQuadraticInterval = {1.0, 0.5, 1./6.};
  typename Basis::Traits<Quadratic, 2>::array_t exactQuadraticQuadrilateral = {1.0, 0.5, 0.5, 1./6., 1./4., 1./6.};
  typename Basis::Traits<Quadratic, 2>::array_t exactQuadraticTriangle = {0.5, 1./6., 1./6., 1./24., 1./24., 1./24.};
  typename Basis::Traits<Quadratic, 3>::array_t exactQuadraticTetrahedron =
    {1./6., 1./24., 1./24., 1./24., 1./120., 1./120., 1./120., 1./120, 1./120., 1./120.};

  // apply translation operator to exact results
  template<Basis::Type type, size_t dim>
  typename Basis::Traits<type, dim>::array_t
  makeTranslatedExact(typename Basis::Traits<type, dim>::array_t const& values,
                      Point<dim> const& point) {
    typename Basis::Traits<type, dim>::array_t tex;
    auto tf = Basis::transfactor<dim>(type, point);
    for (int i = 0; i < tf.size(); i++) {
      tex[i] = 0.;
      for (int j = 0; j < tf.size(); j++) {
        tex[i] += tf[i][j] * values[j];
      }
    }
    return tex;
  }
};

// Basic tests on unit interval, square, cube, wedge, tet
TEST_F(OperatorTest, UnitaryIntervalBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Unitary, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactUnitaryInterval[i], 1.e-12);
}

TEST_F(OperatorTest, LinearIntervalBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Linear, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactLinearInterval[i], 1.e-12);
  ASSERT_NEAR(result[1][0], exactLinearInterval[1], 1.e-12);
}

TEST_F(OperatorTest, QuadraticIntervalBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Quadratic, Interval>>(interval_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactQuadraticInterval[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryQuadrilateralBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Unitary, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactUnitaryQuadrilateral[i], 1.e-12);
}

TEST_F(OperatorTest, LinearQuadrilateralBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Linear, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactLinearQuadrilateral[i], 1.e-12);
}

TEST_F(OperatorTest, QuadraticQuadrilateralBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Quadratic, Quadrilateral>>(quadrilateral_points_, result, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactQuadraticQuadrilateral[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryTriangleBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Unitary, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactUnitaryTriangle[i], 1.e-12);
}

TEST_F(OperatorTest, LinearTriangleBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Linear, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactLinearTriangle[i], 1.e-12);
}

TEST_F(OperatorTest, QuadraticTriangleBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Quadratic, Triangle>>(triangle_points_, result, false);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactQuadraticTriangle[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryHexahedronBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Unitary, Hexahedron>>(hexahedron_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactUnitaryHexahedron[i], 1.e-12);
}

TEST_F(OperatorTest, LinearHexahedronBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Linear, Hexahedron>>(hexahedron_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactLinearHexahedron[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryWedgeBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Unitary, Wedge>>(wedge_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactUnitaryWedge[i], 1.e-12);
}

TEST_F(OperatorTest, LinearWedgeBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Linear, Wedge>>(wedge_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactLinearWedge[i], 1.e-12);
}

TEST_F(OperatorTest, UnitaryTetrahedronBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Unitary, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactUnitaryTetrahedron[i], 1.e-12);
}

TEST_F(OperatorTest, LinearTetrahedronBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Linear, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactLinearTetrahedron[i], 1.e-12);
}

TEST_F(OperatorTest, QuadraticTetrahedronBasic) {
  std::vector<std::vector<double>> result;
  get_result<OP<VolumeIntegral, Quadratic, Tetrahedron>>(tetrahedron_points_, result, false);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < result.size(); i++)
    ASSERT_NEAR(result[i][0], exactQuadraticTetrahedron[i], 1.e-12);
}

// Test transfactor operation on unit interval, square, cube, wedge, tet

TEST_F(OperatorTest, UnitaryIntervalTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  get_result<OP<VolumeIntegral, Unitary, Interval>>(points, result);
  auto tex = makeTranslatedExact<Unitary, 1>(exactUnitaryInterval, shift1d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearIntervalTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  get_result<OP<VolumeIntegral, Linear, Interval>>(points, result);
  auto tex = makeTranslatedExact<Linear, 1>(exactLinearInterval, shift1d);
  ASSERT_EQ(result.size(), 2);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticIntervalTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  get_result<OP<VolumeIntegral, Quadratic, Interval>>(points, result);
  auto tex = makeTranslatedExact<Quadratic, 1>(exactQuadraticInterval, shift1d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryQuadrilateralTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<2>(quadrilateral_points_, shift2d);
  get_result<OP<VolumeIntegral, Unitary, Quadrilateral>>(points, result);
  auto tex = makeTranslatedExact<Unitary, 2>(exactUnitaryQuadrilateral, shift2d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearQuadrilateralTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<2>(quadrilateral_points_, shift2d);
  get_result<OP<VolumeIntegral, Linear, Quadrilateral>>(points, result);
  auto tex = makeTranslatedExact<Linear, 2>(exactLinearQuadrilateral, shift2d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticQuadrilateralTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<2>(quadrilateral_points_, shift2d);
  get_result<OP<VolumeIntegral, Quadratic, Quadrilateral>>(points, result);
  auto tex = makeTranslatedExact<Quadratic, 2>(exactQuadraticQuadrilateral, shift2d);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryTriangleTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<2>(triangle_points_, shift2d);
  get_result<OP<VolumeIntegral, Unitary, Triangle>>(points, result);
  auto tex = makeTranslatedExact<Unitary, 2>(exactUnitaryTriangle, shift2d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearTriangleTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<2>(triangle_points_, shift2d);
  get_result<OP<VolumeIntegral, Linear, Triangle>>(points, result);
  auto tex = makeTranslatedExact<Linear, 2>(exactLinearTriangle, shift2d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticTriangleTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<2>(triangle_points_, shift2d);
  get_result<OP<VolumeIntegral, Quadratic, Triangle>>(points, result);
  auto tex = makeTranslatedExact<Quadratic, 2>(exactQuadraticTriangle, shift2d);
  ASSERT_EQ(result.size(), 6);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryHexahedronTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(hexahedron_points_, shift3d);
  get_result<OP<VolumeIntegral, Unitary, Hexahedron>>(points, result);
  auto tex = makeTranslatedExact<Unitary, 3>(exactUnitaryHexahedron, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearHexahedronTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(hexahedron_points_, shift3d);
  get_result<OP<VolumeIntegral, Linear, Hexahedron>>(points, result);
  auto tex = makeTranslatedExact<Linear, 3>(exactLinearHexahedron, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryWedgeTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(wedge_points_, shift3d);
  get_result<OP<VolumeIntegral, Unitary, Wedge>>(points, result);
  auto tex = makeTranslatedExact<Unitary, 3>(exactUnitaryWedge, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearWedgeTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(wedge_points_, shift3d);
  get_result<OP<VolumeIntegral, Linear, Wedge>>(points, result);
  auto tex = makeTranslatedExact<Linear, 3>(exactLinearWedge, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryTetrahedronTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  get_result<OP<VolumeIntegral, Unitary, Tetrahedron>>(points, result);
  auto tex = makeTranslatedExact<Unitary, 3>(exactUnitaryTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, LinearTetrahedronTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  get_result<OP<VolumeIntegral, Linear, Tetrahedron>>(points, result);
  auto tex = makeTranslatedExact<Linear, 3>(exactLinearTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 4);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, QuadraticTetrahedronTF) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  get_result<OP<VolumeIntegral, Quadratic, Tetrahedron>>(points, result);
  auto tex = makeTranslatedExact<Quadratic, 3>(exactQuadraticTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

// Test correctness of volume integrals under rotations and stretches

TEST_F(OperatorTest, UnitaryQuadrilateralDeform) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<2>(quadrilateral_points_, matrix2);
  get_result<OP<VolumeIntegral, Unitary, Quadrilateral>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryQuadrilateral[0] * determinant2, 1.e-12);
}

TEST_F(OperatorTest, UnitaryTriangleDeform) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<2>(triangle_points_, matrix2);
  get_result<OP<VolumeIntegral, Unitary, Triangle>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryTriangle[0] * determinant2, 1.e-12);
}

TEST_F(OperatorTest, UnitaryHexahedronDeform) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<3>(hexahedron_points_, matrix3);
  get_result<OP<VolumeIntegral, Unitary, Hexahedron>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryHexahedron[0] * determinant3, 1.e-12);
}

TEST_F(OperatorTest, UnitaryWedgeDeform) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<3>(wedge_points_, matrix3);
  get_result<OP<VolumeIntegral, Unitary, Wedge>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryWedge[0] * determinant3, 1.e-12);
}

TEST_F(OperatorTest, UnitaryTetrahedronDeform) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<3>(tetrahedron_points_, matrix3);
  get_result<OP<VolumeIntegral, Unitary, Tetrahedron>>(points, result);
  ASSERT_NEAR(result[0][0], exactUnitaryTetrahedron[0] * determinant3, 1.e-12);
}

// Test dynamic interface.

TEST_F(OperatorTest, QuadraticIntervalTFDynamic) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<1>(interval_points_, shift1d);
  apply<1>(VolumeIntegral, Quadratic, Interval, points, result);
  auto tex = makeTranslatedExact<Quadratic, 1>(exactQuadraticInterval, shift1d);
  ASSERT_EQ(result.size(), 3);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}

TEST_F(OperatorTest, UnitaryQuadrilateralDeformDynamic) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<2>(quadrilateral_points_, matrix2);
  apply<2>(VolumeIntegral, Unitary, Quadrilateral, points, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryQuadrilateral[0] * determinant2, 1.e-12);
}

TEST_F(OperatorTest, UnitaryHexahedronDeformDynamic) {
  std::vector<std::vector<double>> result;
  auto points = deform_points<3>(hexahedron_points_, matrix3);
  apply<3>(VolumeIntegral, Unitary, Hexahedron, points, result);
  ASSERT_EQ(result.size(), 1);
  ASSERT_EQ(result[0].size(), 1);
  ASSERT_NEAR(result[0][0], exactUnitaryHexahedron[0] * determinant3, 1.e-12);
}

TEST_F(OperatorTest, QuadraticTetrahedronTFDynamic) {
  std::vector<std::vector<double>> result;
  auto points = shift_points<3>(tetrahedron_points_, shift3d);
  apply<3>(VolumeIntegral, Quadratic, Tetrahedron, points, result);
  auto tex = makeTranslatedExact<Quadratic, 3>(exactQuadraticTetrahedron, shift3d);
  ASSERT_EQ(result.size(), 10);
  ASSERT_EQ(result[0].size(), 1);
  for (int i = 0; i < tex.size(); i++)
    ASSERT_NEAR(result[i][0], tex[i], 1.e-7 * fabs(tex[i]));
}
