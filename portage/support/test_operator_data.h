/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <cmath>

#include "portage/support/operator.h"
#include "portage/support/portage.h"
#include "wonton/support/Point.h"

using Wonton::Point;
using Portage::Meshfree::Basis::Unitary;
using Portage::Meshfree::Basis::Linear;
using Portage::Meshfree::Basis::Quadratic;
using Portage::Meshfree::Operator::Interval;
using Portage::Meshfree::Operator::Quadrilateral;
using Portage::Meshfree::Operator::Triangle;
using Portage::Meshfree::Operator::Hexahedron;
using Portage::Meshfree::Operator::Tetrahedron;
using Portage::Meshfree::Operator::Wedge;
using Portage::Meshfree::Operator::dimension;

// reference points for different domains

std::vector<Point<1>> interval_points_=
  {{0},{1}};
std::vector<Point<2>> quadrilateral_points_=
  {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
std::vector<Point<2>> triangle_points_=
  {{0, 0}, {1, 0}, {0, 1}};
std::vector<Point<3>> hexahedron_points_=
  {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
std::vector<Point<3>> wedge_points_=
  {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
std::vector<Point<3>> tetrahedron_points_=
  {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};

// accessor for different reference points

template<Portage::Meshfree::Operator::Domain domain>
std::vector<Point<dimension(domain)>>
reference_points()
{};

template<>
std::vector<Point<dimension(Interval)>>
reference_points<Interval>()
{return interval_points_;}

template<>
std::vector<Point<dimension(Quadrilateral)>>
reference_points<Quadrilateral>()
{return quadrilateral_points_;}

template<>
std::vector<Point<dimension(Triangle)>>
reference_points<Triangle>()
{return triangle_points_;}

template<>
std::vector<Point<dimension(Hexahedron)>>
reference_points<Hexahedron>()
{return hexahedron_points_;}

template<>
std::vector<Point<dimension(Wedge)>>
reference_points<Wedge>()
{return wedge_points_;}

template<>
std::vector<Point<dimension(Tetrahedron)>>
reference_points<Tetrahedron>()
{return tetrahedron_points_;}

// function to shift reference points

template<size_t dim>
std::vector<Point<dim>> shift_points(const std::vector<Point<dim>> &points,
				const Point<dim> &shift) {
  auto result(points);
  for (int i=0; i<points.size(); i++) {
    for (int j=0; j<dim; j++) {
      result[i][j] = points[i][j] + shift[j];
    }
  }
  return result;
}

// standard shift std::vectors

Point<1> shift1d={1.e8};
Point<2> shift2d={1.e8,2.e8};
Point<3> shift3d={1.e8,2.e8,3.e8};

// standard deformation std::vectors

std::vector<std::vector<double>> matrix2 =
  {{2.8549186535902855, -4.346198279115005}, {0.14208725143260595, 0.0933338790596824}};
std::vector<std::vector<double>> matrix3 =
  {{3.1611889909865836, -3.1727215693209625, -2.6421056009990864},
   {0.0636728533375156, 0.13338461548842906, -0.0839899523685015},
   {1.4212135017018008, 0.22338659728810717, 1.4321838606591486}};
double determinant2 = 0.884;
double determinant3 = 1.79452;

// function to deform reference points

template<size_t dim>
std::vector<Point<dim>> deform_points(const std::vector<Point<dim>> &points,
				 const std::vector<std::vector<double>> &matrix) {
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

// exact results for integrations of bases over different domains

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

// apply translation operator to exact results

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
