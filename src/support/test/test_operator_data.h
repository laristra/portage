/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <vector>
#include <cmath>
using std::vector;
using std::array;

#include "portage/support/operator.h"
#include "portage/support/Point.h"
using Portage::Point;
using Portage::Meshfree::Basis::Unitary;
using Portage::Meshfree::Basis::Linear;
using Portage::Meshfree::Basis::Quadratic;

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

Point<1> shift1d={1.e8};
Point<2> shift2d={1.e8,2.e8};
Point<3> shift3d={1.e8,2.e8,3.e8};

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
