/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_SUPPORT_OPERATOR_REFERENCES_H_
#define PORTAGE_SUPPORT_OPERATOR_REFERENCES_H_

#include <vector>
#include <cmath>
#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "portage/support/operator.h"

using Wonton::Point;

namespace Portage { namespace Meshfree { namespace reference {

// avoid compiler confusion on 'Operator'
template<oper::Type type,
         basis::Type basis_type,
         oper::Domain domain_type>
using OP = oper::Operator<type, basis_type, domain_type>;

// reference points for different domains
template<oper::Domain domain>
constexpr std::vector<Point<oper::dimension(domain)>> points() {
  switch (domain) {
    case oper::Interval:      return {{0}, {1}};
    case oper::Quadrilateral: return {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    case oper::Triangle:      return {{0, 0}, {1, 0}, {0, 1}};
    case oper::Tetrahedron:   return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    case oper::Wedge:         return {{0, 0, 0}, {1, 0, 0}, {0, 1, 0},
                                      {0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
    case oper::Hexahedron:    return {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                      {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    default:
      throw std::runtime_error("invalid domain");
  }
};

// function to shift reference points
template<int dim>
std::vector<Point<dim>> shift_points(const std::vector<Point<dim>>& points,
                                     const Point<dim>& shift) {
  auto result(points);
  int const num_points = points.size();
  for (int i = 0; i < num_points; i++) {
    for (int j = 0; j < dim; j++) {
      result[i][j] = points[i][j] + shift[j];
    }
  }
  return result;
}

// standard shift vectors
Point<1> const shift1d = {1.e8};
Point<2> const shift2d = {1.e8, 2.e8};
Point<3> const shift3d = {1.e8, 2.e8, 3.e8};

// standard deformation vectors
std::vector<std::vector<double>> const matrix2 =
  {{2.8549186535902855,  -4.346198279115005},
   {0.14208725143260595, 0.0933338790596824}};

std::vector<std::vector<double>> const matrix3 =
  {{3.1611889909865836, -3.1727215693209625, -2.6421056009990864},
   {0.0636728533375156, 0.13338461548842906, -0.0839899523685015},
   {1.4212135017018008, 0.22338659728810717, 1.4321838606591486}};

double const determinant2 = 0.884;
double const determinant3 = 1.79452;

// function to deform reference points
template<int dim>
std::vector<Point<dim>> deform_points(const std::vector<Point<dim>>& points,
                                      const std::vector<std::vector<double>>& matrix) {
  auto result(points);
  int const num_points = points.size();
  for (int i = 0; i < num_points; i++) {
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
template<int dim> using B0 = typename basis::Traits<basis::Unitary, dim>::array_t;
template<int dim> using B1 = typename basis::Traits<basis::Linear, dim>::array_t;
template<int dim> using B2 = typename basis::Traits<basis::Quadratic, dim>::array_t;

B0<1> const unitary_interval        = {1.0};
B0<2> const unitary_quadrilateral   = {1.0};
B0<2> const unitary_triangle        = {0.5};
B0<3> const unitary_hexahedron      = {1.0};
B0<3> const unitary_wedge           = {0.5};
B0<3> const unitary_tetrahedron     = {1./6.};
B1<1> const linear_interval         = {1.0, 0.5};
B1<2> const linear_quadrilateral    = {1.0, 0.5, 0.5};
B1<2> const linear_triangle         = {0.5, 1./6., 1./6.};
B1<3> const linear_hexahedron       = {1.0, 0.5, 0.5, 0.5};
B1<3> const linear_wedge            = {0.5, 1./6., 1./6., 1./4.};
B1<3> const linear_tetrahedron      = {1./6., 1./24., 1./24., 1./24.};
B2<1> const quadratic_interval      = {1.0, 0.5, 1./6.};
B2<2> const quadratic_quadrilateral = {1.0, 0.5, 0.5, 1./6., 1./4., 1./6.};
B2<2> const quadratic_triangle      = {0.5, 1./6., 1./6., 1./24., 1./24., 1./24.};
B2<3> const quadratic_tetrahedron   = {1./6., 1./24., 1./24., 1./24., 1./120.,
                                       1./120., 1./120., 1./120, 1./120., 1./120.};

// apply translation operator to exact results
template<basis::Type type, size_t dim>
typename basis::Traits<type, dim>::array_t
make_translated_exact(typename basis::Traits<type, dim>::array_t const& values,
                      Point<dim> const& point) {
  typename basis::Traits<type, dim>::array_t tex {};
  auto tf = basis::transfactor<dim>(type, point);
  int const num_tf = tf.size();

  for (int i = 0; i < num_tf; i++) {
    tex[i] = 0.;
    for (int j = 0; j < num_tf; j++) {
      tex[i] += tf[i][j] * values[j];
    }
  }
  return tex;
}

}}} // namespace Portage::Meshfree::reference

#endif
