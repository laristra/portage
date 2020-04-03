/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef PORTAGE_SUPPORT_WEIGHT_H_
#define PORTAGE_SUPPORT_WEIGHT_H_

#include <cmath>
#include <array>
#include <cassert>
#include <limits>
#include <vector>

// portage includes
#include "portage/support/portage.h"
#include "wonton/support/Point.h"

namespace Portage { namespace Meshfree { namespace Weight {

//\///////////////////////////////////////////////////////////////////////////
// constants and maths functions
//\///////////////////////////////////////////////////////////////////////////

/**
 * @brief normalization constants for cubic B-spline
 *        (linear, cylindrical, and spherical).
 */
double const normconst[4] = { 2./3., 1./(.7 * M_PI), 1./M_PI, 1./M_PI };

/**
 * @brief Numerical tolerance.
 */
double const epsilon = std::numeric_limits<double>::epsilon();

/**
 * @brief scalar sign function.
 *
 * @param x: given scalar.
 * @return its sign.
 */
inline double sign(double x) { return x > 0. ? 1. : (x < 0. ? -1. : 0.); }

/// scalar step function
inline double unit_step(double x) { return 0.5 * (1. + sign(x)); }

//\///////////////////////////////////////////////////////////////////////////
// various scalar kernels and derivatives
//\///////////////////////////////////////////////////////////////////////////

/**
 * @brief scalar cubic b-spline.
 *
 * @param x: given scalar
 * @return its value at x.
 */
inline double b4(double x) {
  return (std::pow(2. - std::abs(x), 3) * unit_step(2. - std::abs(x))
          * unit_step(-1. + std::abs(x))) / 4.
          + (1. - 1.5 * std::pow(x, 2)
          + 0.75 * std::pow(std::abs(x), 3)) * unit_step(1. - std::abs(x));
}

/**
 * @brief scalar cubic b-spline derivative.
 *
 * @param x: given scalar
 * @return its value at x.
 */
inline double db4(double x) {
  return (-3. * x + (9. * std::pow(x, 2) * sign(x)) / 4.)
         * unit_step(1. - std::abs(x)) - (3. * std::pow(2. - std::abs(x), 2)
         * sign(x) * unit_step(2. - std::abs(x)) * unit_step(-1. + std::abs(x))) / 4.;
}

/**
 * @brief scalar cubic b-spline second derivative.
 *
 * @param x: given scalar
 * @return its value at x.
 */
inline double ddb4(double x) {
  return (-3. + (9. * x * sign(x)) / 2.) * unit_step(1. - std::abs(x))
          + (3. * (2. - std::abs(x)) * unit_step(2. - std::abs(x))
          * unit_step(-1. + std::abs(x))) / 2.;
}

/**
 * @brief scalar cubic b-spline anti-derivative.
 *
 * @param x: given scalar
 * @return its value at x.
 */
inline double ib4(double x) {
  return unit_step(-2. + x) + 0.6666666666666666 *
         ((1.4375 + (0.25 - std::pow(-2. + x, 4) / 4.) / 4.) *
         unit_step(2. - x) * unit_step(-1. + x) +
         (0.75 + x - std::pow(x, 3) / 2. + (3. * std::pow(x, 4)) / 16.) *
         unit_step(1. - x) * unit_step(x) +
         (0.75 + x - std::pow(x, 3) / 2. - (3. * std::pow(x, 4)) / 16.) *
         unit_step(-x) * unit_step(1. + x) +
         (std::pow(2. + x, 4) * unit_step(-1. - x) * unit_step(2. + x)) / 16.);
}

/**
 * @brief scalar left half cubic b-spline.
 *
 * @param x: given scalar
 * @return its value at x.
 */
inline double b4lh(double x) {
  return 0.0625 * std::pow(2. - std::abs(x), 3) * (1. + sign(-1. - x)) *
         (1. + sign(2. + x)) +
         0.5 * (1. - (3. * std::pow(x, 2)) / 2. + (3. * std::pow(std::abs(x), 3)) / 4.) *
         (1. + sign(1. + x)) * unit_step(-x);
}

/**
 * @brief scalar right half cubic b-spline.
 *
 * @param x: given scalar
 * @return its value at x.
 */
inline double b4rh(double x) {
  return 0.0625 * std::pow(2. - std::abs(x), 3) * (1. + sign(2. - x)) *
         (1. + sign(-1. + x)) +
         0.5 * (1. - (3. * std::pow(x, 2)) / 2. + (3. * std::pow(std::abs(x), 3)) / 4.) *
         (1. + sign(1. - x)) * unit_step(x);
}

/**
 * @brief scalar epanechnikov kernel.
 *
 * @param x: given scalar
 * @return its value at x.
 */
inline double epanechnikov(double x) {
  return std::abs(x) >= 2. ? 0. : 0.375 * (1. - x * x * 0.25);
}

/**
 * @brief scalar epanechnikov kernel derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double depanechnikov(double x) { return std::abs(x) >= 2. ? 0. : -.1875 * x; }

/**
 * @brief scalar epanechnikov kernel second derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double ddepanechnikov(double x) { return std::abs(x) >= 2. ? 0. : -.1875; }

/**
 * @brief scalar square kernel.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double square(double x) { return std::abs(x) <= 2. ? 1. : 0.; }

/**
 * @brief scalar square kernel derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double dsquare(double x) { return 0.; }

/**
 * @brief scalar square kernel second derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double ddsquare(double x) { return 0.; }

/**
 * @brief scalar smooth ramp for faceted weight.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double polyramp(double x) {
  return (x < 0. ? 1.
                 : (((1.5 - x) * (1. + sign(1. - x))) * .5
                 + ((2. + (-2. + .5 * x) * x) * (1. + sign(2. - x))
                 * (1. + sign(-1. + x))) * .25) / 1.5);  // normalize to 1 at x=0
}

/**
 * @brief scalar smooth ramp for faceted weight derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double dpolyramp(double x) {
  return (((-1.) * (1. + sign(1. - x))) * .5 +
         ((-2. + x) * (1. + sign(2. - x)) * (1. + sign(-1. + x))) / 4.) / 1.5;
}

/**
 * @brief scalar smooth ramp for faceted weight derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double ddpolyramp(double x) {
  return (((0.) * (1. + sign(1. - x))) * .5 +
         ((1.) * (1. + sign(2. - x)) * (1. + sign(-1. + x))) / 4.) / 1.5;
}

/**
 * @brief inverse square root kernel.
 * 
 * @param x: given scalar.
 * @return its value at x.
 */
inline double invsqrt(double x) {
  double const ax = std::abs(x);
  return 0.5 * (1. + sign(2. - ax)) * ((ax - 2.) * ax + 4.) *
         pow(ax + epsilon, -0.5);
}

/**
 * @brief derivative of inverse square root kernel.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double dinvsqrt(double x) {
  double const sx = 1.; //x > 0 ? 1. : 1.;
  double const ax = std::abs(x);
  return 0.25 * (1. + sign(2. - ax)) * sx * ((3. * ax - 4.) * ax - 4.) *
         std::pow(ax + epsilon, -1.5);
}

/**
 * @brief second derivative of inverse square root kernel.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double ddinvsqrt(double x) {
  double const ax = std::abs(x);
  return 0.125 * (1. + sign(2. - ax)) * ((3. * ax + 4.) * ax + 12.) *
         std::pow(ax + epsilon, -2.5);
}

/**
 * @brief coulomb weight.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double coulomb(double x) { return 1. / (std::abs(x) + epsilon); }

/**
 * @brief coulomb weight derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double dcoulomb(double x) { return -sign(x) / (x * x + epsilon); }

/**
 * @brief coulomb weight second derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double ddcoulomb(double x) { return 2. / (x * x * std::abs(x) + epsilon); }

/**
 * @brief step weight.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double step(double x) { return x <= 2 ? 1. : 0.; }

/**
 * @brief step weight derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double dstep(double x) { return 0.; }

/**
 * @brief step weight second derivative.
 *
 * @param x: given scalar.
 * @return its value at x.
 */
inline double ddstep(double x) { return 0.; }

/**
 * @brief Kernel types.
 */
enum Kernel { B4, SQUARE, EPANECHNIKOV, POLYRAMP, INVSQRT, COULOMB, STEP };

/**
 * @brief General kernel function.
 *
 * @param kern: kernel type.
 * @param x: given scalar.
 * @return related value.
 */
inline double kernel(const Kernel kern, double x) {
  switch (kern) {
    case B4:           return b4(x);
    case SQUARE:       return square(x);
    case EPANECHNIKOV: return epanechnikov(x);
    case POLYRAMP:     return polyramp(x);
    case INVSQRT:      return invsqrt(x);
    case COULOMB:      return coulomb(x);
    case STEP:         return step(x);
    default: return 0.;
  }
}

//\///////////////////////////////////////////////////////////////////////////
// general multi-dimensional evaluation function - what the public uses
//\///////////////////////////////////////////////////////////////////////////

/**
 * @brief Geometry types.
 */
enum Geometry { ELLIPTIC, TENSOR, FACETED };

/**
 * @brief generic elliptically symmetric weight function argument.
 *
 * @tparam dim: spatial dimension.
 * @param x: first point
 * @param y: second point
 * @param h: size metric
 * @return sqrt(sum_i (xi - yi)^2 / hi^2)
 */
template<int dim>
double elliptic(Wonton::Point<dim> const& x,
                Wonton::Point<dim> const& y,
                std::array<double,dim> const& h) {
  double distance = 0.;
  for (int i = 0; i < dim; i++) {
    distance += (x[i] - y[i]) * (x[i] - y[i]) / (h[i] * h[i]);
  }
  return std::sqrt(distance);
}

/**
 * @brief generic tensor weight function arguments.
 *
 * @tparam dim: spatial dimension.
 * @param x: first point
 * @param y: second point
 * @param h: size metric
 * @return [(xi - yi) / hi]_{i=0,dim-1}
 */
template<int dim>
std::array<double,dim> tensor(Wonton::Point<dim> const& x,
                              Wonton::Point<dim> const& y,
                              std::array<double,dim> const& h) {
  std::array<double,dim> result;
  for (int i = 0; i < dim; i++) {
    result[i] = (x[i] - y[i]) / h[i];
  }
  return result;
}

/**
 * @brief evaluation function for elliptic or tensor product weights.
 *
 * @tparam dim: spatial dimension.
 * @param geometry: the geometry to consider.
 * @param kern: the kernel to consider.
 * @param x: first point
 * @param y: second point
 * @param h: size metric
 * @return evaluated kernel value.
 */
template<int dim>
double eval(Geometry const geometry,
            Kernel const kern,
            Wonton::Point<dim> const& x,
            Wonton::Point<dim> const& y,
            std::array<double,dim> const& h) {
  double result = 0.;
  double const norm = kernel(kern, 0.0);
  switch (geometry) {
    case ELLIPTIC: {
      auto const arg = elliptic<dim>(x, y, h);
      result = kernel(kern, arg) / norm;
      break;
    }
    case TENSOR: {
      result = 1.;
      auto const arg = tensor<dim>(x, y, h);
      for (int i = 0; i < dim; i++)
        result *= kernel(kern, arg[i]) / norm;
      break;
    }
    default:
      throw std::runtime_error("invalid weight geometry");
  }
  return result;
}

/**
 * @brief data for specifying a faceted weight.
 *
 * @tparam dim: spatial dimension.
 */
template<int dim>
struct FacetData {
  double smoothing;
  std::array<double,dim> normal;
};

/**
 * @brief faceted weight function
 *
 * @tparam dim: spatial dimension.
 * @param kern: the kernel to use.
 * @param x: first point.
 * @param y: second point.
 * @param facets: list of facets.
 * @param nsides: number of sides.
 * @return evaluated kernel value.
 */
template<int dim>
double faceted(Kernel const kern,
               Wonton::Point<dim> const& x,
               Wonton::Point<dim> const& y,
               std::vector<FacetData<dim>> const& facets,
               size_t nsides) {
  assert(kern == POLYRAMP or kern == STEP);
  double result = 1.;
  for (size_t i = 0; i < nsides; i++) {
    double arg = 0.;
    for (size_t j = 0; j < dim; j++)
      arg += facets[i].normal[j] * (y[j] - x[j]);
    arg /= facets[i].smoothing;
    result *= kernel(kern, arg);
  }
  return result;
}

/**
 * @brief evaluation function for any weight.
 *
 * @tparam dim: spatial dimension.
 * @param geo: the geometry to consider.
 * @param kern: the kernel to consider.
 * @param x: first point.
 * @param y: second point.
 * @param vh: size matrix.
 * @return evaluated kernel value.
 */
template<int dim>
double eval(Geometry const geo,
            Kernel const kern,
            Wonton::Point<dim> const& x,
            Wonton::Point<dim> const& y,
            std::vector<std::vector<double>> const& vh) {
  switch (geo) {
    case TENSOR:
    case ELLIPTIC: {
      std::array<double,dim> h;
      for (size_t i = 0; i < dim; i++)
        h[i] = vh[0][i];
      return eval<dim>(geo, kern, x, y, h);
    }
    case FACETED: {
      int const nsides = vh.size();
      std::vector<FacetData<dim>> facets(nsides);
      for (int i = 0; i < nsides; i++) {
        for (int j = 0; j < dim; j++)
          facets[i].normal[j] = vh[i][j];
        facets[i].smoothing = vh[i][dim];
      }
      return faceted<dim>(kern, x, y, facets, nsides);
    }
    default:
      throw std::runtime_error("invalid weight geometry");
  }
}

}}}  // namespace Portage::Meshfree::Weight

#endif  // PORTAGE_SUPPORT_WEIGHT_H_
