/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef OPERATOR_H_INC_
#define OPERATOR_H_INC_

#include <cassert>
#include <cmath>
#include <array>
#include <vector>
#include <stdexcept>

#include "portage/support/basis.h"
#include "wonton/support/Point.h"

namespace Portage { namespace Meshfree { namespace oper {

using basis::Traits;
using basis::transfactor;

/**
 *
 */
enum Domain { Interval,
              Quadrilateral, Triangle, Circle,
              Hexahedron, Tetrahedron, Wedge, Sphere, LastDomain };

enum Type { VolumeIntegral, SurfaceIntegral, LastOperator };

template<Domain domain> class DomainTraits  { public: static const size_t dimension=0; };
template<> class DomainTraits<Interval>     { public: static const size_t dimension=1; };
template<> class DomainTraits<Quadrilateral>{ public: static const size_t dimension=2; };
template<> class DomainTraits<Triangle>     { public: static const size_t dimension=2; };
template<> class DomainTraits<Circle>       { public: static const size_t dimension=2; };
template<> class DomainTraits<Hexahedron>   { public: static const size_t dimension=3; };
template<> class DomainTraits<Tetrahedron>  { public: static const size_t dimension=3; };
template<> class DomainTraits<Wedge>        { public: static const size_t dimension=3; };
template<> class DomainTraits<Sphere>       { public: static const size_t dimension=3; };

/**
 *
 * @param domain
 * @return
 */
constexpr size_t dimension(Domain domain) {
  switch (domain) {
    case Interval:      return DomainTraits<Interval>::dimension;
    case Quadrilateral: return DomainTraits<Quadrilateral>::dimension;
    case Triangle:      return DomainTraits<Triangle>::dimension;
    case Circle:        return DomainTraits<Circle>::dimension;
    case Hexahedron:    return DomainTraits<Hexahedron>::dimension;
    case Tetrahedron:   return DomainTraits<Tetrahedron>::dimension;
    case Wedge:         return DomainTraits<Wedge>::dimension;
    case Sphere:        return DomainTraits<Sphere>::dimension;
    default: return 0;
  }
}

/**
 *
 * @tparam dim
 * @param points
 * @return
 */
template<int dim>
Domain domain_from_points(std::vector<Wonton::Point<dim>> const& points) {

  int const nb_points = points.size();

  switch(dim) {
    case 1: return Interval;
    case 2:
      switch(nb_points) {
        case 3: return Triangle;
        case 4: return Quadrilateral;
        default: throw std::runtime_error("invalid number of points");
      }
    case 3:
      switch(nb_points) {
        case 4: return Tetrahedron;
        case 6: return Wedge;
        case 8: return Hexahedron;
        default: throw std::runtime_error("invalid number of points");
      }
    default: throw std::runtime_error("invalid dimension");
  }
}

////////////////////////////////////////////////////////////////////////////////
// Template for integral operator base class
////////////////////////////////////////////////////////////////////////////////

template<Type type, basis::Type basis_type, Domain domain_type>
class OperatorBase {
public:
  Type operator_type = type;
  static constexpr basis::Type basis = basis_type;
  static constexpr Domain domain = domain_type;
};

////////////////////////////////////////////////////////////////////////////////
// Template for integral operators
////////////////////////////////////////////////////////////////////////////////

template<Type type, basis::Type basis_type, Domain domain_type>
class Operator: public OperatorBase<type,basis_type,domain_type> {
public:
  static constexpr size_t dim = dimension(domain_type);
  static constexpr size_t operator_size = 0;
  static constexpr size_t basis_size = 0;
  static constexpr size_t point_size = 0;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    throw std::runtime_error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
// Specializations
////////////////////////////////////////////////////////////////////////////////

// 1D
template<>
class Operator<VolumeIntegral, basis::Unitary, Interval>:
  public OperatorBase<VolumeIntegral, basis::Unitary, Interval>
{
public:
  static constexpr size_t dim = dimension(Interval);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Unitary, dim>::function_size;
  static constexpr size_t point_size = 2;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = p[1][0] - p[0][0];
    return result;
  }
};


template<>
class Operator<VolumeIntegral, basis::Linear, Interval>:
  public OperatorBase<VolumeIntegral, basis::Linear, Interval>
{
public:
  static constexpr size_t dim = dimension(Interval);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Linear, dim>::function_size;
  static constexpr size_t point_size = 2;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = p[1][0] - p[0][0];
    result[1][0] = .5*(pow(p[1][0],2.) - pow(p[0][0],2.));
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Quadratic, Interval>:
  public OperatorBase<VolumeIntegral, basis::Quadratic, Interval>
{
public:
  static constexpr size_t dim = dimension(Interval);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Quadratic, dim>::function_size;
  static constexpr size_t point_size = 2;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = p[1][0] - p[0][0];
    result[1][0] = .5*(pow(p[1][0],2.) - pow(p[0][0],2.));
    result[2][0] = (pow(p[1][0],3.) - pow(p[0][0],3.))/6.;
    return result;
  }
};

// 2D Triangle

template<>
class Operator<VolumeIntegral, basis::Unitary, Triangle>:
  public OperatorBase<VolumeIntegral, basis::Unitary, Interval>
{
public:
  static constexpr size_t dim = dimension(Triangle);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Unitary, dim>::function_size;
  static constexpr size_t point_size = 3;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = .5*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]);
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Linear, Triangle>:
  public OperatorBase<VolumeIntegral, basis::Linear, Interval>
{
public:
  static constexpr size_t dim = dimension(Triangle);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Linear, dim>::function_size;
  static constexpr size_t point_size = 3;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = .5*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]);
    result[1][0] = ((p[0][0] + p[1][0] + p[2][0])*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1])+ p[1][0]*p[2][1]))/6.;
    result[2][0] = -((p[0][1] + p[1][1] + p[2][1])*(p[0][1]*(p[1][0] - p[2][0]) + p[1][1]*p[2][0] - p[1][0]*p[2][1] + p[0][0]*(-p[1][1] + p[2][1])))/6.;
    return result;
  }
};

template<>
  class Operator<VolumeIntegral, basis::Quadratic, Triangle>:
    public OperatorBase<VolumeIntegral, basis::Quadratic, Interval>
  {
  public:
    static constexpr size_t dim = dimension(Triangle);
    static constexpr size_t operator_size = 1;
    static constexpr size_t basis_size =
      basis::Traits<basis::Quadratic, dim>::function_size;
    static constexpr size_t point_size = 3;

    using result_t = std::array<std::array<double, operator_size>, basis_size>;
    using points_t = std::array<Wonton::Point<dim>, point_size>;

    static result_t apply(const points_t p) {
      result_t result;
      result[0][0] = .5 * (-(p[1][1] * p[2][0]) + p[0][1] * (-p[1][0] + p[2][0]) + p[0][0] * (p[1][1] - p[2][1]) + p[1][0] * p[2][1]);
      result[1][0] = ((p[0][0] + p[1][0] + p[2][0]) * (-(p[1][1] * p[2][0]) + p[0][1] * (-p[1][0] + p[2][0]) + p[0][0] * (p[1][1] - p[2][1]) + p[1][0] * p[2][1])) / 6.;
      result[2][0] = -((p[0][1] + p[1][1] + p[2][1]) * (p[0][1] * (p[1][0] - p[2][0]) + p[1][1] * p[2][0] - p[1][0] * p[2][1] + p[0][0] * (-p[1][1] + p[2][1]))) / 6.;
      result[3][0] = ((pow(p[0][0], 2) + pow(p[1][0], 2) + pow(p[2][0], 2) + p[1][0] * p[2][0] + p[0][0] * (p[1][0] + p[2][0])) * (-(p[1][1] * p[2][0]) + p[0][1] * (-p[1][0] + p[2][0]) + p[0][0] * (p[1][1] - p[2][1]) + p[1][0] * p[2][1])) / 24.;
      result[4][0] = (-(pow(p[1][1], 2) * pow(p[2][0], 2)) + pow(p[0][1], 2) * (-pow(p[1][0], 2) + pow(p[2][0], 2)) + pow(p[1][0], 2) * pow(p[2][1], 2) - 2 * pow(p[1][1], 2) * p[1][0] * p[2][0] + 2 * pow(p[2][1], 2) * p[1][0] * p[2][0] - 2 * p[0][0] * (-(pow(p[1][1], 2) * p[1][0]) + pow(p[0][1], 2) * (p[1][0] - p[2][0]) + pow(p[2][1], 2) * p[2][0]) + 2 * pow(p[1][0], 2) * p[1][1] * p[2][1] - 2 * pow(p[2][0], 2) * p[1][1] * p[2][1] + pow(p[0][0], 2) * (p[1][1] - p[2][1]) * (2 * p[0][1] + p[1][1] + p[2][1]) + p[0][1] * (-2 * pow(p[1][0], 2) * p[1][1] + 2 * pow(p[2][0], 2) * p[2][1])) / 24.;
      result[5][0] = -((p[0][1] * (p[1][0] - p[2][0]) + p[1][1] * p[2][0] - p[1][0] * p[2][1] + p[0][0] * (-p[1][1] + p[2][1])) * (pow(p[0][1], 2) + pow(p[1][1], 2) + pow(p[2][1], 2) + p[1][1] * p[2][1] + p[0][1] * (p[1][1] + p[2][1]))) / 24.;
      return result;
    }
  };

// 2D Quadrilateral

template<>
class Operator<VolumeIntegral, basis::Unitary, Quadrilateral>:
  public OperatorBase<VolumeIntegral, basis::Unitary, Interval>
{
public:
  static constexpr size_t dim = dimension(Quadrilateral);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Unitary, dim>::function_size;
  static constexpr size_t point_size = 4;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (-(p[1][1] * p[2][0]) + p[1][0] * p[2][1] - p[2][1] * p[3][0] + p[0][1] * (-p[1][0] + p[3][0]) + p[0][0] * (p[1][1] - p[3][1]) + p[2][0] * p[3][1]) / 2;
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Linear, Quadrilateral>:
  public OperatorBase<VolumeIntegral, basis::Linear, Interval>
{
public:
  static constexpr size_t dim = dimension(Quadrilateral);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Linear, dim>::function_size;
  static constexpr size_t point_size = 4;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (-(p[1][1] * p[2][0]) + p[1][0] * p[2][1] - p[2][1] * p[3][0] + p[0][1] * (-p[1][0] + p[3][0]) + p[0][0] * (p[1][1] - p[3][1]) + p[2][0] * p[3][1]) / 2;
    result[1][0] = ((-pow(p[1][0], 2) + pow(p[3][0], 2)) * p[0][1] - pow(p[2][0], 2) * p[1][1] - p[1][0] * p[1][1] * p[2][0] + pow(p[1][0], 2) * p[2][1] - pow(p[3][0], 2) * p[2][1] + p[1][0] * p[2][0] * p[2][1] - p[2][0] * p[2][1] * p[3][0] + pow(p[0][0], 2) * (p[1][1] - p[3][1]) + pow(p[2][0], 2) * p[3][1] + p[2][0] * p[3][0] * p[3][1] + p[0][0] * (p[1][0] * p[1][1] + p[0][1] * (-p[1][0] + p[3][0]) - p[3][0] * p[3][1])) / 6;
    result[2][0] = ((pow(p[1][1], 2) - pow(p[3][1], 2)) * p[0][0] + pow(p[2][1], 2) * p[1][0] - pow(p[1][1], 2) * p[2][0] + pow(p[3][1], 2) * p[2][0] + p[1][0] * p[1][1] * p[2][1] - p[1][1] * p[2][0] * p[2][1] - pow(p[2][1], 2) * p[3][0] + pow(p[0][1], 2) * (-p[1][0] + p[3][0]) + p[2][0] * p[2][1] * p[3][1] - p[2][1] * p[3][0] * p[3][1] + p[0][1] * (-(p[1][0] * p[1][1]) + p[0][0] * (p[1][1] - p[3][1]) + p[3][0] * p[3][1])) / 6;
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Quadratic, Quadrilateral>:
  public OperatorBase<VolumeIntegral, basis::Quadratic, Interval>
{
public:
  static constexpr size_t dim = dimension(Quadrilateral);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Quadratic, dim>::function_size;
  static constexpr size_t point_size = 4;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (-(p[1][1]*p[2][0]) + p[1][0]*p[2][1] - p[2][1]*p[3][0] + p[0][1]*(-p[1][0] + p[3][0]) + p[0][0]*(p[1][1] - p[3][1]) + p[2][0]*p[3][1])/2;
    result[1][0] = ((-pow(p[1][0], 2) + pow(p[3][0], 2))*p[0][1] - pow(p[2][0], 2)*p[1][1] - p[1][0]*p[1][1]*p[2][0] + pow(p[1][0], 2)*p[2][1] - pow(p[3][0], 2)*p[2][1] + p[1][0]*p[2][0]*p[2][1] - p[2][0]*p[2][1]*p[3][0] + pow(p[0][0], 2)*(p[1][1] - p[3][1]) + pow(p[2][0], 2)*p[3][1] + p[2][0]*p[3][0]*p[3][1] + p[0][0]*(p[1][0]*p[1][1] + p[0][1]*(-p[1][0] + p[3][0]) - p[3][0]*p[3][1]))/6;
    result[2][0] = ((pow(p[1][1], 2) - pow(p[3][1], 2))*p[0][0] + pow(p[2][1], 2)*p[1][0] - pow(p[1][1], 2)*p[2][0] + pow(p[3][1], 2)*p[2][0] + p[1][0]*p[1][1]*p[2][1] - p[1][1]*p[2][0]*p[2][1] - pow(p[2][1], 2)*p[3][0] + pow(p[0][1], 2)*(-p[1][0] + p[3][0]) + p[2][0]*p[2][1]*p[3][1] - p[2][1]*p[3][0]*p[3][1] + p[0][1]*(-(p[1][0]*p[1][1]) + p[0][0]*(p[1][1] - p[3][1]) + p[3][0]*p[3][1]))/6;
    result[3][0] = ((-pow(p[1][0], 3) + pow(p[3][0], 3))*p[0][1] - pow(p[2][0], 3)*p[1][1] - pow(p[2][0], 2)*p[1][0]*p[1][1] - pow(p[1][0], 2)*p[1][1]*p[2][0] + pow(p[1][0], 3)*p[2][1] - pow(p[3][0], 3)*p[2][1] + pow(p[2][0], 2)*p[1][0]*p[2][1] + pow(p[1][0], 2)*p[2][0]*p[2][1] - pow(p[3][0], 2)*p[2][0]*p[2][1] - pow(p[2][0], 2)*p[2][1]*p[3][0] + pow(p[0][0], 3)*(p[1][1] - p[3][1]) + pow(p[2][0], 3)*p[3][1] + pow(p[3][0], 2)*p[2][0]*p[3][1] + pow(p[2][0], 2)*p[3][0]*p[3][1] + p[0][0]*((-pow(p[1][0], 2) + pow(p[3][0], 2))*p[0][1] + pow(p[1][0], 2)*p[1][1] - pow(p[3][0], 2)*p[3][1]) + pow(p[0][0], 2)*(p[1][0]*p[1][1] + p[0][1]*(-p[1][0] + p[3][0]) - p[3][0]*p[3][1]))/24;
    result[4][0] = (-(pow(p[1][1], 2)*pow(p[2][0], 2)) + pow(p[1][0], 2)*pow(p[2][1], 2) - pow(p[2][1], 2)*pow(p[3][0], 2) + pow(p[0][1], 2)*(-pow(p[1][0], 2) + pow(p[3][0], 2)) + pow(p[2][0], 2)*pow(p[3][1], 2) - 2*pow(p[1][1], 2)*p[1][0]*p[2][0] + 2*pow(p[2][1], 2)*p[1][0]*p[2][0] + 2*pow(p[1][0], 2)*p[1][1]*p[2][1] - 2*pow(p[2][0], 2)*p[1][1]*p[2][1] - 2*pow(p[2][1], 2)*p[2][0]*p[3][0] + 2*pow(p[3][1], 2)*p[2][0]*p[3][0] - 2*p[0][0]*(-(pow(p[1][1], 2)*p[1][0]) + pow(p[0][1], 2)*(p[1][0] - p[3][0]) + pow(p[3][1], 2)*p[3][0]) + 2*pow(p[2][0], 2)*p[2][1]*p[3][1] - 2*pow(p[3][0], 2)*p[2][1]*p[3][1] + pow(p[0][0], 2)*(p[1][1] - p[3][1])*(2*p[0][1] + p[1][1] + p[3][1]) + p[0][1]*(-2*pow(p[1][0], 2)*p[1][1] + 2*pow(p[3][0], 2)*p[3][1]))/24;
    result[5][0] = ((pow(p[1][1], 3) - pow(p[3][1], 3))*p[0][0] + pow(p[2][1], 3)*p[1][0] + pow(p[2][1], 2)*p[1][0]*p[1][1] - pow(p[1][1], 3)*p[2][0] + pow(p[3][1], 3)*p[2][0] - pow(p[2][1], 2)*p[1][1]*p[2][0] + pow(p[1][1], 2)*p[1][0]*p[2][1] - pow(p[1][1], 2)*p[2][0]*p[2][1] + pow(p[3][1], 2)*p[2][0]*p[2][1] - pow(p[2][1], 3)*p[3][0] - pow(p[3][1], 2)*p[2][1]*p[3][0] + pow(p[0][1], 3)*(-p[1][0] + p[3][0]) + p[0][1]*((pow(p[1][1], 2) - pow(p[3][1], 2))*p[0][0] - pow(p[1][1], 2)*p[1][0] + pow(p[3][1], 2)*p[3][0]) + pow(p[2][1], 2)*p[2][0]*p[3][1] - pow(p[2][1], 2)*p[3][0]*p[3][1] + pow(p[0][1], 2)*(-(p[1][0]*p[1][1]) + p[0][0]*(p[1][1] - p[3][1]) + p[3][0]*p[3][1]))/24;
    return result;
  }
};

// 3D Hexahedron

template<>
class Operator<VolumeIntegral, basis::Unitary, Hexahedron>:
  public OperatorBase<VolumeIntegral, basis::Unitary, Interval>
{
public:
  static constexpr size_t dim = dimension(Hexahedron);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Unitary, dim>::function_size;
  static constexpr size_t point_size = 8;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] - p[0][0]*p[1][1]*p[3][2] + p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] - p[0][0]*p[1][2]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + p[1][2]*p[2][0]*p[5][1] - p[1][0]*p[2][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[0][0]*p[1][1]*p[5][2] - p[1][1]*p[2][0]*p[5][2] + p[1][0]*p[2][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] - p[0][0]*p[4][1]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[1][2]*p[2][1]*p[6][0] + p[1][1]*p[2][2]*p[6][0] - p[2][2]*p[3][1]*p[6][0] + p[2][1]*p[3][2]*p[6][0] + p[1][2]*p[5][1]*p[6][0] + p[2][2]*p[5][1]*p[6][0] - p[4][2]*p[5][1]*p[6][0] - p[1][1]*p[5][2]*p[6][0] - p[2][1]*p[5][2]*p[6][0] + p[4][1]*p[5][2]*p[6][0] + p[1][2]*p[2][0]*p[6][1] - p[1][0]*p[2][2]*p[6][1] + p[2][2]*p[3][0]*p[6][1] - p[2][0]*p[3][2]*p[6][1] - p[1][2]*p[5][0]*p[6][1] - p[2][2]*p[5][0]*p[6][1] + p[4][2]*p[5][0]*p[6][1] + p[1][0]*p[5][2]*p[6][1] + p[2][0]*p[5][2]*p[6][1] - p[4][0]*p[5][2]*p[6][1] - p[1][1]*p[2][0]*p[6][2] + p[1][0]*p[2][1]*p[6][2] - p[2][1]*p[3][0]*p[6][2] + p[2][0]*p[3][1]*p[6][2] + p[1][1]*p[5][0]*p[6][2] + p[2][1]*p[5][0]*p[6][2] - p[4][1]*p[5][0]*p[6][2] - p[1][0]*p[5][1]*p[6][2] - p[2][0]*p[5][1]*p[6][2] + p[4][0]*p[5][1]*p[6][2] - p[2][2]*p[3][1]*p[7][0] + p[2][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] - p[4][2]*p[5][1]*p[7][0] + p[4][1]*p[5][2]*p[7][0] + p[2][2]*p[6][1]*p[7][0] + p[3][2]*p[6][1]*p[7][0] - p[4][2]*p[6][1]*p[7][0] - p[5][2]*p[6][1]*p[7][0] - p[2][1]*p[6][2]*p[7][0] - p[3][1]*p[6][2]*p[7][0] + p[4][1]*p[6][2]*p[7][0] + p[5][1]*p[6][2]*p[7][0] + p[2][2]*p[3][0]*p[7][1] + p[0][0]*p[3][2]*p[7][1] - p[2][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][0]*p[4][2]*p[7][1] + p[4][2]*p[5][0]*p[7][1] - p[4][0]*p[5][2]*p[7][1] - p[2][2]*p[6][0]*p[7][1] - p[3][2]*p[6][0]*p[7][1] + p[4][2]*p[6][0]*p[7][1] + p[5][2]*p[6][0]*p[7][1] + p[2][0]*p[6][2]*p[7][1] + p[3][0]*p[6][2]*p[7][1] - p[4][0]*p[6][2]*p[7][1] - p[5][0]*p[6][2]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][1]*(p[2][0] + p[3][0] - p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - p[3][1] + p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[2][1]*p[3][0]*p[7][2] - p[0][0]*p[3][1]*p[7][2] + p[2][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][0]*p[4][1]*p[7][2] - p[4][1]*p[5][0]*p[7][2] + p[4][0]*p[5][1]*p[7][2] + p[2][1]*p[6][0]*p[7][2] + p[3][1]*p[6][0]*p[7][2] - p[4][1]*p[6][0]*p[7][2] - p[5][1]*p[6][0]*p[7][2] - p[2][0]*p[6][1]*p[7][2] - p[3][0]*p[6][1]*p[7][2] + p[4][0]*p[6][1]*p[7][2] + p[5][0]*p[6][1]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - p[3][0] + p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + p[3][2] - p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2]))/12;
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Linear, Hexahedron>:
  public OperatorBase<VolumeIntegral, basis::Linear, Interval>
{
public:
  static constexpr size_t dim = dimension(Hexahedron);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Linear, dim>::function_size;
  static constexpr size_t point_size = 8;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] - p[0][0]*p[1][1]*p[3][2] + p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] - p[0][0]*p[1][2]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + p[1][2]*p[2][0]*p[5][1] - p[1][0]*p[2][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[0][0]*p[1][1]*p[5][2] - p[1][1]*p[2][0]*p[5][2] + p[1][0]*p[2][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] - p[0][0]*p[4][1]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[1][2]*p[2][1]*p[6][0] + p[1][1]*p[2][2]*p[6][0] - p[2][2]*p[3][1]*p[6][0] + p[2][1]*p[3][2]*p[6][0] + p[1][2]*p[5][1]*p[6][0] + p[2][2]*p[5][1]*p[6][0] - p[4][2]*p[5][1]*p[6][0] - p[1][1]*p[5][2]*p[6][0] - p[2][1]*p[5][2]*p[6][0] + p[4][1]*p[5][2]*p[6][0] + p[1][2]*p[2][0]*p[6][1] - p[1][0]*p[2][2]*p[6][1] + p[2][2]*p[3][0]*p[6][1] - p[2][0]*p[3][2]*p[6][1] - p[1][2]*p[5][0]*p[6][1] - p[2][2]*p[5][0]*p[6][1] + p[4][2]*p[5][0]*p[6][1] + p[1][0]*p[5][2]*p[6][1] + p[2][0]*p[5][2]*p[6][1] - p[4][0]*p[5][2]*p[6][1] - p[1][1]*p[2][0]*p[6][2] + p[1][0]*p[2][1]*p[6][2] - p[2][1]*p[3][0]*p[6][2] + p[2][0]*p[3][1]*p[6][2] + p[1][1]*p[5][0]*p[6][2] + p[2][1]*p[5][0]*p[6][2] - p[4][1]*p[5][0]*p[6][2] - p[1][0]*p[5][1]*p[6][2] - p[2][0]*p[5][1]*p[6][2] + p[4][0]*p[5][1]*p[6][2] - p[2][2]*p[3][1]*p[7][0] + p[2][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] - p[4][2]*p[5][1]*p[7][0] + p[4][1]*p[5][2]*p[7][0] + p[2][2]*p[6][1]*p[7][0] + p[3][2]*p[6][1]*p[7][0] - p[4][2]*p[6][1]*p[7][0] - p[5][2]*p[6][1]*p[7][0] - p[2][1]*p[6][2]*p[7][0] - p[3][1]*p[6][2]*p[7][0] + p[4][1]*p[6][2]*p[7][0] + p[5][1]*p[6][2]*p[7][0] + p[2][2]*p[3][0]*p[7][1] + p[0][0]*p[3][2]*p[7][1] - p[2][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][0]*p[4][2]*p[7][1] + p[4][2]*p[5][0]*p[7][1] - p[4][0]*p[5][2]*p[7][1] - p[2][2]*p[6][0]*p[7][1] - p[3][2]*p[6][0]*p[7][1] + p[4][2]*p[6][0]*p[7][1] + p[5][2]*p[6][0]*p[7][1] + p[2][0]*p[6][2]*p[7][1] + p[3][0]*p[6][2]*p[7][1] - p[4][0]*p[6][2]*p[7][1] - p[5][0]*p[6][2]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][1]*(p[2][0] + p[3][0] - p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - p[3][1] + p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[2][1]*p[3][0]*p[7][2] - p[0][0]*p[3][1]*p[7][2] + p[2][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][0]*p[4][1]*p[7][2] - p[4][1]*p[5][0]*p[7][2] + p[4][0]*p[5][1]*p[7][2] + p[2][1]*p[6][0]*p[7][2] + p[3][1]*p[6][0]*p[7][2] - p[4][1]*p[6][0]*p[7][2] - p[5][1]*p[6][0]*p[7][2] - p[2][0]*p[6][1]*p[7][2] - p[3][0]*p[6][1]*p[7][2] + p[4][0]*p[6][1]*p[7][2] + p[5][0]*p[6][1]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - p[3][0] + p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + p[3][2] - p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2]))/12;
    result[1][0] = (-(pow(p[2][0], 2)*p[0][1]*p[1][2]) - pow(p[3][0], 2)*p[0][1]*p[1][2] + pow(p[4][0], 2)*p[0][1]*p[1][2] + pow(p[5][0], 2)*p[0][1]*p[1][2] - 2*p[0][1]*p[1][0]*p[1][2]*p[2][0] + pow(p[3][0], 2)*p[1][2]*p[2][1] - pow(p[5][0], 2)*p[1][2]*p[2][1] - pow(p[6][0], 2)*p[1][2]*p[2][1] + 2*pow(p[1][0], 2)*p[0][1]*p[2][2] - 2*pow(p[3][0], 2)*p[0][1]*p[2][2] - pow(p[3][0], 2)*p[1][1]*p[2][2] + pow(p[5][0], 2)*p[1][1]*p[2][2] + pow(p[6][0], 2)*p[1][1]*p[2][2] + p[0][1]*p[1][0]*p[2][0]*p[2][2] - p[0][1]*p[1][0]*p[1][2]*p[3][0] - p[0][1]*p[1][2]*p[2][0]*p[3][0] + p[1][0]*p[1][2]*p[2][1]*p[3][0] + 2*p[1][2]*p[2][0]*p[2][1]*p[3][0] - p[1][0]*p[1][1]*p[2][2]*p[3][0] - p[0][1]*p[2][0]*p[2][2]*p[3][0] - 2*p[1][1]*p[2][0]*p[2][2]*p[3][0] - 2*pow(p[2][0], 2)*p[1][2]*p[3][1] - p[1][0]*p[1][2]*p[2][0]*p[3][1] + pow(p[1][0], 2)*p[2][2]*p[3][1] - pow(p[6][0], 2)*p[2][2]*p[3][1] - pow(p[7][0], 2)*p[2][2]*p[3][1] + 2*p[1][0]*p[2][0]*p[2][2]*p[3][1] - p[1][2]*p[2][0]*p[3][0]*p[3][1] + p[1][0]*p[2][2]*p[3][0]*p[3][1] + pow(p[1][0], 2)*p[0][1]*p[3][2] + pow(p[2][0], 2)*p[0][1]*p[3][2] - pow(p[4][0], 2)*p[0][1]*p[3][2] - pow(p[7][0], 2)*p[0][1]*p[3][2] + 2*pow(p[2][0], 2)*p[1][1]*p[3][2] + p[0][1]*p[1][0]*p[2][0]*p[3][2] + p[1][0]*p[1][1]*p[2][0]*p[3][2] - pow(p[1][0], 2)*p[2][1]*p[3][2] + pow(p[6][0], 2)*p[2][1]*p[3][2] + pow(p[7][0], 2)*p[2][1]*p[3][2] - 2*p[1][0]*p[2][0]*p[2][1]*p[3][2] + p[0][1]*p[1][0]*p[3][0]*p[3][2] + 2*p[0][1]*p[2][0]*p[3][0]*p[3][2] + p[1][1]*p[2][0]*p[3][0]*p[3][2] - p[1][0]*p[2][1]*p[3][0]*p[3][2] + p[0][1]*p[1][0]*p[1][2]*p[4][0] - p[0][1]*p[3][0]*p[3][2]*p[4][0] + 2*pow(p[5][0], 2)*p[1][2]*p[4][1] - 2*pow(p[7][0], 2)*p[3][2]*p[4][1] - pow(p[1][0], 2)*p[0][1]*p[4][2] + pow(p[3][0], 2)*p[0][1]*p[4][2] - pow(p[5][0], 2)*p[0][1]*p[4][2] + pow(p[7][0], 2)*p[0][1]*p[4][2] - 2*pow(p[5][0], 2)*p[1][1]*p[4][2] + 2*pow(p[7][0], 2)*p[3][1]*p[4][2] - p[0][1]*p[1][0]*p[4][0]*p[4][2] + p[0][1]*p[3][0]*p[4][0]*p[4][2] + 2*p[0][1]*p[1][0]*p[1][2]*p[5][0] - 2*p[1][0]*p[1][2]*p[2][1]*p[5][0] - p[1][2]*p[2][0]*p[2][1]*p[5][0] + 2*p[1][0]*p[1][1]*p[2][2]*p[5][0] + p[1][1]*p[2][0]*p[2][2]*p[5][0] + p[0][1]*p[1][2]*p[4][0]*p[5][0] + p[1][0]*p[1][2]*p[4][1]*p[5][0] + p[1][2]*p[4][0]*p[4][1]*p[5][0] - p[0][1]*p[1][0]*p[4][2]*p[5][0] - p[1][0]*p[1][1]*p[4][2]*p[5][0] - 2*p[0][1]*p[4][0]*p[4][2]*p[5][0] - p[1][1]*p[4][0]*p[4][2]*p[5][0] + pow(p[2][0], 2)*p[1][2]*p[5][1] - pow(p[4][0], 2)*p[1][2]*p[5][1] + pow(p[6][0], 2)*p[1][2]*p[5][1] + 2*p[1][0]*p[1][2]*p[2][0]*p[5][1] - 2*pow(p[1][0], 2)*p[2][2]*p[5][1] + 2*pow(p[6][0], 2)*p[2][2]*p[5][1] - p[1][0]*p[2][0]*p[2][2]*p[5][1] - p[1][0]*p[1][2]*p[4][0]*p[5][1] + pow(p[1][0], 2)*p[4][2]*p[5][1] - pow(p[6][0], 2)*p[4][2]*p[5][1] - pow(p[7][0], 2)*p[4][2]*p[5][1] + p[1][0]*p[4][0]*p[4][2]*p[5][1] + p[1][2]*p[2][0]*p[5][0]*p[5][1] - p[1][0]*p[2][2]*p[5][0]*p[5][1] - 2*p[1][2]*p[4][0]*p[5][0]*p[5][1] + 2*p[1][0]*p[4][2]*p[5][0]*p[5][1] - 2*pow(p[1][0], 2)*p[0][1]*p[5][2] + 2*pow(p[4][0], 2)*p[0][1]*p[5][2] - pow(p[2][0], 2)*p[1][1]*p[5][2] + pow(p[4][0], 2)*p[1][1]*p[5][2] - pow(p[6][0], 2)*p[1][1]*p[5][2] - 2*p[1][0]*p[1][1]*p[2][0]*p[5][2] + 2*pow(p[1][0], 2)*p[2][1]*p[5][2] - 2*pow(p[6][0], 2)*p[2][1]*p[5][2] + p[1][0]*p[2][0]*p[2][1]*p[5][2] + p[1][0]*p[1][1]*p[4][0]*p[5][2] - pow(p[1][0], 2)*p[4][1]*p[5][2] + pow(p[6][0], 2)*p[4][1]*p[5][2] + pow(p[7][0], 2)*p[4][1]*p[5][2] - p[1][0]*p[4][0]*p[4][1]*p[5][2] - p[0][1]*p[1][0]*p[5][0]*p[5][2] - p[1][1]*p[2][0]*p[5][0]*p[5][2] + p[1][0]*p[2][1]*p[5][0]*p[5][2] + p[0][1]*p[4][0]*p[5][0]*p[5][2] + 2*p[1][1]*p[4][0]*p[5][0]*p[5][2] - 2*p[1][0]*p[4][1]*p[5][0]*p[5][2] - p[1][0]*p[1][2]*p[2][1]*p[6][0] - 2*p[1][2]*p[2][0]*p[2][1]*p[6][0] + p[1][0]*p[1][1]*p[2][2]*p[6][0] + 2*p[1][1]*p[2][0]*p[2][2]*p[6][0] - 2*p[2][0]*p[2][2]*p[3][1]*p[6][0] - p[2][2]*p[3][0]*p[3][1]*p[6][0] + 2*p[2][0]*p[2][1]*p[3][2]*p[6][0] + p[2][1]*p[3][0]*p[3][2]*p[6][0] - p[1][2]*p[2][1]*p[5][0]*p[6][0] + p[1][1]*p[2][2]*p[5][0]*p[6][0] + p[1][0]*p[1][2]*p[5][1]*p[6][0] + p[1][2]*p[2][0]*p[5][1]*p[6][0] + p[2][0]*p[2][2]*p[5][1]*p[6][0] - p[4][0]*p[4][2]*p[5][1]*p[6][0] + 2*p[1][2]*p[5][0]*p[5][1]*p[6][0] + p[2][2]*p[5][0]*p[5][1]*p[6][0] - 2*p[4][2]*p[5][0]*p[5][1]*p[6][0] - p[1][0]*p[1][1]*p[5][2]*p[6][0] - p[1][1]*p[2][0]*p[5][2]*p[6][0] - p[2][0]*p[2][1]*p[5][2]*p[6][0] + p[4][0]*p[4][1]*p[5][2]*p[6][0] - 2*p[1][1]*p[5][0]*p[5][2]*p[6][0] - p[2][1]*p[5][0]*p[5][2]*p[6][0] + 2*p[4][1]*p[5][0]*p[5][2]*p[6][0] + 2*pow(p[2][0], 2)*p[1][2]*p[6][1] - 2*pow(p[5][0], 2)*p[1][2]*p[6][1] + p[1][0]*p[1][2]*p[2][0]*p[6][1] - pow(p[1][0], 2)*p[2][2]*p[6][1] + pow(p[3][0], 2)*p[2][2]*p[6][1] - pow(p[5][0], 2)*p[2][2]*p[6][1] + pow(p[7][0], 2)*p[2][2]*p[6][1] - 2*p[1][0]*p[2][0]*p[2][2]*p[6][1] + 2*p[2][0]*p[2][2]*p[3][0]*p[6][1] - 2*pow(p[2][0], 2)*p[3][2]*p[6][1] + 2*pow(p[7][0], 2)*p[3][2]*p[6][1] - p[2][0]*p[3][0]*p[3][2]*p[6][1] + 2*pow(p[5][0], 2)*p[4][2]*p[6][1] - 2*pow(p[7][0], 2)*p[4][2]*p[6][1] - p[1][0]*p[1][2]*p[5][0]*p[6][1] - p[1][0]*p[2][2]*p[5][0]*p[6][1] - p[2][0]*p[2][2]*p[5][0]*p[6][1] + p[4][0]*p[4][2]*p[5][0]*p[6][1] + pow(p[1][0], 2)*p[5][2]*p[6][1] + pow(p[2][0], 2)*p[5][2]*p[6][1] - pow(p[4][0], 2)*p[5][2]*p[6][1] - pow(p[7][0], 2)*p[5][2]*p[6][1] + p[1][0]*p[2][0]*p[5][2]*p[6][1] + 2*p[1][0]*p[5][0]*p[5][2]*p[6][1] + p[2][0]*p[5][0]*p[5][2]*p[6][1] - 2*p[4][0]*p[5][0]*p[5][2]*p[6][1] + p[1][2]*p[2][0]*p[6][0]*p[6][1] - p[1][0]*p[2][2]*p[6][0]*p[6][1] + p[2][2]*p[3][0]*p[6][0]*p[6][1] - p[2][0]*p[3][2]*p[6][0]*p[6][1] - p[1][2]*p[5][0]*p[6][0]*p[6][1] - 2*p[2][2]*p[5][0]*p[6][0]*p[6][1] + p[4][2]*p[5][0]*p[6][0]*p[6][1] + p[1][0]*p[5][2]*p[6][0]*p[6][1] + 2*p[2][0]*p[5][2]*p[6][0]*p[6][1] - p[4][0]*p[5][2]*p[6][0]*p[6][1] - 2*pow(p[2][0], 2)*p[1][1]*p[6][2] + 2*pow(p[5][0], 2)*p[1][1]*p[6][2] - p[1][0]*p[1][1]*p[2][0]*p[6][2] + pow(p[1][0], 2)*p[2][1]*p[6][2] - pow(p[3][0], 2)*p[2][1]*p[6][2] + pow(p[5][0], 2)*p[2][1]*p[6][2] - pow(p[7][0], 2)*p[2][1]*p[6][2] + 2*p[1][0]*p[2][0]*p[2][1]*p[6][2] - 2*p[2][0]*p[2][1]*p[3][0]*p[6][2] + 2*pow(p[2][0], 2)*p[3][1]*p[6][2] - 2*pow(p[7][0], 2)*p[3][1]*p[6][2] + p[2][0]*p[3][0]*p[3][1]*p[6][2] - 2*pow(p[5][0], 2)*p[4][1]*p[6][2] + 2*pow(p[7][0], 2)*p[4][1]*p[6][2] + p[1][0]*p[1][1]*p[5][0]*p[6][2] + p[1][0]*p[2][1]*p[5][0]*p[6][2] + p[2][0]*p[2][1]*p[5][0]*p[6][2] - p[4][0]*p[4][1]*p[5][0]*p[6][2] - pow(p[1][0], 2)*p[5][1]*p[6][2] - pow(p[2][0], 2)*p[5][1]*p[6][2] + pow(p[4][0], 2)*p[5][1]*p[6][2] + pow(p[7][0], 2)*p[5][1]*p[6][2] - p[1][0]*p[2][0]*p[5][1]*p[6][2] - 2*p[1][0]*p[5][0]*p[5][1]*p[6][2] - p[2][0]*p[5][0]*p[5][1]*p[6][2] + 2*p[4][0]*p[5][0]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[6][0]*p[6][2] + p[1][0]*p[2][1]*p[6][0]*p[6][2] - p[2][1]*p[3][0]*p[6][0]*p[6][2] + p[2][0]*p[3][1]*p[6][0]*p[6][2] + p[1][1]*p[5][0]*p[6][0]*p[6][2] + 2*p[2][1]*p[5][0]*p[6][0]*p[6][2] - p[4][1]*p[5][0]*p[6][0]*p[6][2] - p[1][0]*p[5][1]*p[6][0]*p[6][2] - 2*p[2][0]*p[5][1]*p[6][0]*p[6][2] + p[4][0]*p[5][1]*p[6][0]*p[6][2] - p[2][0]*p[2][2]*p[3][1]*p[7][0] - 2*p[2][2]*p[3][0]*p[3][1]*p[7][0] + p[2][0]*p[2][1]*p[3][2]*p[7][0] - 2*p[0][1]*p[3][0]*p[3][2]*p[7][0] + 2*p[2][1]*p[3][0]*p[3][2]*p[7][0] - p[0][1]*p[3][2]*p[4][0]*p[7][0] - p[3][0]*p[3][2]*p[4][1]*p[7][0] - p[3][2]*p[4][0]*p[4][1]*p[7][0] + p[0][1]*p[3][0]*p[4][2]*p[7][0] + p[3][0]*p[3][1]*p[4][2]*p[7][0] + 2*p[0][1]*p[4][0]*p[4][2]*p[7][0] + p[3][1]*p[4][0]*p[4][2]*p[7][0] - 2*p[4][0]*p[4][2]*p[5][1]*p[7][0] - p[4][2]*p[5][0]*p[5][1]*p[7][0] + 2*p[4][0]*p[4][1]*p[5][2]*p[7][0] + p[4][1]*p[5][0]*p[5][2]*p[7][0] - p[2][2]*p[3][1]*p[6][0]*p[7][0] + p[2][1]*p[3][2]*p[6][0]*p[7][0] - p[4][2]*p[5][1]*p[6][0]*p[7][0] + p[4][1]*p[5][2]*p[6][0]*p[7][0] + p[2][0]*p[2][2]*p[6][1]*p[7][0] + p[2][2]*p[3][0]*p[6][1]*p[7][0] + p[3][0]*p[3][2]*p[6][1]*p[7][0] - p[4][0]*p[4][2]*p[6][1]*p[7][0] - p[4][0]*p[5][2]*p[6][1]*p[7][0] - p[5][0]*p[5][2]*p[6][1]*p[7][0] + 2*p[2][2]*p[6][0]*p[6][1]*p[7][0] + p[3][2]*p[6][0]*p[6][1]*p[7][0] - p[4][2]*p[6][0]*p[6][1]*p[7][0] - 2*p[5][2]*p[6][0]*p[6][1]*p[7][0] - p[2][0]*p[2][1]*p[6][2]*p[7][0] - p[2][1]*p[3][0]*p[6][2]*p[7][0] - p[3][0]*p[3][1]*p[6][2]*p[7][0] + p[4][0]*p[4][1]*p[6][2]*p[7][0] + p[4][0]*p[5][1]*p[6][2]*p[7][0] + p[5][0]*p[5][1]*p[6][2]*p[7][0] - 2*p[2][1]*p[6][0]*p[6][2]*p[7][0] - p[3][1]*p[6][0]*p[6][2]*p[7][0] + p[4][1]*p[6][0]*p[6][2]*p[7][0] + 2*p[5][1]*p[6][0]*p[6][2]*p[7][0] + 2*pow(p[3][0], 2)*p[2][2]*p[7][1] - 2*pow(p[6][0], 2)*p[2][2]*p[7][1] + p[2][0]*p[2][2]*p[3][0]*p[7][1] - pow(p[2][0], 2)*p[3][2]*p[7][1] + pow(p[4][0], 2)*p[3][2]*p[7][1] - pow(p[6][0], 2)*p[3][2]*p[7][1] - 2*p[2][0]*p[3][0]*p[3][2]*p[7][1] + p[3][0]*p[3][2]*p[4][0]*p[7][1] - pow(p[3][0], 2)*p[4][2]*p[7][1] + pow(p[5][0], 2)*p[4][2]*p[7][1] + pow(p[6][0], 2)*p[4][2]*p[7][1] - p[3][0]*p[4][0]*p[4][2]*p[7][1] + 2*p[4][0]*p[4][2]*p[5][0]*p[7][1] - 2*pow(p[4][0], 2)*p[5][2]*p[7][1] + 2*pow(p[6][0], 2)*p[5][2]*p[7][1] - p[4][0]*p[5][0]*p[5][2]*p[7][1] - p[2][0]*p[2][2]*p[6][0]*p[7][1] - p[2][0]*p[3][2]*p[6][0]*p[7][1] - p[3][0]*p[3][2]*p[6][0]*p[7][1] + p[4][0]*p[4][2]*p[6][0]*p[7][1] + p[4][2]*p[5][0]*p[6][0]*p[7][1] + p[5][0]*p[5][2]*p[6][0]*p[7][1] + pow(p[2][0], 2)*p[6][2]*p[7][1] + pow(p[3][0], 2)*p[6][2]*p[7][1] - pow(p[4][0], 2)*p[6][2]*p[7][1] - pow(p[5][0], 2)*p[6][2]*p[7][1] + p[2][0]*p[3][0]*p[6][2]*p[7][1] - p[4][0]*p[5][0]*p[6][2]*p[7][1] + 2*p[2][0]*p[6][0]*p[6][2]*p[7][1] + p[3][0]*p[6][0]*p[6][2]*p[7][1] - p[4][0]*p[6][0]*p[6][2]*p[7][1] - 2*p[5][0]*p[6][0]*p[6][2]*p[7][1] + p[2][2]*p[3][0]*p[7][0]*p[7][1] - p[2][0]*p[3][2]*p[7][0]*p[7][1] + 2*p[3][2]*p[4][0]*p[7][0]*p[7][1] - 2*p[3][0]*p[4][2]*p[7][0]*p[7][1] + p[4][2]*p[5][0]*p[7][0]*p[7][1] - p[4][0]*p[5][2]*p[7][0]*p[7][1] - p[2][2]*p[6][0]*p[7][0]*p[7][1] - 2*p[3][2]*p[6][0]*p[7][0]*p[7][1] + 2*p[4][2]*p[6][0]*p[7][0]*p[7][1] + p[5][2]*p[6][0]*p[7][0]*p[7][1] + p[2][0]*p[6][2]*p[7][0]*p[7][1] + 2*p[3][0]*p[6][2]*p[7][0]*p[7][1] - 2*p[4][0]*p[6][2]*p[7][0]*p[7][1] - p[5][0]*p[6][2]*p[7][0]*p[7][1] + p[0][2]*(2*pow(p[3][0], 2)*p[2][1] + p[2][0]*p[2][1]*p[3][0] - pow(p[2][0], 2)*p[3][1] + pow(p[4][0], 2)*p[3][1] + pow(p[7][0], 2)*p[3][1] - 2*p[2][0]*p[3][0]*p[3][1] + p[3][0]*p[3][1]*p[4][0] - pow(p[3][0], 2)*p[4][1] + pow(p[5][0], 2)*p[4][1] - pow(p[7][0], 2)*p[4][1] - p[3][0]*p[4][0]*p[4][1] + 2*p[4][0]*p[4][1]*p[5][0] + p[1][1]*(pow(p[2][0], 2) + pow(p[3][0], 2) - pow(p[4][0], 2) - pow(p[5][0], 2) + p[2][0]*p[3][0] - p[4][0]*p[5][0]) - 2*pow(p[4][0], 2)*p[5][1] - p[4][0]*p[5][0]*p[5][1] + pow(p[1][0], 2)*(-2*p[2][1] - p[3][1] + p[4][1] + 2*p[5][1]) + p[1][0]*(-(p[3][0]*p[3][1]) - p[2][0]*(p[2][1] + p[3][1]) + p[4][0]*p[4][1] + p[1][1]*(2*p[2][0] + p[3][0] - p[4][0] - 2*p[5][0]) + p[4][1]*p[5][0] + p[5][0]*p[5][1]) + 2*p[3][0]*p[3][1]*p[7][0] + p[3][1]*p[4][0]*p[7][0] - p[3][0]*p[4][1]*p[7][0] - 2*p[4][0]*p[4][1]*p[7][0] - 2*pow(p[3][0], 2)*p[7][1] + 2*pow(p[4][0], 2)*p[7][1] - p[3][0]*p[7][0]*p[7][1] + p[4][0]*p[7][0]*p[7][1]) + 2*pow(p[3][0], 2)*p[0][1]*p[7][2] - 2*pow(p[4][0], 2)*p[0][1]*p[7][2] - 2*pow(p[3][0], 2)*p[2][1]*p[7][2] + 2*pow(p[6][0], 2)*p[2][1]*p[7][2] - p[2][0]*p[2][1]*p[3][0]*p[7][2] + pow(p[2][0], 2)*p[3][1]*p[7][2] - pow(p[4][0], 2)*p[3][1]*p[7][2] + pow(p[6][0], 2)*p[3][1]*p[7][2] + 2*p[2][0]*p[3][0]*p[3][1]*p[7][2] - p[3][0]*p[3][1]*p[4][0]*p[7][2] + pow(p[3][0], 2)*p[4][1]*p[7][2] - pow(p[5][0], 2)*p[4][1]*p[7][2] - pow(p[6][0], 2)*p[4][1]*p[7][2] + p[3][0]*p[4][0]*p[4][1]*p[7][2] - 2*p[4][0]*p[4][1]*p[5][0]*p[7][2] + 2*pow(p[4][0], 2)*p[5][1]*p[7][2] - 2*pow(p[6][0], 2)*p[5][1]*p[7][2] + p[4][0]*p[5][0]*p[5][1]*p[7][2] + p[2][0]*p[2][1]*p[6][0]*p[7][2] + p[2][0]*p[3][1]*p[6][0]*p[7][2] + p[3][0]*p[3][1]*p[6][0]*p[7][2] - p[4][0]*p[4][1]*p[6][0]*p[7][2] - p[4][1]*p[5][0]*p[6][0]*p[7][2] - p[5][0]*p[5][1]*p[6][0]*p[7][2] - pow(p[2][0], 2)*p[6][1]*p[7][2] - pow(p[3][0], 2)*p[6][1]*p[7][2] + pow(p[4][0], 2)*p[6][1]*p[7][2] + pow(p[5][0], 2)*p[6][1]*p[7][2] - p[2][0]*p[3][0]*p[6][1]*p[7][2] + p[4][0]*p[5][0]*p[6][1]*p[7][2] - 2*p[2][0]*p[6][0]*p[6][1]*p[7][2] - p[3][0]*p[6][0]*p[6][1]*p[7][2] + p[4][0]*p[6][0]*p[6][1]*p[7][2] + 2*p[5][0]*p[6][0]*p[6][1]*p[7][2] + p[0][1]*p[3][0]*p[7][0]*p[7][2] - p[2][1]*p[3][0]*p[7][0]*p[7][2] + p[2][0]*p[3][1]*p[7][0]*p[7][2] - p[0][1]*p[4][0]*p[7][0]*p[7][2] - 2*p[3][1]*p[4][0]*p[7][0]*p[7][2] + 2*p[3][0]*p[4][1]*p[7][0]*p[7][2] - p[4][1]*p[5][0]*p[7][0]*p[7][2] + p[4][0]*p[5][1]*p[7][0]*p[7][2] + p[2][1]*p[6][0]*p[7][0]*p[7][2] + 2*p[3][1]*p[6][0]*p[7][0]*p[7][2] - 2*p[4][1]*p[6][0]*p[7][0]*p[7][2] - p[5][1]*p[6][0]*p[7][0]*p[7][2] - p[2][0]*p[6][1]*p[7][0]*p[7][2] - 2*p[3][0]*p[6][1]*p[7][0]*p[7][2] + 2*p[4][0]*p[6][1]*p[7][0]*p[7][2] + p[5][0]*p[6][1]*p[7][0]*p[7][2] + pow(p[0][0], 2)*(p[2][2]*p[3][1] - p[2][1]*p[3][2] + 2*p[3][2]*p[4][1] - 2*p[3][1]*p[4][2] + p[1][2]*(p[2][1] + 2*p[3][1] - 2*p[4][1] - p[5][1]) + p[4][2]*p[5][1] - p[4][1]*p[5][2] + p[1][1]*(-p[2][2] - 2*p[3][2] + 2*p[4][2] + p[5][2]) + p[3][2]*p[7][1] - p[4][2]*p[7][1] - p[3][1]*p[7][2] + p[4][1]*p[7][2]) + p[0][0]*(2*p[1][0]*p[1][2]*p[2][1] + p[1][2]*p[2][0]*p[2][1] - 2*p[1][0]*p[1][1]*p[2][2] - p[1][1]*p[2][0]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[1][0]*p[1][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] + p[2][0]*p[2][2]*p[3][1] + p[1][2]*p[3][0]*p[3][1] + 2*p[2][2]*p[3][0]*p[3][1] - p[1][0]*p[1][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] - p[2][0]*p[2][1]*p[3][2] - p[1][1]*p[3][0]*p[3][2] - 2*p[2][1]*p[3][0]*p[3][2] - p[1][0]*p[1][2]*p[4][1] + p[3][0]*p[3][2]*p[4][1] - p[1][2]*p[4][0]*p[4][1] + p[3][2]*p[4][0]*p[4][1] + p[1][0]*p[1][1]*p[4][2] - p[3][0]*p[3][1]*p[4][2] + p[1][1]*p[4][0]*p[4][2] - p[3][1]*p[4][0]*p[4][2] - 2*p[1][0]*p[1][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + 2*p[4][0]*p[4][2]*p[5][1] - p[1][2]*p[5][0]*p[5][1] + p[4][2]*p[5][0]*p[5][1] + 2*p[1][0]*p[1][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - 2*p[4][0]*p[4][1]*p[5][2] + p[1][1]*p[5][0]*p[5][2] - p[4][1]*p[5][0]*p[5][2] + 2*p[3][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[3][0]*p[4][2]*p[7][1] - 2*p[4][0]*p[4][2]*p[7][1] + p[3][2]*p[7][0]*p[7][1] - p[4][2]*p[7][0]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + 2*p[3][1]*p[4][0] - 2*p[3][0]*p[4][1] + p[1][1]*(p[2][0] + 2*p[3][0] - 2*p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - 2*p[3][1] + 2*p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - 2*p[3][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[3][0]*p[4][1]*p[7][2] + 2*p[4][0]*p[4][1]*p[7][2] - p[3][1]*p[7][0]*p[7][2] + p[4][1]*p[7][0]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - 2*p[3][2]*p[4][0] + 2*p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - 2*p[3][0] + 2*p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + 2*p[3][2] - 2*p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2])))/72;
    result[2][0] = (pow(p[2][1], 2)*p[0][0]*p[1][2] + pow(p[3][1], 2)*p[0][0]*p[1][2] - pow(p[4][1], 2)*p[0][0]*p[1][2] - pow(p[5][1], 2)*p[0][0]*p[1][2] - pow(p[3][1], 2)*p[1][2]*p[2][0] + pow(p[5][1], 2)*p[1][2]*p[2][0] + pow(p[6][1], 2)*p[1][2]*p[2][0] + 2*p[0][0]*p[1][1]*p[1][2]*p[2][1] - 2*pow(p[1][1], 2)*p[0][0]*p[2][2] + 2*pow(p[3][1], 2)*p[0][0]*p[2][2] + pow(p[3][1], 2)*p[1][0]*p[2][2] - pow(p[5][1], 2)*p[1][0]*p[2][2] - pow(p[6][1], 2)*p[1][0]*p[2][2] - p[0][0]*p[1][1]*p[2][1]*p[2][2] + 2*pow(p[2][1], 2)*p[1][2]*p[3][0] + p[1][1]*p[1][2]*p[2][1]*p[3][0] - pow(p[1][1], 2)*p[2][2]*p[3][0] + pow(p[6][1], 2)*p[2][2]*p[3][0] + pow(p[7][1], 2)*p[2][2]*p[3][0] - 2*p[1][1]*p[2][1]*p[2][2]*p[3][0] + p[0][0]*p[1][1]*p[1][2]*p[3][1] - p[1][1]*p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[1][2]*p[2][1]*p[3][1] - 2*p[1][2]*p[2][0]*p[2][1]*p[3][1] + p[1][0]*p[1][1]*p[2][2]*p[3][1] + p[0][0]*p[2][1]*p[2][2]*p[3][1] + 2*p[1][0]*p[2][1]*p[2][2]*p[3][1] + p[1][2]*p[2][1]*p[3][0]*p[3][1] - p[1][1]*p[2][2]*p[3][0]*p[3][1] - pow(p[1][1], 2)*p[0][0]*p[3][2] - pow(p[2][1], 2)*p[0][0]*p[3][2] + pow(p[4][1], 2)*p[0][0]*p[3][2] + pow(p[7][1], 2)*p[0][0]*p[3][2] - 2*pow(p[2][1], 2)*p[1][0]*p[3][2] + pow(p[1][1], 2)*p[2][0]*p[3][2] - pow(p[6][1], 2)*p[2][0]*p[3][2] - pow(p[7][1], 2)*p[2][0]*p[3][2] - p[0][0]*p[1][1]*p[2][1]*p[3][2] - p[1][0]*p[1][1]*p[2][1]*p[3][2] + 2*p[1][1]*p[2][0]*p[2][1]*p[3][2] - p[0][0]*p[1][1]*p[3][1]*p[3][2] + p[1][1]*p[2][0]*p[3][1]*p[3][2] - 2*p[0][0]*p[2][1]*p[3][1]*p[3][2] - p[1][0]*p[2][1]*p[3][1]*p[3][2] - 2*pow(p[5][1], 2)*p[1][2]*p[4][0] + 2*pow(p[7][1], 2)*p[3][2]*p[4][0] - p[0][0]*p[1][1]*p[1][2]*p[4][1] + p[0][0]*p[3][1]*p[3][2]*p[4][1] + pow(p[1][1], 2)*p[0][0]*p[4][2] - pow(p[3][1], 2)*p[0][0]*p[4][2] + pow(p[5][1], 2)*p[0][0]*p[4][2] - pow(p[7][1], 2)*p[0][0]*p[4][2] + 2*pow(p[5][1], 2)*p[1][0]*p[4][2] - 2*pow(p[7][1], 2)*p[3][0]*p[4][2] + p[0][0]*p[1][1]*p[4][1]*p[4][2] - p[0][0]*p[3][1]*p[4][1]*p[4][2] - pow(p[2][1], 2)*p[1][2]*p[5][0] + pow(p[4][1], 2)*p[1][2]*p[5][0] - pow(p[6][1], 2)*p[1][2]*p[5][0] - 2*p[1][1]*p[1][2]*p[2][1]*p[5][0] + 2*pow(p[1][1], 2)*p[2][2]*p[5][0] - 2*pow(p[6][1], 2)*p[2][2]*p[5][0] + p[1][1]*p[2][1]*p[2][2]*p[5][0] + p[1][1]*p[1][2]*p[4][1]*p[5][0] - pow(p[1][1], 2)*p[4][2]*p[5][0] + pow(p[6][1], 2)*p[4][2]*p[5][0] + pow(p[7][1], 2)*p[4][2]*p[5][0] - p[1][1]*p[4][1]*p[4][2]*p[5][0] - 2*p[0][0]*p[1][1]*p[1][2]*p[5][1] + 2*p[1][1]*p[1][2]*p[2][0]*p[5][1] + p[1][2]*p[2][0]*p[2][1]*p[5][1] - 2*p[1][0]*p[1][1]*p[2][2]*p[5][1] - p[1][0]*p[2][1]*p[2][2]*p[5][1] - p[1][1]*p[1][2]*p[4][0]*p[5][1] - p[0][0]*p[1][2]*p[4][1]*p[5][1] - p[1][2]*p[4][0]*p[4][1]*p[5][1] + p[0][0]*p[1][1]*p[4][2]*p[5][1] + p[1][0]*p[1][1]*p[4][2]*p[5][1] + 2*p[0][0]*p[4][1]*p[4][2]*p[5][1] + p[1][0]*p[4][1]*p[4][2]*p[5][1] - p[1][2]*p[2][1]*p[5][0]*p[5][1] + p[1][1]*p[2][2]*p[5][0]*p[5][1] + 2*p[1][2]*p[4][1]*p[5][0]*p[5][1] - 2*p[1][1]*p[4][2]*p[5][0]*p[5][1] + 2*pow(p[1][1], 2)*p[0][0]*p[5][2] - 2*pow(p[4][1], 2)*p[0][0]*p[5][2] + pow(p[2][1], 2)*p[1][0]*p[5][2] - pow(p[4][1], 2)*p[1][0]*p[5][2] + pow(p[6][1], 2)*p[1][0]*p[5][2] - 2*pow(p[1][1], 2)*p[2][0]*p[5][2] + 2*pow(p[6][1], 2)*p[2][0]*p[5][2] + 2*p[1][0]*p[1][1]*p[2][1]*p[5][2] - p[1][1]*p[2][0]*p[2][1]*p[5][2] + pow(p[1][1], 2)*p[4][0]*p[5][2] - pow(p[6][1], 2)*p[4][0]*p[5][2] - pow(p[7][1], 2)*p[4][0]*p[5][2] - p[1][0]*p[1][1]*p[4][1]*p[5][2] + p[1][1]*p[4][0]*p[4][1]*p[5][2] + p[0][0]*p[1][1]*p[5][1]*p[5][2] - p[1][1]*p[2][0]*p[5][1]*p[5][2] + p[1][0]*p[2][1]*p[5][1]*p[5][2] + 2*p[1][1]*p[4][0]*p[5][1]*p[5][2] - p[0][0]*p[4][1]*p[5][1]*p[5][2] - 2*p[1][0]*p[4][1]*p[5][1]*p[5][2] - 2*pow(p[2][1], 2)*p[1][2]*p[6][0] + 2*pow(p[5][1], 2)*p[1][2]*p[6][0] - p[1][1]*p[1][2]*p[2][1]*p[6][0] + pow(p[1][1], 2)*p[2][2]*p[6][0] - pow(p[3][1], 2)*p[2][2]*p[6][0] + pow(p[5][1], 2)*p[2][2]*p[6][0] - pow(p[7][1], 2)*p[2][2]*p[6][0] + 2*p[1][1]*p[2][1]*p[2][2]*p[6][0] - 2*p[2][1]*p[2][2]*p[3][1]*p[6][0] + 2*pow(p[2][1], 2)*p[3][2]*p[6][0] - 2*pow(p[7][1], 2)*p[3][2]*p[6][0] + p[2][1]*p[3][1]*p[3][2]*p[6][0] - 2*pow(p[5][1], 2)*p[4][2]*p[6][0] + 2*pow(p[7][1], 2)*p[4][2]*p[6][0] + p[1][1]*p[1][2]*p[5][1]*p[6][0] + p[1][1]*p[2][2]*p[5][1]*p[6][0] + p[2][1]*p[2][2]*p[5][1]*p[6][0] - p[4][1]*p[4][2]*p[5][1]*p[6][0] - pow(p[1][1], 2)*p[5][2]*p[6][0] - pow(p[2][1], 2)*p[5][2]*p[6][0] + pow(p[4][1], 2)*p[5][2]*p[6][0] + pow(p[7][1], 2)*p[5][2]*p[6][0] - p[1][1]*p[2][1]*p[5][2]*p[6][0] - 2*p[1][1]*p[5][1]*p[5][2]*p[6][0] - p[2][1]*p[5][1]*p[5][2]*p[6][0] + 2*p[4][1]*p[5][1]*p[5][2]*p[6][0] + p[1][1]*p[1][2]*p[2][0]*p[6][1] + 2*p[1][2]*p[2][0]*p[2][1]*p[6][1] - p[1][0]*p[1][1]*p[2][2]*p[6][1] - 2*p[1][0]*p[2][1]*p[2][2]*p[6][1] + 2*p[2][1]*p[2][2]*p[3][0]*p[6][1] + p[2][2]*p[3][0]*p[3][1]*p[6][1] - 2*p[2][0]*p[2][1]*p[3][2]*p[6][1] - p[2][0]*p[3][1]*p[3][2]*p[6][1] - p[1][1]*p[1][2]*p[5][0]*p[6][1] - p[1][2]*p[2][1]*p[5][0]*p[6][1] - p[2][1]*p[2][2]*p[5][0]*p[6][1] + p[4][1]*p[4][2]*p[5][0]*p[6][1] + p[1][2]*p[2][0]*p[5][1]*p[6][1] - p[1][0]*p[2][2]*p[5][1]*p[6][1] - 2*p[1][2]*p[5][0]*p[5][1]*p[6][1] - p[2][2]*p[5][0]*p[5][1]*p[6][1] + 2*p[4][2]*p[5][0]*p[5][1]*p[6][1] + p[1][0]*p[1][1]*p[5][2]*p[6][1] + p[1][0]*p[2][1]*p[5][2]*p[6][1] + p[2][0]*p[2][1]*p[5][2]*p[6][1] - p[4][0]*p[4][1]*p[5][2]*p[6][1] + 2*p[1][0]*p[5][1]*p[5][2]*p[6][1] + p[2][0]*p[5][1]*p[5][2]*p[6][1] - 2*p[4][0]*p[5][1]*p[5][2]*p[6][1] - p[1][2]*p[2][1]*p[6][0]*p[6][1] + p[1][1]*p[2][2]*p[6][0]*p[6][1] - p[2][2]*p[3][1]*p[6][0]*p[6][1] + p[2][1]*p[3][2]*p[6][0]*p[6][1] + p[1][2]*p[5][1]*p[6][0]*p[6][1] + 2*p[2][2]*p[5][1]*p[6][0]*p[6][1] - p[4][2]*p[5][1]*p[6][0]*p[6][1] - p[1][1]*p[5][2]*p[6][0]*p[6][1] - 2*p[2][1]*p[5][2]*p[6][0]*p[6][1] + p[4][1]*p[5][2]*p[6][0]*p[6][1] + 2*pow(p[2][1], 2)*p[1][0]*p[6][2] - 2*pow(p[5][1], 2)*p[1][0]*p[6][2] - pow(p[1][1], 2)*p[2][0]*p[6][2] + pow(p[3][1], 2)*p[2][0]*p[6][2] - pow(p[5][1], 2)*p[2][0]*p[6][2] + pow(p[7][1], 2)*p[2][0]*p[6][2] + p[1][0]*p[1][1]*p[2][1]*p[6][2] - 2*p[1][1]*p[2][0]*p[2][1]*p[6][2] - 2*pow(p[2][1], 2)*p[3][0]*p[6][2] + 2*pow(p[7][1], 2)*p[3][0]*p[6][2] + 2*p[2][0]*p[2][1]*p[3][1]*p[6][2] - p[2][1]*p[3][0]*p[3][1]*p[6][2] + 2*pow(p[5][1], 2)*p[4][0]*p[6][2] - 2*pow(p[7][1], 2)*p[4][0]*p[6][2] + pow(p[1][1], 2)*p[5][0]*p[6][2] + pow(p[2][1], 2)*p[5][0]*p[6][2] - pow(p[4][1], 2)*p[5][0]*p[6][2] - pow(p[7][1], 2)*p[5][0]*p[6][2] + p[1][1]*p[2][1]*p[5][0]*p[6][2] - p[1][0]*p[1][1]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[5][1]*p[6][2] - p[2][0]*p[2][1]*p[5][1]*p[6][2] + p[4][0]*p[4][1]*p[5][1]*p[6][2] + 2*p[1][1]*p[5][0]*p[5][1]*p[6][2] + p[2][1]*p[5][0]*p[5][1]*p[6][2] - 2*p[4][1]*p[5][0]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[6][1]*p[6][2] + p[1][0]*p[2][1]*p[6][1]*p[6][2] - p[2][1]*p[3][0]*p[6][1]*p[6][2] + p[2][0]*p[3][1]*p[6][1]*p[6][2] + p[1][1]*p[5][0]*p[6][1]*p[6][2] + 2*p[2][1]*p[5][0]*p[6][1]*p[6][2] - p[4][1]*p[5][0]*p[6][1]*p[6][2] - p[1][0]*p[5][1]*p[6][1]*p[6][2] - 2*p[2][0]*p[5][1]*p[6][1]*p[6][2] + p[4][0]*p[5][1]*p[6][1]*p[6][2] - 2*pow(p[3][1], 2)*p[2][2]*p[7][0] + 2*pow(p[6][1], 2)*p[2][2]*p[7][0] - p[2][1]*p[2][2]*p[3][1]*p[7][0] + pow(p[2][1], 2)*p[3][2]*p[7][0] - pow(p[4][1], 2)*p[3][2]*p[7][0] + pow(p[6][1], 2)*p[3][2]*p[7][0] + 2*p[2][1]*p[3][1]*p[3][2]*p[7][0] - p[3][1]*p[3][2]*p[4][1]*p[7][0] + pow(p[3][1], 2)*p[4][2]*p[7][0] - pow(p[5][1], 2)*p[4][2]*p[7][0] - pow(p[6][1], 2)*p[4][2]*p[7][0] + p[3][1]*p[4][1]*p[4][2]*p[7][0] - 2*p[4][1]*p[4][2]*p[5][1]*p[7][0] + 2*pow(p[4][1], 2)*p[5][2]*p[7][0] - 2*pow(p[6][1], 2)*p[5][2]*p[7][0] + p[4][1]*p[5][1]*p[5][2]*p[7][0] + p[2][1]*p[2][2]*p[6][1]*p[7][0] + p[2][1]*p[3][2]*p[6][1]*p[7][0] + p[3][1]*p[3][2]*p[6][1]*p[7][0] - p[4][1]*p[4][2]*p[6][1]*p[7][0] - p[4][2]*p[5][1]*p[6][1]*p[7][0] - p[5][1]*p[5][2]*p[6][1]*p[7][0] - pow(p[2][1], 2)*p[6][2]*p[7][0] - pow(p[3][1], 2)*p[6][2]*p[7][0] + pow(p[4][1], 2)*p[6][2]*p[7][0] + pow(p[5][1], 2)*p[6][2]*p[7][0] - p[2][1]*p[3][1]*p[6][2]*p[7][0] + p[4][1]*p[5][1]*p[6][2]*p[7][0] - 2*p[2][1]*p[6][1]*p[6][2]*p[7][0] - p[3][1]*p[6][1]*p[6][2]*p[7][0] + p[4][1]*p[6][1]*p[6][2]*p[7][0] + 2*p[5][1]*p[6][1]*p[6][2]*p[7][0] + p[2][1]*p[2][2]*p[3][0]*p[7][1] + 2*p[2][2]*p[3][0]*p[3][1]*p[7][1] - p[2][0]*p[2][1]*p[3][2]*p[7][1] + 2*p[0][0]*p[3][1]*p[3][2]*p[7][1] - 2*p[2][0]*p[3][1]*p[3][2]*p[7][1] + p[3][1]*p[3][2]*p[4][0]*p[7][1] + p[0][0]*p[3][2]*p[4][1]*p[7][1] + p[3][2]*p[4][0]*p[4][1]*p[7][1] - p[0][0]*p[3][1]*p[4][2]*p[7][1] - p[3][0]*p[3][1]*p[4][2]*p[7][1] - 2*p[0][0]*p[4][1]*p[4][2]*p[7][1] - p[3][0]*p[4][1]*p[4][2]*p[7][1] + 2*p[4][1]*p[4][2]*p[5][0]*p[7][1] + p[4][2]*p[5][0]*p[5][1]*p[7][1] - 2*p[4][0]*p[4][1]*p[5][2]*p[7][1] - p[4][0]*p[5][1]*p[5][2]*p[7][1] - p[2][1]*p[2][2]*p[6][0]*p[7][1] - p[2][2]*p[3][1]*p[6][0]*p[7][1] - p[3][1]*p[3][2]*p[6][0]*p[7][1] + p[4][1]*p[4][2]*p[6][0]*p[7][1] + p[4][1]*p[5][2]*p[6][0]*p[7][1] + p[5][1]*p[5][2]*p[6][0]*p[7][1] + p[2][2]*p[3][0]*p[6][1]*p[7][1] - p[2][0]*p[3][2]*p[6][1]*p[7][1] + p[4][2]*p[5][0]*p[6][1]*p[7][1] - p[4][0]*p[5][2]*p[6][1]*p[7][1] - 2*p[2][2]*p[6][0]*p[6][1]*p[7][1] - p[3][2]*p[6][0]*p[6][1]*p[7][1] + p[4][2]*p[6][0]*p[6][1]*p[7][1] + 2*p[5][2]*p[6][0]*p[6][1]*p[7][1] + p[2][0]*p[2][1]*p[6][2]*p[7][1] + p[2][0]*p[3][1]*p[6][2]*p[7][1] + p[3][0]*p[3][1]*p[6][2]*p[7][1] - p[4][0]*p[4][1]*p[6][2]*p[7][1] - p[4][1]*p[5][0]*p[6][2]*p[7][1] - p[5][0]*p[5][1]*p[6][2]*p[7][1] + 2*p[2][0]*p[6][1]*p[6][2]*p[7][1] + p[3][0]*p[6][1]*p[6][2]*p[7][1] - p[4][0]*p[6][1]*p[6][2]*p[7][1] - 2*p[5][0]*p[6][1]*p[6][2]*p[7][1] - p[2][2]*p[3][1]*p[7][0]*p[7][1] + p[2][1]*p[3][2]*p[7][0]*p[7][1] - 2*p[3][2]*p[4][1]*p[7][0]*p[7][1] + 2*p[3][1]*p[4][2]*p[7][0]*p[7][1] - p[4][2]*p[5][1]*p[7][0]*p[7][1] + p[4][1]*p[5][2]*p[7][0]*p[7][1] + p[2][2]*p[6][1]*p[7][0]*p[7][1] + 2*p[3][2]*p[6][1]*p[7][0]*p[7][1] - 2*p[4][2]*p[6][1]*p[7][0]*p[7][1] - p[5][2]*p[6][1]*p[7][0]*p[7][1] - p[2][1]*p[6][2]*p[7][0]*p[7][1] - 2*p[3][1]*p[6][2]*p[7][0]*p[7][1] + 2*p[4][1]*p[6][2]*p[7][0]*p[7][1] + p[5][1]*p[6][2]*p[7][0]*p[7][1] + p[0][2]*(-2*pow(p[3][1], 2)*p[2][0] + pow(p[2][1], 2)*p[3][0] - pow(p[4][1], 2)*p[3][0] - pow(p[7][1], 2)*p[3][0] - p[2][0]*p[2][1]*p[3][1] + 2*p[2][1]*p[3][0]*p[3][1] + pow(p[3][1], 2)*p[4][0] - pow(p[5][1], 2)*p[4][0] + pow(p[7][1], 2)*p[4][0] - p[3][0]*p[3][1]*p[4][1] + p[3][1]*p[4][0]*p[4][1] + pow(p[1][1], 2)*(2*p[2][0] + p[3][0] - p[4][0] - 2*p[5][0]) + 2*pow(p[4][1], 2)*p[5][0] - 2*p[4][0]*p[4][1]*p[5][1] + p[4][1]*p[5][0]*p[5][1] + p[1][0]*(-pow(p[2][1], 2) - pow(p[3][1], 2) + pow(p[4][1], 2) + pow(p[5][1], 2) - p[2][1]*p[3][1] + p[4][1]*p[5][1]) + p[1][1]*(p[2][0]*p[2][1] + p[2][1]*p[3][0] + p[3][0]*p[3][1] - p[4][0]*p[4][1] - p[4][0]*p[5][1] - p[5][0]*p[5][1] + p[1][0]*(-2*p[2][1] - p[3][1] + p[4][1] + 2*p[5][1])) + 2*pow(p[3][1], 2)*p[7][0] - 2*pow(p[4][1], 2)*p[7][0] - 2*p[3][0]*p[3][1]*p[7][1] + p[3][1]*p[4][0]*p[7][1] - p[3][0]*p[4][1]*p[7][1] + 2*p[4][0]*p[4][1]*p[7][1] + p[3][1]*p[7][0]*p[7][1] - p[4][1]*p[7][0]*p[7][1]) - 2*pow(p[3][1], 2)*p[0][0]*p[7][2] + 2*pow(p[4][1], 2)*p[0][0]*p[7][2] + 2*pow(p[3][1], 2)*p[2][0]*p[7][2] - 2*pow(p[6][1], 2)*p[2][0]*p[7][2] - pow(p[2][1], 2)*p[3][0]*p[7][2] + pow(p[4][1], 2)*p[3][0]*p[7][2] - pow(p[6][1], 2)*p[3][0]*p[7][2] + p[2][0]*p[2][1]*p[3][1]*p[7][2] - 2*p[2][1]*p[3][0]*p[3][1]*p[7][2] - pow(p[3][1], 2)*p[4][0]*p[7][2] + pow(p[5][1], 2)*p[4][0]*p[7][2] + pow(p[6][1], 2)*p[4][0]*p[7][2] + p[3][0]*p[3][1]*p[4][1]*p[7][2] - p[3][1]*p[4][0]*p[4][1]*p[7][2] - 2*pow(p[4][1], 2)*p[5][0]*p[7][2] + 2*pow(p[6][1], 2)*p[5][0]*p[7][2] + 2*p[4][0]*p[4][1]*p[5][1]*p[7][2] - p[4][1]*p[5][0]*p[5][1]*p[7][2] + pow(p[2][1], 2)*p[6][0]*p[7][2] + pow(p[3][1], 2)*p[6][0]*p[7][2] - pow(p[4][1], 2)*p[6][0]*p[7][2] - pow(p[5][1], 2)*p[6][0]*p[7][2] + p[2][1]*p[3][1]*p[6][0]*p[7][2] - p[4][1]*p[5][1]*p[6][0]*p[7][2] - p[2][0]*p[2][1]*p[6][1]*p[7][2] - p[2][1]*p[3][0]*p[6][1]*p[7][2] - p[3][0]*p[3][1]*p[6][1]*p[7][2] + p[4][0]*p[4][1]*p[6][1]*p[7][2] + p[4][0]*p[5][1]*p[6][1]*p[7][2] + p[5][0]*p[5][1]*p[6][1]*p[7][2] + 2*p[2][1]*p[6][0]*p[6][1]*p[7][2] + p[3][1]*p[6][0]*p[6][1]*p[7][2] - p[4][1]*p[6][0]*p[6][1]*p[7][2] - 2*p[5][1]*p[6][0]*p[6][1]*p[7][2] - p[2][1]*p[3][0]*p[7][1]*p[7][2] - p[0][0]*p[3][1]*p[7][1]*p[7][2] + p[2][0]*p[3][1]*p[7][1]*p[7][2] - 2*p[3][1]*p[4][0]*p[7][1]*p[7][2] + p[0][0]*p[4][1]*p[7][1]*p[7][2] + 2*p[3][0]*p[4][1]*p[7][1]*p[7][2] - p[4][1]*p[5][0]*p[7][1]*p[7][2] + p[4][0]*p[5][1]*p[7][1]*p[7][2] + p[2][1]*p[6][0]*p[7][1]*p[7][2] + 2*p[3][1]*p[6][0]*p[7][1]*p[7][2] - 2*p[4][1]*p[6][0]*p[7][1]*p[7][2] - p[5][1]*p[6][0]*p[7][1]*p[7][2] - p[2][0]*p[6][1]*p[7][1]*p[7][2] - 2*p[3][0]*p[6][1]*p[7][1]*p[7][2] + 2*p[4][0]*p[6][1]*p[7][1]*p[7][2] + p[5][0]*p[6][1]*p[7][1]*p[7][2] + pow(p[0][1], 2)*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - 2*p[3][2]*p[4][0] + 2*p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - 2*p[3][0] + 2*p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + 2*p[3][2] - 2*p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2]) + p[0][1]*(p[0][0]*p[1][2]*p[2][1] - p[1][2]*p[2][0]*p[2][1] + p[1][0]*p[2][1]*p[2][2] - p[2][1]*p[2][2]*p[3][0] + 2*p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] - p[1][2]*p[3][0]*p[3][1] - 2*p[2][2]*p[3][0]*p[3][1] - p[0][0]*p[2][1]*p[3][2] + p[2][0]*p[2][1]*p[3][2] + p[1][0]*p[3][1]*p[3][2] + 2*p[2][0]*p[3][1]*p[3][2] - p[3][1]*p[3][2]*p[4][0] - 2*p[0][0]*p[1][2]*p[4][1] + 2*p[0][0]*p[3][2]*p[4][1] + p[1][2]*p[4][0]*p[4][1] - p[3][2]*p[4][0]*p[4][1] - 2*p[0][0]*p[3][1]*p[4][2] + p[3][0]*p[3][1]*p[4][2] - p[1][0]*p[4][1]*p[4][2] + p[3][0]*p[4][1]*p[4][2] + p[1][2]*p[4][1]*p[5][0] - 2*p[4][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][2]*p[5][0]*p[5][1] - p[4][2]*p[5][0]*p[5][1] - p[0][0]*p[4][1]*p[5][2] - p[1][0]*p[4][1]*p[5][2] + 2*p[4][0]*p[4][1]*p[5][2] - p[1][0]*p[5][1]*p[5][2] + p[4][0]*p[5][1]*p[5][2] + p[1][1]*(2*p[1][0]*p[2][2] - p[2][2]*p[3][0] + p[1][0]*p[3][2] + p[2][0]*p[3][2] - p[1][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-2*p[2][0] - p[3][0] + p[4][0] + 2*p[5][0]) - 2*p[1][0]*p[5][2] + p[4][0]*p[5][2] + p[0][0]*(-p[2][2] - 2*p[3][2] + 2*p[4][2] + p[5][2])) - 2*p[3][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] + 2*p[4][1]*p[4][2]*p[7][0] + p[0][0]*p[3][2]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][2]*p[7][0]*p[7][1] + p[4][2]*p[7][0]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + 2*p[3][1]*p[4][0] - 2*p[3][0]*p[4][1] + p[1][1]*(p[2][0] + 2*p[3][0] - 2*p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - 2*p[3][1] + 2*p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[0][0]*p[3][1]*p[7][2] + 2*p[3][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][0]*p[4][1]*p[7][2] - 2*p[4][0]*p[4][1]*p[7][2] + p[3][0]*p[7][1]*p[7][2] - p[4][0]*p[7][1]*p[7][2]))/72;
    result[3][0] = (-(pow(p[2][2], 2)*p[0][0]*p[1][1]) - pow(p[3][2], 2)*p[0][0]*p[1][1] + pow(p[4][2], 2)*p[0][0]*p[1][1] + pow(p[5][2], 2)*p[0][0]*p[1][1] + pow(p[3][2], 2)*p[1][1]*p[2][0] - pow(p[5][2], 2)*p[1][1]*p[2][0] - pow(p[6][2], 2)*p[1][1]*p[2][0] + 2*pow(p[1][2], 2)*p[0][0]*p[2][1] - 2*pow(p[3][2], 2)*p[0][0]*p[2][1] - pow(p[3][2], 2)*p[1][0]*p[2][1] + pow(p[5][2], 2)*p[1][0]*p[2][1] + pow(p[6][2], 2)*p[1][0]*p[2][1] - 2*p[0][0]*p[1][1]*p[1][2]*p[2][2] + p[0][0]*p[1][2]*p[2][1]*p[2][2] - 2*pow(p[2][2], 2)*p[1][1]*p[3][0] + pow(p[1][2], 2)*p[2][1]*p[3][0] - pow(p[6][2], 2)*p[2][1]*p[3][0] - pow(p[7][2], 2)*p[2][1]*p[3][0] - p[1][1]*p[1][2]*p[2][2]*p[3][0] + 2*p[1][2]*p[2][1]*p[2][2]*p[3][0] + pow(p[1][2], 2)*p[0][0]*p[3][1] + pow(p[2][2], 2)*p[0][0]*p[3][1] - pow(p[4][2], 2)*p[0][0]*p[3][1] - pow(p[7][2], 2)*p[0][0]*p[3][1] + 2*pow(p[2][2], 2)*p[1][0]*p[3][1] - pow(p[1][2], 2)*p[2][0]*p[3][1] + pow(p[6][2], 2)*p[2][0]*p[3][1] + pow(p[7][2], 2)*p[2][0]*p[3][1] + p[0][0]*p[1][2]*p[2][2]*p[3][1] + p[1][0]*p[1][2]*p[2][2]*p[3][1] - 2*p[1][2]*p[2][0]*p[2][2]*p[3][1] - p[0][0]*p[1][1]*p[1][2]*p[3][2] + p[1][1]*p[1][2]*p[2][0]*p[3][2] - p[1][0]*p[1][2]*p[2][1]*p[3][2] - p[0][0]*p[1][1]*p[2][2]*p[3][2] + 2*p[1][1]*p[2][0]*p[2][2]*p[3][2] - p[0][0]*p[2][1]*p[2][2]*p[3][2] - 2*p[1][0]*p[2][1]*p[2][2]*p[3][2] + p[1][2]*p[2][1]*p[3][0]*p[3][2] - p[1][1]*p[2][2]*p[3][0]*p[3][2] + p[0][0]*p[1][2]*p[3][1]*p[3][2] - p[1][2]*p[2][0]*p[3][1]*p[3][2] + 2*p[0][0]*p[2][2]*p[3][1]*p[3][2] + p[1][0]*p[2][2]*p[3][1]*p[3][2] + 2*pow(p[5][2], 2)*p[1][1]*p[4][0] - 2*pow(p[7][2], 2)*p[3][1]*p[4][0] - pow(p[1][2], 2)*p[0][0]*p[4][1] + pow(p[3][2], 2)*p[0][0]*p[4][1] - pow(p[5][2], 2)*p[0][0]*p[4][1] + pow(p[7][2], 2)*p[0][0]*p[4][1] - 2*pow(p[5][2], 2)*p[1][0]*p[4][1] + 2*pow(p[7][2], 2)*p[3][0]*p[4][1] + p[0][0]*p[1][1]*p[1][2]*p[4][2] - p[0][0]*p[3][1]*p[3][2]*p[4][2] - p[0][0]*p[1][2]*p[4][1]*p[4][2] + p[0][0]*p[3][2]*p[4][1]*p[4][2] + pow(p[2][2], 2)*p[1][1]*p[5][0] - pow(p[4][2], 2)*p[1][1]*p[5][0] + pow(p[6][2], 2)*p[1][1]*p[5][0] - 2*pow(p[1][2], 2)*p[2][1]*p[5][0] + 2*pow(p[6][2], 2)*p[2][1]*p[5][0] + 2*p[1][1]*p[1][2]*p[2][2]*p[5][0] - p[1][2]*p[2][1]*p[2][2]*p[5][0] + pow(p[1][2], 2)*p[4][1]*p[5][0] - pow(p[6][2], 2)*p[4][1]*p[5][0] - pow(p[7][2], 2)*p[4][1]*p[5][0] - p[1][1]*p[1][2]*p[4][2]*p[5][0] + p[1][2]*p[4][1]*p[4][2]*p[5][0] - 2*pow(p[1][2], 2)*p[0][0]*p[5][1] + 2*pow(p[4][2], 2)*p[0][0]*p[5][1] - pow(p[2][2], 2)*p[1][0]*p[5][1] + pow(p[4][2], 2)*p[1][0]*p[5][1] - pow(p[6][2], 2)*p[1][0]*p[5][1] + 2*pow(p[1][2], 2)*p[2][0]*p[5][1] - 2*pow(p[6][2], 2)*p[2][0]*p[5][1] - 2*p[1][0]*p[1][2]*p[2][2]*p[5][1] + p[1][2]*p[2][0]*p[2][2]*p[5][1] - pow(p[1][2], 2)*p[4][0]*p[5][1] + pow(p[6][2], 2)*p[4][0]*p[5][1] + pow(p[7][2], 2)*p[4][0]*p[5][1] + p[1][0]*p[1][2]*p[4][2]*p[5][1] - p[1][2]*p[4][0]*p[4][2]*p[5][1] + 2*p[0][0]*p[1][1]*p[1][2]*p[5][2] - 2*p[1][1]*p[1][2]*p[2][0]*p[5][2] + 2*p[1][0]*p[1][2]*p[2][1]*p[5][2] - p[1][1]*p[2][0]*p[2][2]*p[5][2] + p[1][0]*p[2][1]*p[2][2]*p[5][2] + p[1][1]*p[1][2]*p[4][0]*p[5][2] - p[0][0]*p[1][2]*p[4][1]*p[5][2] - p[1][0]*p[1][2]*p[4][1]*p[5][2] + p[0][0]*p[1][1]*p[4][2]*p[5][2] + p[1][1]*p[4][0]*p[4][2]*p[5][2] - 2*p[0][0]*p[4][1]*p[4][2]*p[5][2] - p[1][0]*p[4][1]*p[4][2]*p[5][2] - p[1][2]*p[2][1]*p[5][0]*p[5][2] + p[1][1]*p[2][2]*p[5][0]*p[5][2] + 2*p[1][2]*p[4][1]*p[5][0]*p[5][2] - 2*p[1][1]*p[4][2]*p[5][0]*p[5][2] - p[0][0]*p[1][2]*p[5][1]*p[5][2] + p[1][2]*p[2][0]*p[5][1]*p[5][2] - p[1][0]*p[2][2]*p[5][1]*p[5][2] - 2*p[1][2]*p[4][0]*p[5][1]*p[5][2] + p[0][0]*p[4][2]*p[5][1]*p[5][2] + 2*p[1][0]*p[4][2]*p[5][1]*p[5][2] + 2*pow(p[2][2], 2)*p[1][1]*p[6][0] - 2*pow(p[5][2], 2)*p[1][1]*p[6][0] - pow(p[1][2], 2)*p[2][1]*p[6][0] + pow(p[3][2], 2)*p[2][1]*p[6][0] - pow(p[5][2], 2)*p[2][1]*p[6][0] + pow(p[7][2], 2)*p[2][1]*p[6][0] + p[1][1]*p[1][2]*p[2][2]*p[6][0] - 2*p[1][2]*p[2][1]*p[2][2]*p[6][0] - 2*pow(p[2][2], 2)*p[3][1]*p[6][0] + 2*pow(p[7][2], 2)*p[3][1]*p[6][0] + 2*p[2][1]*p[2][2]*p[3][2]*p[6][0] - p[2][2]*p[3][1]*p[3][2]*p[6][0] + 2*pow(p[5][2], 2)*p[4][1]*p[6][0] - 2*pow(p[7][2], 2)*p[4][1]*p[6][0] + pow(p[1][2], 2)*p[5][1]*p[6][0] + pow(p[2][2], 2)*p[5][1]*p[6][0] - pow(p[4][2], 2)*p[5][1]*p[6][0] - pow(p[7][2], 2)*p[5][1]*p[6][0] + p[1][2]*p[2][2]*p[5][1]*p[6][0] - p[1][1]*p[1][2]*p[5][2]*p[6][0] - p[1][2]*p[2][1]*p[5][2]*p[6][0] - p[2][1]*p[2][2]*p[5][2]*p[6][0] + p[4][1]*p[4][2]*p[5][2]*p[6][0] + 2*p[1][2]*p[5][1]*p[5][2]*p[6][0] + p[2][2]*p[5][1]*p[5][2]*p[6][0] - 2*p[4][2]*p[5][1]*p[5][2]*p[6][0] - 2*pow(p[2][2], 2)*p[1][0]*p[6][1] + 2*pow(p[5][2], 2)*p[1][0]*p[6][1] + pow(p[1][2], 2)*p[2][0]*p[6][1] - pow(p[3][2], 2)*p[2][0]*p[6][1] + pow(p[5][2], 2)*p[2][0]*p[6][1] - pow(p[7][2], 2)*p[2][0]*p[6][1] - p[1][0]*p[1][2]*p[2][2]*p[6][1] + 2*p[1][2]*p[2][0]*p[2][2]*p[6][1] + 2*pow(p[2][2], 2)*p[3][0]*p[6][1] - 2*pow(p[7][2], 2)*p[3][0]*p[6][1] - 2*p[2][0]*p[2][2]*p[3][2]*p[6][1] + p[2][2]*p[3][0]*p[3][2]*p[6][1] - 2*pow(p[5][2], 2)*p[4][0]*p[6][1] + 2*pow(p[7][2], 2)*p[4][0]*p[6][1] - pow(p[1][2], 2)*p[5][0]*p[6][1] - pow(p[2][2], 2)*p[5][0]*p[6][1] + pow(p[4][2], 2)*p[5][0]*p[6][1] + pow(p[7][2], 2)*p[5][0]*p[6][1] - p[1][2]*p[2][2]*p[5][0]*p[6][1] + p[1][0]*p[1][2]*p[5][2]*p[6][1] + p[1][2]*p[2][0]*p[5][2]*p[6][1] + p[2][0]*p[2][2]*p[5][2]*p[6][1] - p[4][0]*p[4][2]*p[5][2]*p[6][1] - 2*p[1][2]*p[5][0]*p[5][2]*p[6][1] - p[2][2]*p[5][0]*p[5][2]*p[6][1] + 2*p[4][2]*p[5][0]*p[5][2]*p[6][1] - p[1][1]*p[1][2]*p[2][0]*p[6][2] + p[1][0]*p[1][2]*p[2][1]*p[6][2] - 2*p[1][1]*p[2][0]*p[2][2]*p[6][2] + 2*p[1][0]*p[2][1]*p[2][2]*p[6][2] - 2*p[2][1]*p[2][2]*p[3][0]*p[6][2] + 2*p[2][0]*p[2][2]*p[3][1]*p[6][2] - p[2][1]*p[3][0]*p[3][2]*p[6][2] + p[2][0]*p[3][1]*p[3][2]*p[6][2] + p[1][1]*p[1][2]*p[5][0]*p[6][2] + p[1][1]*p[2][2]*p[5][0]*p[6][2] + p[2][1]*p[2][2]*p[5][0]*p[6][2] - p[4][1]*p[4][2]*p[5][0]*p[6][2] - p[1][0]*p[1][2]*p[5][1]*p[6][2] - p[1][0]*p[2][2]*p[5][1]*p[6][2] - p[2][0]*p[2][2]*p[5][1]*p[6][2] + p[4][0]*p[4][2]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[5][2]*p[6][2] + p[1][0]*p[2][1]*p[5][2]*p[6][2] + 2*p[1][1]*p[5][0]*p[5][2]*p[6][2] + p[2][1]*p[5][0]*p[5][2]*p[6][2] - 2*p[4][1]*p[5][0]*p[5][2]*p[6][2] - 2*p[1][0]*p[5][1]*p[5][2]*p[6][2] - p[2][0]*p[5][1]*p[5][2]*p[6][2] + 2*p[4][0]*p[5][1]*p[5][2]*p[6][2] - p[1][2]*p[2][1]*p[6][0]*p[6][2] + p[1][1]*p[2][2]*p[6][0]*p[6][2] - p[2][2]*p[3][1]*p[6][0]*p[6][2] + p[2][1]*p[3][2]*p[6][0]*p[6][2] + p[1][2]*p[5][1]*p[6][0]*p[6][2] + 2*p[2][2]*p[5][1]*p[6][0]*p[6][2] - p[4][2]*p[5][1]*p[6][0]*p[6][2] - p[1][1]*p[5][2]*p[6][0]*p[6][2] - 2*p[2][1]*p[5][2]*p[6][0]*p[6][2] + p[4][1]*p[5][2]*p[6][0]*p[6][2] + p[1][2]*p[2][0]*p[6][1]*p[6][2] - p[1][0]*p[2][2]*p[6][1]*p[6][2] + p[2][2]*p[3][0]*p[6][1]*p[6][2] - p[2][0]*p[3][2]*p[6][1]*p[6][2] - p[1][2]*p[5][0]*p[6][1]*p[6][2] - 2*p[2][2]*p[5][0]*p[6][1]*p[6][2] + p[4][2]*p[5][0]*p[6][1]*p[6][2] + p[1][0]*p[5][2]*p[6][1]*p[6][2] + 2*p[2][0]*p[5][2]*p[6][1]*p[6][2] - p[4][0]*p[5][2]*p[6][1]*p[6][2] + 2*pow(p[3][2], 2)*p[2][1]*p[7][0] - 2*pow(p[6][2], 2)*p[2][1]*p[7][0] - pow(p[2][2], 2)*p[3][1]*p[7][0] + pow(p[4][2], 2)*p[3][1]*p[7][0] - pow(p[6][2], 2)*p[3][1]*p[7][0] + p[2][1]*p[2][2]*p[3][2]*p[7][0] - 2*p[2][2]*p[3][1]*p[3][2]*p[7][0] - pow(p[3][2], 2)*p[4][1]*p[7][0] + pow(p[5][2], 2)*p[4][1]*p[7][0] + pow(p[6][2], 2)*p[4][1]*p[7][0] + p[3][1]*p[3][2]*p[4][2]*p[7][0] - p[3][2]*p[4][1]*p[4][2]*p[7][0] - 2*pow(p[4][2], 2)*p[5][1]*p[7][0] + 2*pow(p[6][2], 2)*p[5][1]*p[7][0] + 2*p[4][1]*p[4][2]*p[5][2]*p[7][0] - p[4][2]*p[5][1]*p[5][2]*p[7][0] + pow(p[2][2], 2)*p[6][1]*p[7][0] + pow(p[3][2], 2)*p[6][1]*p[7][0] - pow(p[4][2], 2)*p[6][1]*p[7][0] - pow(p[5][2], 2)*p[6][1]*p[7][0] + p[2][2]*p[3][2]*p[6][1]*p[7][0] - p[4][2]*p[5][2]*p[6][1]*p[7][0] - p[2][1]*p[2][2]*p[6][2]*p[7][0] - p[2][2]*p[3][1]*p[6][2]*p[7][0] - p[3][1]*p[3][2]*p[6][2]*p[7][0] + p[4][1]*p[4][2]*p[6][2]*p[7][0] + p[4][1]*p[5][2]*p[6][2]*p[7][0] + p[5][1]*p[5][2]*p[6][2]*p[7][0] + 2*p[2][2]*p[6][1]*p[6][2]*p[7][0] + p[3][2]*p[6][1]*p[6][2]*p[7][0] - p[4][2]*p[6][1]*p[6][2]*p[7][0] - 2*p[5][2]*p[6][1]*p[6][2]*p[7][0] + 2*pow(p[3][2], 2)*p[0][0]*p[7][1] - 2*pow(p[4][2], 2)*p[0][0]*p[7][1] - 2*pow(p[3][2], 2)*p[2][0]*p[7][1] + 2*pow(p[6][2], 2)*p[2][0]*p[7][1] + pow(p[2][2], 2)*p[3][0]*p[7][1] - pow(p[4][2], 2)*p[3][0]*p[7][1] + pow(p[6][2], 2)*p[3][0]*p[7][1] - p[2][0]*p[2][2]*p[3][2]*p[7][1] + 2*p[2][2]*p[3][0]*p[3][2]*p[7][1] + pow(p[3][2], 2)*p[4][0]*p[7][1] - pow(p[5][2], 2)*p[4][0]*p[7][1] - pow(p[6][2], 2)*p[4][0]*p[7][1] - p[3][0]*p[3][2]*p[4][2]*p[7][1] + p[3][2]*p[4][0]*p[4][2]*p[7][1] + 2*pow(p[4][2], 2)*p[5][0]*p[7][1] - 2*pow(p[6][2], 2)*p[5][0]*p[7][1] - 2*p[4][0]*p[4][2]*p[5][2]*p[7][1] + p[4][2]*p[5][0]*p[5][2]*p[7][1] - pow(p[2][2], 2)*p[6][0]*p[7][1] - pow(p[3][2], 2)*p[6][0]*p[7][1] + pow(p[4][2], 2)*p[6][0]*p[7][1] + pow(p[5][2], 2)*p[6][0]*p[7][1] - p[2][2]*p[3][2]*p[6][0]*p[7][1] + p[4][2]*p[5][2]*p[6][0]*p[7][1] + p[2][0]*p[2][2]*p[6][2]*p[7][1] + p[2][2]*p[3][0]*p[6][2]*p[7][1] + p[3][0]*p[3][2]*p[6][2]*p[7][1] - p[4][0]*p[4][2]*p[6][2]*p[7][1] - p[4][0]*p[5][2]*p[6][2]*p[7][1] - p[5][0]*p[5][2]*p[6][2]*p[7][1] - 2*p[2][2]*p[6][0]*p[6][2]*p[7][1] - p[3][2]*p[6][0]*p[6][2]*p[7][1] + p[4][2]*p[6][0]*p[6][2]*p[7][1] + 2*p[5][2]*p[6][0]*p[6][2]*p[7][1] + pow(p[0][2], 2)*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + 2*p[3][1]*p[4][0] - 2*p[3][0]*p[4][1] + p[1][1]*(p[2][0] + 2*p[3][0] - 2*p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - 2*p[3][1] + 2*p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[2][1]*p[2][2]*p[3][0]*p[7][2] + p[2][0]*p[2][2]*p[3][1]*p[7][2] - 2*p[2][1]*p[3][0]*p[3][2]*p[7][2] - 2*p[0][0]*p[3][1]*p[3][2]*p[7][2] + 2*p[2][0]*p[3][1]*p[3][2]*p[7][2] - p[3][1]*p[3][2]*p[4][0]*p[7][2] + p[0][0]*p[3][2]*p[4][1]*p[7][2] + p[3][0]*p[3][2]*p[4][1]*p[7][2] - p[0][0]*p[3][1]*p[4][2]*p[7][2] - p[3][1]*p[4][0]*p[4][2]*p[7][2] + 2*p[0][0]*p[4][1]*p[4][2]*p[7][2] + p[3][0]*p[4][1]*p[4][2]*p[7][2] - 2*p[4][1]*p[4][2]*p[5][0]*p[7][2] + 2*p[4][0]*p[4][2]*p[5][1]*p[7][2] - p[4][1]*p[5][0]*p[5][2]*p[7][2] + p[4][0]*p[5][1]*p[5][2]*p[7][2] + p[2][1]*p[2][2]*p[6][0]*p[7][2] + p[2][1]*p[3][2]*p[6][0]*p[7][2] + p[3][1]*p[3][2]*p[6][0]*p[7][2] - p[4][1]*p[4][2]*p[6][0]*p[7][2] - p[4][2]*p[5][1]*p[6][0]*p[7][2] - p[5][1]*p[5][2]*p[6][0]*p[7][2] - p[2][0]*p[2][2]*p[6][1]*p[7][2] - p[2][0]*p[3][2]*p[6][1]*p[7][2] - p[3][0]*p[3][2]*p[6][1]*p[7][2] + p[4][0]*p[4][2]*p[6][1]*p[7][2] + p[4][2]*p[5][0]*p[6][1]*p[7][2] + p[5][0]*p[5][2]*p[6][1]*p[7][2] - p[2][1]*p[3][0]*p[6][2]*p[7][2] + p[2][0]*p[3][1]*p[6][2]*p[7][2] - p[4][1]*p[5][0]*p[6][2]*p[7][2] + p[4][0]*p[5][1]*p[6][2]*p[7][2] + 2*p[2][1]*p[6][0]*p[6][2]*p[7][2] + p[3][1]*p[6][0]*p[6][2]*p[7][2] - p[4][1]*p[6][0]*p[6][2]*p[7][2] - 2*p[5][1]*p[6][0]*p[6][2]*p[7][2] - 2*p[2][0]*p[6][1]*p[6][2]*p[7][2] - p[3][0]*p[6][1]*p[6][2]*p[7][2] + p[4][0]*p[6][1]*p[6][2]*p[7][2] + 2*p[5][0]*p[6][1]*p[6][2]*p[7][2] - p[2][2]*p[3][1]*p[7][0]*p[7][2] + p[2][1]*p[3][2]*p[7][0]*p[7][2] - 2*p[3][2]*p[4][1]*p[7][0]*p[7][2] + 2*p[3][1]*p[4][2]*p[7][0]*p[7][2] - p[4][2]*p[5][1]*p[7][0]*p[7][2] + p[4][1]*p[5][2]*p[7][0]*p[7][2] + p[2][2]*p[6][1]*p[7][0]*p[7][2] + 2*p[3][2]*p[6][1]*p[7][0]*p[7][2] - 2*p[4][2]*p[6][1]*p[7][0]*p[7][2] - p[5][2]*p[6][1]*p[7][0]*p[7][2] - p[2][1]*p[6][2]*p[7][0]*p[7][2] - 2*p[3][1]*p[6][2]*p[7][0]*p[7][2] + 2*p[4][1]*p[6][2]*p[7][0]*p[7][2] + p[5][1]*p[6][2]*p[7][0]*p[7][2] + p[2][2]*p[3][0]*p[7][1]*p[7][2] + p[0][0]*p[3][2]*p[7][1]*p[7][2] - p[2][0]*p[3][2]*p[7][1]*p[7][2] + 2*p[3][2]*p[4][0]*p[7][1]*p[7][2] - p[0][0]*p[4][2]*p[7][1]*p[7][2] - 2*p[3][0]*p[4][2]*p[7][1]*p[7][2] + p[4][2]*p[5][0]*p[7][1]*p[7][2] - p[4][0]*p[5][2]*p[7][1]*p[7][2] - p[2][2]*p[6][0]*p[7][1]*p[7][2] - 2*p[3][2]*p[6][0]*p[7][1]*p[7][2] + 2*p[4][2]*p[6][0]*p[7][1]*p[7][2] + p[5][2]*p[6][0]*p[7][1]*p[7][2] + p[2][0]*p[6][2]*p[7][1]*p[7][2] + 2*p[3][0]*p[6][2]*p[7][1]*p[7][2] - 2*p[4][0]*p[6][2]*p[7][1]*p[7][2] - p[5][0]*p[6][2]*p[7][1]*p[7][2] + p[0][1]*(2*pow(p[3][2], 2)*p[2][0] - pow(p[2][2], 2)*p[3][0] + pow(p[4][2], 2)*p[3][0] + pow(p[7][2], 2)*p[3][0] + p[2][0]*p[2][2]*p[3][2] - 2*p[2][2]*p[3][0]*p[3][2] - pow(p[3][2], 2)*p[4][0] + pow(p[5][2], 2)*p[4][0] - pow(p[7][2], 2)*p[4][0] + p[3][0]*p[3][2]*p[4][2] - p[3][2]*p[4][0]*p[4][2] - 2*pow(p[4][2], 2)*p[5][0] + pow(p[1][2], 2)*(-2*p[2][0] - p[3][0] + p[4][0] + 2*p[5][0]) + 2*p[4][0]*p[4][2]*p[5][2] - p[4][2]*p[5][0]*p[5][2] + p[1][0]*(pow(p[2][2], 2) + pow(p[3][2], 2) - pow(p[4][2], 2) - pow(p[5][2], 2) + p[2][2]*p[3][2] - p[4][2]*p[5][2]) + p[1][2]*(-(p[2][0]*p[2][2]) - p[2][2]*p[3][0] - p[3][0]*p[3][2] + p[4][0]*p[4][2] + p[1][0]*(2*p[2][2] + p[3][2] - p[4][2] - 2*p[5][2]) + p[4][0]*p[5][2] + p[5][0]*p[5][2]) - 2*pow(p[3][2], 2)*p[7][0] + 2*pow(p[4][2], 2)*p[7][0] + 2*p[3][0]*p[3][2]*p[7][2] - p[3][2]*p[4][0]*p[7][2] + p[3][0]*p[4][2]*p[7][2] - 2*p[4][0]*p[4][2]*p[7][2] - p[3][2]*p[7][0]*p[7][2] + p[4][2]*p[7][0]*p[7][2]) + p[0][2]*(p[0][0]*p[1][2]*p[2][1] - 2*p[1][0]*p[1][2]*p[2][1] - p[1][0]*p[2][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] + p[2][1]*p[2][2]*p[3][0] + 2*p[0][0]*p[1][2]*p[3][1] - p[1][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[2][0]*p[2][2]*p[3][1] - p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] + 2*p[2][1]*p[3][0]*p[3][2] - p[1][0]*p[3][1]*p[3][2] - 2*p[2][0]*p[3][1]*p[3][2] + p[3][1]*p[3][2]*p[4][0] - 2*p[0][0]*p[1][2]*p[4][1] + p[1][0]*p[1][2]*p[4][1] + 2*p[0][0]*p[3][2]*p[4][1] - p[3][0]*p[3][2]*p[4][1] - 2*p[0][0]*p[3][1]*p[4][2] + p[3][1]*p[4][0]*p[4][2] + p[1][0]*p[4][1]*p[4][2] - p[3][0]*p[4][1]*p[4][2] + p[1][2]*p[4][1]*p[5][0] + 2*p[4][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + 2*p[1][0]*p[1][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][0]*p[4][2]*p[5][1] - 2*p[4][0]*p[4][2]*p[5][1] - p[0][0]*p[4][1]*p[5][2] + p[4][1]*p[5][0]*p[5][2] + p[1][0]*p[5][1]*p[5][2] - p[4][0]*p[5][1]*p[5][2] + p[1][1]*(p[2][0]*p[2][2] + p[2][0]*p[3][2] + p[3][0]*p[3][2] - p[4][0]*p[4][2] + p[1][2]*(2*p[2][0] + p[3][0] - p[4][0] - 2*p[5][0]) - p[4][2]*p[5][0] - p[5][0]*p[5][2] + p[0][0]*(-p[2][2] - 2*p[3][2] + 2*p[4][2] + p[5][2])) + 2*p[3][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] - 2*p[4][1]*p[4][2]*p[7][0] + p[0][0]*p[3][2]*p[7][1] - 2*p[3][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][0]*p[4][2]*p[7][1] + 2*p[4][0]*p[4][2]*p[7][1] - p[0][0]*p[3][1]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][1]*p[7][0]*p[7][2] - p[4][1]*p[7][0]*p[7][2] - p[3][0]*p[7][1]*p[7][2] + p[4][0]*p[7][1]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - 2*p[3][2]*p[4][0] + 2*p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - 2*p[3][0] + 2*p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + 2*p[3][2] - 2*p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2])))/72;
    return result;
  }
};

// 3D Wedge

template<>
class Operator<VolumeIntegral, basis::Unitary, Wedge>:
  public OperatorBase<VolumeIntegral, basis::Unitary, Interval>
{
public:
  static constexpr size_t dim = dimension(Hexahedron);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Unitary, dim>::function_size;
  static constexpr size_t point_size = 6;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (2*p[0][0]*p[1][2]*p[2][1] - 2*p[0][0]*p[1][1]*p[2][2] - p[0][0]*p[1][2]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[0][0]*p[1][1]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][2]*p[2][1]*p[4][0] + p[1][1]*p[2][2]*p[4][0] + p[1][2]*p[3][1]*p[4][0] - p[1][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + p[1][2]*p[2][0]*p[4][1] - p[1][0]*p[2][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[1][1]*p[2][0]*p[4][2] + p[1][0]*p[2][1]*p[4][2] + p[1][1]*p[3][0]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] + p[2][2]*p[4][1]*p[5][0] - 2*p[3][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[2][1]*p[4][2]*p[5][0] + 2*p[3][1]*p[4][2]*p[5][0] + p[1][2]*p[2][0]*p[5][1] + p[0][0]*p[2][2]*p[5][1] - p[1][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][0]*p[3][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] - p[2][2]*p[4][0]*p[5][1] + 2*p[3][2]*p[4][0]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[2][0]*p[4][2]*p[5][1] - 2*p[3][0]*p[4][2]*p[5][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][1]*(2*p[2][0] - p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-2*p[2][1] + p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[1][1]*p[2][0]*p[5][2] - p[0][0]*p[2][1]*p[5][2] + p[1][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][0]*p[3][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] + p[2][1]*p[4][0]*p[5][2] - 2*p[3][1]*p[4][0]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[2][0]*p[4][1]*p[5][2] + 2*p[3][0]*p[4][1]*p[5][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-2*p[2][0] + p[3][0] + p[4][0]) + p[1][0]*(2*p[2][2] - p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2]))/12;
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Linear, Wedge>:
  public OperatorBase<VolumeIntegral, basis::Linear, Interval>
{
public:
  static constexpr size_t dim = dimension(Wedge);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Linear, dim>::function_size;
  static constexpr size_t point_size = 6;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (2*p[0][0]*p[1][2]*p[2][1] - 2*p[0][0]*p[1][1]*p[2][2] - p[0][0]*p[1][2]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[0][0]*p[1][1]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][2]*p[2][1]*p[4][0] + p[1][1]*p[2][2]*p[4][0] + p[1][2]*p[3][1]*p[4][0] - p[1][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + p[1][2]*p[2][0]*p[4][1] - p[1][0]*p[2][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[1][1]*p[2][0]*p[4][2] + p[1][0]*p[2][1]*p[4][2] + p[1][1]*p[3][0]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] + p[2][2]*p[4][1]*p[5][0] - 2*p[3][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[2][1]*p[4][2]*p[5][0] + 2*p[3][1]*p[4][2]*p[5][0] + p[1][2]*p[2][0]*p[5][1] + p[0][0]*p[2][2]*p[5][1] - p[1][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][0]*p[3][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] - p[2][2]*p[4][0]*p[5][1] + 2*p[3][2]*p[4][0]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[2][0]*p[4][2]*p[5][1] - 2*p[3][0]*p[4][2]*p[5][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][1]*(2*p[2][0] - p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-2*p[2][1] + p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[1][1]*p[2][0]*p[5][2] - p[0][0]*p[2][1]*p[5][2] + p[1][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][0]*p[3][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] + p[2][1]*p[4][0]*p[5][2] - 2*p[3][1]*p[4][0]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[2][0]*p[4][1]*p[5][2] + 2*p[3][0]*p[4][1]*p[5][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-2*p[2][0] + p[3][0] + p[4][0]) + p[1][0]*(2*p[2][2] - p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2]))/12;
    result[1][0] = (-3*pow(p[2][0], 2)*p[0][1]*p[1][2] + pow(p[3][0], 2)*p[0][1]*p[1][2] + pow(p[4][0], 2)*p[0][1]*p[1][2] - 3*p[0][1]*p[1][0]*p[1][2]*p[2][0] - pow(p[4][0], 2)*p[1][2]*p[2][1] - pow(p[5][0], 2)*p[1][2]*p[2][1] + 3*pow(p[1][0], 2)*p[0][1]*p[2][2] - pow(p[3][0], 2)*p[0][1]*p[2][2] - pow(p[5][0], 2)*p[0][1]*p[2][2] + pow(p[4][0], 2)*p[1][1]*p[2][2] + pow(p[5][0], 2)*p[1][1]*p[2][2] + 3*p[0][1]*p[1][0]*p[2][0]*p[2][2] + p[0][1]*p[1][0]*p[1][2]*p[3][0] - p[0][1]*p[2][0]*p[2][2]*p[3][0] + 2*pow(p[4][0], 2)*p[1][2]*p[3][1] - 2*pow(p[5][0], 2)*p[2][2]*p[3][1] - pow(p[1][0], 2)*p[0][1]*p[3][2] + pow(p[2][0], 2)*p[0][1]*p[3][2] - pow(p[4][0], 2)*p[0][1]*p[3][2] + pow(p[5][0], 2)*p[0][1]*p[3][2] - 2*pow(p[4][0], 2)*p[1][1]*p[3][2] + 2*pow(p[5][0], 2)*p[2][1]*p[3][2] - p[0][1]*p[1][0]*p[3][0]*p[3][2] + p[0][1]*p[2][0]*p[3][0]*p[3][2] + 2*p[0][1]*p[1][0]*p[1][2]*p[4][0] - 2*p[1][0]*p[1][2]*p[2][1]*p[4][0] - p[1][2]*p[2][0]*p[2][1]*p[4][0] + 2*p[1][0]*p[1][1]*p[2][2]*p[4][0] + p[1][1]*p[2][0]*p[2][2]*p[4][0] + p[0][1]*p[1][2]*p[3][0]*p[4][0] + p[1][0]*p[1][2]*p[3][1]*p[4][0] + p[1][2]*p[3][0]*p[3][1]*p[4][0] - p[0][1]*p[1][0]*p[3][2]*p[4][0] - p[1][0]*p[1][1]*p[3][2]*p[4][0] - 2*p[0][1]*p[3][0]*p[3][2]*p[4][0] - p[1][1]*p[3][0]*p[3][2]*p[4][0] + pow(p[2][0], 2)*p[1][2]*p[4][1] - pow(p[3][0], 2)*p[1][2]*p[4][1] + pow(p[5][0], 2)*p[1][2]*p[4][1] + 2*p[1][0]*p[1][2]*p[2][0]*p[4][1] - 2*pow(p[1][0], 2)*p[2][2]*p[4][1] + 2*pow(p[5][0], 2)*p[2][2]*p[4][1] - p[1][0]*p[2][0]*p[2][2]*p[4][1] - p[1][0]*p[1][2]*p[3][0]*p[4][1] + pow(p[1][0], 2)*p[3][2]*p[4][1] - 3*pow(p[5][0], 2)*p[3][2]*p[4][1] + p[1][0]*p[3][0]*p[3][2]*p[4][1] + p[1][2]*p[2][0]*p[4][0]*p[4][1] - p[1][0]*p[2][2]*p[4][0]*p[4][1] - 2*p[1][2]*p[3][0]*p[4][0]*p[4][1] + 2*p[1][0]*p[3][2]*p[4][0]*p[4][1] - 2*pow(p[1][0], 2)*p[0][1]*p[4][2] + 2*pow(p[3][0], 2)*p[0][1]*p[4][2] - pow(p[2][0], 2)*p[1][1]*p[4][2] + pow(p[3][0], 2)*p[1][1]*p[4][2] - pow(p[5][0], 2)*p[1][1]*p[4][2] - 2*p[1][0]*p[1][1]*p[2][0]*p[4][2] + 2*pow(p[1][0], 2)*p[2][1]*p[4][2] - 2*pow(p[5][0], 2)*p[2][1]*p[4][2] + p[1][0]*p[2][0]*p[2][1]*p[4][2] + p[1][0]*p[1][1]*p[3][0]*p[4][2] - pow(p[1][0], 2)*p[3][1]*p[4][2] + 3*pow(p[5][0], 2)*p[3][1]*p[4][2] - p[1][0]*p[3][0]*p[3][1]*p[4][2] - p[0][1]*p[1][0]*p[4][0]*p[4][2] - p[1][1]*p[2][0]*p[4][0]*p[4][2] + p[1][0]*p[2][1]*p[4][0]*p[4][2] + p[0][1]*p[3][0]*p[4][0]*p[4][2] + 2*p[1][1]*p[3][0]*p[4][0]*p[4][2] - 2*p[1][0]*p[3][1]*p[4][0]*p[4][2] - p[1][0]*p[1][2]*p[2][1]*p[5][0] - 2*p[1][2]*p[2][0]*p[2][1]*p[5][0] + p[1][0]*p[1][1]*p[2][2]*p[5][0] - 2*p[0][1]*p[2][0]*p[2][2]*p[5][0] + 2*p[1][1]*p[2][0]*p[2][2]*p[5][0] - p[0][1]*p[2][2]*p[3][0]*p[5][0] - p[2][0]*p[2][2]*p[3][1]*p[5][0] - p[2][2]*p[3][0]*p[3][1]*p[5][0] + p[0][1]*p[2][0]*p[3][2]*p[5][0] + p[2][0]*p[2][1]*p[3][2]*p[5][0] + 2*p[0][1]*p[3][0]*p[3][2]*p[5][0] + p[2][1]*p[3][0]*p[3][2]*p[5][0] - p[1][2]*p[2][1]*p[4][0]*p[5][0] + p[1][1]*p[2][2]*p[4][0]*p[5][0] + p[1][0]*p[1][2]*p[4][1]*p[5][0] + p[1][2]*p[2][0]*p[4][1]*p[5][0] + p[2][0]*p[2][2]*p[4][1]*p[5][0] - 3*p[3][0]*p[3][2]*p[4][1]*p[5][0] + 2*p[1][2]*p[4][0]*p[4][1]*p[5][0] + p[2][2]*p[4][0]*p[4][1]*p[5][0] - 3*p[3][2]*p[4][0]*p[4][1]*p[5][0] - p[1][0]*p[1][1]*p[4][2]*p[5][0] - p[1][1]*p[2][0]*p[4][2]*p[5][0] - p[2][0]*p[2][1]*p[4][2]*p[5][0] + 3*p[3][0]*p[3][1]*p[4][2]*p[5][0] - 2*p[1][1]*p[4][0]*p[4][2]*p[5][0] - p[2][1]*p[4][0]*p[4][2]*p[5][0] + 3*p[3][1]*p[4][0]*p[4][2]*p[5][0] + 2*pow(p[2][0], 2)*p[1][2]*p[5][1] - 2*pow(p[4][0], 2)*p[1][2]*p[5][1] + p[1][0]*p[1][2]*p[2][0]*p[5][1] - pow(p[1][0], 2)*p[2][2]*p[5][1] + pow(p[3][0], 2)*p[2][2]*p[5][1] - pow(p[4][0], 2)*p[2][2]*p[5][1] - 2*p[1][0]*p[2][0]*p[2][2]*p[5][1] + p[2][0]*p[2][2]*p[3][0]*p[5][1] - pow(p[2][0], 2)*p[3][2]*p[5][1] + 3*pow(p[4][0], 2)*p[3][2]*p[5][1] - p[2][0]*p[3][0]*p[3][2]*p[5][1] - p[1][0]*p[1][2]*p[4][0]*p[5][1] - p[1][0]*p[2][2]*p[4][0]*p[5][1] - p[2][0]*p[2][2]*p[4][0]*p[5][1] + 3*p[3][0]*p[3][2]*p[4][0]*p[5][1] + pow(p[1][0], 2)*p[4][2]*p[5][1] + pow(p[2][0], 2)*p[4][2]*p[5][1] - 3*pow(p[3][0], 2)*p[4][2]*p[5][1] + p[1][0]*p[2][0]*p[4][2]*p[5][1] + 2*p[1][0]*p[4][0]*p[4][2]*p[5][1] + p[2][0]*p[4][0]*p[4][2]*p[5][1] - 3*p[3][0]*p[4][0]*p[4][2]*p[5][1] + p[1][2]*p[2][0]*p[5][0]*p[5][1] - p[1][0]*p[2][2]*p[5][0]*p[5][1] + 2*p[2][2]*p[3][0]*p[5][0]*p[5][1] - 2*p[2][0]*p[3][2]*p[5][0]*p[5][1] - p[1][2]*p[4][0]*p[5][0]*p[5][1] - 2*p[2][2]*p[4][0]*p[5][0]*p[5][1] + 3*p[3][2]*p[4][0]*p[5][0]*p[5][1] + p[1][0]*p[4][2]*p[5][0]*p[5][1] + 2*p[2][0]*p[4][2]*p[5][0]*p[5][1] - 3*p[3][0]*p[4][2]*p[5][0]*p[5][1] + p[0][2]*(pow(p[3][0], 2)*p[2][1] + pow(p[5][0], 2)*p[2][1] + p[2][0]*p[2][1]*p[3][0] - pow(p[2][0], 2)*p[3][1] + pow(p[4][0], 2)*p[3][1] - pow(p[5][0], 2)*p[3][1] - p[2][0]*p[3][0]*p[3][1] + 2*p[3][0]*p[3][1]*p[4][0] + p[1][1]*(3*pow(p[2][0], 2) - pow(p[3][0], 2) - pow(p[4][0], 2) - p[3][0]*p[4][0]) - 2*pow(p[3][0], 2)*p[4][1] - p[3][0]*p[4][0]*p[4][1] + pow(p[1][0], 2)*(-3*p[2][1] + p[3][1] + 2*p[4][1]) + p[1][0]*(-3*p[2][0]*p[2][1] + p[3][0]*p[3][1] + p[1][1]*(3*p[2][0] - p[3][0] - 2*p[4][0]) + p[3][1]*p[4][0] + p[4][0]*p[4][1]) + 2*p[2][0]*p[2][1]*p[5][0] + p[2][1]*p[3][0]*p[5][0] - p[2][0]*p[3][1]*p[5][0] - 2*p[3][0]*p[3][1]*p[5][0] - 2*pow(p[2][0], 2)*p[5][1] + 2*pow(p[3][0], 2)*p[5][1] - p[2][0]*p[5][0]*p[5][1] + p[3][0]*p[5][0]*p[5][1]) + 2*pow(p[2][0], 2)*p[0][1]*p[5][2] - 2*pow(p[3][0], 2)*p[0][1]*p[5][2] - 2*pow(p[2][0], 2)*p[1][1]*p[5][2] + 2*pow(p[4][0], 2)*p[1][1]*p[5][2] - p[1][0]*p[1][1]*p[2][0]*p[5][2] + pow(p[1][0], 2)*p[2][1]*p[5][2] - pow(p[3][0], 2)*p[2][1]*p[5][2] + pow(p[4][0], 2)*p[2][1]*p[5][2] + 2*p[1][0]*p[2][0]*p[2][1]*p[5][2] - p[2][0]*p[2][1]*p[3][0]*p[5][2] + pow(p[2][0], 2)*p[3][1]*p[5][2] - 3*pow(p[4][0], 2)*p[3][1]*p[5][2] + p[2][0]*p[3][0]*p[3][1]*p[5][2] + p[1][0]*p[1][1]*p[4][0]*p[5][2] + p[1][0]*p[2][1]*p[4][0]*p[5][2] + p[2][0]*p[2][1]*p[4][0]*p[5][2] - 3*p[3][0]*p[3][1]*p[4][0]*p[5][2] - pow(p[1][0], 2)*p[4][1]*p[5][2] - pow(p[2][0], 2)*p[4][1]*p[5][2] + 3*pow(p[3][0], 2)*p[4][1]*p[5][2] - p[1][0]*p[2][0]*p[4][1]*p[5][2] - 2*p[1][0]*p[4][0]*p[4][1]*p[5][2] - p[2][0]*p[4][0]*p[4][1]*p[5][2] + 3*p[3][0]*p[4][0]*p[4][1]*p[5][2] + p[0][1]*p[2][0]*p[5][0]*p[5][2] - p[1][1]*p[2][0]*p[5][0]*p[5][2] + p[1][0]*p[2][1]*p[5][0]*p[5][2] - p[0][1]*p[3][0]*p[5][0]*p[5][2] - 2*p[2][1]*p[3][0]*p[5][0]*p[5][2] + 2*p[2][0]*p[3][1]*p[5][0]*p[5][2] + p[1][1]*p[4][0]*p[5][0]*p[5][2] + 2*p[2][1]*p[4][0]*p[5][0]*p[5][2] - 3*p[3][1]*p[4][0]*p[5][0]*p[5][2] - p[1][0]*p[4][1]*p[5][0]*p[5][2] - 2*p[2][0]*p[4][1]*p[5][0]*p[5][2] + 3*p[3][0]*p[4][1]*p[5][0]*p[5][2] + pow(p[0][0], 2)*(2*p[2][2]*p[3][1] - 2*p[2][1]*p[3][2] + p[1][2]*(3*p[2][1] - 2*p[3][1] - p[4][1]) + p[3][2]*p[4][1] - p[3][1]*p[4][2] + p[1][1]*(-3*p[2][2] + 2*p[3][2] + p[4][2]) + p[2][2]*p[5][1] - p[3][2]*p[5][1] - p[2][1]*p[5][2] + p[3][1]*p[5][2]) + p[0][0]*(3*p[1][0]*p[1][2]*p[2][1] + 3*p[1][2]*p[2][0]*p[2][1] - 3*p[1][0]*p[1][1]*p[2][2] - 3*p[1][1]*p[2][0]*p[2][2] - p[1][0]*p[1][2]*p[3][1] + p[2][0]*p[2][2]*p[3][1] - p[1][2]*p[3][0]*p[3][1] + p[2][2]*p[3][0]*p[3][1] + p[1][0]*p[1][1]*p[3][2] - p[2][0]*p[2][1]*p[3][2] + p[1][1]*p[3][0]*p[3][2] - p[2][1]*p[3][0]*p[3][2] - 2*p[1][0]*p[1][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[1][0]*p[3][2]*p[4][1] + 2*p[3][0]*p[3][2]*p[4][1] - p[1][2]*p[4][0]*p[4][1] + p[3][2]*p[4][0]*p[4][1] + 2*p[1][0]*p[1][1]*p[4][2] + p[1][1]*p[3][0]*p[4][2] - p[1][0]*p[3][1]*p[4][2] - 2*p[3][0]*p[3][1]*p[4][2] + p[1][1]*p[4][0]*p[4][2] - p[3][1]*p[4][0]*p[4][2] + 2*p[2][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[2][0]*p[3][2]*p[5][1] - 2*p[3][0]*p[3][2]*p[5][1] + p[2][2]*p[5][0]*p[5][1] - p[3][2]*p[5][0]*p[5][1] + p[0][2]*(2*p[2][1]*p[3][0] - 2*p[2][0]*p[3][1] + p[1][1]*(3*p[2][0] - 2*p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-3*p[2][1] + 2*p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - 2*p[2][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[2][0]*p[3][1]*p[5][2] + 2*p[3][0]*p[3][1]*p[5][2] - p[2][1]*p[5][0]*p[5][2] + p[3][1]*p[5][0]*p[5][2] + p[0][1]*(-2*p[2][2]*p[3][0] + 2*p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + 2*p[3][0] + p[4][0]) + p[1][0]*(3*p[2][2] - 2*p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2])))/72;
    result[2][0] = (3*pow(p[2][1], 2)*p[0][0]*p[1][2] - pow(p[3][1], 2)*p[0][0]*p[1][2] - pow(p[4][1], 2)*p[0][0]*p[1][2] + pow(p[4][1], 2)*p[1][2]*p[2][0] + pow(p[5][1], 2)*p[1][2]*p[2][0] + 3*p[0][0]*p[1][1]*p[1][2]*p[2][1] - 3*pow(p[1][1], 2)*p[0][0]*p[2][2] + pow(p[3][1], 2)*p[0][0]*p[2][2] + pow(p[5][1], 2)*p[0][0]*p[2][2] - pow(p[4][1], 2)*p[1][0]*p[2][2] - pow(p[5][1], 2)*p[1][0]*p[2][2] - 3*p[0][0]*p[1][1]*p[2][1]*p[2][2] - 2*pow(p[4][1], 2)*p[1][2]*p[3][0] + 2*pow(p[5][1], 2)*p[2][2]*p[3][0] - p[0][0]*p[1][1]*p[1][2]*p[3][1] + p[0][0]*p[2][1]*p[2][2]*p[3][1] + pow(p[1][1], 2)*p[0][0]*p[3][2] - pow(p[2][1], 2)*p[0][0]*p[3][2] + pow(p[4][1], 2)*p[0][0]*p[3][2] - pow(p[5][1], 2)*p[0][0]*p[3][2] + 2*pow(p[4][1], 2)*p[1][0]*p[3][2] - 2*pow(p[5][1], 2)*p[2][0]*p[3][2] + p[0][0]*p[1][1]*p[3][1]*p[3][2] - p[0][0]*p[2][1]*p[3][1]*p[3][2] - pow(p[2][1], 2)*p[1][2]*p[4][0] + pow(p[3][1], 2)*p[1][2]*p[4][0] - pow(p[5][1], 2)*p[1][2]*p[4][0] - 2*p[1][1]*p[1][2]*p[2][1]*p[4][0] + 2*pow(p[1][1], 2)*p[2][2]*p[4][0] - 2*pow(p[5][1], 2)*p[2][2]*p[4][0] + p[1][1]*p[2][1]*p[2][2]*p[4][0] + p[1][1]*p[1][2]*p[3][1]*p[4][0] - pow(p[1][1], 2)*p[3][2]*p[4][0] + 3*pow(p[5][1], 2)*p[3][2]*p[4][0] - p[1][1]*p[3][1]*p[3][2]*p[4][0] - 2*p[0][0]*p[1][1]*p[1][2]*p[4][1] + 2*p[1][1]*p[1][2]*p[2][0]*p[4][1] + p[1][2]*p[2][0]*p[2][1]*p[4][1] - 2*p[1][0]*p[1][1]*p[2][2]*p[4][1] - p[1][0]*p[2][1]*p[2][2]*p[4][1] - p[1][1]*p[1][2]*p[3][0]*p[4][1] - p[0][0]*p[1][2]*p[3][1]*p[4][1] - p[1][2]*p[3][0]*p[3][1]*p[4][1] + p[0][0]*p[1][1]*p[3][2]*p[4][1] + p[1][0]*p[1][1]*p[3][2]*p[4][1] + 2*p[0][0]*p[3][1]*p[3][2]*p[4][1] + p[1][0]*p[3][1]*p[3][2]*p[4][1] - p[1][2]*p[2][1]*p[4][0]*p[4][1] + p[1][1]*p[2][2]*p[4][0]*p[4][1] + 2*p[1][2]*p[3][1]*p[4][0]*p[4][1] - 2*p[1][1]*p[3][2]*p[4][0]*p[4][1] + 2*pow(p[1][1], 2)*p[0][0]*p[4][2] - 2*pow(p[3][1], 2)*p[0][0]*p[4][2] + pow(p[2][1], 2)*p[1][0]*p[4][2] - pow(p[3][1], 2)*p[1][0]*p[4][2] + pow(p[5][1], 2)*p[1][0]*p[4][2] - 2*pow(p[1][1], 2)*p[2][0]*p[4][2] + 2*pow(p[5][1], 2)*p[2][0]*p[4][2] + 2*p[1][0]*p[1][1]*p[2][1]*p[4][2] - p[1][1]*p[2][0]*p[2][1]*p[4][2] + pow(p[1][1], 2)*p[3][0]*p[4][2] - 3*pow(p[5][1], 2)*p[3][0]*p[4][2] - p[1][0]*p[1][1]*p[3][1]*p[4][2] + p[1][1]*p[3][0]*p[3][1]*p[4][2] + p[0][0]*p[1][1]*p[4][1]*p[4][2] - p[1][1]*p[2][0]*p[4][1]*p[4][2] + p[1][0]*p[2][1]*p[4][1]*p[4][2] + 2*p[1][1]*p[3][0]*p[4][1]*p[4][2] - p[0][0]*p[3][1]*p[4][1]*p[4][2] - 2*p[1][0]*p[3][1]*p[4][1]*p[4][2] - 2*pow(p[2][1], 2)*p[1][2]*p[5][0] + 2*pow(p[4][1], 2)*p[1][2]*p[5][0] - p[1][1]*p[1][2]*p[2][1]*p[5][0] + pow(p[1][1], 2)*p[2][2]*p[5][0] - pow(p[3][1], 2)*p[2][2]*p[5][0] + pow(p[4][1], 2)*p[2][2]*p[5][0] + 2*p[1][1]*p[2][1]*p[2][2]*p[5][0] - p[2][1]*p[2][2]*p[3][1]*p[5][0] + pow(p[2][1], 2)*p[3][2]*p[5][0] - 3*pow(p[4][1], 2)*p[3][2]*p[5][0] + p[2][1]*p[3][1]*p[3][2]*p[5][0] + p[1][1]*p[1][2]*p[4][1]*p[5][0] + p[1][1]*p[2][2]*p[4][1]*p[5][0] + p[2][1]*p[2][2]*p[4][1]*p[5][0] - 3*p[3][1]*p[3][2]*p[4][1]*p[5][0] - pow(p[1][1], 2)*p[4][2]*p[5][0] - pow(p[2][1], 2)*p[4][2]*p[5][0] + 3*pow(p[3][1], 2)*p[4][2]*p[5][0] - p[1][1]*p[2][1]*p[4][2]*p[5][0] - 2*p[1][1]*p[4][1]*p[4][2]*p[5][0] - p[2][1]*p[4][1]*p[4][2]*p[5][0] + 3*p[3][1]*p[4][1]*p[4][2]*p[5][0] + p[1][1]*p[1][2]*p[2][0]*p[5][1] + 2*p[1][2]*p[2][0]*p[2][1]*p[5][1] - p[1][0]*p[1][1]*p[2][2]*p[5][1] + 2*p[0][0]*p[2][1]*p[2][2]*p[5][1] - 2*p[1][0]*p[2][1]*p[2][2]*p[5][1] + p[2][1]*p[2][2]*p[3][0]*p[5][1] + p[0][0]*p[2][2]*p[3][1]*p[5][1] + p[2][2]*p[3][0]*p[3][1]*p[5][1] - p[0][0]*p[2][1]*p[3][2]*p[5][1] - p[2][0]*p[2][1]*p[3][2]*p[5][1] - 2*p[0][0]*p[3][1]*p[3][2]*p[5][1] - p[2][0]*p[3][1]*p[3][2]*p[5][1] - p[1][1]*p[1][2]*p[4][0]*p[5][1] - p[1][2]*p[2][1]*p[4][0]*p[5][1] - p[2][1]*p[2][2]*p[4][0]*p[5][1] + 3*p[3][1]*p[3][2]*p[4][0]*p[5][1] + p[1][2]*p[2][0]*p[4][1]*p[5][1] - p[1][0]*p[2][2]*p[4][1]*p[5][1] - 2*p[1][2]*p[4][0]*p[4][1]*p[5][1] - p[2][2]*p[4][0]*p[4][1]*p[5][1] + 3*p[3][2]*p[4][0]*p[4][1]*p[5][1] + p[1][0]*p[1][1]*p[4][2]*p[5][1] + p[1][0]*p[2][1]*p[4][2]*p[5][1] + p[2][0]*p[2][1]*p[4][2]*p[5][1] - 3*p[3][0]*p[3][1]*p[4][2]*p[5][1] + 2*p[1][0]*p[4][1]*p[4][2]*p[5][1] + p[2][0]*p[4][1]*p[4][2]*p[5][1] - 3*p[3][0]*p[4][1]*p[4][2]*p[5][1] - p[1][2]*p[2][1]*p[5][0]*p[5][1] + p[1][1]*p[2][2]*p[5][0]*p[5][1] - 2*p[2][2]*p[3][1]*p[5][0]*p[5][1] + 2*p[2][1]*p[3][2]*p[5][0]*p[5][1] + p[1][2]*p[4][1]*p[5][0]*p[5][1] + 2*p[2][2]*p[4][1]*p[5][0]*p[5][1] - 3*p[3][2]*p[4][1]*p[5][0]*p[5][1] - p[1][1]*p[4][2]*p[5][0]*p[5][1] - 2*p[2][1]*p[4][2]*p[5][0]*p[5][1] + 3*p[3][1]*p[4][2]*p[5][0]*p[5][1] + p[0][2]*(-(pow(p[3][1], 2)*p[2][0]) - pow(p[5][1], 2)*p[2][0] + pow(p[2][1], 2)*p[3][0] - pow(p[4][1], 2)*p[3][0] + pow(p[5][1], 2)*p[3][0] - p[2][0]*p[2][1]*p[3][1] + p[2][1]*p[3][0]*p[3][1] + pow(p[1][1], 2)*(3*p[2][0] - p[3][0] - 2*p[4][0]) + 2*pow(p[3][1], 2)*p[4][0] - 2*p[3][0]*p[3][1]*p[4][1] + p[3][1]*p[4][0]*p[4][1] + p[1][0]*(-3*pow(p[2][1], 2) + pow(p[3][1], 2) + pow(p[4][1], 2) + p[3][1]*p[4][1]) - p[1][1]*(-3*p[2][0]*p[2][1] + p[3][0]*p[3][1] + p[1][0]*(3*p[2][1] - p[3][1] - 2*p[4][1]) + p[3][0]*p[4][1] + p[4][0]*p[4][1]) + 2*pow(p[2][1], 2)*p[5][0] - 2*pow(p[3][1], 2)*p[5][0] - 2*p[2][0]*p[2][1]*p[5][1] + p[2][1]*p[3][0]*p[5][1] - p[2][0]*p[3][1]*p[5][1] + 2*p[3][0]*p[3][1]*p[5][1] + p[2][1]*p[5][0]*p[5][1] - p[3][1]*p[5][0]*p[5][1]) - 2*pow(p[2][1], 2)*p[0][0]*p[5][2] + 2*pow(p[3][1], 2)*p[0][0]*p[5][2] + 2*pow(p[2][1], 2)*p[1][0]*p[5][2] - 2*pow(p[4][1], 2)*p[1][0]*p[5][2] - pow(p[1][1], 2)*p[2][0]*p[5][2] + pow(p[3][1], 2)*p[2][0]*p[5][2] - pow(p[4][1], 2)*p[2][0]*p[5][2] + p[1][0]*p[1][1]*p[2][1]*p[5][2] - 2*p[1][1]*p[2][0]*p[2][1]*p[5][2] - pow(p[2][1], 2)*p[3][0]*p[5][2] + 3*pow(p[4][1], 2)*p[3][0]*p[5][2] + p[2][0]*p[2][1]*p[3][1]*p[5][2] - p[2][1]*p[3][0]*p[3][1]*p[5][2] + pow(p[1][1], 2)*p[4][0]*p[5][2] + pow(p[2][1], 2)*p[4][0]*p[5][2] - 3*pow(p[3][1], 2)*p[4][0]*p[5][2] + p[1][1]*p[2][1]*p[4][0]*p[5][2] - p[1][0]*p[1][1]*p[4][1]*p[5][2] - p[1][1]*p[2][0]*p[4][1]*p[5][2] - p[2][0]*p[2][1]*p[4][1]*p[5][2] + 3*p[3][0]*p[3][1]*p[4][1]*p[5][2] + 2*p[1][1]*p[4][0]*p[4][1]*p[5][2] + p[2][1]*p[4][0]*p[4][1]*p[5][2] - 3*p[3][1]*p[4][0]*p[4][1]*p[5][2] - p[1][1]*p[2][0]*p[5][1]*p[5][2] - p[0][0]*p[2][1]*p[5][1]*p[5][2] + p[1][0]*p[2][1]*p[5][1]*p[5][2] - 2*p[2][1]*p[3][0]*p[5][1]*p[5][2] + p[0][0]*p[3][1]*p[5][1]*p[5][2] + 2*p[2][0]*p[3][1]*p[5][1]*p[5][2] + p[1][1]*p[4][0]*p[5][1]*p[5][2] + 2*p[2][1]*p[4][0]*p[5][1]*p[5][2] - 3*p[3][1]*p[4][0]*p[5][1]*p[5][2] - p[1][0]*p[4][1]*p[5][1]*p[5][2] - 2*p[2][0]*p[4][1]*p[5][1]*p[5][2] + 3*p[3][0]*p[4][1]*p[5][1]*p[5][2] + pow(p[0][1], 2)*(-2*p[2][2]*p[3][0] + 2*p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + 2*p[3][0] + p[4][0]) + p[1][0]*(3*p[2][2] - 2*p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2]) + p[0][1]*(3*p[0][0]*p[1][2]*p[2][1] - 3*p[1][2]*p[2][0]*p[2][1] + 3*p[1][0]*p[2][1]*p[2][2] - p[2][1]*p[2][2]*p[3][0] - 2*p[0][0]*p[1][2]*p[3][1] + 2*p[0][0]*p[2][2]*p[3][1] + p[1][2]*p[3][0]*p[3][1] - p[2][2]*p[3][0]*p[3][1] - 2*p[0][0]*p[2][1]*p[3][2] + p[2][0]*p[2][1]*p[3][2] - p[1][0]*p[3][1]*p[3][2] + p[2][0]*p[3][1]*p[3][2] + p[1][2]*p[3][1]*p[4][0] - 2*p[3][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][2]*p[4][0]*p[4][1] - p[3][2]*p[4][0]*p[4][1] - p[0][0]*p[3][1]*p[4][2] - p[1][0]*p[3][1]*p[4][2] + 2*p[3][0]*p[3][1]*p[4][2] - p[1][0]*p[4][1]*p[4][2] + p[3][0]*p[4][1]*p[4][2] + p[1][1]*(3*p[1][0]*p[2][2] - p[1][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + p[3][0] + 2*p[4][0]) - 2*p[1][0]*p[4][2] + p[3][0]*p[4][2] + p[0][0]*(-3*p[2][2] + 2*p[3][2] + p[4][2])) - 2*p[2][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] + 2*p[3][1]*p[3][2]*p[5][0] + p[0][0]*p[2][2]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][2]*p[5][0]*p[5][1] + p[3][2]*p[5][0]*p[5][1] + p[0][2]*(2*p[2][1]*p[3][0] - 2*p[2][0]*p[3][1] + p[1][1]*(3*p[2][0] - 2*p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-3*p[2][1] + 2*p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[0][0]*p[2][1]*p[5][2] + 2*p[2][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][0]*p[3][1]*p[5][2] - 2*p[3][0]*p[3][1]*p[5][2] + p[2][0]*p[5][1]*p[5][2] - p[3][0]*p[5][1]*p[5][2]))/72;
    result[3][0] = (-3*pow(p[2][2], 2)*p[0][0]*p[1][1] + pow(p[3][2], 2)*p[0][0]*p[1][1] + pow(p[4][2], 2)*p[0][0]*p[1][1] - pow(p[4][2], 2)*p[1][1]*p[2][0] - pow(p[5][2], 2)*p[1][1]*p[2][0] + 3*pow(p[1][2], 2)*p[0][0]*p[2][1] - pow(p[3][2], 2)*p[0][0]*p[2][1] - pow(p[5][2], 2)*p[0][0]*p[2][1] + pow(p[4][2], 2)*p[1][0]*p[2][1] + pow(p[5][2], 2)*p[1][0]*p[2][1] - 3*p[0][0]*p[1][1]*p[1][2]*p[2][2] + 3*p[0][0]*p[1][2]*p[2][1]*p[2][2] + 2*pow(p[4][2], 2)*p[1][1]*p[3][0] - 2*pow(p[5][2], 2)*p[2][1]*p[3][0] - pow(p[1][2], 2)*p[0][0]*p[3][1] + pow(p[2][2], 2)*p[0][0]*p[3][1] - pow(p[4][2], 2)*p[0][0]*p[3][1] + pow(p[5][2], 2)*p[0][0]*p[3][1] - 2*pow(p[4][2], 2)*p[1][0]*p[3][1] + 2*pow(p[5][2], 2)*p[2][0]*p[3][1] + p[0][0]*p[1][1]*p[1][2]*p[3][2] - p[0][0]*p[2][1]*p[2][2]*p[3][2] - p[0][0]*p[1][2]*p[3][1]*p[3][2] + p[0][0]*p[2][2]*p[3][1]*p[3][2] + pow(p[2][2], 2)*p[1][1]*p[4][0] - pow(p[3][2], 2)*p[1][1]*p[4][0] + pow(p[5][2], 2)*p[1][1]*p[4][0] - 2*pow(p[1][2], 2)*p[2][1]*p[4][0] + 2*pow(p[5][2], 2)*p[2][1]*p[4][0] + 2*p[1][1]*p[1][2]*p[2][2]*p[4][0] - p[1][2]*p[2][1]*p[2][2]*p[4][0] + pow(p[1][2], 2)*p[3][1]*p[4][0] - 3*pow(p[5][2], 2)*p[3][1]*p[4][0] - p[1][1]*p[1][2]*p[3][2]*p[4][0] + p[1][2]*p[3][1]*p[3][2]*p[4][0] - 2*pow(p[1][2], 2)*p[0][0]*p[4][1] + 2*pow(p[3][2], 2)*p[0][0]*p[4][1] - pow(p[2][2], 2)*p[1][0]*p[4][1] + pow(p[3][2], 2)*p[1][0]*p[4][1] - pow(p[5][2], 2)*p[1][0]*p[4][1] + 2*pow(p[1][2], 2)*p[2][0]*p[4][1] - 2*pow(p[5][2], 2)*p[2][0]*p[4][1] - 2*p[1][0]*p[1][2]*p[2][2]*p[4][1] + p[1][2]*p[2][0]*p[2][2]*p[4][1] - pow(p[1][2], 2)*p[3][0]*p[4][1] + 3*pow(p[5][2], 2)*p[3][0]*p[4][1] + p[1][0]*p[1][2]*p[3][2]*p[4][1] - p[1][2]*p[3][0]*p[3][2]*p[4][1] + 2*p[0][0]*p[1][1]*p[1][2]*p[4][2] - 2*p[1][1]*p[1][2]*p[2][0]*p[4][2] + 2*p[1][0]*p[1][2]*p[2][1]*p[4][2] - p[1][1]*p[2][0]*p[2][2]*p[4][2] + p[1][0]*p[2][1]*p[2][2]*p[4][2] + p[1][1]*p[1][2]*p[3][0]*p[4][2] - p[0][0]*p[1][2]*p[3][1]*p[4][2] - p[1][0]*p[1][2]*p[3][1]*p[4][2] + p[0][0]*p[1][1]*p[3][2]*p[4][2] + p[1][1]*p[3][0]*p[3][2]*p[4][2] - 2*p[0][0]*p[3][1]*p[3][2]*p[4][2] - p[1][0]*p[3][1]*p[3][2]*p[4][2] - p[1][2]*p[2][1]*p[4][0]*p[4][2] + p[1][1]*p[2][2]*p[4][0]*p[4][2] + 2*p[1][2]*p[3][1]*p[4][0]*p[4][2] - 2*p[1][1]*p[3][2]*p[4][0]*p[4][2] - p[0][0]*p[1][2]*p[4][1]*p[4][2] + p[1][2]*p[2][0]*p[4][1]*p[4][2] - p[1][0]*p[2][2]*p[4][1]*p[4][2] - 2*p[1][2]*p[3][0]*p[4][1]*p[4][2] + p[0][0]*p[3][2]*p[4][1]*p[4][2] + 2*p[1][0]*p[3][2]*p[4][1]*p[4][2] + 2*pow(p[2][2], 2)*p[1][1]*p[5][0] - 2*pow(p[4][2], 2)*p[1][1]*p[5][0] - pow(p[1][2], 2)*p[2][1]*p[5][0] + pow(p[3][2], 2)*p[2][1]*p[5][0] - pow(p[4][2], 2)*p[2][1]*p[5][0] + p[1][1]*p[1][2]*p[2][2]*p[5][0] - 2*p[1][2]*p[2][1]*p[2][2]*p[5][0] - pow(p[2][2], 2)*p[3][1]*p[5][0] + 3*pow(p[4][2], 2)*p[3][1]*p[5][0] + p[2][1]*p[2][2]*p[3][2]*p[5][0] - p[2][2]*p[3][1]*p[3][2]*p[5][0] + pow(p[1][2], 2)*p[4][1]*p[5][0] + pow(p[2][2], 2)*p[4][1]*p[5][0] - 3*pow(p[3][2], 2)*p[4][1]*p[5][0] + p[1][2]*p[2][2]*p[4][1]*p[5][0] - p[1][1]*p[1][2]*p[4][2]*p[5][0] - p[1][2]*p[2][1]*p[4][2]*p[5][0] - p[2][1]*p[2][2]*p[4][2]*p[5][0] + 3*p[3][1]*p[3][2]*p[4][2]*p[5][0] + 2*p[1][2]*p[4][1]*p[4][2]*p[5][0] + p[2][2]*p[4][1]*p[4][2]*p[5][0] - 3*p[3][2]*p[4][1]*p[4][2]*p[5][0] + 2*pow(p[2][2], 2)*p[0][0]*p[5][1] - 2*pow(p[3][2], 2)*p[0][0]*p[5][1] - 2*pow(p[2][2], 2)*p[1][0]*p[5][1] + 2*pow(p[4][2], 2)*p[1][0]*p[5][1] + pow(p[1][2], 2)*p[2][0]*p[5][1] - pow(p[3][2], 2)*p[2][0]*p[5][1] + pow(p[4][2], 2)*p[2][0]*p[5][1] - p[1][0]*p[1][2]*p[2][2]*p[5][1] + 2*p[1][2]*p[2][0]*p[2][2]*p[5][1] + pow(p[2][2], 2)*p[3][0]*p[5][1] - 3*pow(p[4][2], 2)*p[3][0]*p[5][1] - p[2][0]*p[2][2]*p[3][2]*p[5][1] + p[2][2]*p[3][0]*p[3][2]*p[5][1] - pow(p[1][2], 2)*p[4][0]*p[5][1] - pow(p[2][2], 2)*p[4][0]*p[5][1] + 3*pow(p[3][2], 2)*p[4][0]*p[5][1] - p[1][2]*p[2][2]*p[4][0]*p[5][1] + p[1][0]*p[1][2]*p[4][2]*p[5][1] + p[1][2]*p[2][0]*p[4][2]*p[5][1] + p[2][0]*p[2][2]*p[4][2]*p[5][1] - 3*p[3][0]*p[3][2]*p[4][2]*p[5][1] - 2*p[1][2]*p[4][0]*p[4][2]*p[5][1] - p[2][2]*p[4][0]*p[4][2]*p[5][1] + 3*p[3][2]*p[4][0]*p[4][2]*p[5][1] + pow(p[0][2], 2)*(2*p[2][1]*p[3][0] - 2*p[2][0]*p[3][1] + p[1][1]*(3*p[2][0] - 2*p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-3*p[2][1] + 2*p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[1][1]*p[1][2]*p[2][0]*p[5][2] + p[1][0]*p[1][2]*p[2][1]*p[5][2] - 2*p[1][1]*p[2][0]*p[2][2]*p[5][2] - 2*p[0][0]*p[2][1]*p[2][2]*p[5][2] + 2*p[1][0]*p[2][1]*p[2][2]*p[5][2] - p[2][1]*p[2][2]*p[3][0]*p[5][2] + p[0][0]*p[2][2]*p[3][1]*p[5][2] + p[2][0]*p[2][2]*p[3][1]*p[5][2] - p[0][0]*p[2][1]*p[3][2]*p[5][2] - p[2][1]*p[3][0]*p[3][2]*p[5][2] + 2*p[0][0]*p[3][1]*p[3][2]*p[5][2] + p[2][0]*p[3][1]*p[3][2]*p[5][2] + p[1][1]*p[1][2]*p[4][0]*p[5][2] + p[1][1]*p[2][2]*p[4][0]*p[5][2] + p[2][1]*p[2][2]*p[4][0]*p[5][2] - 3*p[3][1]*p[3][2]*p[4][0]*p[5][2] - p[1][0]*p[1][2]*p[4][1]*p[5][2] - p[1][0]*p[2][2]*p[4][1]*p[5][2] - p[2][0]*p[2][2]*p[4][1]*p[5][2] + 3*p[3][0]*p[3][2]*p[4][1]*p[5][2] - p[1][1]*p[2][0]*p[4][2]*p[5][2] + p[1][0]*p[2][1]*p[4][2]*p[5][2] + 2*p[1][1]*p[4][0]*p[4][2]*p[5][2] + p[2][1]*p[4][0]*p[4][2]*p[5][2] - 3*p[3][1]*p[4][0]*p[4][2]*p[5][2] - 2*p[1][0]*p[4][1]*p[4][2]*p[5][2] - p[2][0]*p[4][1]*p[4][2]*p[5][2] + 3*p[3][0]*p[4][1]*p[4][2]*p[5][2] - p[1][2]*p[2][1]*p[5][0]*p[5][2] + p[1][1]*p[2][2]*p[5][0]*p[5][2] - 2*p[2][2]*p[3][1]*p[5][0]*p[5][2] + 2*p[2][1]*p[3][2]*p[5][0]*p[5][2] + p[1][2]*p[4][1]*p[5][0]*p[5][2] + 2*p[2][2]*p[4][1]*p[5][0]*p[5][2] - 3*p[3][2]*p[4][1]*p[5][0]*p[5][2] - p[1][1]*p[4][2]*p[5][0]*p[5][2] - 2*p[2][1]*p[4][2]*p[5][0]*p[5][2] + 3*p[3][1]*p[4][2]*p[5][0]*p[5][2] + p[1][2]*p[2][0]*p[5][1]*p[5][2] + p[0][0]*p[2][2]*p[5][1]*p[5][2] - p[1][0]*p[2][2]*p[5][1]*p[5][2] + 2*p[2][2]*p[3][0]*p[5][1]*p[5][2] - p[0][0]*p[3][2]*p[5][1]*p[5][2] - 2*p[2][0]*p[3][2]*p[5][1]*p[5][2] - p[1][2]*p[4][0]*p[5][1]*p[5][2] - 2*p[2][2]*p[4][0]*p[5][1]*p[5][2] + 3*p[3][2]*p[4][0]*p[5][1]*p[5][2] + p[1][0]*p[4][2]*p[5][1]*p[5][2] + 2*p[2][0]*p[4][2]*p[5][1]*p[5][2] - 3*p[3][0]*p[4][2]*p[5][1]*p[5][2] + p[0][1]*(pow(p[3][2], 2)*p[2][0] + pow(p[5][2], 2)*p[2][0] - pow(p[2][2], 2)*p[3][0] + pow(p[4][2], 2)*p[3][0] - pow(p[5][2], 2)*p[3][0] + p[2][0]*p[2][2]*p[3][2] - p[2][2]*p[3][0]*p[3][2] - 2*pow(p[3][2], 2)*p[4][0] + pow(p[1][2], 2)*(-3*p[2][0] + p[3][0] + 2*p[4][0]) + 2*p[3][0]*p[3][2]*p[4][2] - p[3][2]*p[4][0]*p[4][2] + p[1][0]*(3*pow(p[2][2], 2) - pow(p[3][2], 2) - pow(p[4][2], 2) - p[3][2]*p[4][2]) + p[1][2]*(-3*p[2][0]*p[2][2] + p[3][0]*p[3][2] + p[1][0]*(3*p[2][2] - p[3][2] - 2*p[4][2]) + p[3][0]*p[4][2] + p[4][0]*p[4][2]) - 2*pow(p[2][2], 2)*p[5][0] + 2*pow(p[3][2], 2)*p[5][0] + 2*p[2][0]*p[2][2]*p[5][2] - p[2][2]*p[3][0]*p[5][2] + p[2][0]*p[3][2]*p[5][2] - 2*p[3][0]*p[3][2]*p[5][2] - p[2][2]*p[5][0]*p[5][2] + p[3][2]*p[5][0]*p[5][2]) + p[0][2]*(3*p[0][0]*p[1][2]*p[2][1] - 3*p[1][0]*p[1][2]*p[2][1] - 3*p[1][0]*p[2][1]*p[2][2] + p[2][1]*p[2][2]*p[3][0] - 2*p[0][0]*p[1][2]*p[3][1] + p[1][0]*p[1][2]*p[3][1] + 2*p[0][0]*p[2][2]*p[3][1] - p[2][0]*p[2][2]*p[3][1] - 2*p[0][0]*p[2][1]*p[3][2] + p[2][1]*p[3][0]*p[3][2] + p[1][0]*p[3][1]*p[3][2] - p[2][0]*p[3][1]*p[3][2] + p[1][2]*p[3][1]*p[4][0] + 2*p[3][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + 2*p[1][0]*p[1][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][0]*p[3][2]*p[4][1] - 2*p[3][0]*p[3][2]*p[4][1] - p[0][0]*p[3][1]*p[4][2] + p[3][1]*p[4][0]*p[4][2] + p[1][0]*p[4][1]*p[4][2] - p[3][0]*p[4][1]*p[4][2] + p[1][1]*(3*p[2][0]*p[2][2] - p[3][0]*p[3][2] + p[1][2]*(3*p[2][0] - p[3][0] - 2*p[4][0]) - p[3][2]*p[4][0] - p[4][0]*p[4][2] + p[0][0]*(-3*p[2][2] + 2*p[3][2] + p[4][2])) + 2*p[2][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] - 2*p[3][1]*p[3][2]*p[5][0] + p[0][0]*p[2][2]*p[5][1] - 2*p[2][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][0]*p[3][2]*p[5][1] + 2*p[3][0]*p[3][2]*p[5][1] - p[0][0]*p[2][1]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][1]*p[5][0]*p[5][2] - p[3][1]*p[5][0]*p[5][2] - p[2][0]*p[5][1]*p[5][2] + p[3][0]*p[5][1]*p[5][2] + p[0][1]*(-2*p[2][2]*p[3][0] + 2*p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + 2*p[3][0] + p[4][0]) + p[1][0]*(3*p[2][2] - 2*p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2])))/72;
    return result;
  }
};

// 3D Tetrahedron

template<>
class Operator<VolumeIntegral, basis::Unitary, Tetrahedron>:
  public OperatorBase<VolumeIntegral, basis::Unitary, Interval>
{
public:
  static constexpr size_t dim = dimension(Tetrahedron);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Unitary, dim>::function_size;
  static constexpr size_t point_size = 4;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2]))/6;
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Linear, Tetrahedron>:
  public OperatorBase<VolumeIntegral, basis::Linear, Interval>
{
public:
  static constexpr size_t dim = dimension(Tetrahedron);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Linear, dim>::function_size;
  static constexpr size_t point_size = 4;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2]))/6;
    result[1][0] = ((p[0][0] + p[1][0] + p[2][0] + p[3][0])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;
    result[2][0] = ((p[0][1] + p[1][1] + p[2][1] + p[3][1])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;
    result[3][0] = ((p[0][2] + p[1][2] + p[2][2] + p[3][2])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;
    return result;
  }
};

template<>
class Operator<VolumeIntegral, basis::Quadratic, Tetrahedron>:
  public OperatorBase<VolumeIntegral, basis::Quadratic, Interval>
{
public:
  static constexpr size_t dim = dimension(Tetrahedron);
  static constexpr size_t operator_size = 1;
  static constexpr size_t basis_size =
    basis::Traits<basis::Quadratic, dim>::function_size;
  static constexpr size_t point_size = 4;

  using result_t = std::array<std::array<double, operator_size>, basis_size>;
  using points_t = std::array<Wonton::Point<dim>, point_size>;

  static result_t apply(const points_t p) {
    result_t result;
    result[0][0] = (p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2]))/6;
    result[1][0] = ((p[0][0] + p[1][0] + p[2][0] + p[3][0])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;
    result[2][0] = ((p[0][1] + p[1][1] + p[2][1] + p[3][1])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;
    result[3][0] = ((p[0][2] + p[1][2] + p[2][2] + p[3][2])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;
    result[4][0] = ((pow(p[0][0], 2) + pow(p[1][0], 2) + pow(p[2][0], 2) + pow(p[3][0], 2) + p[2][0]*p[3][0] + p[1][0]*(p[2][0] + p[3][0]) + p[0][0]*(p[1][0] + p[2][0] + p[3][0]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;
    result[5][0] = -((2*p[1][0]*p[1][1] + p[1][1]*p[2][0] + p[1][0]*p[2][1] + 2*p[2][0]*p[2][1] + p[1][1]*p[3][0] + p[2][1]*p[3][0] + p[0][1]*(p[1][0] + p[2][0] + p[3][0]) + p[1][0]*p[3][1] + p[2][0]*p[3][1] + 2*p[3][0]*p[3][1] + p[0][0]*(2*p[0][1] + p[1][1] + p[2][1] + p[3][1]))*(-(p[0][0]*p[1][2]*p[2][1]) + p[0][0]*p[1][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] - p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] + p[0][2]*(-(p[2][1]*p[3][0]) + p[1][1]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][1] - p[3][1]) + p[2][0]*p[3][1]) - p[0][0]*p[1][1]*p[3][2] + p[1][1]*p[2][0]*p[3][2] + p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] + p[0][1]*(p[1][2]*(p[2][0] - p[3][0]) + p[2][2]*p[3][0] - p[2][0]*p[3][2] + p[1][0]*(-p[2][2] + p[3][2]))))/120;
    result[6][0] = ((2*p[1][0]*p[1][2] + p[1][2]*p[2][0] + p[1][0]*p[2][2] + 2*p[2][0]*p[2][2] + p[1][2]*p[3][0] + p[2][2]*p[3][0] + p[0][2]*(p[1][0] + p[2][0] + p[3][0]) + p[1][0]*p[3][2] + p[2][0]*p[3][2] + 2*p[3][0]*p[3][2] + p[0][0]*(2*p[0][2] + p[1][2] + p[2][2] + p[3][2]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;
    result[7][0] = ((pow(p[0][1], 2) + pow(p[1][1], 2) + pow(p[2][1], 2) + pow(p[3][1], 2) + p[2][1]*p[3][1] + p[1][1]*(p[2][1] + p[3][1]) + p[0][1]*(p[1][1] + p[2][1] + p[3][1]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;
    result[8][0] = ((2*p[1][1]*p[1][2] + p[1][2]*p[2][1] + p[1][1]*p[2][2] + 2*p[2][1]*p[2][2] + p[1][2]*p[3][1] + p[2][2]*p[3][1] + p[0][2]*(p[1][1] + p[2][1] + p[3][1]) + p[1][1]*p[3][2] + p[2][1]*p[3][2] + 2*p[3][1]*p[3][2] + p[0][1]*(2*p[0][2] + p[1][2] + p[2][2] + p[3][2]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;
    result[9][0] = ((pow(p[0][2], 2) + pow(p[1][2], 2) + pow(p[2][2], 2) + pow(p[3][2], 2) + p[2][2]*p[3][2] + p[1][2]*(p[2][2] + p[3][2]) + p[0][2]*(p[1][2] + p[2][2] + p[3][2]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;
    return result;
  }
};

////////////////////////////////////////////////////////////////////////////////
// Dynamic size information
////////////////////////////////////////////////////////////////////////////////
inline std::array<size_t,3> size_info(Type type, basis::Type basis, Domain domain) {

  if (type == VolumeIntegral) {
    switch (basis) {
      case basis::Unitary:
        switch (domain) {
          case Interval: {
            using OP = Operator<VolumeIntegral, basis::Unitary, Interval>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Quadrilateral: {
            using OP = Operator<VolumeIntegral, basis::Unitary, Quadrilateral>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Triangle: {
            using OP = Operator<VolumeIntegral, basis::Unitary, Triangle>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Hexahedron: {
            using OP = Operator<VolumeIntegral, basis::Unitary, Hexahedron>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Wedge: {
            using OP = Operator<VolumeIntegral, basis::Unitary, Wedge>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Tetrahedron: {
            using OP = Operator<VolumeIntegral, basis::Unitary, Tetrahedron>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          default: break;
        }
        break;
      case basis::Linear:
        switch (domain) {
          case Interval: {
            using OP = Operator<VolumeIntegral, basis::Linear, Interval>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Quadrilateral: {
            using OP = Operator<VolumeIntegral, basis::Linear, Quadrilateral>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Triangle: {
            using OP = Operator<VolumeIntegral, basis::Linear, Triangle>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Hexahedron: {
            using OP = Operator<VolumeIntegral, basis::Linear, Hexahedron>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Wedge: {
            using OP = Operator<VolumeIntegral, basis::Linear, Wedge>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Tetrahedron: {
            using OP = Operator<VolumeIntegral, basis::Linear, Tetrahedron>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          default: break;
        }
        break;
      case basis::Quadratic:
        switch (domain) {
          case Interval: {
            using OP = Operator<VolumeIntegral, basis::Quadratic, Interval>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Quadrilateral: {
            using OP = Operator<VolumeIntegral, basis::Quadratic, Quadrilateral>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Triangle: {
            using OP = Operator<VolumeIntegral, basis::Quadratic, Triangle>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          case Tetrahedron: {
            using OP = Operator<VolumeIntegral, basis::Quadratic, Tetrahedron>;
            return { OP::operator_size, OP::basis_size, OP::point_size };
          }
          default: break;
        }
        break;
      default: break;
    }
  }
  return { 0, 0, 0 };
}

////////////////////////////////////////////////////////////////////////////////
// Helper template functions
////////////////////////////////////////////////////////////////////////////////

template<class OP>
inline void copy_points(const std::vector<Wonton::Point<OP::dim>>& points,
                        typename OP::points_t& apts) {
  for (int i = 0; i < OP::point_size; i++) {
    apts[i] = points[i];
  }
}

template<class OP>
inline Wonton::Point<OP::dim> centroid(const typename OP::points_t& apts) {
  Wonton::Point<OP::dim> cent;
  for (int j = 0; j < OP::dim; j++) cent[j] = 0.;
  for (int i = 0; i < OP::point_size; i++)
    for (int j = 0; j < OP::dim; j++)
      cent[j] += apts[i][j];
  for (int j = 0; j < OP::dim; j++)
    cent[j] /= OP::point_size;
  return cent;
}

template<class OP>
inline void shift_points(const Wonton::Point<OP::dim> c, typename OP::points_t& apts) {
  for (int i = 0; i < OP::point_size; i++)
    for (int j = 0; j < OP::dim; j++)
      apts[i][j] = apts[i][j] - c[j];
}

template<class OP>
inline void resize_result(std::vector<std::vector<double>>& result) {
  result.resize(OP::basis_size);
  for (auto iter = result.begin(); iter != result.end(); iter++) {
    iter->resize(OP::operator_size);
  }
}

template<class OP>
inline void copy_result(const typename OP::result_t& ares,
                        std::vector<std::vector<double>>& result) {
  for (int i = 0; i < OP::basis_size; i++) {
    for (int j = 0; j < OP::operator_size; j++) {
      result[i][j] = ares[i][j];
    }
  }
}

template<class OP>
inline void get_result(const std::vector<Wonton::Point<OP::dim>>& points,
                       std::vector<std::vector<double>>& result, const bool center = true) {
  resize_result<OP>(result);
  typename OP::points_t apts;
  copy_points<OP>(points, apts);
  if (center) {
    Wonton::Point<OP::dim> c = centroid<OP>(apts);
    shift_points<OP>(c, apts);
    auto tf = basis::transfactor<OP::dim>(OP::basis, c);
    typename OP::result_t ares = OP::apply(apts);
    for (int i = 0; i < OP::basis_size; i++)
      for (int j = 0; j < OP::operator_size; j++) {
        result[i][j] = 0.;
      }
    for (int j = 0; j < OP::operator_size; j++)
      for (int i = 0; i < OP::basis_size; i++)
        for (int k = 0; k < OP::basis_size; k++) {
          result[i][j] += tf[i][k] * ares[k][j];
        }
  } else {
    typename OP::result_t ares = OP::apply(apts);
    copy_result<OP>(ares, result);
  }
}

template<class OP>
inline std::vector<std::vector<double>>
get_result(const std::vector<Wonton::Point<OP::dim>>& points, bool center = true) {
  std::vector<std::vector<double>> result;
  get_result<OP>(points, result, center);
  return result;
}

////////////////////////////////////////////////////////////////////////////////
// Dynamic Accessor
////////////////////////////////////////////////////////////////////////////////

template<size_t dim>
void apply(const Type type, const basis::Type basis_type,
           const Domain domain_type, const std::vector<Wonton::Point<dim>> &points,
           std::vector<std::vector<double>> &result);

template<>
inline
void apply<1>(const Type type, const basis::Type basis_type,
              const Domain domain_type, const std::vector<Wonton::Point<1>>& points,
              std::vector<std::vector<double>>& result) {
  bool center = true;
  switch (type) {
    case VolumeIntegral:
      switch (domain_type) {
        case Interval:
          switch (basis_type) {
            case basis::Unitary: {
              using OP = Operator<VolumeIntegral, basis::Unitary, Interval>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Linear: {
              using OP = Operator<VolumeIntegral, basis::Linear, Interval>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Quadratic: {
              using OP = Operator<VolumeIntegral, basis::Quadratic, Interval>;
              get_result<OP>(points, result, center);
            }
              break;
            default:
              throw (std::runtime_error("invalid basis"));
          }
          break;
        default:
          throw (std::runtime_error("invalid domain"));
      }
      break;
    default:
      throw (std::runtime_error("invalid operator"));
  }
}

template<>
inline
void apply<2>(const Type type, const basis::Type basis_type,
              const Domain domain_type, const std::vector<Wonton::Point<2>>& points,
              std::vector<std::vector<double>>& result) {
  bool center = true;
  switch (type) {
    case VolumeIntegral:
      switch (domain_type) {
        case Quadrilateral:
          switch (basis_type) {
            case basis::Unitary: {
              using OP = Operator<VolumeIntegral, basis::Unitary, Quadrilateral>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Linear: {
              using OP = Operator<VolumeIntegral, basis::Linear, Quadrilateral>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Quadratic: {
              using OP = Operator<VolumeIntegral, basis::Quadratic, Quadrilateral>;
              get_result<OP>(points, result, center);
            }
              break;
            default:
              throw (std::runtime_error("invalid basis"));
          }
          break;
        case Triangle:
          switch (basis_type) {
            case basis::Unitary: {
              using OP = Operator<VolumeIntegral, basis::Unitary, Triangle>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Linear: {
              using OP = Operator<VolumeIntegral, basis::Linear, Triangle>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Quadratic: {
              using OP = Operator<VolumeIntegral, basis::Quadratic, Triangle>;
              get_result<OP>(points, result, center);
            }
              break;
            default:
              throw (std::runtime_error("invalid basis"));
          }
          break;
        default:
          throw (std::runtime_error("invalid domain"));
      }
      break;
    default:
      throw (std::runtime_error("invalid operator"));
  }
}

template<>
inline
void apply<3>(const Type type, const basis::Type basis_type,
              const Domain domain_type, const std::vector<Wonton::Point<3>>& points,
              std::vector<std::vector<double>>& result) {
  bool center = true;
  switch (type) {
    case VolumeIntegral:
      switch (domain_type) {
        case Hexahedron:
          switch (basis_type) {
            case basis::Unitary: {
              using OP = Operator<VolumeIntegral, basis::Unitary, Hexahedron>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Linear: {
              using OP = Operator<VolumeIntegral, basis::Linear, Hexahedron>;
              get_result<OP>(points, result, center);
            }
              break;
            default:
              throw (std::runtime_error("invalid basis"));
          }
          break;
        case Wedge:
          switch (basis_type) {
            case basis::Unitary: {
              using OP = Operator<VolumeIntegral, basis::Unitary, Wedge>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Linear: {
              using OP = Operator<VolumeIntegral, basis::Linear, Wedge>;
              get_result<OP>(points, result, center);
            }
              break;
            default:
              throw (std::runtime_error("invalid basis"));
          }
          break;
        case Tetrahedron:
          switch (basis_type) {
            case basis::Unitary: {
              using OP = Operator<VolumeIntegral, basis::Unitary, Tetrahedron>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Linear: {
              using OP = Operator<VolumeIntegral, basis::Linear, Tetrahedron>;
              get_result<OP>(points, result, center);
            }
              break;
            case basis::Quadratic: {
              using OP = Operator<VolumeIntegral, basis::Quadratic, Tetrahedron>;
              get_result<OP>(points, result, center);
            }
              break;
            default:
              throw (std::runtime_error("invalid basis"));
          }
          break;
        default:
          throw (std::runtime_error("invalid domain"));
      }
      break;
    default:
      throw (std::runtime_error("invalid operator"));
  }
}

}}}  // namespace Portage::Meshfree::Operator

#endif // OPERATOR_H_INC_
