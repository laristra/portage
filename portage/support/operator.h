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
namespace Portage {
  namespace Meshfree {
    namespace Operator {

      using std::array;
      using std::vector;
      using Basis::Traits;
      using Basis::transfactor;

      enum Domain{
        Interval,         // 1D
        Quadrilateral,    // 2D
        Triangle,
        Circle,
        Hexahedron,       // 3D
        Tetrahedron,
        Wedge,
        Sphere,
        LastDomain
      };

      enum Type {VolumeIntegral, SurfaceIntegral, LastOperator};

      template<Domain domain> class DomainTraits  {public: static const size_t dimension=0;};
      template<> class DomainTraits<Interval>     {public: static const size_t dimension=1;};
      template<> class DomainTraits<Quadrilateral>{public: static const size_t dimension=2;};
      template<> class DomainTraits<Triangle>     {public: static const size_t dimension=2;};
      template<> class DomainTraits<Circle>       {public: static const size_t dimension=2;};
      template<> class DomainTraits<Hexahedron>   {public: static const size_t dimension=3;};
      template<> class DomainTraits<Tetrahedron>  {public: static const size_t dimension=3;};
      template<> class DomainTraits<Wedge>        {public: static const size_t dimension=3;};
      template<> class DomainTraits<Sphere>       {public: static const size_t dimension=3;};

      constexpr size_t dimension(Domain domain) {
        switch (domain) {
        case Interval:      return DomainTraits<Interval>::dimension;      break;
        case Quadrilateral: return DomainTraits<Quadrilateral>::dimension; break;
        case Triangle:      return DomainTraits<Triangle>::dimension;      break;
        case Circle:        return DomainTraits<Circle>::dimension;        break;
        case Hexahedron:    return DomainTraits<Hexahedron>::dimension;    break;
        case Tetrahedron:   return DomainTraits<Tetrahedron>::dimension;   break;
        case Wedge:         return DomainTraits<Wedge>::dimension;         break;
        case Sphere:        return DomainTraits<Sphere>::dimension;        break;
        default: return 0;
        }
      }

      template<size_t dim>
      Domain domain_from_points(std::vector<Wonton::Point<dim>> &points) {
        Domain result;
        switch(dim) {
        case 1:
          result = Interval; break;
        case 2:
          switch(points.size()) {
          case 3:
            result = Triangle; break;
          case 4:
            result = Quadrilateral; break;
          default:
            throw(std::runtime_error("invalid number of points for this dimension"));
          }
          break;
        case 3:
          switch(points.size()) {
          case 4:
            result = Tetrahedron; break;
          case 6:
            result = Wedge; break;
          case 8:
            result = Hexahedron; break;
          default:
            throw(std::runtime_error("invalid number of points for this dimension"));
          }
          break;
        }
        return result;
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Template for integral operator base class
      ////////////////////////////////////////////////////////////////////////////////

      template<Type type, Basis::Type basis_type, Domain domain_type>
        class OperatorBase
      {
      public:
        Type operator_type=type;
        static constexpr Basis::Type basis=basis_type;
        static constexpr Domain domain=domain_type;
      };

      ////////////////////////////////////////////////////////////////////////////////
      // Template for integral operators
      ////////////////////////////////////////////////////////////////////////////////

      template<Type type, Basis::Type basis_type, Domain domain_type>
        class Operator: public OperatorBase<type,basis_type,domain_type>
      {
      public:
        static constexpr size_t dim = dimension(domain_type);
        static constexpr size_t operator_size=0;
        static constexpr size_t basis_size=0;
        static constexpr size_t point_size=0;

        using result_t = array<array<double, operator_size>, basis_size>;
        using points_t = array<Wonton::Point<dim>, point_size>;

        static result_t apply(const points_t p) {
          result_t result;
          throw(std::runtime_error("invalid operator"));
          return result;
        }
      };

      ////////////////////////////////////////////////////////////////////////////////
      // Specializations
      ////////////////////////////////////////////////////////////////////////////////

      // 1D

      template<>
        class Operator<VolumeIntegral, Basis::Unitary, Interval>:
          public OperatorBase<VolumeIntegral, Basis::Unitary, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Interval);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Unitary, dim>::function_size;
          static constexpr size_t point_size=2;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;
            result[0][0] = p[1][0] - p[0][0];
            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Linear, Interval>:
          public OperatorBase<VolumeIntegral, Basis::Linear, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Interval);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Linear, dim>::function_size;
          static constexpr size_t point_size=2;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;
            result[0][0] = p[1][0] - p[0][0];
            result[1][0] = .5*(pow(p[1][0],2.) - pow(p[0][0],2.));
            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Quadratic, Interval>:
          public OperatorBase<VolumeIntegral, Basis::Quadratic, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Interval);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Quadratic, dim>::function_size;
          static constexpr size_t point_size=2;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

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
        class Operator<VolumeIntegral, Basis::Unitary, Triangle>:
          public OperatorBase<VolumeIntegral, Basis::Unitary, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Triangle);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Unitary, dim>::function_size;
          static constexpr size_t point_size=3;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

            result[0][0] =
              .5*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) +
                  p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]);

            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Linear, Triangle>:
          public OperatorBase<VolumeIntegral, Basis::Linear, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Triangle);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Linear, dim>::function_size;
          static constexpr size_t point_size=3;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;
            result[0][0] =
              .5*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]);

            result[1][0] =
              ((p[0][0] + p[1][0] + p[2][0])*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]))/6.;

            result[2][0] =
              -((p[0][1] + p[1][1] + p[2][1])*(p[0][1]*(p[1][0] - p[2][0]) + p[1][1]*p[2][0] - p[1][0]*p[2][1] + p[0][0]*(-p[1][1] + p[2][1])))/6.;
            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Quadratic, Triangle>:
          public OperatorBase<VolumeIntegral, Basis::Quadratic, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Triangle);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Quadratic, dim>::function_size;
          static constexpr size_t point_size=3;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;
            result[0][0] =
              .5*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]);

            result[1][0] =
              ((p[0][0] + p[1][0] + p[2][0])*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]))/6.;

            result[2][0] =
              -((p[0][1] + p[1][1] + p[2][1])*(p[0][1]*(p[1][0] - p[2][0]) + p[1][1]*p[2][0] - p[1][0]*p[2][1] + p[0][0]*(-p[1][1] + p[2][1])))/6.;

            result[3][0] =
              ((pow(p[0][0], 2) + pow(p[1][0], 2) + pow(p[2][0], 2) + p[1][0]*p[2][0] + p[0][0]*(p[1][0] + p[2][0]))*(-(p[1][1]*p[2][0]) + p[0][1]*(-p[1][0] + p[2][0]) + p[0][0]*(p[1][1] - p[2][1]) + p[1][0]*p[2][1]))/24.;

            result[4][0] =
              (-(pow(p[1][1], 2)*pow(p[2][0], 2)) + pow(p[0][1], 2)*(-pow(p[1][0], 2) + pow(p[2][0], 2)) + pow(p[1][0], 2)*pow(p[2][1], 2) - 2*pow(p[1][1], 2)*p[1][0]*p[2][0] + 2*pow(p[2][1], 2)*p[1][0]*p[2][0] - 2*p[0][0]*(-(pow(p[1][1], 2)*p[1][0]) + pow(p[0][1], 2)*(p[1][0] - p[2][0]) + pow(p[2][1], 2)*p[2][0]) + 2*pow(p[1][0], 2)*p[1][1]*p[2][1] - 2*pow(p[2][0], 2)*p[1][1]*p[2][1] + pow(p[0][0], 2)*(p[1][1] - p[2][1])*(2*p[0][1] + p[1][1] + p[2][1]) + p[0][1]*(-2*pow(p[1][0], 2)*p[1][1] + 2*pow(p[2][0], 2)*p[2][1]))/24.;

            result[5][0] =
              -((p[0][1]*(p[1][0] - p[2][0]) + p[1][1]*p[2][0] - p[1][0]*p[2][1] + p[0][0]*(-p[1][1] + p[2][1]))*(pow(p[0][1], 2) + pow(p[1][1], 2) + pow(p[2][1], 2) + p[1][1]*p[2][1] + p[0][1]*(p[1][1] + p[2][1])))/24.;
            return result;
          }
        };

      // 2D Quadrilateral

      template<>
        class Operator<VolumeIntegral, Basis::Unitary, Quadrilateral>:
          public OperatorBase<VolumeIntegral, Basis::Unitary, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Quadrilateral);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Unitary, dim>::function_size;
          static constexpr size_t point_size=4;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

            result[0][0] =
              (-(p[1][1]*p[2][0]) + p[1][0]*p[2][1] - p[2][1]*p[3][0] + p[0][1]*(-p[1][0] + p[3][0]) + p[0][0]*(p[1][1] - p[3][1]) + p[2][0]*p[3][1])/2;

            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Linear, Quadrilateral>:
          public OperatorBase<VolumeIntegral, Basis::Linear, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Quadrilateral);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Linear, dim>::function_size;
          static constexpr size_t point_size=4;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

            result[0][0] =
              (-(p[1][1]*p[2][0]) + p[1][0]*p[2][1] - p[2][1]*p[3][0] + p[0][1]*(-p[1][0] + p[3][0]) + p[0][0]*(p[1][1] - p[3][1]) + p[2][0]*p[3][1])/2;

            result[1][0] =
              ((-pow(p[1][0], 2) + pow(p[3][0], 2))*p[0][1] - pow(p[2][0], 2)*p[1][1] - p[1][0]*p[1][1]*p[2][0] + pow(p[1][0], 2)*p[2][1] - pow(p[3][0], 2)*p[2][1] + p[1][0]*p[2][0]*p[2][1] - p[2][0]*p[2][1]*p[3][0] + pow(p[0][0], 2)*(p[1][1] - p[3][1]) + pow(p[2][0], 2)*p[3][1] + p[2][0]*p[3][0]*p[3][1] + p[0][0]*(p[1][0]*p[1][1] + p[0][1]*(-p[1][0] + p[3][0]) - p[3][0]*p[3][1]))/6;

            result[2][0] =
              ((pow(p[1][1], 2) - pow(p[3][1], 2))*p[0][0] + pow(p[2][1], 2)*p[1][0] - pow(p[1][1], 2)*p[2][0] + pow(p[3][1], 2)*p[2][0] + p[1][0]*p[1][1]*p[2][1] - p[1][1]*p[2][0]*p[2][1] - pow(p[2][1], 2)*p[3][0] + pow(p[0][1], 2)*(-p[1][0] + p[3][0]) + p[2][0]*p[2][1]*p[3][1] - p[2][1]*p[3][0]*p[3][1] + p[0][1]*(-(p[1][0]*p[1][1]) + p[0][0]*(p[1][1] - p[3][1]) + p[3][0]*p[3][1]))/6;

            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Quadratic, Quadrilateral>:
          public OperatorBase<VolumeIntegral, Basis::Quadratic, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Quadrilateral);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Quadratic, dim>::function_size;
          static constexpr size_t point_size=4;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

            result[0][0] =
              (-(p[1][1]*p[2][0]) + p[1][0]*p[2][1] - p[2][1]*p[3][0] + p[0][1]*(-p[1][0] + p[3][0]) + p[0][0]*(p[1][1] - p[3][1]) + p[2][0]*p[3][1])/2;

            result[1][0] =
              ((-pow(p[1][0], 2) + pow(p[3][0], 2))*p[0][1] - pow(p[2][0], 2)*p[1][1] - p[1][0]*p[1][1]*p[2][0] + pow(p[1][0], 2)*p[2][1] - pow(p[3][0], 2)*p[2][1] + p[1][0]*p[2][0]*p[2][1] - p[2][0]*p[2][1]*p[3][0] + pow(p[0][0], 2)*(p[1][1] - p[3][1]) + pow(p[2][0], 2)*p[3][1] + p[2][0]*p[3][0]*p[3][1] + p[0][0]*(p[1][0]*p[1][1] + p[0][1]*(-p[1][0] + p[3][0]) - p[3][0]*p[3][1]))/6;

            result[2][0] =
              ((pow(p[1][1], 2) - pow(p[3][1], 2))*p[0][0] + pow(p[2][1], 2)*p[1][0] - pow(p[1][1], 2)*p[2][0] + pow(p[3][1], 2)*p[2][0] + p[1][0]*p[1][1]*p[2][1] - p[1][1]*p[2][0]*p[2][1] - pow(p[2][1], 2)*p[3][0] + pow(p[0][1], 2)*(-p[1][0] + p[3][0]) + p[2][0]*p[2][1]*p[3][1] - p[2][1]*p[3][0]*p[3][1] + p[0][1]*(-(p[1][0]*p[1][1]) + p[0][0]*(p[1][1] - p[3][1]) + p[3][0]*p[3][1]))/6;

            result[3][0] =
              ((-pow(p[1][0], 3) + pow(p[3][0], 3))*p[0][1] - pow(p[2][0], 3)*p[1][1] - pow(p[2][0], 2)*p[1][0]*p[1][1] - pow(p[1][0], 2)*p[1][1]*p[2][0] + pow(p[1][0], 3)*p[2][1] - pow(p[3][0], 3)*p[2][1] + pow(p[2][0], 2)*p[1][0]*p[2][1] + pow(p[1][0], 2)*p[2][0]*p[2][1] - pow(p[3][0], 2)*p[2][0]*p[2][1] - pow(p[2][0], 2)*p[2][1]*p[3][0] + pow(p[0][0], 3)*(p[1][1] - p[3][1]) + pow(p[2][0], 3)*p[3][1] + pow(p[3][0], 2)*p[2][0]*p[3][1] + pow(p[2][0], 2)*p[3][0]*p[3][1] + p[0][0]*((-pow(p[1][0], 2) + pow(p[3][0], 2))*p[0][1] + pow(p[1][0], 2)*p[1][1] - pow(p[3][0], 2)*p[3][1]) + pow(p[0][0], 2)*(p[1][0]*p[1][1] + p[0][1]*(-p[1][0] + p[3][0]) - p[3][0]*p[3][1]))/24;

            result[4][0] =
              (-(pow(p[1][1], 2)*pow(p[2][0], 2)) + pow(p[1][0], 2)*pow(p[2][1], 2) - pow(p[2][1], 2)*pow(p[3][0], 2) + pow(p[0][1], 2)*(-pow(p[1][0], 2) + pow(p[3][0], 2)) + pow(p[2][0], 2)*pow(p[3][1], 2) - 2*pow(p[1][1], 2)*p[1][0]*p[2][0] + 2*pow(p[2][1], 2)*p[1][0]*p[2][0] + 2*pow(p[1][0], 2)*p[1][1]*p[2][1] - 2*pow(p[2][0], 2)*p[1][1]*p[2][1] - 2*pow(p[2][1], 2)*p[2][0]*p[3][0] + 2*pow(p[3][1], 2)*p[2][0]*p[3][0] - 2*p[0][0]*(-(pow(p[1][1], 2)*p[1][0]) + pow(p[0][1], 2)*(p[1][0] - p[3][0]) + pow(p[3][1], 2)*p[3][0]) + 2*pow(p[2][0], 2)*p[2][1]*p[3][1] - 2*pow(p[3][0], 2)*p[2][1]*p[3][1] + pow(p[0][0], 2)*(p[1][1] - p[3][1])*(2*p[0][1] + p[1][1] + p[3][1]) + p[0][1]*(-2*pow(p[1][0], 2)*p[1][1] + 2*pow(p[3][0], 2)*p[3][1]))/24;

            result[5][0] =
              ((pow(p[1][1], 3) - pow(p[3][1], 3))*p[0][0] + pow(p[2][1], 3)*p[1][0] + pow(p[2][1], 2)*p[1][0]*p[1][1] - pow(p[1][1], 3)*p[2][0] + pow(p[3][1], 3)*p[2][0] - pow(p[2][1], 2)*p[1][1]*p[2][0] + pow(p[1][1], 2)*p[1][0]*p[2][1] - pow(p[1][1], 2)*p[2][0]*p[2][1] + pow(p[3][1], 2)*p[2][0]*p[2][1] - pow(p[2][1], 3)*p[3][0] - pow(p[3][1], 2)*p[2][1]*p[3][0] + pow(p[0][1], 3)*(-p[1][0] + p[3][0]) + p[0][1]*((pow(p[1][1], 2) - pow(p[3][1], 2))*p[0][0] - pow(p[1][1], 2)*p[1][0] + pow(p[3][1], 2)*p[3][0]) + pow(p[2][1], 2)*p[2][0]*p[3][1] - pow(p[2][1], 2)*p[3][0]*p[3][1] + pow(p[0][1], 2)*(-(p[1][0]*p[1][1]) + p[0][0]*(p[1][1] - p[3][1]) + p[3][0]*p[3][1]))/24;

            return result;
          }
        };

      // 3D Hexahedron

      template<>
        class Operator<VolumeIntegral, Basis::Unitary, Hexahedron>:
          public OperatorBase<VolumeIntegral, Basis::Unitary, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Hexahedron);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Unitary, dim>::function_size;
          static constexpr size_t point_size=8;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

            result[0][0] =
              (p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] - p[0][0]*p[1][1]*p[3][2] + p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] - p[0][0]*p[1][2]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + p[1][2]*p[2][0]*p[5][1] - p[1][0]*p[2][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[0][0]*p[1][1]*p[5][2] - p[1][1]*p[2][0]*p[5][2] + p[1][0]*p[2][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] - p[0][0]*p[4][1]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[1][2]*p[2][1]*p[6][0] + p[1][1]*p[2][2]*p[6][0] - p[2][2]*p[3][1]*p[6][0] + p[2][1]*p[3][2]*p[6][0] + p[1][2]*p[5][1]*p[6][0] + p[2][2]*p[5][1]*p[6][0] - p[4][2]*p[5][1]*p[6][0] - p[1][1]*p[5][2]*p[6][0] - p[2][1]*p[5][2]*p[6][0] + p[4][1]*p[5][2]*p[6][0] + p[1][2]*p[2][0]*p[6][1] - p[1][0]*p[2][2]*p[6][1] + p[2][2]*p[3][0]*p[6][1] - p[2][0]*p[3][2]*p[6][1] - p[1][2]*p[5][0]*p[6][1] - p[2][2]*p[5][0]*p[6][1] + p[4][2]*p[5][0]*p[6][1] + p[1][0]*p[5][2]*p[6][1] + p[2][0]*p[5][2]*p[6][1] - p[4][0]*p[5][2]*p[6][1] - p[1][1]*p[2][0]*p[6][2] + p[1][0]*p[2][1]*p[6][2] - p[2][1]*p[3][0]*p[6][2] + p[2][0]*p[3][1]*p[6][2] + p[1][1]*p[5][0]*p[6][2] + p[2][1]*p[5][0]*p[6][2] - p[4][1]*p[5][0]*p[6][2] - p[1][0]*p[5][1]*p[6][2] - p[2][0]*p[5][1]*p[6][2] + p[4][0]*p[5][1]*p[6][2] - p[2][2]*p[3][1]*p[7][0] + p[2][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] - p[4][2]*p[5][1]*p[7][0] + p[4][1]*p[5][2]*p[7][0] + p[2][2]*p[6][1]*p[7][0] + p[3][2]*p[6][1]*p[7][0] - p[4][2]*p[6][1]*p[7][0] - p[5][2]*p[6][1]*p[7][0] - p[2][1]*p[6][2]*p[7][0] - p[3][1]*p[6][2]*p[7][0] + p[4][1]*p[6][2]*p[7][0] + p[5][1]*p[6][2]*p[7][0] + p[2][2]*p[3][0]*p[7][1] + p[0][0]*p[3][2]*p[7][1] - p[2][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][0]*p[4][2]*p[7][1] + p[4][2]*p[5][0]*p[7][1] - p[4][0]*p[5][2]*p[7][1] - p[2][2]*p[6][0]*p[7][1] - p[3][2]*p[6][0]*p[7][1] + p[4][2]*p[6][0]*p[7][1] + p[5][2]*p[6][0]*p[7][1] + p[2][0]*p[6][2]*p[7][1] + p[3][0]*p[6][2]*p[7][1] - p[4][0]*p[6][2]*p[7][1] - p[5][0]*p[6][2]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][1]*(p[2][0] + p[3][0] - p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - p[3][1] + p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[2][1]*p[3][0]*p[7][2] - p[0][0]*p[3][1]*p[7][2] + p[2][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][0]*p[4][1]*p[7][2] - p[4][1]*p[5][0]*p[7][2] + p[4][0]*p[5][1]*p[7][2] + p[2][1]*p[6][0]*p[7][2] + p[3][1]*p[6][0]*p[7][2] - p[4][1]*p[6][0]*p[7][2] - p[5][1]*p[6][0]*p[7][2] - p[2][0]*p[6][1]*p[7][2] - p[3][0]*p[6][1]*p[7][2] + p[4][0]*p[6][1]*p[7][2] + p[5][0]*p[6][1]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - p[3][0] + p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + p[3][2] - p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2]))/12;

            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Linear, Hexahedron>:
          public OperatorBase<VolumeIntegral, Basis::Linear, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Hexahedron);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Linear, dim>::function_size;
          static constexpr size_t point_size=8;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

            result[0][0] =
              (p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] - p[0][0]*p[1][1]*p[3][2] + p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] - p[0][0]*p[1][2]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + p[1][2]*p[2][0]*p[5][1] - p[1][0]*p[2][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[0][0]*p[1][1]*p[5][2] - p[1][1]*p[2][0]*p[5][2] + p[1][0]*p[2][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] - p[0][0]*p[4][1]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[1][2]*p[2][1]*p[6][0] + p[1][1]*p[2][2]*p[6][0] - p[2][2]*p[3][1]*p[6][0] + p[2][1]*p[3][2]*p[6][0] + p[1][2]*p[5][1]*p[6][0] + p[2][2]*p[5][1]*p[6][0] - p[4][2]*p[5][1]*p[6][0] - p[1][1]*p[5][2]*p[6][0] - p[2][1]*p[5][2]*p[6][0] + p[4][1]*p[5][2]*p[6][0] + p[1][2]*p[2][0]*p[6][1] - p[1][0]*p[2][2]*p[6][1] + p[2][2]*p[3][0]*p[6][1] - p[2][0]*p[3][2]*p[6][1] - p[1][2]*p[5][0]*p[6][1] - p[2][2]*p[5][0]*p[6][1] + p[4][2]*p[5][0]*p[6][1] + p[1][0]*p[5][2]*p[6][1] + p[2][0]*p[5][2]*p[6][1] - p[4][0]*p[5][2]*p[6][1] - p[1][1]*p[2][0]*p[6][2] + p[1][0]*p[2][1]*p[6][2] - p[2][1]*p[3][0]*p[6][2] + p[2][0]*p[3][1]*p[6][2] + p[1][1]*p[5][0]*p[6][2] + p[2][1]*p[5][0]*p[6][2] - p[4][1]*p[5][0]*p[6][2] - p[1][0]*p[5][1]*p[6][2] - p[2][0]*p[5][1]*p[6][2] + p[4][0]*p[5][1]*p[6][2] - p[2][2]*p[3][1]*p[7][0] + p[2][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] - p[4][2]*p[5][1]*p[7][0] + p[4][1]*p[5][2]*p[7][0] + p[2][2]*p[6][1]*p[7][0] + p[3][2]*p[6][1]*p[7][0] - p[4][2]*p[6][1]*p[7][0] - p[5][2]*p[6][1]*p[7][0] - p[2][1]*p[6][2]*p[7][0] - p[3][1]*p[6][2]*p[7][0] + p[4][1]*p[6][2]*p[7][0] + p[5][1]*p[6][2]*p[7][0] + p[2][2]*p[3][0]*p[7][1] + p[0][0]*p[3][2]*p[7][1] - p[2][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][0]*p[4][2]*p[7][1] + p[4][2]*p[5][0]*p[7][1] - p[4][0]*p[5][2]*p[7][1] - p[2][2]*p[6][0]*p[7][1] - p[3][2]*p[6][0]*p[7][1] + p[4][2]*p[6][0]*p[7][1] + p[5][2]*p[6][0]*p[7][1] + p[2][0]*p[6][2]*p[7][1] + p[3][0]*p[6][2]*p[7][1] - p[4][0]*p[6][2]*p[7][1] - p[5][0]*p[6][2]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][1]*(p[2][0] + p[3][0] - p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - p[3][1] + p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[2][1]*p[3][0]*p[7][2] - p[0][0]*p[3][1]*p[7][2] + p[2][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][0]*p[4][1]*p[7][2] - p[4][1]*p[5][0]*p[7][2] + p[4][0]*p[5][1]*p[7][2] + p[2][1]*p[6][0]*p[7][2] + p[3][1]*p[6][0]*p[7][2] - p[4][1]*p[6][0]*p[7][2] - p[5][1]*p[6][0]*p[7][2] - p[2][0]*p[6][1]*p[7][2] - p[3][0]*p[6][1]*p[7][2] + p[4][0]*p[6][1]*p[7][2] + p[5][0]*p[6][1]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - p[3][0] + p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + p[3][2] - p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2]))/12;

            result[1][0] =
              (-(pow(p[2][0], 2)*p[0][1]*p[1][2]) - pow(p[3][0], 2)*p[0][1]*p[1][2] + pow(p[4][0], 2)*p[0][1]*p[1][2] + pow(p[5][0], 2)*p[0][1]*p[1][2] - 2*p[0][1]*p[1][0]*p[1][2]*p[2][0] + pow(p[3][0], 2)*p[1][2]*p[2][1] - pow(p[5][0], 2)*p[1][2]*p[2][1] - pow(p[6][0], 2)*p[1][2]*p[2][1] + 2*pow(p[1][0], 2)*p[0][1]*p[2][2] - 2*pow(p[3][0], 2)*p[0][1]*p[2][2] - pow(p[3][0], 2)*p[1][1]*p[2][2] + pow(p[5][0], 2)*p[1][1]*p[2][2] + pow(p[6][0], 2)*p[1][1]*p[2][2] + p[0][1]*p[1][0]*p[2][0]*p[2][2] - p[0][1]*p[1][0]*p[1][2]*p[3][0] - p[0][1]*p[1][2]*p[2][0]*p[3][0] + p[1][0]*p[1][2]*p[2][1]*p[3][0] + 2*p[1][2]*p[2][0]*p[2][1]*p[3][0] - p[1][0]*p[1][1]*p[2][2]*p[3][0] - p[0][1]*p[2][0]*p[2][2]*p[3][0] - 2*p[1][1]*p[2][0]*p[2][2]*p[3][0] - 2*pow(p[2][0], 2)*p[1][2]*p[3][1] - p[1][0]*p[1][2]*p[2][0]*p[3][1] + pow(p[1][0], 2)*p[2][2]*p[3][1] - pow(p[6][0], 2)*p[2][2]*p[3][1] - pow(p[7][0], 2)*p[2][2]*p[3][1] + 2*p[1][0]*p[2][0]*p[2][2]*p[3][1] - p[1][2]*p[2][0]*p[3][0]*p[3][1] + p[1][0]*p[2][2]*p[3][0]*p[3][1] + pow(p[1][0], 2)*p[0][1]*p[3][2] + pow(p[2][0], 2)*p[0][1]*p[3][2] - pow(p[4][0], 2)*p[0][1]*p[3][2] - pow(p[7][0], 2)*p[0][1]*p[3][2] + 2*pow(p[2][0], 2)*p[1][1]*p[3][2] + p[0][1]*p[1][0]*p[2][0]*p[3][2] + p[1][0]*p[1][1]*p[2][0]*p[3][2] - pow(p[1][0], 2)*p[2][1]*p[3][2] + pow(p[6][0], 2)*p[2][1]*p[3][2] + pow(p[7][0], 2)*p[2][1]*p[3][2] - 2*p[1][0]*p[2][0]*p[2][1]*p[3][2] + p[0][1]*p[1][0]*p[3][0]*p[3][2] + 2*p[0][1]*p[2][0]*p[3][0]*p[3][2] + p[1][1]*p[2][0]*p[3][0]*p[3][2] - p[1][0]*p[2][1]*p[3][0]*p[3][2] + p[0][1]*p[1][0]*p[1][2]*p[4][0] - p[0][1]*p[3][0]*p[3][2]*p[4][0] + 2*pow(p[5][0], 2)*p[1][2]*p[4][1] - 2*pow(p[7][0], 2)*p[3][2]*p[4][1] - pow(p[1][0], 2)*p[0][1]*p[4][2] + pow(p[3][0], 2)*p[0][1]*p[4][2] - pow(p[5][0], 2)*p[0][1]*p[4][2] + pow(p[7][0], 2)*p[0][1]*p[4][2] - 2*pow(p[5][0], 2)*p[1][1]*p[4][2] + 2*pow(p[7][0], 2)*p[3][1]*p[4][2] - p[0][1]*p[1][0]*p[4][0]*p[4][2] + p[0][1]*p[3][0]*p[4][0]*p[4][2] + 2*p[0][1]*p[1][0]*p[1][2]*p[5][0] - 2*p[1][0]*p[1][2]*p[2][1]*p[5][0] - p[1][2]*p[2][0]*p[2][1]*p[5][0] + 2*p[1][0]*p[1][1]*p[2][2]*p[5][0] + p[1][1]*p[2][0]*p[2][2]*p[5][0] + p[0][1]*p[1][2]*p[4][0]*p[5][0] + p[1][0]*p[1][2]*p[4][1]*p[5][0] + p[1][2]*p[4][0]*p[4][1]*p[5][0] - p[0][1]*p[1][0]*p[4][2]*p[5][0] - p[1][0]*p[1][1]*p[4][2]*p[5][0] - 2*p[0][1]*p[4][0]*p[4][2]*p[5][0] - p[1][1]*p[4][0]*p[4][2]*p[5][0] + pow(p[2][0], 2)*p[1][2]*p[5][1] - pow(p[4][0], 2)*p[1][2]*p[5][1] + pow(p[6][0], 2)*p[1][2]*p[5][1] + 2*p[1][0]*p[1][2]*p[2][0]*p[5][1] - 2*pow(p[1][0], 2)*p[2][2]*p[5][1] + 2*pow(p[6][0], 2)*p[2][2]*p[5][1] - p[1][0]*p[2][0]*p[2][2]*p[5][1] - p[1][0]*p[1][2]*p[4][0]*p[5][1] + pow(p[1][0], 2)*p[4][2]*p[5][1] - pow(p[6][0], 2)*p[4][2]*p[5][1] - pow(p[7][0], 2)*p[4][2]*p[5][1] + p[1][0]*p[4][0]*p[4][2]*p[5][1] + p[1][2]*p[2][0]*p[5][0]*p[5][1] - p[1][0]*p[2][2]*p[5][0]*p[5][1] - 2*p[1][2]*p[4][0]*p[5][0]*p[5][1] + 2*p[1][0]*p[4][2]*p[5][0]*p[5][1] - 2*pow(p[1][0], 2)*p[0][1]*p[5][2] + 2*pow(p[4][0], 2)*p[0][1]*p[5][2] - pow(p[2][0], 2)*p[1][1]*p[5][2] + pow(p[4][0], 2)*p[1][1]*p[5][2] - pow(p[6][0], 2)*p[1][1]*p[5][2] - 2*p[1][0]*p[1][1]*p[2][0]*p[5][2] + 2*pow(p[1][0], 2)*p[2][1]*p[5][2] - 2*pow(p[6][0], 2)*p[2][1]*p[5][2] + p[1][0]*p[2][0]*p[2][1]*p[5][2] + p[1][0]*p[1][1]*p[4][0]*p[5][2] - pow(p[1][0], 2)*p[4][1]*p[5][2] + pow(p[6][0], 2)*p[4][1]*p[5][2] + pow(p[7][0], 2)*p[4][1]*p[5][2] - p[1][0]*p[4][0]*p[4][1]*p[5][2] - p[0][1]*p[1][0]*p[5][0]*p[5][2] - p[1][1]*p[2][0]*p[5][0]*p[5][2] + p[1][0]*p[2][1]*p[5][0]*p[5][2] + p[0][1]*p[4][0]*p[5][0]*p[5][2] + 2*p[1][1]*p[4][0]*p[5][0]*p[5][2] - 2*p[1][0]*p[4][1]*p[5][0]*p[5][2] - p[1][0]*p[1][2]*p[2][1]*p[6][0] - 2*p[1][2]*p[2][0]*p[2][1]*p[6][0] + p[1][0]*p[1][1]*p[2][2]*p[6][0] + 2*p[1][1]*p[2][0]*p[2][2]*p[6][0] - 2*p[2][0]*p[2][2]*p[3][1]*p[6][0] - p[2][2]*p[3][0]*p[3][1]*p[6][0] + 2*p[2][0]*p[2][1]*p[3][2]*p[6][0] + p[2][1]*p[3][0]*p[3][2]*p[6][0] - p[1][2]*p[2][1]*p[5][0]*p[6][0] + p[1][1]*p[2][2]*p[5][0]*p[6][0] + p[1][0]*p[1][2]*p[5][1]*p[6][0] + p[1][2]*p[2][0]*p[5][1]*p[6][0] + p[2][0]*p[2][2]*p[5][1]*p[6][0] - p[4][0]*p[4][2]*p[5][1]*p[6][0] + 2*p[1][2]*p[5][0]*p[5][1]*p[6][0] + p[2][2]*p[5][0]*p[5][1]*p[6][0] - 2*p[4][2]*p[5][0]*p[5][1]*p[6][0] - p[1][0]*p[1][1]*p[5][2]*p[6][0] - p[1][1]*p[2][0]*p[5][2]*p[6][0] - p[2][0]*p[2][1]*p[5][2]*p[6][0] + p[4][0]*p[4][1]*p[5][2]*p[6][0] - 2*p[1][1]*p[5][0]*p[5][2]*p[6][0] - p[2][1]*p[5][0]*p[5][2]*p[6][0] + 2*p[4][1]*p[5][0]*p[5][2]*p[6][0] + 2*pow(p[2][0], 2)*p[1][2]*p[6][1] - 2*pow(p[5][0], 2)*p[1][2]*p[6][1] + p[1][0]*p[1][2]*p[2][0]*p[6][1] - pow(p[1][0], 2)*p[2][2]*p[6][1] + pow(p[3][0], 2)*p[2][2]*p[6][1] - pow(p[5][0], 2)*p[2][2]*p[6][1] + pow(p[7][0], 2)*p[2][2]*p[6][1] - 2*p[1][0]*p[2][0]*p[2][2]*p[6][1] + 2*p[2][0]*p[2][2]*p[3][0]*p[6][1] - 2*pow(p[2][0], 2)*p[3][2]*p[6][1] + 2*pow(p[7][0], 2)*p[3][2]*p[6][1] - p[2][0]*p[3][0]*p[3][2]*p[6][1] + 2*pow(p[5][0], 2)*p[4][2]*p[6][1] - 2*pow(p[7][0], 2)*p[4][2]*p[6][1] - p[1][0]*p[1][2]*p[5][0]*p[6][1] - p[1][0]*p[2][2]*p[5][0]*p[6][1] - p[2][0]*p[2][2]*p[5][0]*p[6][1] + p[4][0]*p[4][2]*p[5][0]*p[6][1] + pow(p[1][0], 2)*p[5][2]*p[6][1] + pow(p[2][0], 2)*p[5][2]*p[6][1] - pow(p[4][0], 2)*p[5][2]*p[6][1] - pow(p[7][0], 2)*p[5][2]*p[6][1] + p[1][0]*p[2][0]*p[5][2]*p[6][1] + 2*p[1][0]*p[5][0]*p[5][2]*p[6][1] + p[2][0]*p[5][0]*p[5][2]*p[6][1] - 2*p[4][0]*p[5][0]*p[5][2]*p[6][1] + p[1][2]*p[2][0]*p[6][0]*p[6][1] - p[1][0]*p[2][2]*p[6][0]*p[6][1] + p[2][2]*p[3][0]*p[6][0]*p[6][1] - p[2][0]*p[3][2]*p[6][0]*p[6][1] - p[1][2]*p[5][0]*p[6][0]*p[6][1] - 2*p[2][2]*p[5][0]*p[6][0]*p[6][1] + p[4][2]*p[5][0]*p[6][0]*p[6][1] + p[1][0]*p[5][2]*p[6][0]*p[6][1] + 2*p[2][0]*p[5][2]*p[6][0]*p[6][1] - p[4][0]*p[5][2]*p[6][0]*p[6][1] - 2*pow(p[2][0], 2)*p[1][1]*p[6][2] + 2*pow(p[5][0], 2)*p[1][1]*p[6][2] - p[1][0]*p[1][1]*p[2][0]*p[6][2] + pow(p[1][0], 2)*p[2][1]*p[6][2] - pow(p[3][0], 2)*p[2][1]*p[6][2] + pow(p[5][0], 2)*p[2][1]*p[6][2] - pow(p[7][0], 2)*p[2][1]*p[6][2] + 2*p[1][0]*p[2][0]*p[2][1]*p[6][2] - 2*p[2][0]*p[2][1]*p[3][0]*p[6][2] + 2*pow(p[2][0], 2)*p[3][1]*p[6][2] - 2*pow(p[7][0], 2)*p[3][1]*p[6][2] + p[2][0]*p[3][0]*p[3][1]*p[6][2] - 2*pow(p[5][0], 2)*p[4][1]*p[6][2] + 2*pow(p[7][0], 2)*p[4][1]*p[6][2] + p[1][0]*p[1][1]*p[5][0]*p[6][2] + p[1][0]*p[2][1]*p[5][0]*p[6][2] + p[2][0]*p[2][1]*p[5][0]*p[6][2] - p[4][0]*p[4][1]*p[5][0]*p[6][2] - pow(p[1][0], 2)*p[5][1]*p[6][2] - pow(p[2][0], 2)*p[5][1]*p[6][2] + pow(p[4][0], 2)*p[5][1]*p[6][2] + pow(p[7][0], 2)*p[5][1]*p[6][2] - p[1][0]*p[2][0]*p[5][1]*p[6][2] - 2*p[1][0]*p[5][0]*p[5][1]*p[6][2] - p[2][0]*p[5][0]*p[5][1]*p[6][2] + 2*p[4][0]*p[5][0]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[6][0]*p[6][2] + p[1][0]*p[2][1]*p[6][0]*p[6][2] - p[2][1]*p[3][0]*p[6][0]*p[6][2] + p[2][0]*p[3][1]*p[6][0]*p[6][2] + p[1][1]*p[5][0]*p[6][0]*p[6][2] + 2*p[2][1]*p[5][0]*p[6][0]*p[6][2] - p[4][1]*p[5][0]*p[6][0]*p[6][2] - p[1][0]*p[5][1]*p[6][0]*p[6][2] - 2*p[2][0]*p[5][1]*p[6][0]*p[6][2] + p[4][0]*p[5][1]*p[6][0]*p[6][2] - p[2][0]*p[2][2]*p[3][1]*p[7][0] - 2*p[2][2]*p[3][0]*p[3][1]*p[7][0] + p[2][0]*p[2][1]*p[3][2]*p[7][0] - 2*p[0][1]*p[3][0]*p[3][2]*p[7][0] + 2*p[2][1]*p[3][0]*p[3][2]*p[7][0] - p[0][1]*p[3][2]*p[4][0]*p[7][0] - p[3][0]*p[3][2]*p[4][1]*p[7][0] - p[3][2]*p[4][0]*p[4][1]*p[7][0] + p[0][1]*p[3][0]*p[4][2]*p[7][0] + p[3][0]*p[3][1]*p[4][2]*p[7][0] + 2*p[0][1]*p[4][0]*p[4][2]*p[7][0] + p[3][1]*p[4][0]*p[4][2]*p[7][0] - 2*p[4][0]*p[4][2]*p[5][1]*p[7][0] - p[4][2]*p[5][0]*p[5][1]*p[7][0] + 2*p[4][0]*p[4][1]*p[5][2]*p[7][0] + p[4][1]*p[5][0]*p[5][2]*p[7][0] - p[2][2]*p[3][1]*p[6][0]*p[7][0] + p[2][1]*p[3][2]*p[6][0]*p[7][0] - p[4][2]*p[5][1]*p[6][0]*p[7][0] + p[4][1]*p[5][2]*p[6][0]*p[7][0] + p[2][0]*p[2][2]*p[6][1]*p[7][0] + p[2][2]*p[3][0]*p[6][1]*p[7][0] + p[3][0]*p[3][2]*p[6][1]*p[7][0] - p[4][0]*p[4][2]*p[6][1]*p[7][0] - p[4][0]*p[5][2]*p[6][1]*p[7][0] - p[5][0]*p[5][2]*p[6][1]*p[7][0] + 2*p[2][2]*p[6][0]*p[6][1]*p[7][0] + p[3][2]*p[6][0]*p[6][1]*p[7][0] - p[4][2]*p[6][0]*p[6][1]*p[7][0] - 2*p[5][2]*p[6][0]*p[6][1]*p[7][0] - p[2][0]*p[2][1]*p[6][2]*p[7][0] - p[2][1]*p[3][0]*p[6][2]*p[7][0] - p[3][0]*p[3][1]*p[6][2]*p[7][0] + p[4][0]*p[4][1]*p[6][2]*p[7][0] + p[4][0]*p[5][1]*p[6][2]*p[7][0] + p[5][0]*p[5][1]*p[6][2]*p[7][0] - 2*p[2][1]*p[6][0]*p[6][2]*p[7][0] - p[3][1]*p[6][0]*p[6][2]*p[7][0] + p[4][1]*p[6][0]*p[6][2]*p[7][0] + 2*p[5][1]*p[6][0]*p[6][2]*p[7][0] + 2*pow(p[3][0], 2)*p[2][2]*p[7][1] - 2*pow(p[6][0], 2)*p[2][2]*p[7][1] + p[2][0]*p[2][2]*p[3][0]*p[7][1] - pow(p[2][0], 2)*p[3][2]*p[7][1] + pow(p[4][0], 2)*p[3][2]*p[7][1] - pow(p[6][0], 2)*p[3][2]*p[7][1] - 2*p[2][0]*p[3][0]*p[3][2]*p[7][1] + p[3][0]*p[3][2]*p[4][0]*p[7][1] - pow(p[3][0], 2)*p[4][2]*p[7][1] + pow(p[5][0], 2)*p[4][2]*p[7][1] + pow(p[6][0], 2)*p[4][2]*p[7][1] - p[3][0]*p[4][0]*p[4][2]*p[7][1] + 2*p[4][0]*p[4][2]*p[5][0]*p[7][1] - 2*pow(p[4][0], 2)*p[5][2]*p[7][1] + 2*pow(p[6][0], 2)*p[5][2]*p[7][1] - p[4][0]*p[5][0]*p[5][2]*p[7][1] - p[2][0]*p[2][2]*p[6][0]*p[7][1] - p[2][0]*p[3][2]*p[6][0]*p[7][1] - p[3][0]*p[3][2]*p[6][0]*p[7][1] + p[4][0]*p[4][2]*p[6][0]*p[7][1] + p[4][2]*p[5][0]*p[6][0]*p[7][1] + p[5][0]*p[5][2]*p[6][0]*p[7][1] + pow(p[2][0], 2)*p[6][2]*p[7][1] + pow(p[3][0], 2)*p[6][2]*p[7][1] - pow(p[4][0], 2)*p[6][2]*p[7][1] - pow(p[5][0], 2)*p[6][2]*p[7][1] + p[2][0]*p[3][0]*p[6][2]*p[7][1] - p[4][0]*p[5][0]*p[6][2]*p[7][1] + 2*p[2][0]*p[6][0]*p[6][2]*p[7][1] + p[3][0]*p[6][0]*p[6][2]*p[7][1] - p[4][0]*p[6][0]*p[6][2]*p[7][1] - 2*p[5][0]*p[6][0]*p[6][2]*p[7][1] + p[2][2]*p[3][0]*p[7][0]*p[7][1] - p[2][0]*p[3][2]*p[7][0]*p[7][1] + 2*p[3][2]*p[4][0]*p[7][0]*p[7][1] - 2*p[3][0]*p[4][2]*p[7][0]*p[7][1] + p[4][2]*p[5][0]*p[7][0]*p[7][1] - p[4][0]*p[5][2]*p[7][0]*p[7][1] - p[2][2]*p[6][0]*p[7][0]*p[7][1] - 2*p[3][2]*p[6][0]*p[7][0]*p[7][1] + 2*p[4][2]*p[6][0]*p[7][0]*p[7][1] + p[5][2]*p[6][0]*p[7][0]*p[7][1] + p[2][0]*p[6][2]*p[7][0]*p[7][1] + 2*p[3][0]*p[6][2]*p[7][0]*p[7][1] - 2*p[4][0]*p[6][2]*p[7][0]*p[7][1] - p[5][0]*p[6][2]*p[7][0]*p[7][1] + p[0][2]*(2*pow(p[3][0], 2)*p[2][1] + p[2][0]*p[2][1]*p[3][0] - pow(p[2][0], 2)*p[3][1] + pow(p[4][0], 2)*p[3][1] + pow(p[7][0], 2)*p[3][1] - 2*p[2][0]*p[3][0]*p[3][1] + p[3][0]*p[3][1]*p[4][0] - pow(p[3][0], 2)*p[4][1] + pow(p[5][0], 2)*p[4][1] - pow(p[7][0], 2)*p[4][1] - p[3][0]*p[4][0]*p[4][1] + 2*p[4][0]*p[4][1]*p[5][0] + p[1][1]*(pow(p[2][0], 2) + pow(p[3][0], 2) - pow(p[4][0], 2) - pow(p[5][0], 2) + p[2][0]*p[3][0] - p[4][0]*p[5][0]) - 2*pow(p[4][0], 2)*p[5][1] - p[4][0]*p[5][0]*p[5][1] + pow(p[1][0], 2)*(-2*p[2][1] - p[3][1] + p[4][1] + 2*p[5][1]) + p[1][0]*(-(p[3][0]*p[3][1]) - p[2][0]*(p[2][1] + p[3][1]) + p[4][0]*p[4][1] + p[1][1]*(2*p[2][0] + p[3][0] - p[4][0] - 2*p[5][0]) + p[4][1]*p[5][0] + p[5][0]*p[5][1]) + 2*p[3][0]*p[3][1]*p[7][0] + p[3][1]*p[4][0]*p[7][0] - p[3][0]*p[4][1]*p[7][0] - 2*p[4][0]*p[4][1]*p[7][0] - 2*pow(p[3][0], 2)*p[7][1] + 2*pow(p[4][0], 2)*p[7][1] - p[3][0]*p[7][0]*p[7][1] + p[4][0]*p[7][0]*p[7][1]) + 2*pow(p[3][0], 2)*p[0][1]*p[7][2] - 2*pow(p[4][0], 2)*p[0][1]*p[7][2] - 2*pow(p[3][0], 2)*p[2][1]*p[7][2] + 2*pow(p[6][0], 2)*p[2][1]*p[7][2] - p[2][0]*p[2][1]*p[3][0]*p[7][2] + pow(p[2][0], 2)*p[3][1]*p[7][2] - pow(p[4][0], 2)*p[3][1]*p[7][2] + pow(p[6][0], 2)*p[3][1]*p[7][2] + 2*p[2][0]*p[3][0]*p[3][1]*p[7][2] - p[3][0]*p[3][1]*p[4][0]*p[7][2] + pow(p[3][0], 2)*p[4][1]*p[7][2] - pow(p[5][0], 2)*p[4][1]*p[7][2] - pow(p[6][0], 2)*p[4][1]*p[7][2] + p[3][0]*p[4][0]*p[4][1]*p[7][2] - 2*p[4][0]*p[4][1]*p[5][0]*p[7][2] + 2*pow(p[4][0], 2)*p[5][1]*p[7][2] - 2*pow(p[6][0], 2)*p[5][1]*p[7][2] + p[4][0]*p[5][0]*p[5][1]*p[7][2] + p[2][0]*p[2][1]*p[6][0]*p[7][2] + p[2][0]*p[3][1]*p[6][0]*p[7][2] + p[3][0]*p[3][1]*p[6][0]*p[7][2] - p[4][0]*p[4][1]*p[6][0]*p[7][2] - p[4][1]*p[5][0]*p[6][0]*p[7][2] - p[5][0]*p[5][1]*p[6][0]*p[7][2] - pow(p[2][0], 2)*p[6][1]*p[7][2] - pow(p[3][0], 2)*p[6][1]*p[7][2] + pow(p[4][0], 2)*p[6][1]*p[7][2] + pow(p[5][0], 2)*p[6][1]*p[7][2] - p[2][0]*p[3][0]*p[6][1]*p[7][2] + p[4][0]*p[5][0]*p[6][1]*p[7][2] - 2*p[2][0]*p[6][0]*p[6][1]*p[7][2] - p[3][0]*p[6][0]*p[6][1]*p[7][2] + p[4][0]*p[6][0]*p[6][1]*p[7][2] + 2*p[5][0]*p[6][0]*p[6][1]*p[7][2] + p[0][1]*p[3][0]*p[7][0]*p[7][2] - p[2][1]*p[3][0]*p[7][0]*p[7][2] + p[2][0]*p[3][1]*p[7][0]*p[7][2] - p[0][1]*p[4][0]*p[7][0]*p[7][2] - 2*p[3][1]*p[4][0]*p[7][0]*p[7][2] + 2*p[3][0]*p[4][1]*p[7][0]*p[7][2] - p[4][1]*p[5][0]*p[7][0]*p[7][2] + p[4][0]*p[5][1]*p[7][0]*p[7][2] + p[2][1]*p[6][0]*p[7][0]*p[7][2] + 2*p[3][1]*p[6][0]*p[7][0]*p[7][2] - 2*p[4][1]*p[6][0]*p[7][0]*p[7][2] - p[5][1]*p[6][0]*p[7][0]*p[7][2] - p[2][0]*p[6][1]*p[7][0]*p[7][2] - 2*p[3][0]*p[6][1]*p[7][0]*p[7][2] + 2*p[4][0]*p[6][1]*p[7][0]*p[7][2] + p[5][0]*p[6][1]*p[7][0]*p[7][2] + pow(p[0][0], 2)*(p[2][2]*p[3][1] - p[2][1]*p[3][2] + 2*p[3][2]*p[4][1] - 2*p[3][1]*p[4][2] + p[1][2]*(p[2][1] + 2*p[3][1] - 2*p[4][1] - p[5][1]) + p[4][2]*p[5][1] - p[4][1]*p[5][2] + p[1][1]*(-p[2][2] - 2*p[3][2] + 2*p[4][2] + p[5][2]) + p[3][2]*p[7][1] - p[4][2]*p[7][1] - p[3][1]*p[7][2] + p[4][1]*p[7][2]) + p[0][0]*(2*p[1][0]*p[1][2]*p[2][1] + p[1][2]*p[2][0]*p[2][1] - 2*p[1][0]*p[1][1]*p[2][2] - p[1][1]*p[2][0]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[1][0]*p[1][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] + p[2][0]*p[2][2]*p[3][1] + p[1][2]*p[3][0]*p[3][1] + 2*p[2][2]*p[3][0]*p[3][1] - p[1][0]*p[1][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] - p[2][0]*p[2][1]*p[3][2] - p[1][1]*p[3][0]*p[3][2] - 2*p[2][1]*p[3][0]*p[3][2] - p[1][0]*p[1][2]*p[4][1] + p[3][0]*p[3][2]*p[4][1] - p[1][2]*p[4][0]*p[4][1] + p[3][2]*p[4][0]*p[4][1] + p[1][0]*p[1][1]*p[4][2] - p[3][0]*p[3][1]*p[4][2] + p[1][1]*p[4][0]*p[4][2] - p[3][1]*p[4][0]*p[4][2] - 2*p[1][0]*p[1][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + 2*p[4][0]*p[4][2]*p[5][1] - p[1][2]*p[5][0]*p[5][1] + p[4][2]*p[5][0]*p[5][1] + 2*p[1][0]*p[1][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - 2*p[4][0]*p[4][1]*p[5][2] + p[1][1]*p[5][0]*p[5][2] - p[4][1]*p[5][0]*p[5][2] + 2*p[3][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[3][0]*p[4][2]*p[7][1] - 2*p[4][0]*p[4][2]*p[7][1] + p[3][2]*p[7][0]*p[7][1] - p[4][2]*p[7][0]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + 2*p[3][1]*p[4][0] - 2*p[3][0]*p[4][1] + p[1][1]*(p[2][0] + 2*p[3][0] - 2*p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - 2*p[3][1] + 2*p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - 2*p[3][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[3][0]*p[4][1]*p[7][2] + 2*p[4][0]*p[4][1]*p[7][2] - p[3][1]*p[7][0]*p[7][2] + p[4][1]*p[7][0]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - 2*p[3][2]*p[4][0] + 2*p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - 2*p[3][0] + 2*p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + 2*p[3][2] - 2*p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2])))/72;

            result[2][0] =
              (pow(p[2][1], 2)*p[0][0]*p[1][2] + pow(p[3][1], 2)*p[0][0]*p[1][2] - pow(p[4][1], 2)*p[0][0]*p[1][2] - pow(p[5][1], 2)*p[0][0]*p[1][2] - pow(p[3][1], 2)*p[1][2]*p[2][0] + pow(p[5][1], 2)*p[1][2]*p[2][0] + pow(p[6][1], 2)*p[1][2]*p[2][0] + 2*p[0][0]*p[1][1]*p[1][2]*p[2][1] - 2*pow(p[1][1], 2)*p[0][0]*p[2][2] + 2*pow(p[3][1], 2)*p[0][0]*p[2][2] + pow(p[3][1], 2)*p[1][0]*p[2][2] - pow(p[5][1], 2)*p[1][0]*p[2][2] - pow(p[6][1], 2)*p[1][0]*p[2][2] - p[0][0]*p[1][1]*p[2][1]*p[2][2] + 2*pow(p[2][1], 2)*p[1][2]*p[3][0] + p[1][1]*p[1][2]*p[2][1]*p[3][0] - pow(p[1][1], 2)*p[2][2]*p[3][0] + pow(p[6][1], 2)*p[2][2]*p[3][0] + pow(p[7][1], 2)*p[2][2]*p[3][0] - 2*p[1][1]*p[2][1]*p[2][2]*p[3][0] + p[0][0]*p[1][1]*p[1][2]*p[3][1] - p[1][1]*p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[1][2]*p[2][1]*p[3][1] - 2*p[1][2]*p[2][0]*p[2][1]*p[3][1] + p[1][0]*p[1][1]*p[2][2]*p[3][1] + p[0][0]*p[2][1]*p[2][2]*p[3][1] + 2*p[1][0]*p[2][1]*p[2][2]*p[3][1] + p[1][2]*p[2][1]*p[3][0]*p[3][1] - p[1][1]*p[2][2]*p[3][0]*p[3][1] - pow(p[1][1], 2)*p[0][0]*p[3][2] - pow(p[2][1], 2)*p[0][0]*p[3][2] + pow(p[4][1], 2)*p[0][0]*p[3][2] + pow(p[7][1], 2)*p[0][0]*p[3][2] - 2*pow(p[2][1], 2)*p[1][0]*p[3][2] + pow(p[1][1], 2)*p[2][0]*p[3][2] - pow(p[6][1], 2)*p[2][0]*p[3][2] - pow(p[7][1], 2)*p[2][0]*p[3][2] - p[0][0]*p[1][1]*p[2][1]*p[3][2] - p[1][0]*p[1][1]*p[2][1]*p[3][2] + 2*p[1][1]*p[2][0]*p[2][1]*p[3][2] - p[0][0]*p[1][1]*p[3][1]*p[3][2] + p[1][1]*p[2][0]*p[3][1]*p[3][2] - 2*p[0][0]*p[2][1]*p[3][1]*p[3][2] - p[1][0]*p[2][1]*p[3][1]*p[3][2] - 2*pow(p[5][1], 2)*p[1][2]*p[4][0] + 2*pow(p[7][1], 2)*p[3][2]*p[4][0] - p[0][0]*p[1][1]*p[1][2]*p[4][1] + p[0][0]*p[3][1]*p[3][2]*p[4][1] + pow(p[1][1], 2)*p[0][0]*p[4][2] - pow(p[3][1], 2)*p[0][0]*p[4][2] + pow(p[5][1], 2)*p[0][0]*p[4][2] - pow(p[7][1], 2)*p[0][0]*p[4][2] + 2*pow(p[5][1], 2)*p[1][0]*p[4][2] - 2*pow(p[7][1], 2)*p[3][0]*p[4][2] + p[0][0]*p[1][1]*p[4][1]*p[4][2] - p[0][0]*p[3][1]*p[4][1]*p[4][2] - pow(p[2][1], 2)*p[1][2]*p[5][0] + pow(p[4][1], 2)*p[1][2]*p[5][0] - pow(p[6][1], 2)*p[1][2]*p[5][0] - 2*p[1][1]*p[1][2]*p[2][1]*p[5][0] + 2*pow(p[1][1], 2)*p[2][2]*p[5][0] - 2*pow(p[6][1], 2)*p[2][2]*p[5][0] + p[1][1]*p[2][1]*p[2][2]*p[5][0] + p[1][1]*p[1][2]*p[4][1]*p[5][0] - pow(p[1][1], 2)*p[4][2]*p[5][0] + pow(p[6][1], 2)*p[4][2]*p[5][0] + pow(p[7][1], 2)*p[4][2]*p[5][0] - p[1][1]*p[4][1]*p[4][2]*p[5][0] - 2*p[0][0]*p[1][1]*p[1][2]*p[5][1] + 2*p[1][1]*p[1][2]*p[2][0]*p[5][1] + p[1][2]*p[2][0]*p[2][1]*p[5][1] - 2*p[1][0]*p[1][1]*p[2][2]*p[5][1] - p[1][0]*p[2][1]*p[2][2]*p[5][1] - p[1][1]*p[1][2]*p[4][0]*p[5][1] - p[0][0]*p[1][2]*p[4][1]*p[5][1] - p[1][2]*p[4][0]*p[4][1]*p[5][1] + p[0][0]*p[1][1]*p[4][2]*p[5][1] + p[1][0]*p[1][1]*p[4][2]*p[5][1] + 2*p[0][0]*p[4][1]*p[4][2]*p[5][1] + p[1][0]*p[4][1]*p[4][2]*p[5][1] - p[1][2]*p[2][1]*p[5][0]*p[5][1] + p[1][1]*p[2][2]*p[5][0]*p[5][1] + 2*p[1][2]*p[4][1]*p[5][0]*p[5][1] - 2*p[1][1]*p[4][2]*p[5][0]*p[5][1] + 2*pow(p[1][1], 2)*p[0][0]*p[5][2] - 2*pow(p[4][1], 2)*p[0][0]*p[5][2] + pow(p[2][1], 2)*p[1][0]*p[5][2] - pow(p[4][1], 2)*p[1][0]*p[5][2] + pow(p[6][1], 2)*p[1][0]*p[5][2] - 2*pow(p[1][1], 2)*p[2][0]*p[5][2] + 2*pow(p[6][1], 2)*p[2][0]*p[5][2] + 2*p[1][0]*p[1][1]*p[2][1]*p[5][2] - p[1][1]*p[2][0]*p[2][1]*p[5][2] + pow(p[1][1], 2)*p[4][0]*p[5][2] - pow(p[6][1], 2)*p[4][0]*p[5][2] - pow(p[7][1], 2)*p[4][0]*p[5][2] - p[1][0]*p[1][1]*p[4][1]*p[5][2] + p[1][1]*p[4][0]*p[4][1]*p[5][2] + p[0][0]*p[1][1]*p[5][1]*p[5][2] - p[1][1]*p[2][0]*p[5][1]*p[5][2] + p[1][0]*p[2][1]*p[5][1]*p[5][2] + 2*p[1][1]*p[4][0]*p[5][1]*p[5][2] - p[0][0]*p[4][1]*p[5][1]*p[5][2] - 2*p[1][0]*p[4][1]*p[5][1]*p[5][2] - 2*pow(p[2][1], 2)*p[1][2]*p[6][0] + 2*pow(p[5][1], 2)*p[1][2]*p[6][0] - p[1][1]*p[1][2]*p[2][1]*p[6][0] + pow(p[1][1], 2)*p[2][2]*p[6][0] - pow(p[3][1], 2)*p[2][2]*p[6][0] + pow(p[5][1], 2)*p[2][2]*p[6][0] - pow(p[7][1], 2)*p[2][2]*p[6][0] + 2*p[1][1]*p[2][1]*p[2][2]*p[6][0] - 2*p[2][1]*p[2][2]*p[3][1]*p[6][0] + 2*pow(p[2][1], 2)*p[3][2]*p[6][0] - 2*pow(p[7][1], 2)*p[3][2]*p[6][0] + p[2][1]*p[3][1]*p[3][2]*p[6][0] - 2*pow(p[5][1], 2)*p[4][2]*p[6][0] + 2*pow(p[7][1], 2)*p[4][2]*p[6][0] + p[1][1]*p[1][2]*p[5][1]*p[6][0] + p[1][1]*p[2][2]*p[5][1]*p[6][0] + p[2][1]*p[2][2]*p[5][1]*p[6][0] - p[4][1]*p[4][2]*p[5][1]*p[6][0] - pow(p[1][1], 2)*p[5][2]*p[6][0] - pow(p[2][1], 2)*p[5][2]*p[6][0] + pow(p[4][1], 2)*p[5][2]*p[6][0] + pow(p[7][1], 2)*p[5][2]*p[6][0] - p[1][1]*p[2][1]*p[5][2]*p[6][0] - 2*p[1][1]*p[5][1]*p[5][2]*p[6][0] - p[2][1]*p[5][1]*p[5][2]*p[6][0] + 2*p[4][1]*p[5][1]*p[5][2]*p[6][0] + p[1][1]*p[1][2]*p[2][0]*p[6][1] + 2*p[1][2]*p[2][0]*p[2][1]*p[6][1] - p[1][0]*p[1][1]*p[2][2]*p[6][1] - 2*p[1][0]*p[2][1]*p[2][2]*p[6][1] + 2*p[2][1]*p[2][2]*p[3][0]*p[6][1] + p[2][2]*p[3][0]*p[3][1]*p[6][1] - 2*p[2][0]*p[2][1]*p[3][2]*p[6][1] - p[2][0]*p[3][1]*p[3][2]*p[6][1] - p[1][1]*p[1][2]*p[5][0]*p[6][1] - p[1][2]*p[2][1]*p[5][0]*p[6][1] - p[2][1]*p[2][2]*p[5][0]*p[6][1] + p[4][1]*p[4][2]*p[5][0]*p[6][1] + p[1][2]*p[2][0]*p[5][1]*p[6][1] - p[1][0]*p[2][2]*p[5][1]*p[6][1] - 2*p[1][2]*p[5][0]*p[5][1]*p[6][1] - p[2][2]*p[5][0]*p[5][1]*p[6][1] + 2*p[4][2]*p[5][0]*p[5][1]*p[6][1] + p[1][0]*p[1][1]*p[5][2]*p[6][1] + p[1][0]*p[2][1]*p[5][2]*p[6][1] + p[2][0]*p[2][1]*p[5][2]*p[6][1] - p[4][0]*p[4][1]*p[5][2]*p[6][1] + 2*p[1][0]*p[5][1]*p[5][2]*p[6][1] + p[2][0]*p[5][1]*p[5][2]*p[6][1] - 2*p[4][0]*p[5][1]*p[5][2]*p[6][1] - p[1][2]*p[2][1]*p[6][0]*p[6][1] + p[1][1]*p[2][2]*p[6][0]*p[6][1] - p[2][2]*p[3][1]*p[6][0]*p[6][1] + p[2][1]*p[3][2]*p[6][0]*p[6][1] + p[1][2]*p[5][1]*p[6][0]*p[6][1] + 2*p[2][2]*p[5][1]*p[6][0]*p[6][1] - p[4][2]*p[5][1]*p[6][0]*p[6][1] - p[1][1]*p[5][2]*p[6][0]*p[6][1] - 2*p[2][1]*p[5][2]*p[6][0]*p[6][1] + p[4][1]*p[5][2]*p[6][0]*p[6][1] + 2*pow(p[2][1], 2)*p[1][0]*p[6][2] - 2*pow(p[5][1], 2)*p[1][0]*p[6][2] - pow(p[1][1], 2)*p[2][0]*p[6][2] + pow(p[3][1], 2)*p[2][0]*p[6][2] - pow(p[5][1], 2)*p[2][0]*p[6][2] + pow(p[7][1], 2)*p[2][0]*p[6][2] + p[1][0]*p[1][1]*p[2][1]*p[6][2] - 2*p[1][1]*p[2][0]*p[2][1]*p[6][2] - 2*pow(p[2][1], 2)*p[3][0]*p[6][2] + 2*pow(p[7][1], 2)*p[3][0]*p[6][2] + 2*p[2][0]*p[2][1]*p[3][1]*p[6][2] - p[2][1]*p[3][0]*p[3][1]*p[6][2] + 2*pow(p[5][1], 2)*p[4][0]*p[6][2] - 2*pow(p[7][1], 2)*p[4][0]*p[6][2] + pow(p[1][1], 2)*p[5][0]*p[6][2] + pow(p[2][1], 2)*p[5][0]*p[6][2] - pow(p[4][1], 2)*p[5][0]*p[6][2] - pow(p[7][1], 2)*p[5][0]*p[6][2] + p[1][1]*p[2][1]*p[5][0]*p[6][2] - p[1][0]*p[1][1]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[5][1]*p[6][2] - p[2][0]*p[2][1]*p[5][1]*p[6][2] + p[4][0]*p[4][1]*p[5][1]*p[6][2] + 2*p[1][1]*p[5][0]*p[5][1]*p[6][2] + p[2][1]*p[5][0]*p[5][1]*p[6][2] - 2*p[4][1]*p[5][0]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[6][1]*p[6][2] + p[1][0]*p[2][1]*p[6][1]*p[6][2] - p[2][1]*p[3][0]*p[6][1]*p[6][2] + p[2][0]*p[3][1]*p[6][1]*p[6][2] + p[1][1]*p[5][0]*p[6][1]*p[6][2] + 2*p[2][1]*p[5][0]*p[6][1]*p[6][2] - p[4][1]*p[5][0]*p[6][1]*p[6][2] - p[1][0]*p[5][1]*p[6][1]*p[6][2] - 2*p[2][0]*p[5][1]*p[6][1]*p[6][2] + p[4][0]*p[5][1]*p[6][1]*p[6][2] - 2*pow(p[3][1], 2)*p[2][2]*p[7][0] + 2*pow(p[6][1], 2)*p[2][2]*p[7][0] - p[2][1]*p[2][2]*p[3][1]*p[7][0] + pow(p[2][1], 2)*p[3][2]*p[7][0] - pow(p[4][1], 2)*p[3][2]*p[7][0] + pow(p[6][1], 2)*p[3][2]*p[7][0] + 2*p[2][1]*p[3][1]*p[3][2]*p[7][0] - p[3][1]*p[3][2]*p[4][1]*p[7][0] + pow(p[3][1], 2)*p[4][2]*p[7][0] - pow(p[5][1], 2)*p[4][2]*p[7][0] - pow(p[6][1], 2)*p[4][2]*p[7][0] + p[3][1]*p[4][1]*p[4][2]*p[7][0] - 2*p[4][1]*p[4][2]*p[5][1]*p[7][0] + 2*pow(p[4][1], 2)*p[5][2]*p[7][0] - 2*pow(p[6][1], 2)*p[5][2]*p[7][0] + p[4][1]*p[5][1]*p[5][2]*p[7][0] + p[2][1]*p[2][2]*p[6][1]*p[7][0] + p[2][1]*p[3][2]*p[6][1]*p[7][0] + p[3][1]*p[3][2]*p[6][1]*p[7][0] - p[4][1]*p[4][2]*p[6][1]*p[7][0] - p[4][2]*p[5][1]*p[6][1]*p[7][0] - p[5][1]*p[5][2]*p[6][1]*p[7][0] - pow(p[2][1], 2)*p[6][2]*p[7][0] - pow(p[3][1], 2)*p[6][2]*p[7][0] + pow(p[4][1], 2)*p[6][2]*p[7][0] + pow(p[5][1], 2)*p[6][2]*p[7][0] - p[2][1]*p[3][1]*p[6][2]*p[7][0] + p[4][1]*p[5][1]*p[6][2]*p[7][0] - 2*p[2][1]*p[6][1]*p[6][2]*p[7][0] - p[3][1]*p[6][1]*p[6][2]*p[7][0] + p[4][1]*p[6][1]*p[6][2]*p[7][0] + 2*p[5][1]*p[6][1]*p[6][2]*p[7][0] + p[2][1]*p[2][2]*p[3][0]*p[7][1] + 2*p[2][2]*p[3][0]*p[3][1]*p[7][1] - p[2][0]*p[2][1]*p[3][2]*p[7][1] + 2*p[0][0]*p[3][1]*p[3][2]*p[7][1] - 2*p[2][0]*p[3][1]*p[3][2]*p[7][1] + p[3][1]*p[3][2]*p[4][0]*p[7][1] + p[0][0]*p[3][2]*p[4][1]*p[7][1] + p[3][2]*p[4][0]*p[4][1]*p[7][1] - p[0][0]*p[3][1]*p[4][2]*p[7][1] - p[3][0]*p[3][1]*p[4][2]*p[7][1] - 2*p[0][0]*p[4][1]*p[4][2]*p[7][1] - p[3][0]*p[4][1]*p[4][2]*p[7][1] + 2*p[4][1]*p[4][2]*p[5][0]*p[7][1] + p[4][2]*p[5][0]*p[5][1]*p[7][1] - 2*p[4][0]*p[4][1]*p[5][2]*p[7][1] - p[4][0]*p[5][1]*p[5][2]*p[7][1] - p[2][1]*p[2][2]*p[6][0]*p[7][1] - p[2][2]*p[3][1]*p[6][0]*p[7][1] - p[3][1]*p[3][2]*p[6][0]*p[7][1] + p[4][1]*p[4][2]*p[6][0]*p[7][1] + p[4][1]*p[5][2]*p[6][0]*p[7][1] + p[5][1]*p[5][2]*p[6][0]*p[7][1] + p[2][2]*p[3][0]*p[6][1]*p[7][1] - p[2][0]*p[3][2]*p[6][1]*p[7][1] + p[4][2]*p[5][0]*p[6][1]*p[7][1] - p[4][0]*p[5][2]*p[6][1]*p[7][1] - 2*p[2][2]*p[6][0]*p[6][1]*p[7][1] - p[3][2]*p[6][0]*p[6][1]*p[7][1] + p[4][2]*p[6][0]*p[6][1]*p[7][1] + 2*p[5][2]*p[6][0]*p[6][1]*p[7][1] + p[2][0]*p[2][1]*p[6][2]*p[7][1] + p[2][0]*p[3][1]*p[6][2]*p[7][1] + p[3][0]*p[3][1]*p[6][2]*p[7][1] - p[4][0]*p[4][1]*p[6][2]*p[7][1] - p[4][1]*p[5][0]*p[6][2]*p[7][1] - p[5][0]*p[5][1]*p[6][2]*p[7][1] + 2*p[2][0]*p[6][1]*p[6][2]*p[7][1] + p[3][0]*p[6][1]*p[6][2]*p[7][1] - p[4][0]*p[6][1]*p[6][2]*p[7][1] - 2*p[5][0]*p[6][1]*p[6][2]*p[7][1] - p[2][2]*p[3][1]*p[7][0]*p[7][1] + p[2][1]*p[3][2]*p[7][0]*p[7][1] - 2*p[3][2]*p[4][1]*p[7][0]*p[7][1] + 2*p[3][1]*p[4][2]*p[7][0]*p[7][1] - p[4][2]*p[5][1]*p[7][0]*p[7][1] + p[4][1]*p[5][2]*p[7][0]*p[7][1] + p[2][2]*p[6][1]*p[7][0]*p[7][1] + 2*p[3][2]*p[6][1]*p[7][0]*p[7][1] - 2*p[4][2]*p[6][1]*p[7][0]*p[7][1] - p[5][2]*p[6][1]*p[7][0]*p[7][1] - p[2][1]*p[6][2]*p[7][0]*p[7][1] - 2*p[3][1]*p[6][2]*p[7][0]*p[7][1] + 2*p[4][1]*p[6][2]*p[7][0]*p[7][1] + p[5][1]*p[6][2]*p[7][0]*p[7][1] + p[0][2]*(-2*pow(p[3][1], 2)*p[2][0] + pow(p[2][1], 2)*p[3][0] - pow(p[4][1], 2)*p[3][0] - pow(p[7][1], 2)*p[3][0] - p[2][0]*p[2][1]*p[3][1] + 2*p[2][1]*p[3][0]*p[3][1] + pow(p[3][1], 2)*p[4][0] - pow(p[5][1], 2)*p[4][0] + pow(p[7][1], 2)*p[4][0] - p[3][0]*p[3][1]*p[4][1] + p[3][1]*p[4][0]*p[4][1] + pow(p[1][1], 2)*(2*p[2][0] + p[3][0] - p[4][0] - 2*p[5][0]) + 2*pow(p[4][1], 2)*p[5][0] - 2*p[4][0]*p[4][1]*p[5][1] + p[4][1]*p[5][0]*p[5][1] + p[1][0]*(-pow(p[2][1], 2) - pow(p[3][1], 2) + pow(p[4][1], 2) + pow(p[5][1], 2) - p[2][1]*p[3][1] + p[4][1]*p[5][1]) + p[1][1]*(p[2][0]*p[2][1] + p[2][1]*p[3][0] + p[3][0]*p[3][1] - p[4][0]*p[4][1] - p[4][0]*p[5][1] - p[5][0]*p[5][1] + p[1][0]*(-2*p[2][1] - p[3][1] + p[4][1] + 2*p[5][1])) + 2*pow(p[3][1], 2)*p[7][0] - 2*pow(p[4][1], 2)*p[7][0] - 2*p[3][0]*p[3][1]*p[7][1] + p[3][1]*p[4][0]*p[7][1] - p[3][0]*p[4][1]*p[7][1] + 2*p[4][0]*p[4][1]*p[7][1] + p[3][1]*p[7][0]*p[7][1] - p[4][1]*p[7][0]*p[7][1]) - 2*pow(p[3][1], 2)*p[0][0]*p[7][2] + 2*pow(p[4][1], 2)*p[0][0]*p[7][2] + 2*pow(p[3][1], 2)*p[2][0]*p[7][2] - 2*pow(p[6][1], 2)*p[2][0]*p[7][2] - pow(p[2][1], 2)*p[3][0]*p[7][2] + pow(p[4][1], 2)*p[3][0]*p[7][2] - pow(p[6][1], 2)*p[3][0]*p[7][2] + p[2][0]*p[2][1]*p[3][1]*p[7][2] - 2*p[2][1]*p[3][0]*p[3][1]*p[7][2] - pow(p[3][1], 2)*p[4][0]*p[7][2] + pow(p[5][1], 2)*p[4][0]*p[7][2] + pow(p[6][1], 2)*p[4][0]*p[7][2] + p[3][0]*p[3][1]*p[4][1]*p[7][2] - p[3][1]*p[4][0]*p[4][1]*p[7][2] - 2*pow(p[4][1], 2)*p[5][0]*p[7][2] + 2*pow(p[6][1], 2)*p[5][0]*p[7][2] + 2*p[4][0]*p[4][1]*p[5][1]*p[7][2] - p[4][1]*p[5][0]*p[5][1]*p[7][2] + pow(p[2][1], 2)*p[6][0]*p[7][2] + pow(p[3][1], 2)*p[6][0]*p[7][2] - pow(p[4][1], 2)*p[6][0]*p[7][2] - pow(p[5][1], 2)*p[6][0]*p[7][2] + p[2][1]*p[3][1]*p[6][0]*p[7][2] - p[4][1]*p[5][1]*p[6][0]*p[7][2] - p[2][0]*p[2][1]*p[6][1]*p[7][2] - p[2][1]*p[3][0]*p[6][1]*p[7][2] - p[3][0]*p[3][1]*p[6][1]*p[7][2] + p[4][0]*p[4][1]*p[6][1]*p[7][2] + p[4][0]*p[5][1]*p[6][1]*p[7][2] + p[5][0]*p[5][1]*p[6][1]*p[7][2] + 2*p[2][1]*p[6][0]*p[6][1]*p[7][2] + p[3][1]*p[6][0]*p[6][1]*p[7][2] - p[4][1]*p[6][0]*p[6][1]*p[7][2] - 2*p[5][1]*p[6][0]*p[6][1]*p[7][2] - p[2][1]*p[3][0]*p[7][1]*p[7][2] - p[0][0]*p[3][1]*p[7][1]*p[7][2] + p[2][0]*p[3][1]*p[7][1]*p[7][2] - 2*p[3][1]*p[4][0]*p[7][1]*p[7][2] + p[0][0]*p[4][1]*p[7][1]*p[7][2] + 2*p[3][0]*p[4][1]*p[7][1]*p[7][2] - p[4][1]*p[5][0]*p[7][1]*p[7][2] + p[4][0]*p[5][1]*p[7][1]*p[7][2] + p[2][1]*p[6][0]*p[7][1]*p[7][2] + 2*p[3][1]*p[6][0]*p[7][1]*p[7][2] - 2*p[4][1]*p[6][0]*p[7][1]*p[7][2] - p[5][1]*p[6][0]*p[7][1]*p[7][2] - p[2][0]*p[6][1]*p[7][1]*p[7][2] - 2*p[3][0]*p[6][1]*p[7][1]*p[7][2] + 2*p[4][0]*p[6][1]*p[7][1]*p[7][2] + p[5][0]*p[6][1]*p[7][1]*p[7][2] + pow(p[0][1], 2)*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - 2*p[3][2]*p[4][0] + 2*p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - 2*p[3][0] + 2*p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + 2*p[3][2] - 2*p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2]) + p[0][1]*(p[0][0]*p[1][2]*p[2][1] - p[1][2]*p[2][0]*p[2][1] + p[1][0]*p[2][1]*p[2][2] - p[2][1]*p[2][2]*p[3][0] + 2*p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] - p[1][2]*p[3][0]*p[3][1] - 2*p[2][2]*p[3][0]*p[3][1] - p[0][0]*p[2][1]*p[3][2] + p[2][0]*p[2][1]*p[3][2] + p[1][0]*p[3][1]*p[3][2] + 2*p[2][0]*p[3][1]*p[3][2] - p[3][1]*p[3][2]*p[4][0] - 2*p[0][0]*p[1][2]*p[4][1] + 2*p[0][0]*p[3][2]*p[4][1] + p[1][2]*p[4][0]*p[4][1] - p[3][2]*p[4][0]*p[4][1] - 2*p[0][0]*p[3][1]*p[4][2] + p[3][0]*p[3][1]*p[4][2] - p[1][0]*p[4][1]*p[4][2] + p[3][0]*p[4][1]*p[4][2] + p[1][2]*p[4][1]*p[5][0] - 2*p[4][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][2]*p[5][0]*p[5][1] - p[4][2]*p[5][0]*p[5][1] - p[0][0]*p[4][1]*p[5][2] - p[1][0]*p[4][1]*p[5][2] + 2*p[4][0]*p[4][1]*p[5][2] - p[1][0]*p[5][1]*p[5][2] + p[4][0]*p[5][1]*p[5][2] + p[1][1]*(2*p[1][0]*p[2][2] - p[2][2]*p[3][0] + p[1][0]*p[3][2] + p[2][0]*p[3][2] - p[1][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-2*p[2][0] - p[3][0] + p[4][0] + 2*p[5][0]) - 2*p[1][0]*p[5][2] + p[4][0]*p[5][2] + p[0][0]*(-p[2][2] - 2*p[3][2] + 2*p[4][2] + p[5][2])) - 2*p[3][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] + 2*p[4][1]*p[4][2]*p[7][0] + p[0][0]*p[3][2]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][2]*p[7][0]*p[7][1] + p[4][2]*p[7][0]*p[7][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + 2*p[3][1]*p[4][0] - 2*p[3][0]*p[4][1] + p[1][1]*(p[2][0] + 2*p[3][0] - 2*p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - 2*p[3][1] + 2*p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[0][0]*p[3][1]*p[7][2] + 2*p[3][0]*p[3][1]*p[7][2] - p[3][1]*p[4][0]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][0]*p[4][1]*p[7][2] - 2*p[4][0]*p[4][1]*p[7][2] + p[3][0]*p[7][1]*p[7][2] - p[4][0]*p[7][1]*p[7][2]))/72;

            result[3][0] =
              (-(pow(p[2][2], 2)*p[0][0]*p[1][1]) - pow(p[3][2], 2)*p[0][0]*p[1][1] + pow(p[4][2], 2)*p[0][0]*p[1][1] + pow(p[5][2], 2)*p[0][0]*p[1][1] + pow(p[3][2], 2)*p[1][1]*p[2][0] - pow(p[5][2], 2)*p[1][1]*p[2][0] - pow(p[6][2], 2)*p[1][1]*p[2][0] + 2*pow(p[1][2], 2)*p[0][0]*p[2][1] - 2*pow(p[3][2], 2)*p[0][0]*p[2][1] - pow(p[3][2], 2)*p[1][0]*p[2][1] + pow(p[5][2], 2)*p[1][0]*p[2][1] + pow(p[6][2], 2)*p[1][0]*p[2][1] - 2*p[0][0]*p[1][1]*p[1][2]*p[2][2] + p[0][0]*p[1][2]*p[2][1]*p[2][2] - 2*pow(p[2][2], 2)*p[1][1]*p[3][0] + pow(p[1][2], 2)*p[2][1]*p[3][0] - pow(p[6][2], 2)*p[2][1]*p[3][0] - pow(p[7][2], 2)*p[2][1]*p[3][0] - p[1][1]*p[1][2]*p[2][2]*p[3][0] + 2*p[1][2]*p[2][1]*p[2][2]*p[3][0] + pow(p[1][2], 2)*p[0][0]*p[3][1] + pow(p[2][2], 2)*p[0][0]*p[3][1] - pow(p[4][2], 2)*p[0][0]*p[3][1] - pow(p[7][2], 2)*p[0][0]*p[3][1] + 2*pow(p[2][2], 2)*p[1][0]*p[3][1] - pow(p[1][2], 2)*p[2][0]*p[3][1] + pow(p[6][2], 2)*p[2][0]*p[3][1] + pow(p[7][2], 2)*p[2][0]*p[3][1] + p[0][0]*p[1][2]*p[2][2]*p[3][1] + p[1][0]*p[1][2]*p[2][2]*p[3][1] - 2*p[1][2]*p[2][0]*p[2][2]*p[3][1] - p[0][0]*p[1][1]*p[1][2]*p[3][2] + p[1][1]*p[1][2]*p[2][0]*p[3][2] - p[1][0]*p[1][2]*p[2][1]*p[3][2] - p[0][0]*p[1][1]*p[2][2]*p[3][2] + 2*p[1][1]*p[2][0]*p[2][2]*p[3][2] - p[0][0]*p[2][1]*p[2][2]*p[3][2] - 2*p[1][0]*p[2][1]*p[2][2]*p[3][2] + p[1][2]*p[2][1]*p[3][0]*p[3][2] - p[1][1]*p[2][2]*p[3][0]*p[3][2] + p[0][0]*p[1][2]*p[3][1]*p[3][2] - p[1][2]*p[2][0]*p[3][1]*p[3][2] + 2*p[0][0]*p[2][2]*p[3][1]*p[3][2] + p[1][0]*p[2][2]*p[3][1]*p[3][2] + 2*pow(p[5][2], 2)*p[1][1]*p[4][0] - 2*pow(p[7][2], 2)*p[3][1]*p[4][0] - pow(p[1][2], 2)*p[0][0]*p[4][1] + pow(p[3][2], 2)*p[0][0]*p[4][1] - pow(p[5][2], 2)*p[0][0]*p[4][1] + pow(p[7][2], 2)*p[0][0]*p[4][1] - 2*pow(p[5][2], 2)*p[1][0]*p[4][1] + 2*pow(p[7][2], 2)*p[3][0]*p[4][1] + p[0][0]*p[1][1]*p[1][2]*p[4][2] - p[0][0]*p[3][1]*p[3][2]*p[4][2] - p[0][0]*p[1][2]*p[4][1]*p[4][2] + p[0][0]*p[3][2]*p[4][1]*p[4][2] + pow(p[2][2], 2)*p[1][1]*p[5][0] - pow(p[4][2], 2)*p[1][1]*p[5][0] + pow(p[6][2], 2)*p[1][1]*p[5][0] - 2*pow(p[1][2], 2)*p[2][1]*p[5][0] + 2*pow(p[6][2], 2)*p[2][1]*p[5][0] + 2*p[1][1]*p[1][2]*p[2][2]*p[5][0] - p[1][2]*p[2][1]*p[2][2]*p[5][0] + pow(p[1][2], 2)*p[4][1]*p[5][0] - pow(p[6][2], 2)*p[4][1]*p[5][0] - pow(p[7][2], 2)*p[4][1]*p[5][0] - p[1][1]*p[1][2]*p[4][2]*p[5][0] + p[1][2]*p[4][1]*p[4][2]*p[5][0] - 2*pow(p[1][2], 2)*p[0][0]*p[5][1] + 2*pow(p[4][2], 2)*p[0][0]*p[5][1] - pow(p[2][2], 2)*p[1][0]*p[5][1] + pow(p[4][2], 2)*p[1][0]*p[5][1] - pow(p[6][2], 2)*p[1][0]*p[5][1] + 2*pow(p[1][2], 2)*p[2][0]*p[5][1] - 2*pow(p[6][2], 2)*p[2][0]*p[5][1] - 2*p[1][0]*p[1][2]*p[2][2]*p[5][1] + p[1][2]*p[2][0]*p[2][2]*p[5][1] - pow(p[1][2], 2)*p[4][0]*p[5][1] + pow(p[6][2], 2)*p[4][0]*p[5][1] + pow(p[7][2], 2)*p[4][0]*p[5][1] + p[1][0]*p[1][2]*p[4][2]*p[5][1] - p[1][2]*p[4][0]*p[4][2]*p[5][1] + 2*p[0][0]*p[1][1]*p[1][2]*p[5][2] - 2*p[1][1]*p[1][2]*p[2][0]*p[5][2] + 2*p[1][0]*p[1][2]*p[2][1]*p[5][2] - p[1][1]*p[2][0]*p[2][2]*p[5][2] + p[1][0]*p[2][1]*p[2][2]*p[5][2] + p[1][1]*p[1][2]*p[4][0]*p[5][2] - p[0][0]*p[1][2]*p[4][1]*p[5][2] - p[1][0]*p[1][2]*p[4][1]*p[5][2] + p[0][0]*p[1][1]*p[4][2]*p[5][2] + p[1][1]*p[4][0]*p[4][2]*p[5][2] - 2*p[0][0]*p[4][1]*p[4][2]*p[5][2] - p[1][0]*p[4][1]*p[4][2]*p[5][2] - p[1][2]*p[2][1]*p[5][0]*p[5][2] + p[1][1]*p[2][2]*p[5][0]*p[5][2] + 2*p[1][2]*p[4][1]*p[5][0]*p[5][2] - 2*p[1][1]*p[4][2]*p[5][0]*p[5][2] - p[0][0]*p[1][2]*p[5][1]*p[5][2] + p[1][2]*p[2][0]*p[5][1]*p[5][2] - p[1][0]*p[2][2]*p[5][1]*p[5][2] - 2*p[1][2]*p[4][0]*p[5][1]*p[5][2] + p[0][0]*p[4][2]*p[5][1]*p[5][2] + 2*p[1][0]*p[4][2]*p[5][1]*p[5][2] + 2*pow(p[2][2], 2)*p[1][1]*p[6][0] - 2*pow(p[5][2], 2)*p[1][1]*p[6][0] - pow(p[1][2], 2)*p[2][1]*p[6][0] + pow(p[3][2], 2)*p[2][1]*p[6][0] - pow(p[5][2], 2)*p[2][1]*p[6][0] + pow(p[7][2], 2)*p[2][1]*p[6][0] + p[1][1]*p[1][2]*p[2][2]*p[6][0] - 2*p[1][2]*p[2][1]*p[2][2]*p[6][0] - 2*pow(p[2][2], 2)*p[3][1]*p[6][0] + 2*pow(p[7][2], 2)*p[3][1]*p[6][0] + 2*p[2][1]*p[2][2]*p[3][2]*p[6][0] - p[2][2]*p[3][1]*p[3][2]*p[6][0] + 2*pow(p[5][2], 2)*p[4][1]*p[6][0] - 2*pow(p[7][2], 2)*p[4][1]*p[6][0] + pow(p[1][2], 2)*p[5][1]*p[6][0] + pow(p[2][2], 2)*p[5][1]*p[6][0] - pow(p[4][2], 2)*p[5][1]*p[6][0] - pow(p[7][2], 2)*p[5][1]*p[6][0] + p[1][2]*p[2][2]*p[5][1]*p[6][0] - p[1][1]*p[1][2]*p[5][2]*p[6][0] - p[1][2]*p[2][1]*p[5][2]*p[6][0] - p[2][1]*p[2][2]*p[5][2]*p[6][0] + p[4][1]*p[4][2]*p[5][2]*p[6][0] + 2*p[1][2]*p[5][1]*p[5][2]*p[6][0] + p[2][2]*p[5][1]*p[5][2]*p[6][0] - 2*p[4][2]*p[5][1]*p[5][2]*p[6][0] - 2*pow(p[2][2], 2)*p[1][0]*p[6][1] + 2*pow(p[5][2], 2)*p[1][0]*p[6][1] + pow(p[1][2], 2)*p[2][0]*p[6][1] - pow(p[3][2], 2)*p[2][0]*p[6][1] + pow(p[5][2], 2)*p[2][0]*p[6][1] - pow(p[7][2], 2)*p[2][0]*p[6][1] - p[1][0]*p[1][2]*p[2][2]*p[6][1] + 2*p[1][2]*p[2][0]*p[2][2]*p[6][1] + 2*pow(p[2][2], 2)*p[3][0]*p[6][1] - 2*pow(p[7][2], 2)*p[3][0]*p[6][1] - 2*p[2][0]*p[2][2]*p[3][2]*p[6][1] + p[2][2]*p[3][0]*p[3][2]*p[6][1] - 2*pow(p[5][2], 2)*p[4][0]*p[6][1] + 2*pow(p[7][2], 2)*p[4][0]*p[6][1] - pow(p[1][2], 2)*p[5][0]*p[6][1] - pow(p[2][2], 2)*p[5][0]*p[6][1] + pow(p[4][2], 2)*p[5][0]*p[6][1] + pow(p[7][2], 2)*p[5][0]*p[6][1] - p[1][2]*p[2][2]*p[5][0]*p[6][1] + p[1][0]*p[1][2]*p[5][2]*p[6][1] + p[1][2]*p[2][0]*p[5][2]*p[6][1] + p[2][0]*p[2][2]*p[5][2]*p[6][1] - p[4][0]*p[4][2]*p[5][2]*p[6][1] - 2*p[1][2]*p[5][0]*p[5][2]*p[6][1] - p[2][2]*p[5][0]*p[5][2]*p[6][1] + 2*p[4][2]*p[5][0]*p[5][2]*p[6][1] - p[1][1]*p[1][2]*p[2][0]*p[6][2] + p[1][0]*p[1][2]*p[2][1]*p[6][2] - 2*p[1][1]*p[2][0]*p[2][2]*p[6][2] + 2*p[1][0]*p[2][1]*p[2][2]*p[6][2] - 2*p[2][1]*p[2][2]*p[3][0]*p[6][2] + 2*p[2][0]*p[2][2]*p[3][1]*p[6][2] - p[2][1]*p[3][0]*p[3][2]*p[6][2] + p[2][0]*p[3][1]*p[3][2]*p[6][2] + p[1][1]*p[1][2]*p[5][0]*p[6][2] + p[1][1]*p[2][2]*p[5][0]*p[6][2] + p[2][1]*p[2][2]*p[5][0]*p[6][2] - p[4][1]*p[4][2]*p[5][0]*p[6][2] - p[1][0]*p[1][2]*p[5][1]*p[6][2] - p[1][0]*p[2][2]*p[5][1]*p[6][2] - p[2][0]*p[2][2]*p[5][1]*p[6][2] + p[4][0]*p[4][2]*p[5][1]*p[6][2] - p[1][1]*p[2][0]*p[5][2]*p[6][2] + p[1][0]*p[2][1]*p[5][2]*p[6][2] + 2*p[1][1]*p[5][0]*p[5][2]*p[6][2] + p[2][1]*p[5][0]*p[5][2]*p[6][2] - 2*p[4][1]*p[5][0]*p[5][2]*p[6][2] - 2*p[1][0]*p[5][1]*p[5][2]*p[6][2] - p[2][0]*p[5][1]*p[5][2]*p[6][2] + 2*p[4][0]*p[5][1]*p[5][2]*p[6][2] - p[1][2]*p[2][1]*p[6][0]*p[6][2] + p[1][1]*p[2][2]*p[6][0]*p[6][2] - p[2][2]*p[3][1]*p[6][0]*p[6][2] + p[2][1]*p[3][2]*p[6][0]*p[6][2] + p[1][2]*p[5][1]*p[6][0]*p[6][2] + 2*p[2][2]*p[5][1]*p[6][0]*p[6][2] - p[4][2]*p[5][1]*p[6][0]*p[6][2] - p[1][1]*p[5][2]*p[6][0]*p[6][2] - 2*p[2][1]*p[5][2]*p[6][0]*p[6][2] + p[4][1]*p[5][2]*p[6][0]*p[6][2] + p[1][2]*p[2][0]*p[6][1]*p[6][2] - p[1][0]*p[2][2]*p[6][1]*p[6][2] + p[2][2]*p[3][0]*p[6][1]*p[6][2] - p[2][0]*p[3][2]*p[6][1]*p[6][2] - p[1][2]*p[5][0]*p[6][1]*p[6][2] - 2*p[2][2]*p[5][0]*p[6][1]*p[6][2] + p[4][2]*p[5][0]*p[6][1]*p[6][2] + p[1][0]*p[5][2]*p[6][1]*p[6][2] + 2*p[2][0]*p[5][2]*p[6][1]*p[6][2] - p[4][0]*p[5][2]*p[6][1]*p[6][2] + 2*pow(p[3][2], 2)*p[2][1]*p[7][0] - 2*pow(p[6][2], 2)*p[2][1]*p[7][0] - pow(p[2][2], 2)*p[3][1]*p[7][0] + pow(p[4][2], 2)*p[3][1]*p[7][0] - pow(p[6][2], 2)*p[3][1]*p[7][0] + p[2][1]*p[2][2]*p[3][2]*p[7][0] - 2*p[2][2]*p[3][1]*p[3][2]*p[7][0] - pow(p[3][2], 2)*p[4][1]*p[7][0] + pow(p[5][2], 2)*p[4][1]*p[7][0] + pow(p[6][2], 2)*p[4][1]*p[7][0] + p[3][1]*p[3][2]*p[4][2]*p[7][0] - p[3][2]*p[4][1]*p[4][2]*p[7][0] - 2*pow(p[4][2], 2)*p[5][1]*p[7][0] + 2*pow(p[6][2], 2)*p[5][1]*p[7][0] + 2*p[4][1]*p[4][2]*p[5][2]*p[7][0] - p[4][2]*p[5][1]*p[5][2]*p[7][0] + pow(p[2][2], 2)*p[6][1]*p[7][0] + pow(p[3][2], 2)*p[6][1]*p[7][0] - pow(p[4][2], 2)*p[6][1]*p[7][0] - pow(p[5][2], 2)*p[6][1]*p[7][0] + p[2][2]*p[3][2]*p[6][1]*p[7][0] - p[4][2]*p[5][2]*p[6][1]*p[7][0] - p[2][1]*p[2][2]*p[6][2]*p[7][0] - p[2][2]*p[3][1]*p[6][2]*p[7][0] - p[3][1]*p[3][2]*p[6][2]*p[7][0] + p[4][1]*p[4][2]*p[6][2]*p[7][0] + p[4][1]*p[5][2]*p[6][2]*p[7][0] + p[5][1]*p[5][2]*p[6][2]*p[7][0] + 2*p[2][2]*p[6][1]*p[6][2]*p[7][0] + p[3][2]*p[6][1]*p[6][2]*p[7][0] - p[4][2]*p[6][1]*p[6][2]*p[7][0] - 2*p[5][2]*p[6][1]*p[6][2]*p[7][0] + 2*pow(p[3][2], 2)*p[0][0]*p[7][1] - 2*pow(p[4][2], 2)*p[0][0]*p[7][1] - 2*pow(p[3][2], 2)*p[2][0]*p[7][1] + 2*pow(p[6][2], 2)*p[2][0]*p[7][1] + pow(p[2][2], 2)*p[3][0]*p[7][1] - pow(p[4][2], 2)*p[3][0]*p[7][1] + pow(p[6][2], 2)*p[3][0]*p[7][1] - p[2][0]*p[2][2]*p[3][2]*p[7][1] + 2*p[2][2]*p[3][0]*p[3][2]*p[7][1] + pow(p[3][2], 2)*p[4][0]*p[7][1] - pow(p[5][2], 2)*p[4][0]*p[7][1] - pow(p[6][2], 2)*p[4][0]*p[7][1] - p[3][0]*p[3][2]*p[4][2]*p[7][1] + p[3][2]*p[4][0]*p[4][2]*p[7][1] + 2*pow(p[4][2], 2)*p[5][0]*p[7][1] - 2*pow(p[6][2], 2)*p[5][0]*p[7][1] - 2*p[4][0]*p[4][2]*p[5][2]*p[7][1] + p[4][2]*p[5][0]*p[5][2]*p[7][1] - pow(p[2][2], 2)*p[6][0]*p[7][1] - pow(p[3][2], 2)*p[6][0]*p[7][1] + pow(p[4][2], 2)*p[6][0]*p[7][1] + pow(p[5][2], 2)*p[6][0]*p[7][1] - p[2][2]*p[3][2]*p[6][0]*p[7][1] + p[4][2]*p[5][2]*p[6][0]*p[7][1] + p[2][0]*p[2][2]*p[6][2]*p[7][1] + p[2][2]*p[3][0]*p[6][2]*p[7][1] + p[3][0]*p[3][2]*p[6][2]*p[7][1] - p[4][0]*p[4][2]*p[6][2]*p[7][1] - p[4][0]*p[5][2]*p[6][2]*p[7][1] - p[5][0]*p[5][2]*p[6][2]*p[7][1] - 2*p[2][2]*p[6][0]*p[6][2]*p[7][1] - p[3][2]*p[6][0]*p[6][2]*p[7][1] + p[4][2]*p[6][0]*p[6][2]*p[7][1] + 2*p[5][2]*p[6][0]*p[6][2]*p[7][1] + pow(p[0][2], 2)*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + 2*p[3][1]*p[4][0] - 2*p[3][0]*p[4][1] + p[1][1]*(p[2][0] + 2*p[3][0] - 2*p[4][0] - p[5][0]) + p[4][1]*p[5][0] - p[4][0]*p[5][1] + p[1][0]*(-p[2][1] - 2*p[3][1] + 2*p[4][1] + p[5][1]) + p[3][1]*p[7][0] - p[4][1]*p[7][0] - p[3][0]*p[7][1] + p[4][0]*p[7][1]) - p[2][1]*p[2][2]*p[3][0]*p[7][2] + p[2][0]*p[2][2]*p[3][1]*p[7][2] - 2*p[2][1]*p[3][0]*p[3][2]*p[7][2] - 2*p[0][0]*p[3][1]*p[3][2]*p[7][2] + 2*p[2][0]*p[3][1]*p[3][2]*p[7][2] - p[3][1]*p[3][2]*p[4][0]*p[7][2] + p[0][0]*p[3][2]*p[4][1]*p[7][2] + p[3][0]*p[3][2]*p[4][1]*p[7][2] - p[0][0]*p[3][1]*p[4][2]*p[7][2] - p[3][1]*p[4][0]*p[4][2]*p[7][2] + 2*p[0][0]*p[4][1]*p[4][2]*p[7][2] + p[3][0]*p[4][1]*p[4][2]*p[7][2] - 2*p[4][1]*p[4][2]*p[5][0]*p[7][2] + 2*p[4][0]*p[4][2]*p[5][1]*p[7][2] - p[4][1]*p[5][0]*p[5][2]*p[7][2] + p[4][0]*p[5][1]*p[5][2]*p[7][2] + p[2][1]*p[2][2]*p[6][0]*p[7][2] + p[2][1]*p[3][2]*p[6][0]*p[7][2] + p[3][1]*p[3][2]*p[6][0]*p[7][2] - p[4][1]*p[4][2]*p[6][0]*p[7][2] - p[4][2]*p[5][1]*p[6][0]*p[7][2] - p[5][1]*p[5][2]*p[6][0]*p[7][2] - p[2][0]*p[2][2]*p[6][1]*p[7][2] - p[2][0]*p[3][2]*p[6][1]*p[7][2] - p[3][0]*p[3][2]*p[6][1]*p[7][2] + p[4][0]*p[4][2]*p[6][1]*p[7][2] + p[4][2]*p[5][0]*p[6][1]*p[7][2] + p[5][0]*p[5][2]*p[6][1]*p[7][2] - p[2][1]*p[3][0]*p[6][2]*p[7][2] + p[2][0]*p[3][1]*p[6][2]*p[7][2] - p[4][1]*p[5][0]*p[6][2]*p[7][2] + p[4][0]*p[5][1]*p[6][2]*p[7][2] + 2*p[2][1]*p[6][0]*p[6][2]*p[7][2] + p[3][1]*p[6][0]*p[6][2]*p[7][2] - p[4][1]*p[6][0]*p[6][2]*p[7][2] - 2*p[5][1]*p[6][0]*p[6][2]*p[7][2] - 2*p[2][0]*p[6][1]*p[6][2]*p[7][2] - p[3][0]*p[6][1]*p[6][2]*p[7][2] + p[4][0]*p[6][1]*p[6][2]*p[7][2] + 2*p[5][0]*p[6][1]*p[6][2]*p[7][2] - p[2][2]*p[3][1]*p[7][0]*p[7][2] + p[2][1]*p[3][2]*p[7][0]*p[7][2] - 2*p[3][2]*p[4][1]*p[7][0]*p[7][2] + 2*p[3][1]*p[4][2]*p[7][0]*p[7][2] - p[4][2]*p[5][1]*p[7][0]*p[7][2] + p[4][1]*p[5][2]*p[7][0]*p[7][2] + p[2][2]*p[6][1]*p[7][0]*p[7][2] + 2*p[3][2]*p[6][1]*p[7][0]*p[7][2] - 2*p[4][2]*p[6][1]*p[7][0]*p[7][2] - p[5][2]*p[6][1]*p[7][0]*p[7][2] - p[2][1]*p[6][2]*p[7][0]*p[7][2] - 2*p[3][1]*p[6][2]*p[7][0]*p[7][2] + 2*p[4][1]*p[6][2]*p[7][0]*p[7][2] + p[5][1]*p[6][2]*p[7][0]*p[7][2] + p[2][2]*p[3][0]*p[7][1]*p[7][2] + p[0][0]*p[3][2]*p[7][1]*p[7][2] - p[2][0]*p[3][2]*p[7][1]*p[7][2] + 2*p[3][2]*p[4][0]*p[7][1]*p[7][2] - p[0][0]*p[4][2]*p[7][1]*p[7][2] - 2*p[3][0]*p[4][2]*p[7][1]*p[7][2] + p[4][2]*p[5][0]*p[7][1]*p[7][2] - p[4][0]*p[5][2]*p[7][1]*p[7][2] - p[2][2]*p[6][0]*p[7][1]*p[7][2] - 2*p[3][2]*p[6][0]*p[7][1]*p[7][2] + 2*p[4][2]*p[6][0]*p[7][1]*p[7][2] + p[5][2]*p[6][0]*p[7][1]*p[7][2] + p[2][0]*p[6][2]*p[7][1]*p[7][2] + 2*p[3][0]*p[6][2]*p[7][1]*p[7][2] - 2*p[4][0]*p[6][2]*p[7][1]*p[7][2] - p[5][0]*p[6][2]*p[7][1]*p[7][2] + p[0][1]*(2*pow(p[3][2], 2)*p[2][0] - pow(p[2][2], 2)*p[3][0] + pow(p[4][2], 2)*p[3][0] + pow(p[7][2], 2)*p[3][0] + p[2][0]*p[2][2]*p[3][2] - 2*p[2][2]*p[3][0]*p[3][2] - pow(p[3][2], 2)*p[4][0] + pow(p[5][2], 2)*p[4][0] - pow(p[7][2], 2)*p[4][0] + p[3][0]*p[3][2]*p[4][2] - p[3][2]*p[4][0]*p[4][2] - 2*pow(p[4][2], 2)*p[5][0] + pow(p[1][2], 2)*(-2*p[2][0] - p[3][0] + p[4][0] + 2*p[5][0]) + 2*p[4][0]*p[4][2]*p[5][2] - p[4][2]*p[5][0]*p[5][2] + p[1][0]*(pow(p[2][2], 2) + pow(p[3][2], 2) - pow(p[4][2], 2) - pow(p[5][2], 2) + p[2][2]*p[3][2] - p[4][2]*p[5][2]) + p[1][2]*(-(p[2][0]*p[2][2]) - p[2][2]*p[3][0] - p[3][0]*p[3][2] + p[4][0]*p[4][2] + p[1][0]*(2*p[2][2] + p[3][2] - p[4][2] - 2*p[5][2]) + p[4][0]*p[5][2] + p[5][0]*p[5][2]) - 2*pow(p[3][2], 2)*p[7][0] + 2*pow(p[4][2], 2)*p[7][0] + 2*p[3][0]*p[3][2]*p[7][2] - p[3][2]*p[4][0]*p[7][2] + p[3][0]*p[4][2]*p[7][2] - 2*p[4][0]*p[4][2]*p[7][2] - p[3][2]*p[7][0]*p[7][2] + p[4][2]*p[7][0]*p[7][2]) + p[0][2]*(p[0][0]*p[1][2]*p[2][1] - 2*p[1][0]*p[1][2]*p[2][1] - p[1][0]*p[2][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] + p[2][1]*p[2][2]*p[3][0] + 2*p[0][0]*p[1][2]*p[3][1] - p[1][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[2][0]*p[2][2]*p[3][1] - p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] + 2*p[2][1]*p[3][0]*p[3][2] - p[1][0]*p[3][1]*p[3][2] - 2*p[2][0]*p[3][1]*p[3][2] + p[3][1]*p[3][2]*p[4][0] - 2*p[0][0]*p[1][2]*p[4][1] + p[1][0]*p[1][2]*p[4][1] + 2*p[0][0]*p[3][2]*p[4][1] - p[3][0]*p[3][2]*p[4][1] - 2*p[0][0]*p[3][1]*p[4][2] + p[3][1]*p[4][0]*p[4][2] + p[1][0]*p[4][1]*p[4][2] - p[3][0]*p[4][1]*p[4][2] + p[1][2]*p[4][1]*p[5][0] + 2*p[4][1]*p[4][2]*p[5][0] - p[0][0]*p[1][2]*p[5][1] + 2*p[1][0]*p[1][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] + p[0][0]*p[4][2]*p[5][1] + p[1][0]*p[4][2]*p[5][1] - 2*p[4][0]*p[4][2]*p[5][1] - p[0][0]*p[4][1]*p[5][2] + p[4][1]*p[5][0]*p[5][2] + p[1][0]*p[5][1]*p[5][2] - p[4][0]*p[5][1]*p[5][2] + p[1][1]*(p[2][0]*p[2][2] + p[2][0]*p[3][2] + p[3][0]*p[3][2] - p[4][0]*p[4][2] + p[1][2]*(2*p[2][0] + p[3][0] - p[4][0] - 2*p[5][0]) - p[4][2]*p[5][0] - p[5][0]*p[5][2] + p[0][0]*(-p[2][2] - 2*p[3][2] + 2*p[4][2] + p[5][2])) + 2*p[3][1]*p[3][2]*p[7][0] - p[3][2]*p[4][1]*p[7][0] + p[3][1]*p[4][2]*p[7][0] - 2*p[4][1]*p[4][2]*p[7][0] + p[0][0]*p[3][2]*p[7][1] - 2*p[3][0]*p[3][2]*p[7][1] + p[3][2]*p[4][0]*p[7][1] - p[0][0]*p[4][2]*p[7][1] - p[3][0]*p[4][2]*p[7][1] + 2*p[4][0]*p[4][2]*p[7][1] - p[0][0]*p[3][1]*p[7][2] + p[0][0]*p[4][1]*p[7][2] + p[3][1]*p[7][0]*p[7][2] - p[4][1]*p[7][0]*p[7][2] - p[3][0]*p[7][1]*p[7][2] + p[4][0]*p[7][1]*p[7][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - 2*p[3][2]*p[4][0] + 2*p[3][0]*p[4][2] - p[4][2]*p[5][0] + p[1][2]*(-p[2][0] - 2*p[3][0] + 2*p[4][0] + p[5][0]) + p[1][0]*(p[2][2] + 2*p[3][2] - 2*p[4][2] - p[5][2]) + p[4][0]*p[5][2] - p[3][2]*p[7][0] + p[4][2]*p[7][0] + p[3][0]*p[7][2] - p[4][0]*p[7][2])))/72;

            return result;
          }
        };

      // 3D Wedge

      template<>
        class Operator<VolumeIntegral, Basis::Unitary, Wedge>:
          public OperatorBase<VolumeIntegral, Basis::Unitary, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Hexahedron);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Unitary, dim>::function_size;
          static constexpr size_t point_size=6;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

          result[0][0]=
(2*p[0][0]*p[1][2]*p[2][1] - 2*p[0][0]*p[1][1]*p[2][2] - p[0][0]*p[1][2]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[0][0]*p[1][1]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][2]*p[2][1]*p[4][0] + p[1][1]*p[2][2]*p[4][0] + p[1][2]*p[3][1]*p[4][0] - p[1][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + p[1][2]*p[2][0]*p[4][1] - p[1][0]*p[2][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[1][1]*p[2][0]*p[4][2] + p[1][0]*p[2][1]*p[4][2] + p[1][1]*p[3][0]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] + p[2][2]*p[4][1]*p[5][0] - 2*p[3][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[2][1]*p[4][2]*p[5][0] + 2*p[3][1]*p[4][2]*p[5][0] + p[1][2]*p[2][0]*p[5][1] + p[0][0]*p[2][2]*p[5][1] - p[1][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][0]*p[3][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] - p[2][2]*p[4][0]*p[5][1] + 2*p[3][2]*p[4][0]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[2][0]*p[4][2]*p[5][1] - 2*p[3][0]*p[4][2]*p[5][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][1]*(2*p[2][0] - p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-2*p[2][1] + p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[1][1]*p[2][0]*p[5][2] - p[0][0]*p[2][1]*p[5][2] + p[1][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][0]*p[3][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] + p[2][1]*p[4][0]*p[5][2] - 2*p[3][1]*p[4][0]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[2][0]*p[4][1]*p[5][2] + 2*p[3][0]*p[4][1]*p[5][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-2*p[2][0] + p[3][0] + p[4][0]) + p[1][0]*(2*p[2][2] - p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2]))/12;

            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Linear, Wedge>:
          public OperatorBase<VolumeIntegral, Basis::Linear, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Wedge);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Linear, dim>::function_size;
          static constexpr size_t point_size=6;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

          result[0][0]=
(2*p[0][0]*p[1][2]*p[2][1] - 2*p[0][0]*p[1][1]*p[2][2] - p[0][0]*p[1][2]*p[3][1] + p[0][0]*p[2][2]*p[3][1] + p[0][0]*p[1][1]*p[3][2] - p[0][0]*p[2][1]*p[3][2] - p[1][2]*p[2][1]*p[4][0] + p[1][1]*p[2][2]*p[4][0] + p[1][2]*p[3][1]*p[4][0] - p[1][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + p[1][2]*p[2][0]*p[4][1] - p[1][0]*p[2][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][0]*p[3][2]*p[4][1] + p[0][0]*p[1][1]*p[4][2] - p[1][1]*p[2][0]*p[4][2] + p[1][0]*p[2][1]*p[4][2] + p[1][1]*p[3][0]*p[4][2] - p[0][0]*p[3][1]*p[4][2] - p[1][0]*p[3][1]*p[4][2] - p[1][2]*p[2][1]*p[5][0] + p[1][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] + p[1][2]*p[4][1]*p[5][0] + p[2][2]*p[4][1]*p[5][0] - 2*p[3][2]*p[4][1]*p[5][0] - p[1][1]*p[4][2]*p[5][0] - p[2][1]*p[4][2]*p[5][0] + 2*p[3][1]*p[4][2]*p[5][0] + p[1][2]*p[2][0]*p[5][1] + p[0][0]*p[2][2]*p[5][1] - p[1][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][0]*p[3][2]*p[5][1] - p[1][2]*p[4][0]*p[5][1] - p[2][2]*p[4][0]*p[5][1] + 2*p[3][2]*p[4][0]*p[5][1] + p[1][0]*p[4][2]*p[5][1] + p[2][0]*p[4][2]*p[5][1] - 2*p[3][0]*p[4][2]*p[5][1] + p[0][2]*(p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][1]*(2*p[2][0] - p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-2*p[2][1] + p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[1][1]*p[2][0]*p[5][2] - p[0][0]*p[2][1]*p[5][2] + p[1][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][0]*p[3][1]*p[5][2] + p[1][1]*p[4][0]*p[5][2] + p[2][1]*p[4][0]*p[5][2] - 2*p[3][1]*p[4][0]*p[5][2] - p[1][0]*p[4][1]*p[5][2] - p[2][0]*p[4][1]*p[5][2] + 2*p[3][0]*p[4][1]*p[5][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-2*p[2][0] + p[3][0] + p[4][0]) + p[1][0]*(2*p[2][2] - p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2]))/12;

          result[1][0]=
(-3*pow(p[2][0], 2)*p[0][1]*p[1][2] + pow(p[3][0], 2)*p[0][1]*p[1][2] + pow(p[4][0], 2)*p[0][1]*p[1][2] - 3*p[0][1]*p[1][0]*p[1][2]*p[2][0] - pow(p[4][0], 2)*p[1][2]*p[2][1] - pow(p[5][0], 2)*p[1][2]*p[2][1] + 3*pow(p[1][0], 2)*p[0][1]*p[2][2] - pow(p[3][0], 2)*p[0][1]*p[2][2] - pow(p[5][0], 2)*p[0][1]*p[2][2] + pow(p[4][0], 2)*p[1][1]*p[2][2] + pow(p[5][0], 2)*p[1][1]*p[2][2] + 3*p[0][1]*p[1][0]*p[2][0]*p[2][2] + p[0][1]*p[1][0]*p[1][2]*p[3][0] - p[0][1]*p[2][0]*p[2][2]*p[3][0] + 2*pow(p[4][0], 2)*p[1][2]*p[3][1] - 2*pow(p[5][0], 2)*p[2][2]*p[3][1] - pow(p[1][0], 2)*p[0][1]*p[3][2] + pow(p[2][0], 2)*p[0][1]*p[3][2] - pow(p[4][0], 2)*p[0][1]*p[3][2] + pow(p[5][0], 2)*p[0][1]*p[3][2] - 2*pow(p[4][0], 2)*p[1][1]*p[3][2] + 2*pow(p[5][0], 2)*p[2][1]*p[3][2] - p[0][1]*p[1][0]*p[3][0]*p[3][2] + p[0][1]*p[2][0]*p[3][0]*p[3][2] + 2*p[0][1]*p[1][0]*p[1][2]*p[4][0] - 2*p[1][0]*p[1][2]*p[2][1]*p[4][0] - p[1][2]*p[2][0]*p[2][1]*p[4][0] + 2*p[1][0]*p[1][1]*p[2][2]*p[4][0] + p[1][1]*p[2][0]*p[2][2]*p[4][0] + p[0][1]*p[1][2]*p[3][0]*p[4][0] + p[1][0]*p[1][2]*p[3][1]*p[4][0] + p[1][2]*p[3][0]*p[3][1]*p[4][0] - p[0][1]*p[1][0]*p[3][2]*p[4][0] - p[1][0]*p[1][1]*p[3][2]*p[4][0] - 2*p[0][1]*p[3][0]*p[3][2]*p[4][0] - p[1][1]*p[3][0]*p[3][2]*p[4][0] + pow(p[2][0], 2)*p[1][2]*p[4][1] - pow(p[3][0], 2)*p[1][2]*p[4][1] + pow(p[5][0], 2)*p[1][2]*p[4][1] + 2*p[1][0]*p[1][2]*p[2][0]*p[4][1] - 2*pow(p[1][0], 2)*p[2][2]*p[4][1] + 2*pow(p[5][0], 2)*p[2][2]*p[4][1] - p[1][0]*p[2][0]*p[2][2]*p[4][1] - p[1][0]*p[1][2]*p[3][0]*p[4][1] + pow(p[1][0], 2)*p[3][2]*p[4][1] - 3*pow(p[5][0], 2)*p[3][2]*p[4][1] + p[1][0]*p[3][0]*p[3][2]*p[4][1] + p[1][2]*p[2][0]*p[4][0]*p[4][1] - p[1][0]*p[2][2]*p[4][0]*p[4][1] - 2*p[1][2]*p[3][0]*p[4][0]*p[4][1] + 2*p[1][0]*p[3][2]*p[4][0]*p[4][1] - 2*pow(p[1][0], 2)*p[0][1]*p[4][2] + 2*pow(p[3][0], 2)*p[0][1]*p[4][2] - pow(p[2][0], 2)*p[1][1]*p[4][2] + pow(p[3][0], 2)*p[1][1]*p[4][2] - pow(p[5][0], 2)*p[1][1]*p[4][2] - 2*p[1][0]*p[1][1]*p[2][0]*p[4][2] + 2*pow(p[1][0], 2)*p[2][1]*p[4][2] - 2*pow(p[5][0], 2)*p[2][1]*p[4][2] + p[1][0]*p[2][0]*p[2][1]*p[4][2] + p[1][0]*p[1][1]*p[3][0]*p[4][2] - pow(p[1][0], 2)*p[3][1]*p[4][2] + 3*pow(p[5][0], 2)*p[3][1]*p[4][2] - p[1][0]*p[3][0]*p[3][1]*p[4][2] - p[0][1]*p[1][0]*p[4][0]*p[4][2] - p[1][1]*p[2][0]*p[4][0]*p[4][2] + p[1][0]*p[2][1]*p[4][0]*p[4][2] + p[0][1]*p[3][0]*p[4][0]*p[4][2] + 2*p[1][1]*p[3][0]*p[4][0]*p[4][2] - 2*p[1][0]*p[3][1]*p[4][0]*p[4][2] - p[1][0]*p[1][2]*p[2][1]*p[5][0] - 2*p[1][2]*p[2][0]*p[2][1]*p[5][0] + p[1][0]*p[1][1]*p[2][2]*p[5][0] - 2*p[0][1]*p[2][0]*p[2][2]*p[5][0] + 2*p[1][1]*p[2][0]*p[2][2]*p[5][0] - p[0][1]*p[2][2]*p[3][0]*p[5][0] - p[2][0]*p[2][2]*p[3][1]*p[5][0] - p[2][2]*p[3][0]*p[3][1]*p[5][0] + p[0][1]*p[2][0]*p[3][2]*p[5][0] + p[2][0]*p[2][1]*p[3][2]*p[5][0] + 2*p[0][1]*p[3][0]*p[3][2]*p[5][0] + p[2][1]*p[3][0]*p[3][2]*p[5][0] - p[1][2]*p[2][1]*p[4][0]*p[5][0] + p[1][1]*p[2][2]*p[4][0]*p[5][0] + p[1][0]*p[1][2]*p[4][1]*p[5][0] + p[1][2]*p[2][0]*p[4][1]*p[5][0] + p[2][0]*p[2][2]*p[4][1]*p[5][0] - 3*p[3][0]*p[3][2]*p[4][1]*p[5][0] + 2*p[1][2]*p[4][0]*p[4][1]*p[5][0] + p[2][2]*p[4][0]*p[4][1]*p[5][0] - 3*p[3][2]*p[4][0]*p[4][1]*p[5][0] - p[1][0]*p[1][1]*p[4][2]*p[5][0] - p[1][1]*p[2][0]*p[4][2]*p[5][0] - p[2][0]*p[2][1]*p[4][2]*p[5][0] + 3*p[3][0]*p[3][1]*p[4][2]*p[5][0] - 2*p[1][1]*p[4][0]*p[4][2]*p[5][0] - p[2][1]*p[4][0]*p[4][2]*p[5][0] + 3*p[3][1]*p[4][0]*p[4][2]*p[5][0] + 2*pow(p[2][0], 2)*p[1][2]*p[5][1] - 2*pow(p[4][0], 2)*p[1][2]*p[5][1] + p[1][0]*p[1][2]*p[2][0]*p[5][1] - pow(p[1][0], 2)*p[2][2]*p[5][1] + pow(p[3][0], 2)*p[2][2]*p[5][1] - pow(p[4][0], 2)*p[2][2]*p[5][1] - 2*p[1][0]*p[2][0]*p[2][2]*p[5][1] + p[2][0]*p[2][2]*p[3][0]*p[5][1] - pow(p[2][0], 2)*p[3][2]*p[5][1] + 3*pow(p[4][0], 2)*p[3][2]*p[5][1] - p[2][0]*p[3][0]*p[3][2]*p[5][1] - p[1][0]*p[1][2]*p[4][0]*p[5][1] - p[1][0]*p[2][2]*p[4][0]*p[5][1] - p[2][0]*p[2][2]*p[4][0]*p[5][1] + 3*p[3][0]*p[3][2]*p[4][0]*p[5][1] + pow(p[1][0], 2)*p[4][2]*p[5][1] + pow(p[2][0], 2)*p[4][2]*p[5][1] - 3*pow(p[3][0], 2)*p[4][2]*p[5][1] + p[1][0]*p[2][0]*p[4][2]*p[5][1] + 2*p[1][0]*p[4][0]*p[4][2]*p[5][1] + p[2][0]*p[4][0]*p[4][2]*p[5][1] - 3*p[3][0]*p[4][0]*p[4][2]*p[5][1] + p[1][2]*p[2][0]*p[5][0]*p[5][1] - p[1][0]*p[2][2]*p[5][0]*p[5][1] + 2*p[2][2]*p[3][0]*p[5][0]*p[5][1] - 2*p[2][0]*p[3][2]*p[5][0]*p[5][1] - p[1][2]*p[4][0]*p[5][0]*p[5][1] - 2*p[2][2]*p[4][0]*p[5][0]*p[5][1] + 3*p[3][2]*p[4][0]*p[5][0]*p[5][1] + p[1][0]*p[4][2]*p[5][0]*p[5][1] + 2*p[2][0]*p[4][2]*p[5][0]*p[5][1] - 3*p[3][0]*p[4][2]*p[5][0]*p[5][1] + p[0][2]*(pow(p[3][0], 2)*p[2][1] + pow(p[5][0], 2)*p[2][1] + p[2][0]*p[2][1]*p[3][0] - pow(p[2][0], 2)*p[3][1] + pow(p[4][0], 2)*p[3][1] - pow(p[5][0], 2)*p[3][1] - p[2][0]*p[3][0]*p[3][1] + 2*p[3][0]*p[3][1]*p[4][0] + p[1][1]*(3*pow(p[2][0], 2) - pow(p[3][0], 2) - pow(p[4][0], 2) - p[3][0]*p[4][0]) - 2*pow(p[3][0], 2)*p[4][1] - p[3][0]*p[4][0]*p[4][1] + pow(p[1][0], 2)*(-3*p[2][1] + p[3][1] + 2*p[4][1]) + p[1][0]*(-3*p[2][0]*p[2][1] + p[3][0]*p[3][1] + p[1][1]*(3*p[2][0] - p[3][0] - 2*p[4][0]) + p[3][1]*p[4][0] + p[4][0]*p[4][1]) + 2*p[2][0]*p[2][1]*p[5][0] + p[2][1]*p[3][0]*p[5][0] - p[2][0]*p[3][1]*p[5][0] - 2*p[3][0]*p[3][1]*p[5][0] - 2*pow(p[2][0], 2)*p[5][1] + 2*pow(p[3][0], 2)*p[5][1] - p[2][0]*p[5][0]*p[5][1] + p[3][0]*p[5][0]*p[5][1]) + 2*pow(p[2][0], 2)*p[0][1]*p[5][2] - 2*pow(p[3][0], 2)*p[0][1]*p[5][2] - 2*pow(p[2][0], 2)*p[1][1]*p[5][2] + 2*pow(p[4][0], 2)*p[1][1]*p[5][2] - p[1][0]*p[1][1]*p[2][0]*p[5][2] + pow(p[1][0], 2)*p[2][1]*p[5][2] - pow(p[3][0], 2)*p[2][1]*p[5][2] + pow(p[4][0], 2)*p[2][1]*p[5][2] + 2*p[1][0]*p[2][0]*p[2][1]*p[5][2] - p[2][0]*p[2][1]*p[3][0]*p[5][2] + pow(p[2][0], 2)*p[3][1]*p[5][2] - 3*pow(p[4][0], 2)*p[3][1]*p[5][2] + p[2][0]*p[3][0]*p[3][1]*p[5][2] + p[1][0]*p[1][1]*p[4][0]*p[5][2] + p[1][0]*p[2][1]*p[4][0]*p[5][2] + p[2][0]*p[2][1]*p[4][0]*p[5][2] - 3*p[3][0]*p[3][1]*p[4][0]*p[5][2] - pow(p[1][0], 2)*p[4][1]*p[5][2] - pow(p[2][0], 2)*p[4][1]*p[5][2] + 3*pow(p[3][0], 2)*p[4][1]*p[5][2] - p[1][0]*p[2][0]*p[4][1]*p[5][2] - 2*p[1][0]*p[4][0]*p[4][1]*p[5][2] - p[2][0]*p[4][0]*p[4][1]*p[5][2] + 3*p[3][0]*p[4][0]*p[4][1]*p[5][2] + p[0][1]*p[2][0]*p[5][0]*p[5][2] - p[1][1]*p[2][0]*p[5][0]*p[5][2] + p[1][0]*p[2][1]*p[5][0]*p[5][2] - p[0][1]*p[3][0]*p[5][0]*p[5][2] - 2*p[2][1]*p[3][0]*p[5][0]*p[5][2] + 2*p[2][0]*p[3][1]*p[5][0]*p[5][2] + p[1][1]*p[4][0]*p[5][0]*p[5][2] + 2*p[2][1]*p[4][0]*p[5][0]*p[5][2] - 3*p[3][1]*p[4][0]*p[5][0]*p[5][2] - p[1][0]*p[4][1]*p[5][0]*p[5][2] - 2*p[2][0]*p[4][1]*p[5][0]*p[5][2] + 3*p[3][0]*p[4][1]*p[5][0]*p[5][2] + pow(p[0][0], 2)*(2*p[2][2]*p[3][1] - 2*p[2][1]*p[3][2] + p[1][2]*(3*p[2][1] - 2*p[3][1] - p[4][1]) + p[3][2]*p[4][1] - p[3][1]*p[4][2] + p[1][1]*(-3*p[2][2] + 2*p[3][2] + p[4][2]) + p[2][2]*p[5][1] - p[3][2]*p[5][1] - p[2][1]*p[5][2] + p[3][1]*p[5][2]) + p[0][0]*(3*p[1][0]*p[1][2]*p[2][1] + 3*p[1][2]*p[2][0]*p[2][1] - 3*p[1][0]*p[1][1]*p[2][2] - 3*p[1][1]*p[2][0]*p[2][2] - p[1][0]*p[1][2]*p[3][1] + p[2][0]*p[2][2]*p[3][1] - p[1][2]*p[3][0]*p[3][1] + p[2][2]*p[3][0]*p[3][1] + p[1][0]*p[1][1]*p[3][2] - p[2][0]*p[2][1]*p[3][2] + p[1][1]*p[3][0]*p[3][2] - p[2][1]*p[3][0]*p[3][2] - 2*p[1][0]*p[1][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[1][0]*p[3][2]*p[4][1] + 2*p[3][0]*p[3][2]*p[4][1] - p[1][2]*p[4][0]*p[4][1] + p[3][2]*p[4][0]*p[4][1] + 2*p[1][0]*p[1][1]*p[4][2] + p[1][1]*p[3][0]*p[4][2] - p[1][0]*p[3][1]*p[4][2] - 2*p[3][0]*p[3][1]*p[4][2] + p[1][1]*p[4][0]*p[4][2] - p[3][1]*p[4][0]*p[4][2] + 2*p[2][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[2][0]*p[3][2]*p[5][1] - 2*p[3][0]*p[3][2]*p[5][1] + p[2][2]*p[5][0]*p[5][1] - p[3][2]*p[5][0]*p[5][1] + p[0][2]*(2*p[2][1]*p[3][0] - 2*p[2][0]*p[3][1] + p[1][1]*(3*p[2][0] - 2*p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-3*p[2][1] + 2*p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - 2*p[2][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[2][0]*p[3][1]*p[5][2] + 2*p[3][0]*p[3][1]*p[5][2] - p[2][1]*p[5][0]*p[5][2] + p[3][1]*p[5][0]*p[5][2] + p[0][1]*(-2*p[2][2]*p[3][0] + 2*p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + 2*p[3][0] + p[4][0]) + p[1][0]*(3*p[2][2] - 2*p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2])))/72;

result[2][0]=
(3*pow(p[2][1], 2)*p[0][0]*p[1][2] - pow(p[3][1], 2)*p[0][0]*p[1][2] - pow(p[4][1], 2)*p[0][0]*p[1][2] + pow(p[4][1], 2)*p[1][2]*p[2][0] + pow(p[5][1], 2)*p[1][2]*p[2][0] + 3*p[0][0]*p[1][1]*p[1][2]*p[2][1] - 3*pow(p[1][1], 2)*p[0][0]*p[2][2] + pow(p[3][1], 2)*p[0][0]*p[2][2] + pow(p[5][1], 2)*p[0][0]*p[2][2] - pow(p[4][1], 2)*p[1][0]*p[2][2] - pow(p[5][1], 2)*p[1][0]*p[2][2] - 3*p[0][0]*p[1][1]*p[2][1]*p[2][2] - 2*pow(p[4][1], 2)*p[1][2]*p[3][0] + 2*pow(p[5][1], 2)*p[2][2]*p[3][0] - p[0][0]*p[1][1]*p[1][2]*p[3][1] + p[0][0]*p[2][1]*p[2][2]*p[3][1] + pow(p[1][1], 2)*p[0][0]*p[3][2] - pow(p[2][1], 2)*p[0][0]*p[3][2] + pow(p[4][1], 2)*p[0][0]*p[3][2] - pow(p[5][1], 2)*p[0][0]*p[3][2] + 2*pow(p[4][1], 2)*p[1][0]*p[3][2] - 2*pow(p[5][1], 2)*p[2][0]*p[3][2] + p[0][0]*p[1][1]*p[3][1]*p[3][2] - p[0][0]*p[2][1]*p[3][1]*p[3][2] - pow(p[2][1], 2)*p[1][2]*p[4][0] + pow(p[3][1], 2)*p[1][2]*p[4][0] - pow(p[5][1], 2)*p[1][2]*p[4][0] - 2*p[1][1]*p[1][2]*p[2][1]*p[4][0] + 2*pow(p[1][1], 2)*p[2][2]*p[4][0] - 2*pow(p[5][1], 2)*p[2][2]*p[4][0] + p[1][1]*p[2][1]*p[2][2]*p[4][0] + p[1][1]*p[1][2]*p[3][1]*p[4][0] - pow(p[1][1], 2)*p[3][2]*p[4][0] + 3*pow(p[5][1], 2)*p[3][2]*p[4][0] - p[1][1]*p[3][1]*p[3][2]*p[4][0] - 2*p[0][0]*p[1][1]*p[1][2]*p[4][1] + 2*p[1][1]*p[1][2]*p[2][0]*p[4][1] + p[1][2]*p[2][0]*p[2][1]*p[4][1] - 2*p[1][0]*p[1][1]*p[2][2]*p[4][1] - p[1][0]*p[2][1]*p[2][2]*p[4][1] - p[1][1]*p[1][2]*p[3][0]*p[4][1] - p[0][0]*p[1][2]*p[3][1]*p[4][1] - p[1][2]*p[3][0]*p[3][1]*p[4][1] + p[0][0]*p[1][1]*p[3][2]*p[4][1] + p[1][0]*p[1][1]*p[3][2]*p[4][1] + 2*p[0][0]*p[3][1]*p[3][2]*p[4][1] + p[1][0]*p[3][1]*p[3][2]*p[4][1] - p[1][2]*p[2][1]*p[4][0]*p[4][1] + p[1][1]*p[2][2]*p[4][0]*p[4][1] + 2*p[1][2]*p[3][1]*p[4][0]*p[4][1] - 2*p[1][1]*p[3][2]*p[4][0]*p[4][1] + 2*pow(p[1][1], 2)*p[0][0]*p[4][2] - 2*pow(p[3][1], 2)*p[0][0]*p[4][2] + pow(p[2][1], 2)*p[1][0]*p[4][2] - pow(p[3][1], 2)*p[1][0]*p[4][2] + pow(p[5][1], 2)*p[1][0]*p[4][2] - 2*pow(p[1][1], 2)*p[2][0]*p[4][2] + 2*pow(p[5][1], 2)*p[2][0]*p[4][2] + 2*p[1][0]*p[1][1]*p[2][1]*p[4][2] - p[1][1]*p[2][0]*p[2][1]*p[4][2] + pow(p[1][1], 2)*p[3][0]*p[4][2] - 3*pow(p[5][1], 2)*p[3][0]*p[4][2] - p[1][0]*p[1][1]*p[3][1]*p[4][2] + p[1][1]*p[3][0]*p[3][1]*p[4][2] + p[0][0]*p[1][1]*p[4][1]*p[4][2] - p[1][1]*p[2][0]*p[4][1]*p[4][2] + p[1][0]*p[2][1]*p[4][1]*p[4][2] + 2*p[1][1]*p[3][0]*p[4][1]*p[4][2] - p[0][0]*p[3][1]*p[4][1]*p[4][2] - 2*p[1][0]*p[3][1]*p[4][1]*p[4][2] - 2*pow(p[2][1], 2)*p[1][2]*p[5][0] + 2*pow(p[4][1], 2)*p[1][2]*p[5][0] - p[1][1]*p[1][2]*p[2][1]*p[5][0] + pow(p[1][1], 2)*p[2][2]*p[5][0] - pow(p[3][1], 2)*p[2][2]*p[5][0] + pow(p[4][1], 2)*p[2][2]*p[5][0] + 2*p[1][1]*p[2][1]*p[2][2]*p[5][0] - p[2][1]*p[2][2]*p[3][1]*p[5][0] + pow(p[2][1], 2)*p[3][2]*p[5][0] - 3*pow(p[4][1], 2)*p[3][2]*p[5][0] + p[2][1]*p[3][1]*p[3][2]*p[5][0] + p[1][1]*p[1][2]*p[4][1]*p[5][0] + p[1][1]*p[2][2]*p[4][1]*p[5][0] + p[2][1]*p[2][2]*p[4][1]*p[5][0] - 3*p[3][1]*p[3][2]*p[4][1]*p[5][0] - pow(p[1][1], 2)*p[4][2]*p[5][0] - pow(p[2][1], 2)*p[4][2]*p[5][0] + 3*pow(p[3][1], 2)*p[4][2]*p[5][0] - p[1][1]*p[2][1]*p[4][2]*p[5][0] - 2*p[1][1]*p[4][1]*p[4][2]*p[5][0] - p[2][1]*p[4][1]*p[4][2]*p[5][0] + 3*p[3][1]*p[4][1]*p[4][2]*p[5][0] + p[1][1]*p[1][2]*p[2][0]*p[5][1] + 2*p[1][2]*p[2][0]*p[2][1]*p[5][1] - p[1][0]*p[1][1]*p[2][2]*p[5][1] + 2*p[0][0]*p[2][1]*p[2][2]*p[5][1] - 2*p[1][0]*p[2][1]*p[2][2]*p[5][1] + p[2][1]*p[2][2]*p[3][0]*p[5][1] + p[0][0]*p[2][2]*p[3][1]*p[5][1] + p[2][2]*p[3][0]*p[3][1]*p[5][1] - p[0][0]*p[2][1]*p[3][2]*p[5][1] - p[2][0]*p[2][1]*p[3][2]*p[5][1] - 2*p[0][0]*p[3][1]*p[3][2]*p[5][1] - p[2][0]*p[3][1]*p[3][2]*p[5][1] - p[1][1]*p[1][2]*p[4][0]*p[5][1] - p[1][2]*p[2][1]*p[4][0]*p[5][1] - p[2][1]*p[2][2]*p[4][0]*p[5][1] + 3*p[3][1]*p[3][2]*p[4][0]*p[5][1] + p[1][2]*p[2][0]*p[4][1]*p[5][1] - p[1][0]*p[2][2]*p[4][1]*p[5][1] - 2*p[1][2]*p[4][0]*p[4][1]*p[5][1] - p[2][2]*p[4][0]*p[4][1]*p[5][1] + 3*p[3][2]*p[4][0]*p[4][1]*p[5][1] + p[1][0]*p[1][1]*p[4][2]*p[5][1] + p[1][0]*p[2][1]*p[4][2]*p[5][1] + p[2][0]*p[2][1]*p[4][2]*p[5][1] - 3*p[3][0]*p[3][1]*p[4][2]*p[5][1] + 2*p[1][0]*p[4][1]*p[4][2]*p[5][1] + p[2][0]*p[4][1]*p[4][2]*p[5][1] - 3*p[3][0]*p[4][1]*p[4][2]*p[5][1] - p[1][2]*p[2][1]*p[5][0]*p[5][1] + p[1][1]*p[2][2]*p[5][0]*p[5][1] - 2*p[2][2]*p[3][1]*p[5][0]*p[5][1] + 2*p[2][1]*p[3][2]*p[5][0]*p[5][1] + p[1][2]*p[4][1]*p[5][0]*p[5][1] + 2*p[2][2]*p[4][1]*p[5][0]*p[5][1] - 3*p[3][2]*p[4][1]*p[5][0]*p[5][1] - p[1][1]*p[4][2]*p[5][0]*p[5][1] - 2*p[2][1]*p[4][2]*p[5][0]*p[5][1] + 3*p[3][1]*p[4][2]*p[5][0]*p[5][1] + p[0][2]*(-(pow(p[3][1], 2)*p[2][0]) - pow(p[5][1], 2)*p[2][0] + pow(p[2][1], 2)*p[3][0] - pow(p[4][1], 2)*p[3][0] + pow(p[5][1], 2)*p[3][0] - p[2][0]*p[2][1]*p[3][1] + p[2][1]*p[3][0]*p[3][1] + pow(p[1][1], 2)*(3*p[2][0] - p[3][0] - 2*p[4][0]) + 2*pow(p[3][1], 2)*p[4][0] - 2*p[3][0]*p[3][1]*p[4][1] + p[3][1]*p[4][0]*p[4][1] + p[1][0]*(-3*pow(p[2][1], 2) + pow(p[3][1], 2) + pow(p[4][1], 2) + p[3][1]*p[4][1]) - p[1][1]*(-3*p[2][0]*p[2][1] + p[3][0]*p[3][1] + p[1][0]*(3*p[2][1] - p[3][1] - 2*p[4][1]) + p[3][0]*p[4][1] + p[4][0]*p[4][1]) + 2*pow(p[2][1], 2)*p[5][0] - 2*pow(p[3][1], 2)*p[5][0] - 2*p[2][0]*p[2][1]*p[5][1] + p[2][1]*p[3][0]*p[5][1] - p[2][0]*p[3][1]*p[5][1] + 2*p[3][0]*p[3][1]*p[5][1] + p[2][1]*p[5][0]*p[5][1] - p[3][1]*p[5][0]*p[5][1]) - 2*pow(p[2][1], 2)*p[0][0]*p[5][2] + 2*pow(p[3][1], 2)*p[0][0]*p[5][2] + 2*pow(p[2][1], 2)*p[1][0]*p[5][2] - 2*pow(p[4][1], 2)*p[1][0]*p[5][2] - pow(p[1][1], 2)*p[2][0]*p[5][2] + pow(p[3][1], 2)*p[2][0]*p[5][2] - pow(p[4][1], 2)*p[2][0]*p[5][2] + p[1][0]*p[1][1]*p[2][1]*p[5][2] - 2*p[1][1]*p[2][0]*p[2][1]*p[5][2] - pow(p[2][1], 2)*p[3][0]*p[5][2] + 3*pow(p[4][1], 2)*p[3][0]*p[5][2] + p[2][0]*p[2][1]*p[3][1]*p[5][2] - p[2][1]*p[3][0]*p[3][1]*p[5][2] + pow(p[1][1], 2)*p[4][0]*p[5][2] + pow(p[2][1], 2)*p[4][0]*p[5][2] - 3*pow(p[3][1], 2)*p[4][0]*p[5][2] + p[1][1]*p[2][1]*p[4][0]*p[5][2] - p[1][0]*p[1][1]*p[4][1]*p[5][2] - p[1][1]*p[2][0]*p[4][1]*p[5][2] - p[2][0]*p[2][1]*p[4][1]*p[5][2] + 3*p[3][0]*p[3][1]*p[4][1]*p[5][2] + 2*p[1][1]*p[4][0]*p[4][1]*p[5][2] + p[2][1]*p[4][0]*p[4][1]*p[5][2] - 3*p[3][1]*p[4][0]*p[4][1]*p[5][2] - p[1][1]*p[2][0]*p[5][1]*p[5][2] - p[0][0]*p[2][1]*p[5][1]*p[5][2] + p[1][0]*p[2][1]*p[5][1]*p[5][2] - 2*p[2][1]*p[3][0]*p[5][1]*p[5][2] + p[0][0]*p[3][1]*p[5][1]*p[5][2] + 2*p[2][0]*p[3][1]*p[5][1]*p[5][2] + p[1][1]*p[4][0]*p[5][1]*p[5][2] + 2*p[2][1]*p[4][0]*p[5][1]*p[5][2] - 3*p[3][1]*p[4][0]*p[5][1]*p[5][2] - p[1][0]*p[4][1]*p[5][1]*p[5][2] - 2*p[2][0]*p[4][1]*p[5][1]*p[5][2] + 3*p[3][0]*p[4][1]*p[5][1]*p[5][2] + pow(p[0][1], 2)*(-2*p[2][2]*p[3][0] + 2*p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + 2*p[3][0] + p[4][0]) + p[1][0]*(3*p[2][2] - 2*p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2]) + p[0][1]*(3*p[0][0]*p[1][2]*p[2][1] - 3*p[1][2]*p[2][0]*p[2][1] + 3*p[1][0]*p[2][1]*p[2][2] - p[2][1]*p[2][2]*p[3][0] - 2*p[0][0]*p[1][2]*p[3][1] + 2*p[0][0]*p[2][2]*p[3][1] + p[1][2]*p[3][0]*p[3][1] - p[2][2]*p[3][0]*p[3][1] - 2*p[0][0]*p[2][1]*p[3][2] + p[2][0]*p[2][1]*p[3][2] - p[1][0]*p[3][1]*p[3][2] + p[2][0]*p[3][1]*p[3][2] + p[1][2]*p[3][1]*p[4][0] - 2*p[3][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][2]*p[4][0]*p[4][1] - p[3][2]*p[4][0]*p[4][1] - p[0][0]*p[3][1]*p[4][2] - p[1][0]*p[3][1]*p[4][2] + 2*p[3][0]*p[3][1]*p[4][2] - p[1][0]*p[4][1]*p[4][2] + p[3][0]*p[4][1]*p[4][2] + p[1][1]*(3*p[1][0]*p[2][2] - p[1][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + p[3][0] + 2*p[4][0]) - 2*p[1][0]*p[4][2] + p[3][0]*p[4][2] + p[0][0]*(-3*p[2][2] + 2*p[3][2] + p[4][2])) - 2*p[2][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] + 2*p[3][1]*p[3][2]*p[5][0] + p[0][0]*p[2][2]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][2]*p[5][0]*p[5][1] + p[3][2]*p[5][0]*p[5][1] + p[0][2]*(2*p[2][1]*p[3][0] - 2*p[2][0]*p[3][1] + p[1][1]*(3*p[2][0] - 2*p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-3*p[2][1] + 2*p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[0][0]*p[2][1]*p[5][2] + 2*p[2][0]*p[2][1]*p[5][2] - p[2][1]*p[3][0]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][0]*p[3][1]*p[5][2] - 2*p[3][0]*p[3][1]*p[5][2] + p[2][0]*p[5][1]*p[5][2] - p[3][0]*p[5][1]*p[5][2]))/72;

result[3][0]=
(-3*pow(p[2][2], 2)*p[0][0]*p[1][1] + pow(p[3][2], 2)*p[0][0]*p[1][1] + pow(p[4][2], 2)*p[0][0]*p[1][1] - pow(p[4][2], 2)*p[1][1]*p[2][0] - pow(p[5][2], 2)*p[1][1]*p[2][0] + 3*pow(p[1][2], 2)*p[0][0]*p[2][1] - pow(p[3][2], 2)*p[0][0]*p[2][1] - pow(p[5][2], 2)*p[0][0]*p[2][1] + pow(p[4][2], 2)*p[1][0]*p[2][1] + pow(p[5][2], 2)*p[1][0]*p[2][1] - 3*p[0][0]*p[1][1]*p[1][2]*p[2][2] + 3*p[0][0]*p[1][2]*p[2][1]*p[2][2] + 2*pow(p[4][2], 2)*p[1][1]*p[3][0] - 2*pow(p[5][2], 2)*p[2][1]*p[3][0] - pow(p[1][2], 2)*p[0][0]*p[3][1] + pow(p[2][2], 2)*p[0][0]*p[3][1] - pow(p[4][2], 2)*p[0][0]*p[3][1] + pow(p[5][2], 2)*p[0][0]*p[3][1] - 2*pow(p[4][2], 2)*p[1][0]*p[3][1] + 2*pow(p[5][2], 2)*p[2][0]*p[3][1] + p[0][0]*p[1][1]*p[1][2]*p[3][2] - p[0][0]*p[2][1]*p[2][2]*p[3][2] - p[0][0]*p[1][2]*p[3][1]*p[3][2] + p[0][0]*p[2][2]*p[3][1]*p[3][2] + pow(p[2][2], 2)*p[1][1]*p[4][0] - pow(p[3][2], 2)*p[1][1]*p[4][0] + pow(p[5][2], 2)*p[1][1]*p[4][0] - 2*pow(p[1][2], 2)*p[2][1]*p[4][0] + 2*pow(p[5][2], 2)*p[2][1]*p[4][0] + 2*p[1][1]*p[1][2]*p[2][2]*p[4][0] - p[1][2]*p[2][1]*p[2][2]*p[4][0] + pow(p[1][2], 2)*p[3][1]*p[4][0] - 3*pow(p[5][2], 2)*p[3][1]*p[4][0] - p[1][1]*p[1][2]*p[3][2]*p[4][0] + p[1][2]*p[3][1]*p[3][2]*p[4][0] - 2*pow(p[1][2], 2)*p[0][0]*p[4][1] + 2*pow(p[3][2], 2)*p[0][0]*p[4][1] - pow(p[2][2], 2)*p[1][0]*p[4][1] + pow(p[3][2], 2)*p[1][0]*p[4][1] - pow(p[5][2], 2)*p[1][0]*p[4][1] + 2*pow(p[1][2], 2)*p[2][0]*p[4][1] - 2*pow(p[5][2], 2)*p[2][0]*p[4][1] - 2*p[1][0]*p[1][2]*p[2][2]*p[4][1] + p[1][2]*p[2][0]*p[2][2]*p[4][1] - pow(p[1][2], 2)*p[3][0]*p[4][1] + 3*pow(p[5][2], 2)*p[3][0]*p[4][1] + p[1][0]*p[1][2]*p[3][2]*p[4][1] - p[1][2]*p[3][0]*p[3][2]*p[4][1] + 2*p[0][0]*p[1][1]*p[1][2]*p[4][2] - 2*p[1][1]*p[1][2]*p[2][0]*p[4][2] + 2*p[1][0]*p[1][2]*p[2][1]*p[4][2] - p[1][1]*p[2][0]*p[2][2]*p[4][2] + p[1][0]*p[2][1]*p[2][2]*p[4][2] + p[1][1]*p[1][2]*p[3][0]*p[4][2] - p[0][0]*p[1][2]*p[3][1]*p[4][2] - p[1][0]*p[1][2]*p[3][1]*p[4][2] + p[0][0]*p[1][1]*p[3][2]*p[4][2] + p[1][1]*p[3][0]*p[3][2]*p[4][2] - 2*p[0][0]*p[3][1]*p[3][2]*p[4][2] - p[1][0]*p[3][1]*p[3][2]*p[4][2] - p[1][2]*p[2][1]*p[4][0]*p[4][2] + p[1][1]*p[2][2]*p[4][0]*p[4][2] + 2*p[1][2]*p[3][1]*p[4][0]*p[4][2] - 2*p[1][1]*p[3][2]*p[4][0]*p[4][2] - p[0][0]*p[1][2]*p[4][1]*p[4][2] + p[1][2]*p[2][0]*p[4][1]*p[4][2] - p[1][0]*p[2][2]*p[4][1]*p[4][2] - 2*p[1][2]*p[3][0]*p[4][1]*p[4][2] + p[0][0]*p[3][2]*p[4][1]*p[4][2] + 2*p[1][0]*p[3][2]*p[4][1]*p[4][2] + 2*pow(p[2][2], 2)*p[1][1]*p[5][0] - 2*pow(p[4][2], 2)*p[1][1]*p[5][0] - pow(p[1][2], 2)*p[2][1]*p[5][0] + pow(p[3][2], 2)*p[2][1]*p[5][0] - pow(p[4][2], 2)*p[2][1]*p[5][0] + p[1][1]*p[1][2]*p[2][2]*p[5][0] - 2*p[1][2]*p[2][1]*p[2][2]*p[5][0] - pow(p[2][2], 2)*p[3][1]*p[5][0] + 3*pow(p[4][2], 2)*p[3][1]*p[5][0] + p[2][1]*p[2][2]*p[3][2]*p[5][0] - p[2][2]*p[3][1]*p[3][2]*p[5][0] + pow(p[1][2], 2)*p[4][1]*p[5][0] + pow(p[2][2], 2)*p[4][1]*p[5][0] - 3*pow(p[3][2], 2)*p[4][1]*p[5][0] + p[1][2]*p[2][2]*p[4][1]*p[5][0] - p[1][1]*p[1][2]*p[4][2]*p[5][0] - p[1][2]*p[2][1]*p[4][2]*p[5][0] - p[2][1]*p[2][2]*p[4][2]*p[5][0] + 3*p[3][1]*p[3][2]*p[4][2]*p[5][0] + 2*p[1][2]*p[4][1]*p[4][2]*p[5][0] + p[2][2]*p[4][1]*p[4][2]*p[5][0] - 3*p[3][2]*p[4][1]*p[4][2]*p[5][0] + 2*pow(p[2][2], 2)*p[0][0]*p[5][1] - 2*pow(p[3][2], 2)*p[0][0]*p[5][1] - 2*pow(p[2][2], 2)*p[1][0]*p[5][1] + 2*pow(p[4][2], 2)*p[1][0]*p[5][1] + pow(p[1][2], 2)*p[2][0]*p[5][1] - pow(p[3][2], 2)*p[2][0]*p[5][1] + pow(p[4][2], 2)*p[2][0]*p[5][1] - p[1][0]*p[1][2]*p[2][2]*p[5][1] + 2*p[1][2]*p[2][0]*p[2][2]*p[5][1] + pow(p[2][2], 2)*p[3][0]*p[5][1] - 3*pow(p[4][2], 2)*p[3][0]*p[5][1] - p[2][0]*p[2][2]*p[3][2]*p[5][1] + p[2][2]*p[3][0]*p[3][2]*p[5][1] - pow(p[1][2], 2)*p[4][0]*p[5][1] - pow(p[2][2], 2)*p[4][0]*p[5][1] + 3*pow(p[3][2], 2)*p[4][0]*p[5][1] - p[1][2]*p[2][2]*p[4][0]*p[5][1] + p[1][0]*p[1][2]*p[4][2]*p[5][1] + p[1][2]*p[2][0]*p[4][2]*p[5][1] + p[2][0]*p[2][2]*p[4][2]*p[5][1] - 3*p[3][0]*p[3][2]*p[4][2]*p[5][1] - 2*p[1][2]*p[4][0]*p[4][2]*p[5][1] - p[2][2]*p[4][0]*p[4][2]*p[5][1] + 3*p[3][2]*p[4][0]*p[4][2]*p[5][1] + pow(p[0][2], 2)*(2*p[2][1]*p[3][0] - 2*p[2][0]*p[3][1] + p[1][1]*(3*p[2][0] - 2*p[3][0] - p[4][0]) + p[3][1]*p[4][0] - p[3][0]*p[4][1] + p[1][0]*(-3*p[2][1] + 2*p[3][1] + p[4][1]) + p[2][1]*p[5][0] - p[3][1]*p[5][0] - p[2][0]*p[5][1] + p[3][0]*p[5][1]) - p[1][1]*p[1][2]*p[2][0]*p[5][2] + p[1][0]*p[1][2]*p[2][1]*p[5][2] - 2*p[1][1]*p[2][0]*p[2][2]*p[5][2] - 2*p[0][0]*p[2][1]*p[2][2]*p[5][2] + 2*p[1][0]*p[2][1]*p[2][2]*p[5][2] - p[2][1]*p[2][2]*p[3][0]*p[5][2] + p[0][0]*p[2][2]*p[3][1]*p[5][2] + p[2][0]*p[2][2]*p[3][1]*p[5][2] - p[0][0]*p[2][1]*p[3][2]*p[5][2] - p[2][1]*p[3][0]*p[3][2]*p[5][2] + 2*p[0][0]*p[3][1]*p[3][2]*p[5][2] + p[2][0]*p[3][1]*p[3][2]*p[5][2] + p[1][1]*p[1][2]*p[4][0]*p[5][2] + p[1][1]*p[2][2]*p[4][0]*p[5][2] + p[2][1]*p[2][2]*p[4][0]*p[5][2] - 3*p[3][1]*p[3][2]*p[4][0]*p[5][2] - p[1][0]*p[1][2]*p[4][1]*p[5][2] - p[1][0]*p[2][2]*p[4][1]*p[5][2] - p[2][0]*p[2][2]*p[4][1]*p[5][2] + 3*p[3][0]*p[3][2]*p[4][1]*p[5][2] - p[1][1]*p[2][0]*p[4][2]*p[5][2] + p[1][0]*p[2][1]*p[4][2]*p[5][2] + 2*p[1][1]*p[4][0]*p[4][2]*p[5][2] + p[2][1]*p[4][0]*p[4][2]*p[5][2] - 3*p[3][1]*p[4][0]*p[4][2]*p[5][2] - 2*p[1][0]*p[4][1]*p[4][2]*p[5][2] - p[2][0]*p[4][1]*p[4][2]*p[5][2] + 3*p[3][0]*p[4][1]*p[4][2]*p[5][2] - p[1][2]*p[2][1]*p[5][0]*p[5][2] + p[1][1]*p[2][2]*p[5][0]*p[5][2] - 2*p[2][2]*p[3][1]*p[5][0]*p[5][2] + 2*p[2][1]*p[3][2]*p[5][0]*p[5][2] + p[1][2]*p[4][1]*p[5][0]*p[5][2] + 2*p[2][2]*p[4][1]*p[5][0]*p[5][2] - 3*p[3][2]*p[4][1]*p[5][0]*p[5][2] - p[1][1]*p[4][2]*p[5][0]*p[5][2] - 2*p[2][1]*p[4][2]*p[5][0]*p[5][2] + 3*p[3][1]*p[4][2]*p[5][0]*p[5][2] + p[1][2]*p[2][0]*p[5][1]*p[5][2] + p[0][0]*p[2][2]*p[5][1]*p[5][2] - p[1][0]*p[2][2]*p[5][1]*p[5][2] + 2*p[2][2]*p[3][0]*p[5][1]*p[5][2] - p[0][0]*p[3][2]*p[5][1]*p[5][2] - 2*p[2][0]*p[3][2]*p[5][1]*p[5][2] - p[1][2]*p[4][0]*p[5][1]*p[5][2] - 2*p[2][2]*p[4][0]*p[5][1]*p[5][2] + 3*p[3][2]*p[4][0]*p[5][1]*p[5][2] + p[1][0]*p[4][2]*p[5][1]*p[5][2] + 2*p[2][0]*p[4][2]*p[5][1]*p[5][2] - 3*p[3][0]*p[4][2]*p[5][1]*p[5][2] + p[0][1]*(pow(p[3][2], 2)*p[2][0] + pow(p[5][2], 2)*p[2][0] - pow(p[2][2], 2)*p[3][0] + pow(p[4][2], 2)*p[3][0] - pow(p[5][2], 2)*p[3][0] + p[2][0]*p[2][2]*p[3][2] - p[2][2]*p[3][0]*p[3][2] - 2*pow(p[3][2], 2)*p[4][0] + pow(p[1][2], 2)*(-3*p[2][0] + p[3][0] + 2*p[4][0]) + 2*p[3][0]*p[3][2]*p[4][2] - p[3][2]*p[4][0]*p[4][2] + p[1][0]*(3*pow(p[2][2], 2) - pow(p[3][2], 2) - pow(p[4][2], 2) - p[3][2]*p[4][2]) + p[1][2]*(-3*p[2][0]*p[2][2] + p[3][0]*p[3][2] + p[1][0]*(3*p[2][2] - p[3][2] - 2*p[4][2]) + p[3][0]*p[4][2] + p[4][0]*p[4][2]) - 2*pow(p[2][2], 2)*p[5][0] + 2*pow(p[3][2], 2)*p[5][0] + 2*p[2][0]*p[2][2]*p[5][2] - p[2][2]*p[3][0]*p[5][2] + p[2][0]*p[3][2]*p[5][2] - 2*p[3][0]*p[3][2]*p[5][2] - p[2][2]*p[5][0]*p[5][2] + p[3][2]*p[5][0]*p[5][2]) + p[0][2]*(3*p[0][0]*p[1][2]*p[2][1] - 3*p[1][0]*p[1][2]*p[2][1] - 3*p[1][0]*p[2][1]*p[2][2] + p[2][1]*p[2][2]*p[3][0] - 2*p[0][0]*p[1][2]*p[3][1] + p[1][0]*p[1][2]*p[3][1] + 2*p[0][0]*p[2][2]*p[3][1] - p[2][0]*p[2][2]*p[3][1] - 2*p[0][0]*p[2][1]*p[3][2] + p[2][1]*p[3][0]*p[3][2] + p[1][0]*p[3][1]*p[3][2] - p[2][0]*p[3][1]*p[3][2] + p[1][2]*p[3][1]*p[4][0] + 2*p[3][1]*p[3][2]*p[4][0] - p[0][0]*p[1][2]*p[4][1] + 2*p[1][0]*p[1][2]*p[4][1] - p[1][2]*p[3][0]*p[4][1] + p[0][0]*p[3][2]*p[4][1] + p[1][0]*p[3][2]*p[4][1] - 2*p[3][0]*p[3][2]*p[4][1] - p[0][0]*p[3][1]*p[4][2] + p[3][1]*p[4][0]*p[4][2] + p[1][0]*p[4][1]*p[4][2] - p[3][0]*p[4][1]*p[4][2] + p[1][1]*(3*p[2][0]*p[2][2] - p[3][0]*p[3][2] + p[1][2]*(3*p[2][0] - p[3][0] - 2*p[4][0]) - p[3][2]*p[4][0] - p[4][0]*p[4][2] + p[0][0]*(-3*p[2][2] + 2*p[3][2] + p[4][2])) + 2*p[2][1]*p[2][2]*p[5][0] - p[2][2]*p[3][1]*p[5][0] + p[2][1]*p[3][2]*p[5][0] - 2*p[3][1]*p[3][2]*p[5][0] + p[0][0]*p[2][2]*p[5][1] - 2*p[2][0]*p[2][2]*p[5][1] + p[2][2]*p[3][0]*p[5][1] - p[0][0]*p[3][2]*p[5][1] - p[2][0]*p[3][2]*p[5][1] + 2*p[3][0]*p[3][2]*p[5][1] - p[0][0]*p[2][1]*p[5][2] + p[0][0]*p[3][1]*p[5][2] + p[2][1]*p[5][0]*p[5][2] - p[3][1]*p[5][0]*p[5][2] - p[2][0]*p[5][1]*p[5][2] + p[3][0]*p[5][1]*p[5][2] + p[0][1]*(-2*p[2][2]*p[3][0] + 2*p[2][0]*p[3][2] - p[3][2]*p[4][0] + p[1][2]*(-3*p[2][0] + 2*p[3][0] + p[4][0]) + p[1][0]*(3*p[2][2] - 2*p[3][2] - p[4][2]) + p[3][0]*p[4][2] - p[2][2]*p[5][0] + p[3][2]*p[5][0] + p[2][0]*p[5][2] - p[3][0]*p[5][2])))/72;

            return result;
          }
        };

      // 3D Tetrahedron

      template<>
        class Operator<VolumeIntegral, Basis::Unitary, Tetrahedron>:
          public OperatorBase<VolumeIntegral, Basis::Unitary, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Tetrahedron);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Unitary, dim>::function_size;
          static constexpr size_t point_size=4;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

result[0][0]=
(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2]))/6;

            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Linear, Tetrahedron>:
          public OperatorBase<VolumeIntegral, Basis::Linear, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Tetrahedron);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Linear, dim>::function_size;
          static constexpr size_t point_size=4;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

result[0][0]=
(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2]))/6;

result[1][0]=
((p[0][0] + p[1][0] + p[2][0] + p[3][0])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;

result[2][0]=
((p[0][1] + p[1][1] + p[2][1] + p[3][1])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;

result[3][0]=
((p[0][2] + p[1][2] + p[2][2] + p[3][2])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;


            return result;
          }
        };

      template<>
        class Operator<VolumeIntegral, Basis::Quadratic, Tetrahedron>:
          public OperatorBase<VolumeIntegral, Basis::Quadratic, Interval>
        {
        public:
          static constexpr size_t dim = dimension(Tetrahedron);
          static constexpr size_t operator_size=1;
          static constexpr size_t basis_size=
            Basis::Traits<Basis::Quadratic, dim>::function_size;
          static constexpr size_t point_size=4;

          using result_t = array<array<double, operator_size>, basis_size>;
          using points_t = array<Wonton::Point<dim>, point_size>;

          static result_t apply(const points_t p) {
            result_t result;

result[0][0]=
(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2]))/6;

result[1][0]=
((p[0][0] + p[1][0] + p[2][0] + p[3][0])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;

result[2][0]=
((p[0][1] + p[1][1] + p[2][1] + p[3][1])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;

result[3][0]=
((p[0][2] + p[1][2] + p[2][2] + p[3][2])*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/24;

result[4][0]=
((pow(p[0][0], 2) + pow(p[1][0], 2) + pow(p[2][0], 2) + pow(p[3][0], 2) + p[2][0]*p[3][0] + p[1][0]*(p[2][0] + p[3][0]) + p[0][0]*(p[1][0] + p[2][0] + p[3][0]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;

result[5][0]=
-((2*p[1][0]*p[1][1] + p[1][1]*p[2][0] + p[1][0]*p[2][1] + 2*p[2][0]*p[2][1] + p[1][1]*p[3][0] + p[2][1]*p[3][0] + p[0][1]*(p[1][0] + p[2][0] + p[3][0]) + p[1][0]*p[3][1] + p[2][0]*p[3][1] + 2*p[3][0]*p[3][1] + p[0][0]*(2*p[0][1] + p[1][1] + p[2][1] + p[3][1]))*(-(p[0][0]*p[1][2]*p[2][1]) + p[0][0]*p[1][1]*p[2][2] + p[1][2]*p[2][1]*p[3][0] - p[1][1]*p[2][2]*p[3][0] + p[0][0]*p[1][2]*p[3][1] - p[1][2]*p[2][0]*p[3][1] - p[0][0]*p[2][2]*p[3][1] + p[1][0]*p[2][2]*p[3][1] + p[0][2]*(-(p[2][1]*p[3][0]) + p[1][1]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][1] - p[3][1]) + p[2][0]*p[3][1]) - p[0][0]*p[1][1]*p[3][2] + p[1][1]*p[2][0]*p[3][2] + p[0][0]*p[2][1]*p[3][2] - p[1][0]*p[2][1]*p[3][2] + p[0][1]*(p[1][2]*(p[2][0] - p[3][0]) + p[2][2]*p[3][0] - p[2][0]*p[3][2] + p[1][0]*(-p[2][2] + p[3][2]))))/120;

result[6][0]=
((2*p[1][0]*p[1][2] + p[1][2]*p[2][0] + p[1][0]*p[2][2] + 2*p[2][0]*p[2][2] + p[1][2]*p[3][0] + p[2][2]*p[3][0] + p[0][2]*(p[1][0] + p[2][0] + p[3][0]) + p[1][0]*p[3][2] + p[2][0]*p[3][2] + 2*p[3][0]*p[3][2] + p[0][0]*(2*p[0][2] + p[1][2] + p[2][2] + p[3][2]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;

result[7][0]=
((pow(p[0][1], 2) + pow(p[1][1], 2) + pow(p[2][1], 2) + pow(p[3][1], 2) + p[2][1]*p[3][1] + p[1][1]*(p[2][1] + p[3][1]) + p[0][1]*(p[1][1] + p[2][1] + p[3][1]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;

result[8][0]=
((2*p[1][1]*p[1][2] + p[1][2]*p[2][1] + p[1][1]*p[2][2] + 2*p[2][1]*p[2][2] + p[1][2]*p[3][1] + p[2][2]*p[3][1] + p[0][2]*(p[1][1] + p[2][1] + p[3][1]) + p[1][1]*p[3][2] + p[2][1]*p[3][2] + 2*p[3][1]*p[3][2] + p[0][1]*(2*p[0][2] + p[1][2] + p[2][2] + p[3][2]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;

result[9][0]=
((pow(p[0][2], 2) + pow(p[1][2], 2) + pow(p[2][2], 2) + pow(p[3][2], 2) + p[2][2]*p[3][2] + p[1][2]*(p[2][2] + p[3][2]) + p[0][2]*(p[1][2] + p[2][2] + p[3][2]))*(p[0][0]*p[1][2]*p[2][1] - p[0][0]*p[1][1]*p[2][2] - p[1][2]*p[2][1]*p[3][0] + p[1][1]*p[2][2]*p[3][0] - p[0][0]*p[1][2]*p[3][1] + p[1][2]*p[2][0]*p[3][1] + p[0][0]*p[2][2]*p[3][1] - p[1][0]*p[2][2]*p[3][1] + p[0][2]*(p[1][1]*(p[2][0] - p[3][0]) + p[2][1]*p[3][0] - p[2][0]*p[3][1] + p[1][0]*(-p[2][1] + p[3][1])) + p[0][0]*p[1][1]*p[3][2] - p[1][1]*p[2][0]*p[3][2] - p[0][0]*p[2][1]*p[3][2] + p[1][0]*p[2][1]*p[3][2] + p[0][1]*(-(p[2][2]*p[3][0]) + p[1][2]*(-p[2][0] + p[3][0]) + p[1][0]*(p[2][2] - p[3][2]) + p[2][0]*p[3][2])))/120;

            return result;
          }
        };

      ////////////////////////////////////////////////////////////////////////////////
      // Dynamic size information
      ////////////////////////////////////////////////////////////////////////////////
         inline
          array<size_t, 3> size_info(Type otype, Basis::Type btype, Domain domain) {
            array<size_t, 3> r;
            if (otype == VolumeIntegral) {
              switch(btype) {
              case Basis::Unitary:
                switch (domain) {
                case Interval:
                  {using OP=Operator<VolumeIntegral,Basis::Unitary,Interval>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Quadrilateral:
                  {using OP=Operator<VolumeIntegral,Basis::Unitary,Quadrilateral>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Triangle:
                  {using OP=Operator<VolumeIntegral,Basis::Unitary,Triangle>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Hexahedron:
                  {using OP=Operator<VolumeIntegral,Basis::Unitary,Hexahedron>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Wedge:
                  {using OP=Operator<VolumeIntegral,Basis::Unitary,Wedge>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Tetrahedron:
                  {using OP=Operator<VolumeIntegral,Basis::Unitary,Tetrahedron>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                }
                break;
              case Basis::Linear:
                switch (domain) {
                case Interval:
                  {using OP=Operator<VolumeIntegral,Basis::Linear,Interval>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Quadrilateral:
                  {using OP=Operator<VolumeIntegral,Basis::Linear,Quadrilateral>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Triangle:
                  {using OP=Operator<VolumeIntegral,Basis::Linear,Triangle>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Hexahedron:
                  {using OP=Operator<VolumeIntegral,Basis::Linear,Hexahedron>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Wedge:
                  {using OP=Operator<VolumeIntegral,Basis::Linear,Wedge>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Tetrahedron:
                  {using OP=Operator<VolumeIntegral,Basis::Linear,Tetrahedron>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                }
                break;
              case Basis::Quadratic:
                switch (domain) {
                case Interval:
                  {using OP=Operator<VolumeIntegral,Basis::Quadratic,Interval>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Quadrilateral:
                  {using OP=Operator<VolumeIntegral,Basis::Quadratic,Quadrilateral>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Triangle:
                  {using OP=Operator<VolumeIntegral,Basis::Quadratic,Triangle>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                case Tetrahedron:
                  {using OP=Operator<VolumeIntegral,Basis::Quadratic,Tetrahedron>;
                    r[0]=OP::operator_size; r[1]=OP::basis_size; r[2]=OP::point_size; break;}
                }
                break;
              }
            }
            return r;
          }

      ////////////////////////////////////////////////////////////////////////////////
      // Helper template functions
      ////////////////////////////////////////////////////////////////////////////////

      template<class OP>
        inline void copy_points(const vector<Wonton::Point<OP::dim>> &points,
                                typename OP::points_t &apts)
        {
          for (int i=0; i<OP::point_size; i++) {
            apts[i] = points[i];
          }
        }

      template<class OP>
        inline Wonton::Point<OP::dim> centroid(const typename OP::points_t &apts)
        {
          Wonton::Point<OP::dim> cent;
          for (int j=0; j<OP::dim; j++) cent[j]=0.;
          for (int i=0; i<OP::point_size; i++)
            for (int j=0; j<OP::dim; j++)
              cent[j] += apts[i][j];
          for (int j=0; j<OP::dim; j++)
            cent[j] /= OP::point_size;
          return cent;
        }

      template<class OP>
        inline void shift_points(const Wonton::Point<OP::dim> c, typename OP::points_t &apts)
        {
          for (int i=0; i<OP::point_size; i++)
            for (int j=0; j<OP::dim; j++)
              apts[i][j] = apts[i][j] - c[j];
        }

      template<class OP>
        inline void resize_result(vector<vector<double>> &result)
        {
          result.resize(OP::basis_size);
          for (auto iter=result.begin(); iter!=result.end(); iter++) {
            iter->resize(OP::operator_size);
          }
        }

      template<class OP>
        inline void copy_result(const typename OP::result_t &ares,
                                vector<vector<double>> &result)
        {
          for (int i=0; i<OP::basis_size; i++) {
            for (int j=0; j<OP::operator_size; j++) {
              result[i][j] = ares[i][j];
            }
          }
        }

      template<class OP>
        inline void get_result(const vector<Wonton::Point<OP::dim>> &points,
                               vector<vector<double>> &result, const bool center=true)
        {
          resize_result<OP>(result);
          typename OP::points_t apts; copy_points<OP>(points, apts);
          if (center) {
            Wonton::Point<OP::dim> c = centroid<OP>(apts);
            shift_points<OP>(c, apts);
            auto tf = Basis::transfactor<OP::dim>(OP::basis, c);
            typename OP::result_t ares = OP::apply(apts);
            for (int i=0; i<OP::basis_size; i++)
              for (int j=0; j<OP::operator_size; j++) {
                result[i][j]=0.;
              }
            for (int j=0; j<OP::operator_size; j++)
              for (int i=0; i<OP::basis_size; i++)
                for (int k=0; k<OP::basis_size; k++) {
                  result[i][j] += tf[i][k]*ares[k][j];
                }
          } else {
            typename OP::result_t ares = OP::apply(apts);
            copy_result<OP>(ares, result);
          }
        }

      ////////////////////////////////////////////////////////////////////////////////
      // Dynamic Accessor
      ////////////////////////////////////////////////////////////////////////////////

      template<size_t dim>
      void apply(const Type type, const Basis::Type basis_type,
                 const Domain domain_type, const vector<Wonton::Point<dim>> &points,
                 vector<vector<double>> &result);

      template<>
      inline
      void apply<1>(const Type type, const Basis::Type basis_type,
                    const Domain domain_type, const vector<Wonton::Point<1>> &points,
                    vector<vector<double>> &result)
      {
        bool center = true;
        switch(type) {
        case VolumeIntegral:
          switch(domain_type) {
          case Interval:
            switch(basis_type) {
            case Basis::Unitary:
              {using OP=Operator<VolumeIntegral, Basis::Unitary, Interval>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Linear:
              {using OP=Operator<VolumeIntegral, Basis::Linear, Interval>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Quadratic:
              {using OP=Operator<VolumeIntegral, Basis::Quadratic, Interval>;
                get_result<OP>(points, result, center);}
              break;
            default:
              throw(std::runtime_error("invalid basis"));
            }
            break;
          default:
            throw(std::runtime_error("invalid domain"));
          }
          break;
        default:
          throw(std::runtime_error("invalid operator"));
        }
      }

      template<>
      inline
      void apply<2>(const Type type, const Basis::Type basis_type,
                    const Domain domain_type, const vector<Wonton::Point<2>> &points,
                    vector<vector<double>> &result)
      {
        bool center = true;
        switch(type) {
        case VolumeIntegral:
          switch(domain_type) {
          case Quadrilateral:
            switch(basis_type) {
            case Basis::Unitary:
              {using OP=Operator<VolumeIntegral, Basis::Unitary, Quadrilateral>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Linear:
              {using OP=Operator<VolumeIntegral, Basis::Linear, Quadrilateral>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Quadratic:
              {using OP=Operator<VolumeIntegral, Basis::Quadratic, Quadrilateral>;
                get_result<OP>(points, result, center);}
              break;
            default:
              throw(std::runtime_error("invalid basis"));
            }
            break;
          case Triangle:
            switch(basis_type) {
            case Basis::Unitary:
              {using OP=Operator<VolumeIntegral, Basis::Unitary, Triangle>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Linear:
              {using OP=Operator<VolumeIntegral, Basis::Linear, Triangle>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Quadratic:
              {using OP=Operator<VolumeIntegral, Basis::Quadratic, Triangle>;
                get_result<OP>(points, result, center);}
              break;
            default:
              throw(std::runtime_error("invalid basis"));
            }
            break;
          default:
            throw(std::runtime_error("invalid domain"));
          }
          break;
        default:
          throw(std::runtime_error("invalid operator"));
        }
      }

      template<>
      inline
      void apply<3>(const Type type, const Basis::Type basis_type,
                    const Domain domain_type, const vector<Wonton::Point<3>> &points,
                    vector<vector<double>> &result)
      {
        bool center = true;
        switch(type) {
        case VolumeIntegral:
          switch(domain_type) {
          case Hexahedron:
            switch (basis_type) {
            case Basis::Unitary:
              {using OP=Operator<VolumeIntegral, Basis::Unitary, Hexahedron>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Linear:
              {using OP=Operator<VolumeIntegral, Basis::Linear, Hexahedron>;
                get_result<OP>(points, result, center);}
              break;
            default:
              throw(std::runtime_error("invalid basis"));
            }
            break;
          case Wedge:
            switch (basis_type) {
            case Basis::Unitary:
              {using OP=Operator<VolumeIntegral, Basis::Unitary, Wedge>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Linear:
              {using OP=Operator<VolumeIntegral, Basis::Linear, Wedge>;
                get_result<OP>(points, result, center);}
              break;
            default:
              throw(std::runtime_error("invalid basis"));
            }
            break;
          case Tetrahedron:
            switch (basis_type) {
            case Basis::Unitary:
              {using OP=Operator<VolumeIntegral, Basis::Unitary, Tetrahedron>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Linear:
              {using OP=Operator<VolumeIntegral, Basis::Linear, Tetrahedron>;
                get_result<OP>(points, result, center);}
              break;
            case Basis::Quadratic:
              {using OP=Operator<VolumeIntegral, Basis::Quadratic, Tetrahedron>;
                get_result<OP>(points, result, center);}
              break;
            default:
              throw(std::runtime_error("invalid basis"));
            }
            break;
          default:
            throw(std::runtime_error("invalid domain"));
          }
          break;
        default:
          throw(std::runtime_error("invalid operator"));
        }
      }
    }  // namespace Operator
  }  // namespace Meshfree
}  // namespace Portage

#endif // OPERATOR_H_INC_
