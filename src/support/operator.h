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

#include "portage/support/Point.h"
#include "portage/support/basis.h"

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

      template<Domain vol> class DomainTraits   {public: static const size_t dimension=0;};
      template<> class DomainTraits<Interval>     {public: static const size_t dimension=1;};
      template<> class DomainTraits<Quadrilateral>{public: static const size_t dimension=2;};
      template<> class DomainTraits<Triangle>     {public: static const size_t dimension=2;};
      template<> class DomainTraits<Circle>       {public: static const size_t dimension=2;};
      template<> class DomainTraits<Hexahedron>   {public: static const size_t dimension=3;};
      template<> class DomainTraits<Tetrahedron>  {public: static const size_t dimension=3;};
      template<> class DomainTraits<Wedge>        {public: static const size_t dimension=3;};
      template<> class DomainTraits<Sphere>       {public: static const size_t dimension=3;};

      constexpr size_t dimension(Domain vol) {
        switch (vol) {
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

      ////////////////////////////////////////////////////////////////////////////////
      // Template for integral operator base class
      ////////////////////////////////////////////////////////////////////////////////

      template<Type type, Basis::Type basis_type, Domain vol_type>
	class OperatorBase
      {
      public:
	Type operator_type=type;
	static constexpr Basis::Type basis=basis_type;
	static constexpr Domain vol=vol_type;
      };

      ////////////////////////////////////////////////////////////////////////////////
      // Template for integral operators
      ////////////////////////////////////////////////////////////////////////////////

      template<Type type, Basis::Type basis_type, Domain vol_type>
	class Operator: public OperatorBase<type,basis_type,vol_type>
      {
      public:
	static constexpr size_t dim = dimension(vol_type);
	static constexpr size_t operator_size=0;
	static constexpr size_t basis_size=0;
	static constexpr size_t point_size=0;

	using result_t = array<array<double, operator_size>, basis_size>;
	using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
	  using points_t = array<Point<dim>, point_size>;

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
      // Helper template functions
      ////////////////////////////////////////////////////////////////////////////////

      template<class OP>
	inline void copy_points(const vector<Point<OP::dim>> &points, 
				typename OP::points_t &apts) 
	{
	  for (int i=0; i<OP::point_size; i++) {
	    apts[i] = points[i];
	  }
	}

      template<class OP>
	inline Point<OP::dim> centroid(const typename OP::points_t &apts)
	{
	  Point<OP::dim> cent;
	  for (int j=0; j<OP::dim; j++) cent[j]=0.;
	  for (int i=0; i<OP::point_size; i++)
	    for (int j=0; j<OP::dim; j++)
	      cent[j] += apts[i][j];
	  for (int j=0; j<OP::dim; j++)
	    cent[j] /= OP::point_size;
	  return cent;
	}

      template<class OP>
	inline void shift_points(const Point<OP::dim> c, typename OP::points_t &apts)
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
	inline void get_result(const vector<Point<OP::dim>> &points, 
			       vector<vector<double>> &result, bool center=true) 
	{
	  resize_result<OP>(result);
	  typename OP::points_t apts; copy_points<OP>(points, apts);
	  if (center) {
	    Point<OP::dim> c = centroid<OP>(apts);
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
      void apply(Type type, Basis::Type basis_type, Domain vol_type, 
		 vector<Point<dim>> &points, vector<vector<double>> &result) 
      {
	switch(type) {
	case VolumeIntegral:
	  switch(vol_type) {
	  case Interval:
	    switch(basis_type) {
	    case Basis::Unitary:
	      {using OP=Operator<VolumeIntegral, Basis::Unitary, Interval>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Linear:
	      {using OP=Operator<VolumeIntegral, Basis::Linear, Interval>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Quadratic:
	      {using OP=Operator<VolumeIntegral, Basis::Quadratic, Interval>;
		get_result<OP>(points, result);}
	      break;
	    default:
	      throw(std::runtime_error("invalid basis"));
	    }
	    break;
	  case Quadrilateral:
	    switch(basis_type) {
	    case Basis::Unitary:
	      {using OP=Operator<VolumeIntegral, Basis::Unitary, Quadrilateral>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Linear:
	      {using OP=Operator<VolumeIntegral, Basis::Linear, Quadrilateral>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Quadratic:
	      {using OP=Operator<VolumeIntegral, Basis::Quadratic, Quadrilateral>;
		get_result<OP>(points, result);}
	      break;
	    default:
	      throw(std::runtime_error("invalid basis"));
	    }
	    break;
	  case Triangle:
	    switch(basis_type) {
	    case Basis::Unitary:
	      {using OP=Operator<VolumeIntegral, Basis::Unitary, Triangle>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Linear:
	      {using OP=Operator<VolumeIntegral, Basis::Linear, Triangle>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Quadratic:
	      {using OP=Operator<VolumeIntegral, Basis::Quadratic, Triangle>;
		get_result<OP>(points, result);}
	      break;
	    default:
	      throw(std::runtime_error("invalid basis"));
	    }
	    break;
	  case Hexahedron:
	    switch (basis_type) {
	    case Basis::Unitary:
	      {using OP=Operator<VolumeIntegral, Basis::Unitary, Hexahedron>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Linear:
	      {using OP=Operator<VolumeIntegral, Basis::Linear, Hexahedron>;
		get_result<OP>(points, result);}
	      break;
	    default:
	      throw(std::runtime_error("invalid basis"));
	    }
	    break;
	  case Wedge:
	    switch (basis_type) {
	    case Basis::Unitary:
	      {using OP=Operator<VolumeIntegral, Basis::Unitary, Wedge>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Linear:
	      {using OP=Operator<VolumeIntegral, Basis::Linear, Wedge>;
		get_result<OP>(points, result);}
	      break;
	    default:
	      throw(std::runtime_error("invalid basis"));
	    }
	    break;
	  case Tetrahedron:
	    switch (basis_type) {
	    case Basis::Unitary:
	      {using OP=Operator<VolumeIntegral, Basis::Unitary, Tetrahedron>;
		get_result<OP>(points, result);}
	      break;
	    case Basis::Linear:
	      {using OP=Operator<VolumeIntegral, Basis::Linear, Tetrahedron>;
		get_result<OP>(points, result);}
	      break;
	    default:
	      throw(std::runtime_error("invalid basis"));
	    }
	    break;
	  default:
	    throw(std::runtime_error("invalid volume"));
	  }
	  break;
	default:
	  throw(std::runtime_error("invalid operator"));
	}
      }
    }
  }
}

#endif // OPERATOR_H_INC_
