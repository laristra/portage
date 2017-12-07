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

      enum Geometry{
	Interval,         // 1D
	Quadrilateral,    // 2D
	Triangle,     
	Circle,       
	Hexahedron,       // 3D    
	Tetrahedron,  
	TrianglePrism, 
	Sphere, 
	LastGeometry
      };

      enum Type {VolumeIntegral, SurfaceIntegral, LastOperator};

      template<Geometry geo> class GeometryTraits   {public: static const size_t dimension=0;};
      template<> class GeometryTraits<Interval>     {public: static const size_t dimension=1;};
      template<> class GeometryTraits<Quadrilateral>{public: static const size_t dimension=2;};
      template<> class GeometryTraits<Triangle>     {public: static const size_t dimension=2;};
      template<> class GeometryTraits<Circle>       {public: static const size_t dimension=2;};
      template<> class GeometryTraits<Hexahedron>   {public: static const size_t dimension=3;};
      template<> class GeometryTraits<Tetrahedron>  {public: static const size_t dimension=3;};
      template<> class GeometryTraits<TrianglePrism>{public: static const size_t dimension=3;};
      template<> class GeometryTraits<Sphere>       {public: static const size_t dimension=3;};

      constexpr size_t dimension(Geometry geo) {
        switch (geo) {
          case Interval:      return GeometryTraits<Interval>::dimension;      break;
          case Quadrilateral: return GeometryTraits<Quadrilateral>::dimension; break;
          case Triangle:      return GeometryTraits<Triangle>::dimension;      break;
          case Circle:        return GeometryTraits<Circle>::dimension;        break;
          case Hexahedron:    return GeometryTraits<Hexahedron>::dimension;    break;
          case Tetrahedron:   return GeometryTraits<Tetrahedron>::dimension;   break;
          case TrianglePrism: return GeometryTraits<TrianglePrism>::dimension; break;
          case Sphere:        return GeometryTraits<Sphere>::dimension;        break;
          default: return 0; 
        }
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Template for integral operators
      ////////////////////////////////////////////////////////////////////////////////

      template<Type type, Basis::Type basis_type, Geometry geo_type>
          class Operator
      {
     public:
	static constexpr size_t dim = dimension(geo_type);
	static constexpr size_t operator_size=0;
	static constexpr size_t basis_size=0;
	static constexpr size_t point_size=0;

	using result_t = array<array<double, operator_size>, basis_size>;
	using points_t = array<Point<dim>, point_size>;

	static result_t apply(const points_t points) {
	  result_t result;
	  throw(std::runtime_error("invalid operator"));
	  return result;
	}
      };

      ////////////////////////////////////////////////////////////////////////////////
      // Specializations 
      ////////////////////////////////////////////////////////////////////////////////

      template<>
          class Operator<VolumeIntegral, Basis::Unitary, Interval>
      {
     public:
	static constexpr size_t dim = dimension(Interval);
	static constexpr size_t operator_size=1;
	static constexpr size_t basis_size=
            Basis::Traits<Basis::Unitary, dim>::function_size;
	static constexpr size_t point_size=2;

	using result_t = array<array<double, operator_size>, basis_size>;
	using points_t = array<Point<dim>, point_size>;

	static result_t apply(const points_t points) {
	  result_t result;
	  result[0][0] = points[1][0] - points[0][0];
	  return result;
	}
      };

      ////////////////////////////////////////////////////////////////////////////////
      // Helper template functions
      ////////////////////////////////////////////////////////////////////////////////

      template<class OP>
          inline void resize_result(vector<vector<double>> &result) 
      {
	result.resize(OP::operator_size); 
	for (auto iter=result.begin(); iter!=result.end(); iter++) {
	  iter->resize(OP::basis_size);
	}
      }

      template<class OP>
          inline void copy_points(const vector<Point<OP::dim>> &points, 
                                  typename OP::points_t &apts) 
      {
	for (int i=0; i<OP::point_size; i++) {
	  apts[i] = points[i];
	}
      }

      template<class OP>
          inline void copy_result(const typename OP::result_t &ares, 
                                  vector<vector<double>> &result) 
      {
	resize_result<OP>(result);
	for (int i=0; i<OP::basis_size; i++) {
	  for (int j=0; j<OP::operator_size; j++) {
	    result[i][j] = ares[i][j];
	  }
	}
      }

      template<class OP>
          inline void get_result(const vector<Point<OP::dim>> &points, 
                                 vector<vector<double>> &result) 
      {
	typename OP::points_t apts; copy_points<OP>(points, apts);
	typename OP::result_t ares = OP::apply(apts);
	copy_result<OP>(ares, result);
      }

      ////////////////////////////////////////////////////////////////////////////////
      // Dynamic Accessor
      ////////////////////////////////////////////////////////////////////////////////

      template<size_t dim>
          void apply(Type type, Basis::Type basis_type, Geometry geo_type, 
                     vector<Point<dim>> &points, vector<vector<double>> &result) 
      {
        assert(dim == dimension(geo_type));
	switch(dim) {
          case 1:
            switch(geo_type) {
              case Interval:
                switch(basis_type) {
                  case Basis::Unitary:
                    switch(type) {
                      case VolumeIntegral:
                        using OP=Operator<VolumeIntegral, Basis::Unitary, Interval>;
                        get_result<OP>(points, result);
                        break;
                      default:
			throw(std::runtime_error("invalid operator"));
                    }
		    break;
                  default:
		    throw(std::runtime_error("invalid basis"));
                }
		break;
              default:
		throw(std::runtime_error("invalid geometry"));
            }
	    break;
          default:
	    throw(std::runtime_error("invalid dimension"));
	}
      }
    }
  }
}

#endif // OPERATOR_H_INC_
