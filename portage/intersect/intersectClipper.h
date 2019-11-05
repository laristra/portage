/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_INTERSECT_INTERSECTCLIPPER_H_
#define PORTAGE_INTERSECT_INTERSECTCLIPPER_H_

#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>

// portage includes
#include "portage/intersect/clipper.hpp"
#include "portage/support/portage.h"
#include "wonton/support/Point.h"

namespace Portage {

/*!
  @brief Return the area and moment of the polygon.
  @param[in] poly A vector of a pair of (x,y) coordinates of the nodes making up the polygon.
  @returns std::vector<double>--area, mx, my
*/
inline
std::vector<double> areaAndMomentPolygon(const std::vector<Wonton::Point<2>> poly){
  double area = 0;
  double cx = 0;
  double cy = 0;
  std::vector <double> ret;
  for(int i=0;i<poly.size()-1;i++){
    double a = (poly[i][0]*poly[i+1][1]-poly[i+1][0]*poly[i][1]);
    area+=a;
    cx+= (poly[i][0]+poly[i+1][0])*a;
    cy+= (poly[i][1]+poly[i+1][1])*a;
  }
  int lastIndex = poly.size()-1;
  //close the polygon
  double a = poly[lastIndex][0]*poly[0][1] - poly[0][0]*poly[lastIndex][1];
  area+= a;
  cx+= (poly[lastIndex][0]+poly[0][0])*a;
  cy+= (poly[lastIndex][1]+poly[0][1])*a;
  ret.emplace_back(.5*area);
  ret.emplace_back(cx/6.);
  ret.emplace_back(cy/6.);
  return ret;
}

/*!
  @class IntersectClipper "intersectClipper.h"
  @brief 2-D intersection algorithm for arbitrary convex and non-convex polyhedra
  @tparam SourceMeshType The mesh type of the input mesh.
  @tparam TargetMeshType The mesh type of the output mesh.

  The intersect class is templated on MeshWrapper type.  You must provide a method to convert
  the template cells to an IntersectClipper::Poly.
*/

template <typename SourceMeshType, typename TargetMeshType=SourceMeshType> class IntersectClipper
{

public:

/// Alias for a collection of Points.
typedef std::vector<Wonton::Point<2>> Poly;
/// Alias to provide volume and centroid
typedef std::pair<double, Wonton::Point<2>> Moment;

/// Constructor taking a source mesh @c s and a target mesh @c t.
IntersectClipper(const SourceMeshType &s, const TargetMeshType &t):sourceMeshWrapper(s), targetMeshWrapper(t){}

/*!
  @brief Intersect two cells and return the first two moments.
  @param[in] cellA first cell index to intersect
  @param[in] cellB second cell index to intersect
  @return list of moments; ret[0] == 0th moment; ret[1] == first moment
*/
std::vector<std::vector<double> > operator() (const int cellA, const int cellB) const {
  Poly polyA, polyB;
  sourceMeshWrapper.cell_get_coordinates(cellA, &polyA);
  targetMeshWrapper.cell_get_coordinates(cellB, &polyB);
  double max_size_poly = 0;
  max_size_poly = IntersectClipper::updateMaxSize(polyA, max_size_poly);
  max_size_poly = IntersectClipper::updateMaxSize(polyB, max_size_poly);

  const ClipperLib::Path intPolyA = IntersectClipper::convertPoly2int(polyA, max_size_poly);
  const ClipperLib::Path intPolyB = IntersectClipper::convertPoly2int(polyB, max_size_poly);

  //Make clipper aware of polyA and polyB
  ClipperLib::Clipper clpr;
  clpr.AddPath(intPolyA, ClipperLib::ptSubject, true);
  clpr.AddPath(intPolyB, ClipperLib::ptClip, true);

  //Compute the intersection of all paths which have been added to clipper (in this case
  //polyA and polyB; use the Even-Odd winding rule in case of self-intersecting polygons
  //For more information on the winding rules see:
  //http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Types/PolyFillType.htm
  ClipperLib::Paths solution;
  clpr.Execute(ClipperLib::ctIntersection, solution,
               ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd);

  std::vector<std::vector<Wonton::Point<2>>> intersectionList;
  for(auto const &i: solution){
    std::vector<Wonton::Point<2>> poly;
    for(auto const &j: i){
      //Build list of intersections and convert back to doubles
      poly.emplace_back(integer2real(j.X, max_size_poly), integer2real(j.Y, max_size_poly));
    }
    intersectionList.emplace_back(poly);
  }
  std::vector<std::vector<double>> moments;
  for(auto const &i: intersectionList){
    moments.emplace_back(areaAndMomentPolygon(i));
  }
  return moments;
}

/// Default constructor.
IntersectClipper() {}

/// Copy constructor (disabled)
IntersectClipper(const IntersectClipper &) = delete;

/// Assignment operator (disabled)
IntersectClipper & operator = (const IntersectClipper &) = delete;


private:

//We must use the same max size for all the polygons, so the number we are looking for is the maximum value in the set--all the X and Y values will be converted using this value
static double updateMaxSize( const std::vector<Wonton::Point<2>> poly, double max_size_poly){
  for(auto const &i: poly){
    double m = std::max(std::abs(i[0]), std::abs(i[1]));
    if (m > max_size_poly) max_size_poly = m;
  }
  return max_size_poly;
}

/*!
  @brief Convert a double to an integer.

  @note Could do this by multiplying by a power of 10 large enough to preserve the
  number of digits in the original double; however it is more precise to
  multiply by a power of 2.

  @param[in] a number to be converted
  @param[in] max_size
*/
static long real2integer(double a, const double max_size){
  int exp;

  // We want to move the decimal point of the floating point number so that the last digit given is before the decimal point.
  // The alogrithm is to find the exponent (in the floating point representation) of the largest number in the polygons (either X or Y coordinate).
  //Multiply by the power of two which is number of digits in mantissa (in this case, 53) - exponent
  //lrint, then rounds to an integer

  //Compute the exponent
  frexp(max_size, &exp);
  //Size of long must be >= size double
  return lrint(ldexp(a, DBL_MANT_DIG-exp));
}

/*!
  @brief Convert integer back to double given the max_size of the numbers within the problem
  @param[in] a number to convert
  @param[in] max_size
*/
static double integer2real(long a, const double max_size){
  int exp;
  frexp(max_size, &exp);
  return ldexp(a, exp-DBL_MANT_DIG);
}

//Convert an entire polygon (specifiied as a std::vector<Wonton::Point>) to a std::vector<IntPoint>
static std::vector<ClipperLib::IntPoint> convertPoly2int(std::vector<Wonton::Point<2>> poly, double max_size_poly){
  std::vector<ClipperLib::IntPoint> intpoly(poly.size());
  std::transform(poly.begin(), poly.end(), intpoly.begin(), [max_size_poly](Wonton::Point<2> point){
      return ClipperLib::IntPoint(real2integer(point[0], max_size_poly), real2integer(point[1], max_size_poly));});
  return intpoly;
}

private:
const SourceMeshType &sourceMeshWrapper;
const TargetMeshType &targetMeshWrapper;

};  // class IntersectClipper

}  // namespace Portage

#endif  // PORTAGE_INTERSECT_INTERSECTCLIPPER_H_
