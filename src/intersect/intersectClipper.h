/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef INTERSECT_CLIPPER_H
#define INTERSECT_CLIPPER_H

#include "search.h"
#include "clipper.hpp"
#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>
/**\brief Return the area and moment of the polygon.
 * poly [input] 
 * \returns std::vector<double>--area, mx, my
 **/
std::vector<double> areaAndMomentPolygon(const std::vector<std::pair<double,double> > poly){
    double area = 0;
    double cx = 0;
    double cy = 0;
    std::vector <double> ret;
    for (int i=0;i<poly.size()-1;i++){
        double a = (poly[i].first*poly[i+1].second-poly[i+1].first*poly[i].second);  
        area+=a;
        cx+= (poly[i].first+poly[i+1].first)*a;
        cy+= (poly[i].second+poly[i+1].second)*a;
    }
    int lastIndex = poly.size()-1;
    //close the polygon
    double a = poly[lastIndex].first*poly[0].second - poly[0].first*poly[lastIndex].second;
    area+= a;
    cx+= (poly[lastIndex].first+poly[0].first)*a;
    cy+= (poly[lastIndex].second+poly[0].second)*a;
    ret.emplace_back(.5*area);
    ret.emplace_back(cx/6.);
    ret.emplace_back(cy/6.);
    return ret;
}

/*!
 * \class IntersectClipper <typename C> 2-D intersection algorithm for arbitrary convex and non-convex polyhedra
 * \brief The intersect class is templated on MeshWrapper type.  You must provide a method to convert the template cells to an IntersectClipper::Poly.  
 */

template <typename SourceMeshType, typename TargetMeshType=SourceMeshType> class IntersectClipper
{

public:
    typedef std::pair<double, double> Point; 
    typedef std::vector<Point> Poly; 
    //Provide volume and centroid
    typedef std::pair<double, Point> Moment;

    IntersectClipper(const SourceMeshType &s, const TargetMeshType &t):sourceMeshWrapper(s), targetMeshWrapper(t){}

    /*! \brief Intersect two cells and return the first two moments.
     * \param[in] cellA first cell index to intersect
     * \param[in] cellB second cell index to intersect
     * \return list of moments; ret[0] == 0th moment; ret[1] == first moment
     */

    std::vector<std::vector<double> > operator() (const int cellA, const int cellB) const {      
        Poly polyA = sourceMeshWrapper.cellToXY(cellA);
        Poly polyB = targetMeshWrapper.cellToXY(cellB);
        double max_size_poly = 0;
        max_size_poly = IntersectClipper::updateMaxSize(polyA, max_size_poly);
        max_size_poly = IntersectClipper::updateMaxSize(polyB, max_size_poly);

        ClipperLib::Path intPolyA = IntersectClipper::convertPoly2int(polyA, max_size_poly);
        ClipperLib::Path intPolyB = IntersectClipper::convertPoly2int(polyB, max_size_poly);

        //Make clipper aware of polyA and polyB
        ClipperLib::Clipper clpr;
        clpr.AddPath(intPolyA, ClipperLib::ptSubject, true);
        clpr.AddPath(intPolyB, ClipperLib::ptClip, true);

        //Compute the intersection of all paths which have been added to clipper (in this case
        //polyA and polyB; use the Even-Odd winding rule in case of self-intersecting polygons
        //For more information on the winding rules see:
        //http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Types/PolyFillType.htm
        ClipperLib::Paths solution;
        clpr.Execute(ClipperLib::ctIntersection, solution, ClipperLib::pftEvenOdd, ClipperLib::pftEvenOdd);

        std::vector< std::vector<std::pair<double, double>> > intersectionList;
        for(int i=0;i<solution.size();i++){
            std::vector<std::pair<double, double> > poly;
            for(int j=0;j<solution[i].size();j++){
                //Build list of intersections and convert back to doubles
                poly.emplace_back(std::make_pair(integer2real(solution[i][j].X, max_size_poly), integer2real(solution[i][j].Y, max_size_poly)));
            }
            intersectionList.emplace_back(poly);
        }
        std::vector<std::vector<double>> moments;
        for (int i=0;i<intersectionList.size();i++){
            moments.emplace_back(areaAndMomentPolygon(intersectionList[i]));
        }    
        return moments;
    }

    IntersectClipper() {}

    //! Copy constructor (disabled)
    IntersectClipper(const IntersectClipper &) = delete;

    //! Assignment operator (disabled)
    IntersectClipper & operator = (const IntersectClipper &) = delete;


private:

    //We must use the same max size for all the polygons, so the number we are looking for is the maximum value in the set--all the X and Y values will be converted using this value
    static double updateMaxSize( const std::vector<std::pair<double, double> > poly, double max_size_poly){
        for(auto i = poly.begin(), e=poly.end();i!=e;++i){
            double m = std::max(std::abs(i->first), std::abs(i->second));
            if (m > max_size_poly) max_size_poly = m;
        }
        return max_size_poly;
    }

    /**
       \poly Convert a double to an integer 
       \note Could do this by multiplying by a power of 10 large enough to preserve the 
       number of digits in the original double; however it is more precise to 
       multiply by a power of 2 
       \param a [in] - number to be converted
       \param max_size [in] 
    **/
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

    /** 
        \brief Convert integer back to double given the max_size of the numbers within the problem
        \param a [in] - number to convert
        \param max_size -
    **/
    static double integer2real(long a, const double max_size){
        int exp;
        frexp(max_size, &exp);  
        return ldexp(a, exp-DBL_MANT_DIG);
    }

    //Convert an entire polygon (specifiied as a std::vector<std::pair<double, double>) to a std::vector<IntPoint>
    static std::vector<ClipperLib::IntPoint> convertPoly2int(  std::vector<std::pair<double, double> > poly, double max_size_poly){
        std::vector<ClipperLib::IntPoint> intpoly(poly.size());
        std::transform(poly.begin(), poly.end(), intpoly.begin(), [max_size_poly](std::pair<double, double> point){
                return ClipperLib::IntPoint(real2integer(point.first, max_size_poly), real2integer(point.second, max_size_poly));}); 
        return intpoly;
    }
  
private:
    const SourceMeshType &sourceMeshWrapper;
    const TargetMeshType &targetMeshWrapper;

}; // class IntersectClipper

#endif // INTERSECT_CLIPPER_H

