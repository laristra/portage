/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SEARCH_SIMPLE_H
#define SEARCH_SIMPLE_H

#include <vector>
#include <algorithm>

#include "portage/support/Point.h"

/// file-local namespace
namespace search_simple {

  /*!
    @brief Given a list of 2d coordinates, find the coordinates of the
    bounding box.
    @param[in] cell_coord List of coordinates as (x,y) pairs of points
    to bound.
    @param[in,out] xlow Minimum @a x of bounding box.
    @param[in,out] xhigh Maximum @a x of bounding box.
    @param[in,out] ylow Minimum @a y of bounding box.
    @param[in,out] yhigh Maximum @a y of bounding box.
   */
void getBoundingBox(
    const std::vector<Portage::Point<2>> &cell_coord,
    double* xlow, double* xhigh,
    double* ylow, double* yhigh)
{
    const double big = 1.e99;
    double xl = big;  double xh = -big;
    double yl = big;  double yh = -big;

    for (const auto& p : cell_coord) {
        const double x = p[0];
        const double y = p[1];
        xl = std::min(xl, x);  xh = std::max(xh, x);
        yl = std::min(yl, y);  yh = std::max(yh, y);
    }

    *xlow = xl;  *xhigh = xh;
    *ylow = yl;  *yhigh = yh;

} // getBoundingBox

} // search_simple


namespace Portage {

  /*!
    @class SearchSimple "search_simple.h"
    @brief A simple, crude search algorithm that utilizes bounding boxes
    in 2d.
    @tparam SourceMeshType The mesh type of the input mesh.
    @tparam TargetMeshType The mesh type of the output mesh.

    This search is only valid for 2d meshes.
   */
template <typename SourceMeshType, typename TargetMeshType>
class SearchSimple {
  public:

    //! Default constructor (disabled)
    SearchSimple() = delete;
    
    // Constructor with Meshes
    /*!
      @brief Builds the search structure for finding intersection.
      @param[in] source_mesh_wrapper Pointer to a mesh wrapper for
      getting the source mesh info.
      @param[in] target_mesh_wrapper Pointer to a mesh wrapper for
      getting the target mesh info.
      
      Constructor for search structure for finding cells from a source
      mesh that overlap the target mesh.
    */
    SearchSimple(const SourceMeshType & source_mesh, 
                 const TargetMeshType & target_mesh)
            : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {

        const int numCells = sourceMesh_.num_owned_cells() + 
                sourceMesh_.num_ghost_cells();
        xlow_.reserve(numCells);
        xhigh_.reserve(numCells);
        ylow_.reserve(numCells);
        yhigh_.reserve(numCells);
        
        // find bounding boxes for all cells
        for (int c = 0; c < numCells; ++c) {
            std::vector<Point<2>> cell_coord;
            sourceMesh_.cell_get_coordinates(c, &cell_coord);
            search_simple::getBoundingBox(cell_coord,
                           &xlow_[c], &xhigh_[c],
                           &ylow_[c], &yhigh_[c]);
        }
    } // SearchSimple::SearchSimple

    //! Copy constructor (disabled)
    SearchSimple(const SearchSimple &) = delete;

    //! Assignment operator (disabled)
    SearchSimple & operator = (const SearchSimple &) = delete;

    //! Destructor
    ~SearchSimple() = default;

    /*!
      @brief Find the source mesh cells potentially overlapping a given
      target cell.
      @param[in] cellId The index of the cell in the target mesh for
      which we wish to find the candidate overlapping cells in the
      source mesh.
      @param[in,out] candidates Pointer to a vector of potential candidate
      cells in the source mesh.
    */
    void operator() (const int cellId, std::vector<int> *candidates) const;

  private:

    // Aggregate data members
    const SourceMeshType & sourceMesh_;
    const TargetMeshType & targetMesh_;
    std::vector<double> xlow_;
    std::vector<double> xhigh_;
    std::vector<double> ylow_;
    std::vector<double> yhigh_;

}; // class SearchSimple



template<typename SourceMeshType, typename TargetMeshType>
void SearchSimple<SourceMeshType, TargetMeshType>::
operator() (const int cellId, std::vector<int> *candidates)
const {
    // find bounding box for target cell
    std::vector<Point<2>> cell_coord;
    targetMesh_.cell_get_coordinates(cellId, &cell_coord);
    double txlow, txhigh, tylow, tyhigh;
    search_simple::getBoundingBox(cell_coord, &txlow, &txhigh, &tylow, &tyhigh);
    
    // now see which sourceMesh cells have bounding boxes overlapping
    // with target cell
    // do a naive linear search
    const int numCells = sourceMesh_.num_owned_cells() +
        sourceMesh_.num_ghost_cells();
    for (int c = 0; c < numCells; ++c) {
        if (std::max(txlow, xlow_[c]) < std::min(txhigh, xhigh_[c]) &&
                std::max(tylow, ylow_[c]) < std::min(tyhigh, yhigh_[c])) {
            candidates->push_back(c);
        }
    }

} // SearchSimple::operator()




} // namespace Portage

#endif // SEARCH_SIMPLE_H

