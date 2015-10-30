/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SEARCH_SIMPLE_H
#define SEARCH_SIMPLE_H

/*!
    \class SearchSimple search_simple.h
    \brief SearchSimple provides...
 */

#include <vector>



namespace { // unnamed

void getBoundingBox(
    const std::vector<std::pair<double,double>> &cell_coord,
    double* xlow, double* xhigh,
    double* ylow, double* yhigh)
{
    const double big = 1.e99;
    double xl = big;  double xh = -big;
    double yl = big;  double yh = -big;

    for (int n = 0; n < cell_coord.size(); ++n) {
        std::pair<double,double> p;
        p = cell_coord[n];
        double x = p.first;
        double y = p.second;
        xl = std::min(xl, x);  xh = std::max(xh, x);
        yl = std::min(yl, y);  yh = std::max(yh, y);
    }

    *xlow = xl;  *xhigh = xh;
    *ylow = yl;  *yhigh = yh;

} // getBoundingBox

} // namespace


namespace Portage {

template <typename SourceMeshType, typename TargetMeshType>
class SearchSimple {
  public:

    //! Default constructor (disabled)
    SearchSimple() = delete;
    
    //! Constructor with Meshes
    /*!
      \brief Builds the search structure for finding intersection
      
      \param source_mesh_wrapper   pointer to wrapper for getting the source mesh info
      \param target_mesh_wrapper   pointer to wrapper for getting the target mesh info 
      
      Constructor for search structure for finding cells from a source
      mesh that overlap the target mesh

    */

    SearchSimple(const SourceMeshType & source_mesh, 
                 const TargetMeshType & target_mesh)
            : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {

        int numCells = sourceMesh_.num_owned_cells() + 
                sourceMesh_.num_ghost_cells();
        xlow_.reserve(numCells);
        xhigh_.reserve(numCells);
        ylow_.reserve(numCells);
        yhigh_.reserve(numCells);
        
        // find bounding boxes for all cells
        for (int c = 0; c < numCells; ++c) {
            std::vector<std::pair<double,double>> cell_coord;
            sourceMesh_.cell_get_coordinates(c, &cell_coord);
            getBoundingBox(cell_coord,
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
      \brief returns source mesh cells potentially overlapping given target cell

      \param cellId the index of the target cell that I pass in...
      \param candidates pointer to vector of potential candidate cells in sourceMesh
    */
    void search(const int cellId, std::vector<int> *candidates) const;

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
void SearchSimple<SourceMeshType,TargetMeshType>::
search(const int cellId, std::vector<int> *candidates)
const {
    // find bounding box for target cell
    std::vector<std::pair<double,double>> cell_coord;
    targetMesh_.cell_get_coordinates(cellId, &cell_coord);
    double txlow, txhigh, tylow, tyhigh;
    getBoundingBox(cell_coord,
                                   &txlow, &txhigh, &tylow, &tyhigh);
    
    // now see which sourceMesh cells have bounding boxes overlapping
    // with target cell
    // do a naive linear search
    int numCells = sourceMesh_.num_owned_cells() + sourceMesh_.num_ghost_cells();
    for (int c = 0; c < numCells; ++c) {
        if (std::max(txlow, xlow_[c]) < std::min(txhigh, xhigh_[c]) &&
                std::max(tylow, ylow_[c]) < std::min(tyhigh, yhigh_[c])) {
            candidates->push_back(c);
        }
    }

} // SearchSimple::search




} // namespace Portage

#endif // SEARCH_SIMPLE_H

/*--------------------------------------------------------------------------~-*
 * Formatting options for Emacs and vim.
 *
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:nil
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *--------------------------------------------------------------------------~-*/
