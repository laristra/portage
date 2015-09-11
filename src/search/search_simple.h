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

template<typename MeshWrapper>
void getBoundingBox(
    const MeshWrapper & mesh,
    const std::vector<int>& nodes,
    double* xlow, double* xhigh,
    double* ylow, double* yhigh)
{
    const double big = 1.e99;
    double xl = big;  double xh = -big;
    double yl = big;  double yh = -big;

    for (int n = 0; n < nodes.size(); ++n) {
        std::pair<double,double> p;
        mesh.node_get_coordinates(nodes[n],&p);
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

template <typename SourceMeshWrapper, typename TargetMeshWrapper>
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

    SearchSimple(const SourceMeshWrapper & source_mesh, 
                 const TargetMeshWrapper & target_mesh);

    //! Copy constructor (disabled)
    SearchSimple(const SearchSimple &) = delete;

    //! Assignment operator (disabled)
    SearchSimple & operator = (const SearchSimple &) = delete;

    //! Destructor
    ~SearchSimple();

    /*!
      \brief returns source mesh cells potentially overlapping given target cell

      \param cellId the index of the target cell that I pass in...
      \param candidates pointer to vector of potential candidate cells in sourceMesh
    */
    void search(const int cellId, std::vector<int> *candidates) const;

  private:

    // Aggregate data members
    const SourceMeshWrapper & sourceMesh_;
    const TargetMeshWrapper & targetMesh_;
    double* xlow_;
    double* xhigh_;
    double* ylow_;
    double* yhigh_;

}; // class SearchSimple




template<typename SourceMeshWrapper, typename TargetMeshWrapper>
SearchSimple<SourceMeshWrapper,TargetMeshWrapper>::
SearchSimple(const SourceMeshWrapper & source_mesh,
             const TargetMeshWrapper & target_mesh)
        : sourceMesh_(source_mesh), targetMesh_(target_mesh)  {

    int numCells = sourceMesh_.num_owned_cells() + sourceMesh_.num_ghost_cells();
    xlow_  = new double[numCells];
    xhigh_ = new double[numCells];
    ylow_  = new double[numCells];
    yhigh_ = new double[numCells];

    // find bounding boxes for all cells
    for (int c = 0; c < numCells; ++c) {
        std::vector<int> nodes;
        sourceMesh_.cell_get_nodes(c,&nodes);
        getBoundingBox<SourceMeshWrapper>(sourceMesh_, nodes,
                &xlow_[c], &xhigh_[c], &ylow_[c], &yhigh_[c]);
    }
} // SearchSimple::SearchSimple


template<typename SourceMeshWrapper, typename TargetMeshWrapper>
SearchSimple<SourceMeshWrapper,TargetMeshWrapper>::~SearchSimple()
{
    delete [] xlow_;
    delete [] xhigh_;
    delete [] ylow_;
    delete [] yhigh_;
} // SearchSimple::~SearchSimple


template<typename SourceMeshWrapper, typename TargetMeshWrapper>
void SearchSimple<SourceMeshWrapper,TargetMeshWrapper>::
search(const int cellId, std::vector<int> *candidates)
const {
    // find bounding box for target cell
    std::vector<int> nodes;
    targetMesh_.cell_get_nodes(cellId,&nodes);
    double txlow, txhigh, tylow, tyhigh;
    getBoundingBox<TargetMeshWrapper>(targetMesh_, nodes, &txlow, &txhigh, &tylow, &tyhigh);
    
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
