/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "search_simple.h"

#include "Mesh.hh"

#include <cstdio>
#include <algorithm>


namespace { // unnamed
void getBoundingBox(
    const Jali::Mesh* mesh,
    const Jali::Entity_ID_List& nodes,
    double* xlow, double* xhigh,
    double* ylow, double* yhigh)
{
    const double big = 1.e99;
    double xl = big;  double xh = -big;
    double yl = big;  double yh = -big;

    for (int n = 0; n < nodes.size(); ++n) {
        JaliGeometry::Point p;
        mesh->node_get_coordinates(nodes[n], &p);
        double x = p.x();
        double y = p.y();
        xl = std::min(xl, x);  xh = std::max(xh, x);
        yl = std::min(yl, y);  yh = std::max(yh, y);
    }

    *xlow = xl;  *xhigh = xh;
    *ylow = yl;  *yhigh = yh;

} // getBoundingBox

} // namespace


namespace Portage {

SearchSimple::SearchSimple(const Jali::Mesh* sourceMesh,
        const Jali::Mesh* targetMesh)
    : sourceMesh_(sourceMesh), targetMesh_(targetMesh)
{
    int numCells = sourceMesh_->num_entities(Jali::CELL, Jali::ALL);
    xlow_  = new double[numCells];
    xhigh_ = new double[numCells];
    ylow_  = new double[numCells];
    yhigh_ = new double[numCells];

    // find bounding boxes for all cells
    for (int c = 0; c < numCells; ++c) {
        Jali::Entity_ID_List nodes;
        sourceMesh_->cell_get_nodes(c, &nodes);
        getBoundingBox(sourceMesh_, nodes,
                &xlow_[c], &xhigh_[c], &ylow_[c], &yhigh_[c]);
    }
} // SearchSimple::SearchSimple


SearchSimple::~SearchSimple()
{
    delete [] xlow_;
    delete [] xhigh_;
    delete [] ylow_;
    delete [] yhigh_;
} // SearchSimple::~SearchSimple


void SearchSimple::search(const Jali::Entity_ID cellId, Jali::Entity_ID_List* candidates)
const
{
    // find bounding box for target cell
    Jali::Entity_ID_List nodes;
    targetMesh_->cell_get_nodes(cellId, &nodes);
    double txlow, txhigh, tylow, tyhigh;
    getBoundingBox(targetMesh_, nodes, &txlow, &txhigh, &tylow, &tyhigh);

    // now see which sourceMesh cells have bounding boxes overlapping
    // with target cell
    // do a naive linear search
    int numCells = sourceMesh_->num_entities(Jali::CELL, Jali::ALL);
    for (int c = 0; c < numCells; ++c) {
        if (std::max(txlow, xlow_[c]) < std::min(txhigh, xhigh_[c]) &&
                std::max(tylow, ylow_[c]) < std::min(tyhigh, yhigh_[c])) {
            candidates->push_back(c);
        }
    }

} // SearchSimple::search

} // namespace Portage

/*-------------------------------------------------------------------------~--*
 * Formatting options for Emacs and vim.
 *
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:nil
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *-------------------------------------------------------------------------~--*/
