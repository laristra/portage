/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef SEARCH_SIMPLE_H
#define SEARCH_SIMPLE_H

#include "Mesh.hh"

/*!
    \class SearchSimple search_simple.h
    \brief SearchSimple provides...
 */

namespace Portage {

    class SearchSimple
{
public:

    //! Default constructor (disabled)
    SearchSimple() = delete;

    //! Constructor with Meshes
    /*!
      \brief Builds the search structure
      
      \param sourceMesh pointer to the mesh we would like to remap
      \param targetMesh pointer to the mesh to which we would like to remap
      
      This method does something useful...
    */
    SearchSimple(const Jali::Mesh* sourceMesh, 
            const Jali::Mesh* targetMesh);

    //! Copy constructor (disabled)
    SearchSimple(const SearchSimple &) = delete;

    //! Assignment operator (disabled)
    SearchSimple & operator = (const SearchSimple &) = delete;

    //! Destructor
    ~SearchSimple();

    /*!
      \brief This method does...

      \param cellId the index of the target cell that I pass in...
      \param candidates pointer to vector of potential candidate cells in sourceMesh

      This method does something useful...
    */
    void search(const Jali::Entity_ID cellId, Jali::Entity_ID_List* candidates) const;

private:

    // Aggregate data members
    const Jali::Mesh* sourceMesh_;
    const Jali::Mesh* targetMesh_;
    double* xlow_;
    double* xhigh_;
    double* ylow_;
    double* yhigh_;

}; // class SearchSimple

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
