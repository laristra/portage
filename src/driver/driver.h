/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef DRIVER_H
#define DRIVER_H

#include "Mesh.hh"

/*!
    \class Driver driver.h
    \brief Driver provides the API to mapping from one mesh to another.
 */

namespace NGC {
namespace Remap {

class Driver
{
public:

    //! Default constructor creates meshes for testing remap
    Driver() {}

    //! Constructor uses established meshes for remap
    Driver(const Jali::Mesh& inputMesh, Jali::Mesh& targetMesh) 
      : inputMesh_(&inputMesh),
        targetMesh_(&targetMesh) {}

    //! Copy constructor (disabled)
    Driver(const Driver &) = delete;

    //! Assignment operator (disabled)
    Driver & operator = (const Driver &) = delete;

    //! Destructor
     ~Driver() {}

    /*!
        \brief This method calls the search, intersect, and remap routines
	needed to map one mesh to another.
     */
    void run();

    Jali::Mesh* targetMesh() { return targetMesh_; }

private:

    // Aggregate data members
    const Jali::Mesh* inputMesh_;
    Jali::Mesh* targetMesh_;

}; // class Driver

} // namespace Remap
} // namespace NGC

#endif // DRIVER_H

/*--------------------------------------------------------------------------~-*
 * Formatting options for Emacs and vim.
 *
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *--------------------------------------------------------------------------~-*/
