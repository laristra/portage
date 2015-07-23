/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef REMAP_H
#define REMAP_H

#include "intersect.h"

/*!
    \class Remap remap.h
    \brief Remap provides...
 */
class Remap
{
public:

    //! Default constructor
    Remap(Intersect* isect) {}

    //! Copy constructor (disabled)
    Remap(const Remap &) = delete;

    //! Assignment operator (disabled)
    Remap & operator = (const Remap &) = delete;

    //! Destructor
     ~Remap() {}

    /*!
        \brief This method is...
     */
    void remap();

private:

    // Aggregate data members
//    double val_;

}; // class Remap

#endif // REMAP_H

/*--------------------------------------------------------------------------~-*
 * Formatting options for Emacs and vim.
 *
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *--------------------------------------------------------------------------~-*/
