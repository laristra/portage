/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef INTERSECT_H
#define INTERSECT_H

/*!
    \class Intersect intersect.h
    \brief Intersect provides...
 */
class Intersect
{
public:

    //! Default constructor
    Intersect() {}

    //! Copy constructor (disabled)
    Intersect(const Intersect &) = delete;

    //! Assignment operator (disabled)
    Intersect & operator = (const Intersect &) = delete;

    //! Destructor
     ~Intersect() {}

    /*!
        \brief This method is...
     */
    void intersect();

private:

    // Aggregate data members
//    double val_;

}; // class Intersect

#endif // INTERSECT_H

/*--------------------------------------------------------------------------~-*
 * Formatting options for Emacs and vim.
 *
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *--------------------------------------------------------------------------~-*/
