/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef DRIVER_H
#define DRIVER_H

/*!
    \class Driver driver.h
    \brief Driver provides...
 */
class Driver
{
public:

    //! Default constructor
    Driver() {}

    //! Copy constructor (disabled)
    Driver(const Driver &) = delete;

    //! Assignment operator (disabled)
    Driver & operator = (const Driver &) = delete;

    //! Destructor
     ~Driver() {}

    /*!
        \brief This method is ...
     */
    void run();

private:

    // Aggregate data members
//    double val_;

}; // class Driver

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
