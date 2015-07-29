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

#include "Mesh.hh"
#include "state.h"

namespace Portage {

class Remap
{
public:

    //! Default constructor

  Remap(Jali::Mesh const & sourceMesh, State const & sourceState, 
        Jali::Mesh const & targetMesh, State & targetState) :
      sourceMesh_(sourceMesh), sourceState_(sourceState), 
      targetMesh_(targetMesh), targetState_(targetState) {}

    //! Copy constructor (disabled)
    Remap(const Remap &) = delete;

    //! Assignment operator (disabled)
    Remap & operator = (const Remap &) = delete;

    //! Destructor
     ~Remap() {}

    /*!
        \brief This method is...
     */
  void remap(std::vector<std::string> const & remap_var_names);

private:

  Jali::Mesh const & sourceMesh_;
  Jali::Mesh const & targetMesh_;
  State const & sourceState_;
  State & targetState_;

}; // class Remap

} // namespace Portage

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
