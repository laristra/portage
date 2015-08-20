/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef REMAP_H
#define REMAP_H

#include "portage/intersect/intersect.h"

/*!
    \class Remap remap.h
    \brief Remap provides...
 */

#include "Mesh.hh"
#include "portage/state/state.h"

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
    // double remap(std::string const & remap_var_name, Jali::Entity_ID cellId, Jali::Entity_ID_List candidates, std::vector<float> moments);
    double remap(std::string const & remap_var_name, Jali::Entity_ID cellId, Jali::Entity_ID_List candidates, std::vector<std::vector<std::vector<double> > > moments) const;

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
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *--------------------------------------------------------------------------~-*/
