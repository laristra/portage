/*---------------------------------------------------------------------------~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef INTERSECT_H
#define INTERSECT_H

#include "search.h"

/*!
    \class Intersect intersect.h
    \brief Intersect provides...
 */
#include "Mesh.hh"

namespace Portage {

class Intersect
{
 public:

  //! Default constructor
  Intersect(Jali::Mesh const & sourceMesh, Jali::Mesh const & targetMesh) :
      sourceMesh_(sourceMesh), targetMesh_(targetMesh) {}

  //! Copy constructor (disabled)
  Intersect(const Intersect &) = delete;
  
  //! Assignment operator (disabled)
  Intersect & operator = (const Intersect &) = delete;
  
  //! Destructor
  ~Intersect() {}
  
  /*!
    \brief This method is...
  */
  void intersect(Jali::Entity_ID cellId, Jali::Entity_ID_List* candidates, std::vector<float>* moments);
  
 private:

  Jali::Mesh const & sourceMesh_;
  Jali::Mesh const & targetMesh_;

}; // class Intersect

} // namespace Portage

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
