/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef REMAP_1STORDER_H
#define REMAP_1STORDER_H

/*!
  \class Remap_1stOrder remap_1st_order.h
  \brief Remap_1stOrder does a 1st order remap of scalars
*/

#include "Mesh.hh"
#include "portage/state/state.h"

namespace Portage {

class Remap_1stOrder {
 public:
  
  Remap_1stOrder(Jali::Mesh const & sourceMesh, State const & sourceState) :
  sourceMesh_(sourceMesh), sourceState_(sourceState)
  {}


  //! Copy constructor (disabled)
  Remap_1stOrder(const Remap_1stOrder &) = delete;
  
  //! Assignment operator (disabled)
  Remap_1stOrder & operator = (const Remap_1stOrder &) = delete;

  //! Destructor
  ~Remap_1stOrder() {}

  
  // Remap functor

  double operator() (Jali::Entity_ID const target_cell, 
                     std::string const & remap_var_name,
                     Jali::Entity_ID_List const & source_cells,
                     std::vector<double> const & weights) const;

 private:

  Jali::Mesh const & sourceMesh_;
  State const & sourceState_;

};

} // namespace portage
                     
                     
#endif
