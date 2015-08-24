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
  
  Remap_1stOrder(Jali::Mesh const & sourceMesh, State const & sourceState,
                 std::string remap_var_name, Jali::Entity_kind on_what) :
      sourceMesh_(sourceMesh), sourceState_(sourceState), 
      remap_var_name_(remap_var_name), on_what_(on_what),
      source_var_ptr_(NULL)
  {
    State::const_iterator it;
    it = sourceState_.find(remap_var_name_, on_what_);
    if (it == sourceState_.cend()) {
      std::cerr << "ERROR: Remap_1stOrder::Remap_1stOrder(...) - Remap variable not found in source mesh" << std::endl;
      exit(-1);
    }

    StateVector const & source_var_ref = *it;
    source_var_ptr_ = &(source_var_ref);
  }


  //! Copy constructor (disabled)
  Remap_1stOrder(const Remap_1stOrder &) = delete;
  
  //! Assignment operator (disabled)
  Remap_1stOrder & operator = (const Remap_1stOrder &) = delete;

  //! Destructor
  ~Remap_1stOrder() {}

  
  // Remap functor - Need to make a pair from the vector of source
  // cells that contribute to the the target cell value and the
  // contribution weights associated with each source cell

  double operator() (std::pair< std::vector<int> const &, std::vector<double> const &>) const;

 private:

  Jali::Mesh const & sourceMesh_;
  State const & sourceState_;
  std::string const remap_var_name_;
  Jali::Entity_kind const on_what_;
  StateVector const * source_var_ptr_;

};

} // namespace portage
                     
                     
#endif
