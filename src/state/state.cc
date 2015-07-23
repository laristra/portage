/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "state.h"

namespace NGC {
namespace Remap {

// add a vector to the state object

void State::add(StateVector const & vector) {
  
  if (find(vector.name(),vector.on_what()) == cend()) {
    // a search of the state vectors by name and kind of entity turned up
    // empty, so add the vector to the list; if not, warn about duplicate
    // state data
    
    state_vectors_.push_back(vector);
  }
  else {      
    // found a state vector by same name living on the same kind of entity
    
    std::cerr << "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
  }
  
}


// find a state vector by name and what type of entity it is on

State::const_iterator State::find(std::string const name, Jali::Entity_kind const on_what) const {
  
  std::vector<StateVector>::const_iterator it = state_vectors_.cbegin();
  while (it != state_vectors_.cend()) {
    StateVector const & vector = *it;
    if (vector.name() == name && vector.on_what() == on_what)
      break;
    else
      ++it;
  }
  return it;
  
}

} // namespace Remap
} // namespace NGC
