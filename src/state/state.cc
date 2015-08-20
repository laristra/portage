/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "state.h"

namespace Portage {

std::ostream & operator<<(std::ostream & os, State const & s) {
  State::const_iterator it = s.cbegin();
  while (it != s.cend()) {
    StateVector const & vec = *it;
    os << vec << std::endl;
    ++it;
  }
}

// add a vector to the state object

State::reference State::add(std::string const name, 
                            Jali::Entity_kind const on_what,
                            double const * const data) {

  iterator it = find(name,on_what);
  if (it == end()) {
    // a search of the state vectors by name and kind of entity turned up
    // empty, so add the vector to the list; if not, warn about duplicate
    // state data
    
    // does this syntax cause copying of data from temporary StateVector
    // to the object that is created in state_vectors_ ? If it does, we
    // can use the C++11 emplace to construct in place

    state_vectors_.push_back(StateVector(name,on_what,mymesh_,data));
    
    // push back may cause reallocation of the vector so the iterator
    // may not be valid. Use [] operator to get reference to vector

    int nvec = state_vectors_.size();
    return state_vectors_[nvec-1];
  }
  else {      
    // found a state vector by same name living on the same kind of entity
    
    std::cerr << "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
    return *it;
  }
  
}


// add a vector to the state object from an already defined vector -
// discouraged as it involves a copying of data

State::reference State::add(StateVector const & vector) {

  iterator it = find(vector.name(),vector.on_what());
  if (it == end()) {
    // a search of the state vectors by name and kind of entity turned up
    // empty, so add the vector to the list
    
    if (mymesh_ != vector.mesh()) {
      // the input vector is defined on a different mesh? copy the
      // vector data onto a vector defined on mymesh and then add

      StateVector vector_copy(vector.name(),vector.on_what(),mymesh_,
                              &(vector[0]));
      state_vectors_.push_back(vector_copy);
    }
    else {

      state_vectors_.push_back(vector);
      
    }

    // push back may cause reallocation of the vector so the iterator
    // may not be valid. Use [] operator to get reference to vector

    int nvec = state_vectors_.size();
    return state_vectors_[nvec-1];
  }
  else {      
    // found a state vector by same name living on the same kind of entity
    
    std::cerr << "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
    return *it;
  }
  
}


// find a state vector by name and what type of entity it is on

State::iterator State::find(std::string const name, 
                            Jali::Entity_kind const on_what) {
  
  iterator it = state_vectors_.begin();
  while (it != state_vectors_.end()) {
    StateVector const & vector = *it;
    if (vector.name() == name && vector.on_what() == on_what)
      break;
    else
      ++it;
  }
  return it;
  
}


State::const_iterator State::find(std::string const name, 
                                  Jali::Entity_kind const on_what) const {
  
  const_iterator it = state_vectors_.cbegin();
  while (it != state_vectors_.cend()) {
    StateVector const & vector = *it;
    if (vector.name() == name && vector.on_what() == on_what)
      break;
    else
      ++it;
  }
  return it;
  
}
  
} // namespace Portage
