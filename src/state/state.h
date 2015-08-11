/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef NGC_REMAP_STATE_H
#define NGC_REMAP_STATE_H

/*!
  \class State state.h
  \brief State is a class that stores all of the state data associated 
  with a mesh
*/
#include <iostream>
#include "Mesh.hh"

#include "state_vector.h"


namespace Portage {
    
class State {
 public:
      
  // Constructors
  
  State(Jali::Mesh const * const mesh) : mymesh_(mesh) {}
      
  // Copy constructor (disabled)
  
  State(const State &) = delete;
  
  // Assignment operator (disabled)
  
  State & operator=(const State &) = delete; 
  
  // Destructor
  
  ~State() {}
  
  // Add state vector - returns 1 if successfully added, 0 otherwise
  
  int add(std::string const name, Jali::Entity_kind const on_what, 
          double const * const data);
  int add(StateVector const & vector);  // discouraged - copies data
  
  
  // Iterators for going through all the state vectors
  
  // How can we have specialized iterators for each entity kind (NODE,
  // EDGE, FACE, CELL)
  
  typedef std::vector<StateVector>::iterator iterator;
  typedef std::vector<StateVector>::const_iterator const_iterator;
  
  iterator begin() { return state_vectors_.begin(); };
  iterator end() { return state_vectors_.end(); };
  const_iterator cbegin() const { return state_vectors_.begin(); }
  const_iterator cend() const { return state_vectors_.end(); }
  
  // references to state vectors and the [] operator
  
  typedef StateVector& reference;
  typedef StateVector const& const_reference;
  reference operator[](int i) { return state_vectors_[i]; }
  const_reference operator[](int i) const { return state_vectors_[i]; }
  
  int size() const {return state_vectors_.size();}
  
  
  // Find state vector by name and what type of entity it is on.
  
  const_iterator find(std::string const name, Jali::Entity_kind const on_what) const;
  

  //  friend std::ostream& operator<<(std::ostream os, State& s) const;
  
 private:
  
  Jali::Mesh const * const mymesh_;
  std::vector<StateVector> state_vectors_; // May want to have separate ones
                                           // for each kind of entity

};



///// ONCE CINCH IS ABLE TO BUILD TESTS WITH LIBRARY FILES AND MULTIPLE SOURCE FILES, WE CAN PUT THESE ROUTINES INTO state.cc

inline
std::ostream & operator<<(std::ostream & os, State const & s) {
  State::const_iterator it = s.cbegin();
  while (it != s.cend()) {
    StateVector const & vec = *it;
    os << vec << std::endl;
    ++it;
  }
}

// add a vector to the state object

inline
int State::add(std::string const name, Jali::Entity_kind const on_what,
               double const * const data) {

  if (find(name,on_what) == cend()) {
    // a search of the state vectors by name and kind of entity turned up
    // empty, so add the vector to the list; if not, warn about duplicate
    // state data
    
    // does this syntax cause copying of data from temporary StateVector
    // to the object that is created in state_vectors_ ? If it does, we
    // can use the C++11 emplace to construct in place

    state_vectors_.push_back(StateVector(name,on_what,mymesh_,data));
    return 1;
  }
  else {      
    // found a state vector by same name living on the same kind of entity
    
    std::cerr << "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
    return 0;
  }
  
  return 0;
}


// add a vector to the state object from an already defined vector -
// discouraged as it involves a copying of data

inline
int State::add(StateVector const & vector) {

  if (mymesh_ != vector.mesh()) {
    std::cerr << "State and StateVector not defined on the same mesh" << std::endl;
    return 0;
  }
  
  if (find(vector.name(),vector.on_what()) == cend()) {
    // a search of the state vectors by name and kind of entity turned up
    // empty, so add the vector to the list; if not, warn about duplicate
    // state data
    
    state_vectors_.push_back(vector);
    return 1;
  }
  else {      
    // found a state vector by same name living on the same kind of entity
    
    std::cerr << "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
    return 0;
  }
  
  return 0;
}


// find a state vector by name and what type of entity it is on

inline
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
  

} // namespace Portage
  
#endif
