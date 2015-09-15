/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_H_
#define JALI_STATE_H_

/*!
  \class State jali_state.h
  \brief State is a class that stores all of the state data associated 
  with a mesh
*/
#include <iostream>
#include "Mesh.hh"    // Jali mesh header

#include "jali_state_vector.h"  // Jali-based state vector


namespace Jali {
    
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
  
  iterator find(std::string const name, Jali::Entity_kind const on_what);
  const_iterator find(std::string const name, Jali::Entity_kind const on_what) const;
  

  // Add state vector - returns reference to added StateVector
  
  reference add(std::string const name, Jali::Entity_kind const on_what, 
                double const * const data);
  reference add(StateVector const & vector);  // discouraged - copies data
  
  
  //  friend std::ostream& operator<<(std::ostream os, State& s) const;
  
 private:
  
  Jali::Mesh const * const mymesh_;
  std::vector<StateVector> state_vectors_; // May want to have separate ones
                                           // for each kind of entity

};


std::ostream & operator<<(std::ostream & os, State const & s);


} // namespace Jali
  
#endif  // JALI_STATE_H_
