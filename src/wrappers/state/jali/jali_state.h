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

namespace std
{
    typedef decltype(nullptr) nullptr_t;
}

#include <boost/iterator/permutation_iterator.hpp>

#include "Mesh.hh"    // Jali mesh header

#include "jali_state_vector.h"  // Jali-based state vector


namespace Jali {

#define NUM_ENTITY_TYPES 5

class State
{
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

  typedef std::vector<std::shared_ptr<BaseStateVector>>::iterator iterator;
  typedef std::vector<std::shared_ptr<BaseStateVector>>::const_iterator const_iterator;
  typedef boost::permutation_iterator<std::vector<std::shared_ptr<BaseStateVector>>::iterator, std::vector<int>::iterator> permutation_type;

  iterator begin() { return state_vectors_.begin(); };
  iterator end() { return state_vectors_.end(); };
  const_iterator cbegin() const { return state_vectors_.begin(); }
  const_iterator cend() const { return state_vectors_.end(); }

  // Iterators for specific entity types

  permutation_type entity_begin(Jali::Entity_kind entityAssociation) { return boost::make_permutation_iterator(state_vectors_.begin(), entity_indexes[entityAssociation].begin()); }
  permutation_type entity_end(Jali::Entity_kind entityAssociation) { return boost::make_permutation_iterator(state_vectors_.begin(), entity_indexes[entityAssociation].end()); }

  // References to state vectors and the [] operator
  
  typedef std::shared_ptr<BaseStateVector> pointer;
  typedef const std::shared_ptr<BaseStateVector> const_pointer;
  pointer operator[](int i) { return state_vectors_[i]; }
  const_pointer operator[](int i) const { return state_vectors_[i]; }
  
  int size() const {return state_vectors_.size();}


  // Find state vector by name and what type of entity it is on.

  iterator find(std::string const name, Jali::Entity_kind const on_what) {

    iterator it = state_vectors_.begin();
    while (it != state_vectors_.end()) {
      BaseStateVector const & vector = *(*it);
      if (vector.name() == name && vector.on_what() == on_what)
        break;
      else
        ++it;
    }
    return it;

  }

  // Find state vector by name and what type of entity it is on - const version

  const_iterator find(std::string const name, Jali::Entity_kind const on_what) const {

    const_iterator it = state_vectors_.cbegin();
    while (it != state_vectors_.cend()) {
      BaseStateVector const & vector = *(*it);
      if (vector.name() == name && vector.on_what() == on_what)
        break;
      else
        ++it;
    }
    return it;

  }


  // Add state vector - returns reference to added StateVector

  template <class T>
  StateVector<T> & add(std::string const name, Jali::Entity_kind const on_what, T* data) {
  
    iterator it = find(name,on_what);
    if (it == end()) {
      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list; if not, warn about duplicate
      // state data
    
      // does this syntax cause copying of data from temporary StateVector
      // to the object that is created in state_vectors_ ? If it does, we
      // can use the C++11 emplace to construct in place

      std::shared_ptr<StateVector<T>> vector(new StateVector<T>(name, on_what, mymesh_, data));
      state_vectors_.push_back(vector);
    
      entity_indexes[on_what].push_back(state_vectors_.size()-1);

      // push back may cause reallocation of the vector so the iterator
      // may not be valid. Use [] operator to get reference to vector

      return (*vector);
    }
    else {      
      // found a state vector by same name living on the same kind of entity
    
      std::cerr << "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
      return (*(std::static_pointer_cast<StateVector<T>>(*it)));
    }

  };

  // Add state vector - discouraged, copies data

  template <class T>
  StateVector<T> & add(StateVector<T> & vector) {

    iterator it = find(vector.name(),vector.on_what());
    if (it == end()) {
      
      // a search of the state vectors by name and kind of entity turned up
      // empty, so add the vector to the list
      if (mymesh_ != vector.mesh()) {
      
        // the input vector is defined on a different mesh? copy the
        // vector data onto a vector defined on mymesh and then add
        std::shared_ptr<StateVector<T>> vector_copy(new StateVector<T>(vector.name(),vector.on_what(),mymesh_,&(vector[0])));
        state_vectors_.push_back(vector_copy);
      }
      else {
        std::shared_ptr<StateVector<T>> vector_copy(new StateVector<T>(vector));
        state_vectors_.push_back(vector_copy);
      }

      entity_indexes[vector.on_what()].push_back(state_vectors_.size()-1);

      // push back may cause reallocation of the vector so the iterator
      // may not be valid. Use [] operator to get reference to vector

      int nvec = state_vectors_.size();
      return (*(std::static_pointer_cast<StateVector<T>>(state_vectors_[nvec-1])));
    }
    else {      
      // found a state vector by same name living on the same kind of entity
    
      std::cerr << "Attempted to add duplicate state vector. Ignoring\n" << std::endl;
      return vector;
    }

  };

 private:
  
  Jali::Mesh const * const mymesh_;
  std::vector<std::shared_ptr<BaseStateVector>> state_vectors_;
  std::vector<int> entity_indexes[NUM_ENTITY_TYPES];

};


std::ostream & operator<<(std::ostream & os, State const & s);

} // namespace Jali
  
#endif  // JALI_STATE_H_
