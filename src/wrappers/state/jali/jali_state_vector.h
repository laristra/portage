/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_VECTOR_H_
#define JALI_STATE_VECTOR_H_

/*!
  \class StateVector jali_state_vector.h
  \brief StateVector provides a mechanism to store state data for mesh entities
*/

#include <iostream>
#include <vector>

#include "Mesh.hh"    // jali mesh header

namespace Jali {

// Provides some limited functionality of a std::vector while adding
// some additional meta-data like the mesh associated with this data.
// Cannot inherit from std::vector because we want this to work only
// type 'double'

class StateVector {
 public:
  
  // Constructors
  
  StateVector(std::string const name, 
              Entity_kind const on_what, 
              Mesh const * const mesh,
              double const * const data) : 
      myname_(name), on_what_(on_what), mymesh_(mesh) {
    
    int num = mymesh_->num_entities(on_what,ALL);
    mydata_.resize(num);
    std::copy(data, data+num, mydata_.begin());
    
  }

  // Copy constructor

  StateVector(StateVector const & in_vector) : 
      myname_(in_vector.myname_), on_what_(in_vector.on_what_),
      mymesh_(in_vector.mymesh_), mydata_(in_vector.mydata_) {
  }

  // Assignment operator

  StateVector & operator=(StateVector const & in_vector) {
    myname_ = in_vector.myname_;
    on_what_ = in_vector.on_what_;
    mymesh_ = in_vector.mymesh_;
    mydata_ = in_vector.mydata_;
  }


  // Destructor
  
  ~StateVector() {}


  // Query Metadata
  
  std::string name() const { return myname_; }
  Jali::Entity_kind on_what() const { return on_what_; }
  Jali::Mesh const * mesh() const { return mymesh_; }
  
  // Subset of std::vector functionality. We can add others as needed
  
  typedef std::vector<double>::iterator iterator;
  typedef std::vector<double>::const_iterator const_iterator;
  
  iterator begin() { return mydata_.begin(); };
  iterator end() { return mydata_.end(); };
  const_iterator cbegin() const { return mydata_.begin(); }
  const_iterator cend() const { return mydata_.end(); }
  
  typedef double& reference;
  typedef double const& const_reference;
  reference operator[](int i) { return mydata_[i]; }
  const_reference operator[](int i) const { return mydata_[i]; }
  
  int size() const {return mydata_.size();}
  void resize(size_t n, double val=0.0) { mydata_.resize(n,val);}
  
  void clear() { mydata_.clear(); }
  
 private:
  
  std::string myname_;
  Entity_kind on_what_;
  Mesh const * mymesh_;
  std::vector<double> mydata_;
  
};


inline
std::ostream & operator<<(std::ostream & os, StateVector const & sv) {
  int size = sv.size();
  os << std::endl;
  os << "Vector \"" << sv.name() << "\" on entity kind " << sv.on_what() << " : " << std::endl;
  os << size << " elements " << std::endl;
  for (int i = 0; i < size; i++)
    os << sv[i] << std::endl;
  os << std::endl;
}


} // namespace Jali
  

#endif  // JALI_STATE_VECTOR_H_
