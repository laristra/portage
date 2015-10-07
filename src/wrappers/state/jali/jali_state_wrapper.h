/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_WRAPPER_H_
#define JALI_STATE_WRAPPER_H_

#include "Mesh.hh"       // Jali mesh declarations
#include "jali_state.h"  // Jali-based state manager declarations

/*!
  \class Jali_State_Wrapper  jali_state_wrapper.h
  \brief Jali_State_Wrapper provides access to data stored in Jali_State
*/

class Jali_State_Wrapper {
 public:

  //! Constructor
  Jali_State_Wrapper(Jali::State & jali_state) : jali_state_(jali_state) {}
  
  //! Copy constructor - not a deep copy
  Jali_State_Wrapper(Jali_State_Wrapper & state) : jali_state_(state.jali_state_) {}
  
  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Jali_State_Wrapper & operator=(Jali_State_Wrapper const &) = delete;
  
  //! Empty destructor
  ~Jali_State_Wrapper() {};
  
  //! Get pointer to scalar data

  template <class T>
  void get_data(int const on_what, std::string var_name, T** data) const {
  
    std::shared_ptr<Jali::BaseStateVector> vector = *(jali_state_.find(var_name, (Jali::Entity_kind) on_what));
    if (vector != 0) (*data) = ((T*)(vector->getData()));
  
  }

 private:

  Jali::State & jali_state_;

}; // Jali_State_Wrapper

#endif // JALI_STATE_WRAPPER_H_
