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
  Jali_State_Wrapper(Jali::State & jali_state) :
      jali_state_(jali_state) 
  {}

  //! Copy constructor - not a deep copy
  Jali_State_Wrapper(Jali_State_Wrapper & state) :
      jali_state_(state.jali_state_)
  {}

  //! Assignment operator (disabled) - don't know how to implement (RVG)
  Jali_State_Wrapper & operator=(Jali_State_Wrapper const &) = delete;

  //! Empty destructor
  ~Jali_State_Wrapper() {};

  //! Get pointer to scalar data

  void get_data(int const on_what, std::string var_name, double **data) const {

    Jali::State::iterator it;
    it = jali_state_.find(var_name, (Jali::Entity_kind) on_what);
    if (it == jali_state_.cend()) {
      std::cerr << "ERROR: Jali_State_Wrapper::get_data - Variable " <<
          var_name << "living on entity type " << on_what << 
          "not found on mesh" << std::endl;
      exit(-1);
    }

    Jali::StateVector & statevec = *it;
    *data = &(statevec[0]);

  } // get_data

 private:

  Jali::State & jali_state_;

}; // Jali_State_Wrapper


#endif // JALI_STATE_WRAPPER_H_
