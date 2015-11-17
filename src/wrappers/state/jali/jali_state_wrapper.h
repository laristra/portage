/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_WRAPPER_H_
#define JALI_STATE_WRAPPER_H_

#include "Mesh.hh"       // Jali mesh declarations
#include "JaliState.h"  // Jali-based state manager declarations

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
    if (vector != 0) (*data) = ((T*)(vector->get_data()));
  
  }

  int get_entity(const std::string var_name) const {

    std::shared_ptr<Jali::BaseStateVector> vector = *(jali_state_.find(var_name, Jali::ANY_KIND));
    if (vector != 0) return vector->on_what();
    return -1;

  }

  const std::type_info& get_type(const std::string var_name) const {
    
    std::shared_ptr<Jali::BaseStateVector> vector = *(jali_state_.find(var_name,Jali::ANY_KIND));
    if (vector != 0) return vector->get_type();
    return typeid(0);

  }

  //! Iterator on vector names
  std::vector<std::string>::iterator names_begin() { return jali_state_.names_begin(); }
  std::vector<std::string>::iterator names_end()   { return jali_state_.names_end();   }

  //! Iterator on vector names of specific entity type
  typedef Jali::State::string_permutation string_permutation;
  string_permutation names_entity_begin(int const on_what) { return jali_state_.names_entity_begin((Jali::Entity_kind)on_what); }
  string_permutation names_entity_end(int const on_what)   { return jali_state_.names_entity_end((Jali::Entity_kind)on_what); }

 private:

  Jali::State & jali_state_;

}; // Jali_State_Wrapper

#endif // JALI_STATE_WRAPPER_H_
