/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef JALI_STATE_WRAPPER_H_
#define JALI_STATE_WRAPPER_H_

#include "portage/support/portage.h"

#include "Mesh.hh"       // Jali mesh declarations
#include "JaliState.h"  // Jali-based state manager declarations

/*!
  @file jali_state_wrapper.h
  @brief Wrapper for interfacing with the Jali state manager
 */

namespace Portage {

/*!
  @class Jali_State_Wrapper "jali_state_wrapper.h"
  @brief Provides access to data stored in Jali_State
*/
class Jali_State_Wrapper {
 public:

  /*!
    @brief Constructor of Jali_State_Wrapper
    @param[in] jali_state A reference to a Jali::State instance 
   */
  Jali_State_Wrapper(Jali::State & jali_state) : jali_state_(jali_state) {}
  
  /*!
    @brief Copy constructor of Jali_State_Wrapper - not a deep copy
    @param[in] state A reference to another Jali_State_Wrapper instance
   */
  Jali_State_Wrapper(Jali_State_Wrapper & state) : jali_state_(state.jali_state_) {}
  
  /*!
    @brief Assignment operator (disabled) - don't know how to implement (RVG)
   */
  Jali_State_Wrapper & operator=(Jali_State_Wrapper const &) = delete;
  
  /*!
    @brief Empty destructor
   */
  ~Jali_State_Wrapper() {};
  
   
  /*!
    @brief Initialize fields from mesh file
   */
  void init_from_mesh() { jali_state_.init_from_mesh(); }

  /*!
    @brief Export fields to mesh file
   */
  void export_to_mesh() {jali_state_.export_to_mesh(); }

  /*!
    @brief Get pointer to scalar data
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of data
   */
  template <class T>
  void get_data(const Entity_kind on_what, const std::string var_name, T** const data) const {
  
    std::shared_ptr<Jali::BaseStateVector> vector = 
        *(jali_state_.find(var_name, (Jali::Entity_kind) on_what));
    if (vector != 0) (*data) = ((T*)(vector->get_data()));
  
  }

  /*!
    @brief Get the entity type on which the given field is defined
    @param[in] var_name The string name of the data field
    @return The Entity_kind enum for the entity type on which the field is defined
   */
  Entity_kind get_entity(const std::string var_name) const {

    std::shared_ptr<Jali::BaseStateVector> vector = 
        *(jali_state_.find(var_name, Jali::ANY_KIND));
    if (vector != 0) return (Portage::Entity_kind) vector->on_what();

    return Portage::UNKNOWN_KIND;

  }

  /*!
    @brief Get the data type of the given field
    @param[in] var_name The string name of the data field
    @return A reference to the type_info struct for the field's data type
   */ 
  const std::type_info& get_type(const std::string var_name) const {
    
    std::shared_ptr<Jali::BaseStateVector> vector = 
        *(jali_state_.find(var_name,Jali::ANY_KIND));
    if (vector != 0) return vector->get_type();
    return typeid(0);

  }

  /*!
    @brief Begin iterator on vector names
    @return Begin iterator on vector of strings
   */
  std::vector<std::string>::iterator names_begin() const { 
    return jali_state_.names_begin(); 
  }

  /*!
    @brief End iterator on vector names
    @return End iterator on vector of strings
   */
  std::vector<std::string>::iterator names_end() const { 
    return jali_state_.names_end(); 
  }

  /*!
    @brief Typedef for permutation iterator on vector of strings
   */
  typedef Jali::State::string_permutation string_permutation;

  /*!
    @brief Begin iterator on vector names of specific entity type
    @param[in] on_what The desired entity type
    @return Permutation iterator to start of string vector
   */
  string_permutation names_entity_begin(Entity_kind const on_what) const { 
    return jali_state_.names_entity_begin((Jali::Entity_kind)on_what); 
  }

  /*!
    @brief End iterator on vector of names of specific entity type
    @param[in] on_what The desired entity type
   */
  string_permutation names_entity_end(Entity_kind const on_what) const { 
    return jali_state_.names_entity_end((Jali::Entity_kind)on_what); 
  }

 private:

  Jali::State & jali_state_;

}; // Jali_State_Wrapper

} // namespace Portage

#endif // JALI_STATE_WRAPPER_H_
