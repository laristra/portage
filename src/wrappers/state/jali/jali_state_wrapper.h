/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



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
  void get_data(const Entity_kind on_what, const std::string var_name, T **data) {
  
    Jali::State::const_iterator it =
        jali_state_.find<T, Jali::Mesh>(var_name, jali_state_.mesh(),
                                        (Jali::Entity_kind) on_what);
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::BaseStateVector> vector = *it;
      if (vector) {
        (*data) = ((T *)(vector->get_raw_data()));
        return;
      }
    }

    std::cerr << "Could not find state variable " << var_name << "\n";
    (*data) = nullptr;
  }

  /*!
    @brief Get pointer to const scalar data
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of const data
   */
  template <class T>
  void get_data(const Entity_kind on_what, const std::string var_name,
                T const **data) const {
  
    Jali::State::const_iterator it =
        jali_state_.find<T, Jali::Mesh>(var_name, jali_state_.mesh(),
                                        (Jali::Entity_kind) on_what);
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::BaseStateVector> vector = *it;
      if (vector) {
        (*data) = ((T const *)(vector->get_raw_data()));
        return;
      }
    }

    std::cerr << "Could not find state variable " << var_name << "\n";
    (*data) = nullptr;
  }


  /*!
    @brief Get the entity type on which the given field is defined
    @param[in] var_name The string name of the data field
    @return The Entity_kind enum for the entity type on which the field is defined

    @todo  THIS ASSUMES ONLY DOUBLE VECTORS - WE HAVE TO ACCOUNT FOR OTHER TYPES
           OR WE HAVE TO GENERALIZE THE FIND FUNCTION!!!
   */
  Entity_kind get_entity(const std::string var_name) const {

    // ****** CHANGE WHEN JALISTATE find FUNCTION IS FIXED TO NOT NEED
    // ****** THE DATA TYPE

    Jali::State::const_iterator it =
        jali_state_.find<double, Jali::Mesh>(var_name, jali_state_.mesh());
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::BaseStateVector> vector = *it;
      if (vector)
        return (Portage::Entity_kind) vector->entity_kind();
    }

    std::cerr << "Could not find state variable3 " << var_name << "\n";
    return Portage::UNKNOWN_KIND;
  }


  /*!
  @brief Get the data size for the given field
  @param[in] on_what  The entity type on which the data field is defined
  @param[in] var_name The string name of the data field
  @return The data size for the field with the given name on the given entity type
 
  @todo  THIS ASSUMES ONLY DOUBLE VECTORS - WE HAVE TO ACCOUNT FOR OTHER TYPES
         OR WE HAVE TO GENERALIZE THE FIND FUNCTION!!!
   */
  int get_data_size(const Entity_kind on_what, const std::string var_name) const {

    // ****** CHANGE WHEN JALISTATE find FUNCTION IS FIXED TO NOT NEED
    // ****** THE DATA TYPE

    Jali::State::const_iterator it =
        jali_state_.find<double, Jali::Mesh>(var_name, jali_state_.mesh(),
                                        (Jali::Entity_kind) on_what);
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::BaseStateVector> vector = *it;
      if (vector)
        return (vector->size());
    }

    std::cerr << "Could not find state variable4 " << var_name << "\n";
    return 0;
  }


#if 0
  /*!
    @brief Get the data type of the given field
    @param[in] var_name The string name of the data field
    @return A reference to the type_info struct for the field's data type
   */ 
  const std::type_info& get_type(const std::string var_name) const {
    
    Jali::State::const_iterator it =
        jali_state_.find(var_name, jali_state_.mesh());
    if (it != jali_state_.cend()) {
      std::shared_ptr<Jali::BaseStateVector> vector = *it;
      if (vector)
        return vector->get_type();
    }

    std::cerr << "Could not find state variable " << var_name << "\n";
    return typeid(0);
  }
#endif

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
