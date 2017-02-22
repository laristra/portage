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



#ifndef FLAT_STATE_WRAPPER_H_
#define FLAT_STATE_WRAPPER_H_

#include <map>
#include <memory>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

/*!
  @file flat_state_wrapper.h
  @brief Wrapper for interfacing with the Flat state manager
 */

namespace Portage {

/*!
  @class Flat_State_Wrapper "flat_state_wrapper.h"
  @brief Stores state data in a flat representation
         
         Currently all fields must be of the same type
*/
template <class T=double>
class Flat_State_Wrapper {
 public:

  /*!
    @brief Constructor of Flat_State_Wrapper
   */
  Flat_State_Wrapper() { };

  template <class State_Wrapper>
  void initialize(State_Wrapper const & input,
                  std::vector<std::string> var_names) 
  {
    for (unsigned int i=0; i<var_names.size(); i++)
    {
      std::string varname = var_names[i];  // get_data wants const string
      Entity_kind entity = input.get_entity(varname);
      int dataSize = input.get_data_size(entity, varname);
      T const *data;
        
      input.get_data(entity, varname, &data);
      std::shared_ptr<std::vector<T>> field = std::make_shared<std::vector<T>>();
      field->resize(dataSize);
      std::copy(data, data+dataSize, field->begin());
      state_.push_back(field);
      entities_.push_back(entity);
      name_map_[varname] = state_.size() - 1;
    }
  }
  
  
  /*!
    @brief Assignment operator (disabled) - don't know how to implement (RVG)
   */
  Flat_State_Wrapper & operator=(Flat_State_Wrapper const &) = delete;
  
  /*!
    @brief Empty destructor
   */
  ~Flat_State_Wrapper() {};
  
   
  /*!
    @brief Get pointer to scalar data
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of data
   */
  template <class D>
  void get_data(const Entity_kind on_what, const std::string var_name, D** data) {
  
    std::map<std::string, int>::const_iterator iter = name_map_.find(var_name);
    if (iter != name_map_.end())
      (*data) = (D *)(&((*(state_[iter->second]))[0]));
  }

  /*!
    @brief Get pointer to const scalar data
    @param[in] on_what The entity type on which to get the data
    @param[in] var_name The string name of the data field
    @param[in,out] data A pointer to an array of const data
   */
  template <class D>
  void get_data(const Entity_kind on_what, const std::string var_name, D const **data) const {
  
    std::map<std::string, int>::const_iterator iter = name_map_.find(var_name);
    if (iter != name_map_.end())
      (*data) = (D const *)(&((*(state_[iter->second]))[0]));
  }

  /*!
    @brief Get the entity type
    @param[in] index The index of the data field
    @return The Entity_kind enum for the entity type on which the field is defined
   */
  Entity_kind get_entity(const int index) const
  {
    return entities_[index];
  }

  /*!
    @brief Get the data vector
  */
  std::shared_ptr<std::vector<T>> get_vector(int index) 
  {
    return state_[index]; 
  }

  /*!
    @brief Get gradients
  */
  std::shared_ptr<std::vector<Portage::Point3>> get_gradients(int index)
  {
    return gradients_[index];
  }

  /*! 
    @brief Get the number of data vectors
  */
  int get_num_vectors() { return state_.size(); }

  /*!
    @brief Add a gradient field
  */
  void add_gradients(std::shared_ptr<std::vector<Portage::Point3>> new_grad)
  {
    if (new_grad->size() <= 0) return;
    gradients_.push_back(new_grad); 
  }

  /*! 
    @brief Get field dimension
  */
  int get_field_dim(int index)
  {
    return 1;
  }

  /*!
    @brief Get the number of gradient vectors
  */
  int get_num_gradients() { return gradients_.size(); }

private:
  std::vector<std::shared_ptr<std::vector<T>>> state_;
  std::vector<Entity_kind> entities_;
  std::map<std::string, int> name_map_;
  std::vector<std::shared_ptr<std::vector<Portage::Point3>>> gradients_;


}; // Flat_State_Wrapper

} // namespace Portage

#endif // FLAT_STATE_WRAPPER_H_
