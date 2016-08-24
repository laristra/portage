/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef FLAT_STATE_WRAPPER_H_
#define FLAT_STATE_WRAPPER_H_

#include <map>

#include "portage/support/portage.h"

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
  template <class State_Wrapper>
  Flat_State_Wrapper(State_Wrapper &input, std::vector<std::string> var_names) 
  {
    for (unsigned int i=0; i<var_names.size(); i++)
    {
      Entity_kind entity = input.get_entity(var_names[i]);
      int dataSize = input.get_data_size(entity, var_names[i]);
      T* data;
      input.get_data(entity, var_names[i], &data);
      std::vector<T> field;
      field.resize(dataSize);
      std::copy(data, data+dataSize, field.begin());
      state_.push_back(field);
      name_map_[var_names[i]] = state_.size() - 1;
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
  void get_data(const Entity_kind on_what, const std::string var_name, D** const data) const {
  
    std::map<std::string, int>::const_iterator iter = name_map_.find(var_name);
    if (iter != name_map_.end())
      (*data) = (D*)(&(state_[iter->second][0]));
  }

  /*!
    @brief Get the data vector
  */
  std::vector<T>& get_vector(int index) 
  {
    if (index < state_.size())
      return state_[index]; 
    else
      return gradients_[index - state_.size()];
  }

  int get_field_dim(int index)
  {
    if (index < state_.size()) return 1;
    return 3;
  }

  /*! 
    @brief Get the number of data vectors
  */
  int get_num_vectors() { return state_.size(); }

  /*!
    @brief Add a gradient field
  */
  void add_gradient(std::vector<std::vector<T>>& new_grad)
  {
    if (new_grad.size() <= 0) return;
    std::vector<T> new_vector(new_grad.size()*new_grad[0].size());
    for (unsigned int i=0; i<new_grad.size(); i++)
      std::copy(new_grad[i].begin(), new_grad[i].end(), new_vector.begin() + new_grad[i].size()*i);
    gradients_.push_back(new_vector);
  }

  int get_num_gradients() { return gradients_.size(); }

  //std::vector<T>& get_gradient(int index)
  //{
  //  return gradients_[index];
  //}

private:
  std::vector<std::vector<T>> state_;
  std::map<std::string, int> name_map_;
  std::vector<std::vector<T>> gradients_;


}; // Flat_State_Wrapper

} // namespace Portage

#endif // FLAT_STATE_WRAPPER_H_
