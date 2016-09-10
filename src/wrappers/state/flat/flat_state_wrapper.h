/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#ifndef FLAT_STATE_WRAPPER_H_
#define FLAT_STATE_WRAPPER_H_

#include <map>
#include <memory>

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
      std::shared_ptr<std::vector<T>> field = std::make_shared<std::vector<T>>();
      field->resize(dataSize);
      std::copy(data, data+dataSize, field->begin());
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
      (*data) = (D*)(&((*(state_[iter->second]))[0]));
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
  std::map<std::string, int> name_map_;
  std::vector<std::shared_ptr<std::vector<Portage::Point3>>> gradients_;


}; // Flat_State_Wrapper

} // namespace Portage

#endif // FLAT_STATE_WRAPPER_H_
