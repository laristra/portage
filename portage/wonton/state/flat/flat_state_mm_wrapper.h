/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef FLAT_STATE_MM_WRAPPER_H_
#define FLAT_STATE_MM_WRAPPER_H_

#include "portage/state/state_manager.h"


/*!
  @file flat_state_mm_wrapper.h
  @brief Definitions for a Flat State Wrapper that can do distributed MM remap
*/
namespace Wonton {

using namespace Portage;
  /*!
    @class Flat_State_Wrapper flat_state_wrapper.h
    @brief A state manager wrapper that allows redistribution of data across
    nodes.
   */
template <class MeshWrapper>
class Flat_State_Wrapper: public StateManager<MeshWrapper> {

 public:
 
 
  /*!
    @brief Constructor for the state wrapper
	  @param[in] mesh							mesh wrapper    
	  @param[in] names        		optional, map from material names to material id
	  @param[in] material_cells   optional, map from material id to vector of cells    

	  Constructor that takes the meshwrapper for the underlying mesh and two optional
	  map arguments: the map from material name to id, and the map from material 
	  id to the cells containing that material
   */
  Flat_State_Wrapper(const MeshWrapper& mesh, 
 			std::unordered_map<std::string,int> names={},
 			std::unordered_map<int,std::vector<int>> material_cells={}
 			) :StateManager<MeshWrapper>(mesh,names,material_cells) { }


  /// Assignment operator (disabled).
  Flat_State_Wrapper & operator=(Flat_State_Wrapper const &) = delete;


  /// Destructor.
  ~Flat_State_Wrapper() { }


  /*!
   * @brief Initialize the state wrapper with another state wrapper and a list of names
   * @param[in] input another state wrapper, which need not be for Flat_State.
   * @param[in] var_names a list of state names to initialize
   *
   * Entities and sizes associated with the given name will be obtained from the input state wrapper.
   *
   * A name can be re-used with a different entity, but a name-entity combination
   * must be unique.
   *
   * A name-entity combination must not introduce a new size for that entity if
   * it has previously been encountered.
   *
   * All existing internal data is forgotten.
   */
  template <class State_Wrapper>
  void initialize(State_Wrapper const & input,
		  std::vector<std::string> var_names)
  {
  
	  clear(); // forget everything

	  for (std::string varname : var_names)
	  {
	  
		  // get the entity kind
		  Entity_kind entity = input.get_entity(varname);

		  // get pointer to data for state from input state wrapper
		  // DO WE NEED TO DO DATA INTROSPECTION TO FIND THE TYPE HERE (DOUBLE)
		  double const* data;
		  input.mesh_get_data(entity, varname, &data);
		  
		  // Note, this get_data_size is of the input wrapper, which is defined
		  // not the base class state manager get_data_size which is not implemented
		  size_t dataSize = input.get_data_size(entity, varname);

	    // create a uni state vector
			std::shared_ptr<StateVectorUni<double>> pv = std::make_shared<StateVectorUni<double>> (
				varname, data, data+dataSize, entity);
				
	    // add to database
	    StateManager<MeshWrapper>::add(pv);

	  }
  }

  
  /*!
    @brief Get the number of data vectors
    @return		The number of state vectors
  */
  size_t get_num_vectors() { return StateManager<MeshWrapper>::state_vectors_.size(); }


  /*!
    @brief Get the data vector
    I am doing a return by reference of the data here. In the flat state manager
    as it existed before this, the state manager returned data vector. In the
    new StateManager, we return a shared pointer to a StateVectorBase that has
    metadata as well. Getting at the data requires one more access that is not
    by shared pointer but by array reference as currently implemented. This breaks
    the current MMDriver code, so I need to change MMDriver as well
  */
  std::vector<double>& get_vector(std::string field_name)
  {
  	std::shared_ptr<StateVectorBase> p = StateManager<MeshWrapper>::get(field_name);
  	return std::static_pointer_cast<StateVectorUni<double>>(p)->get_data();
  	
    
  }
  
  
  /*!
    @brief Get field stride
  */
  size_t get_field_stride(std::string field_name)
  {
  	// FIX
    return 1;
  }

  /*!
    @brief Get the entity type on which the given field is defined
    @param[in] index The index of the data field
    @return The Entity_kind enum for the entity type on which the field is defined
   */
  Entity_kind get_entity(std::string field_name) 
  {
    return StateManager<MeshWrapper>::get(field_name)->get_kind(); 
  }
 
 private:


  void clear()
  {
  	StateManager<MeshWrapper>::state_vectors_.clear();
  	StateManager<MeshWrapper>::material_ids_.clear();
  	StateManager<MeshWrapper>::material_names_.clear();
  	StateManager<MeshWrapper>::material_cells_.clear();
  	StateManager<MeshWrapper>::cell_materials_.clear();
  }
  
};  // class Flat_State_Wrapper

}  // namespace Wonton 

#endif  // FLAT_STATE_MM_WRAPPER_H_
