/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_STATE_STATE_MANAGER_H_
#define SRC_STATE_STATE_MANAGER_H_

#include <string>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <stdexcept>
#include <iterator>
#include <algorithm>

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"
#include "portage/state/state_vector_uni.h"
#include "portage/state/state_vector_multi.h"

namespace Portage {  

template <class MeshWrapper>
class StateManager {

	public:

	 		
		/*!
		  @brief Constructor 
		  @param[in] mesh							mesh wrapper    
		  @param[in] names        		optional, map from material names to material id
		  @param[in] material_cells   optional, map from material id to vector of cells    

		  Constructor that takes the meshwrapper for the underlying mesh and two optional
		  map arguments: the map from material name to id, and the map from material 
		  id to the cells containing that material
		*/
 		StateManager (const MeshWrapper& mesh, 
 			std::unordered_map<std::string,int> names={},
 			std::unordered_map<int,std::vector<int>> material_cells={}
 			):mesh_(mesh){
 				add_material_names(names);
 				add_material_cells(material_cells);
 		}


		//////////////////////////////////////////////////
		// state manager level API functions and material names
		//////////////////////////////////////////////////

		
		/*!
		  @brief Add the names of the materials. 
		  @param[in] names		map from material names to material id

		  Add the names of the materials. The map from material name to id.
		*/
		void add_material_names(std::unordered_map<std::string,int> names){
			material_ids_=names;
			for (auto& kv: names) material_names_[kv.second]=kv.first;				
		}

		
		/*!
		  @brief Return the names registered by the state manager. 
			@return		vector of strings of material names
			
		  Return the names registered by the state manager.
		*/
		std::vector<std::string> const get_state_keys(){
			std::vector<std::string> keys;
			for (auto& kv: state_vectors_){
				keys.push_back(kv.first);
			} 
			return keys;
		}


		/*!
		  @brief Return the name of a material from its material id. 
		  @param[in] m		material id
		  @return					string of the material name

		  Return the name of a material from its material id.
		*/
		std::string material_name(int m) const {return material_names_.at(m); } 


		/*!
		  @brief Return the material id from its name.
		  @param[in] name		material name
		  @return						material id

		  Return the material id from its name.
		*/
		int get_material_id(std::string const& name) const {return material_ids_.at(name);}


		/*!
		  @brief Add the material cells.
		  @param[in] cells		map from material id to vector of cell ids

		  Add the material cells. The argument is a map from material id to a vector
		  of cells for that material.
		*/
		void add_material_cells(std::unordered_map<int,std::vector<int>> cells){
		
			
			// make sure all materials are known by index 
			for (auto &kv : cells){
				if ( material_names_.find(kv.first)==material_names_.end()){
					throw std::runtime_error("Tried to add cells for an unknown material: "
						+ std::to_string(kv.first));
				}
				
				// create the inverse map
				for (int cell : kv.second){
				
					// add the material to the cell
					cell_materials_[cell].insert(kv.first);
				}
			}
			
			// set the 
			material_cells_=cells;
		}


		/*!
		  @brief Return the material cells.
		  @return		map from material id to vector of cell ids

		  Return the material cells. The return value is a map from material id to a vector
		  of cells for that material.
		*/
		std::unordered_map<int,std::vector<int>> const& get_material_cells() const {
			return material_cells_;
		}


		/*!
		  @brief Return the number of cells for each material.
		  @return		map from material id to the number of cells with this material

		  Return the number of cells for each material. The return value is a map 
		  from material id to the integer number of cells.
		*/
		std::unordered_map<int,int> const get_material_shape() const {
		
				std::unordered_map<int,int> shape;
				for (auto& kv: material_cells_) shape[kv.first]=kv.second.size();

				return shape;
		}

			
		/*!
		  @brief Return the number of materials in the problem.
		  @return		the integer number of materials

		  Return the number of materials in the problem.
		*/
		int num_materials() const {return material_names_.size();}


		
		//////////////////////////////////////////////////
		// state vector metadata
		//////////////////////////////////////////////////
 
 		
		
		/*!
		  @brief Return the entity kind for a state vector.
		  @param[in] name		string of the material name
		  @return						entity kind (CELL, NODE) of the state vector

		  Return the entity kind for a state vector. The input parameter is the name
		  of the state vector, and the return type is a Portage::Entity_kind
		*/
		Entity_kind get_entity(std::string const& name) const {
			return state_vectors_.at(name)->get_kind();
		}
  	
  	  	
		/*!
		  @brief Return the field type of a state vector.
		  @param[in] kind		Portage::Entity_kind of the state vector (not needed)
		  	but required by the signature of mmdrive
		  @param[in] name		string of the material name
		  @return 					field type (UNIMATERIAL, MULTIMATERIAL) of the state vector

		  Return the field type of a state value. The kind is not required as the names
		  are unique. The original intent was that there could be e.g. a pressure defined
		  on the cells and nodes but with the same same name "pressure". We moved past
		  this, so if the name is unique, the kind is not required. the arguments are
		  entity kind and the name of the state vector, and the return type is the field
		  type.
		*/
  	Field_type field_type(Entity_kind kind, std::string const& name) const{
  		return state_vectors_.at(name)->get_type();
  	}


		/*!
		  @brief Return the data size .
		  
		  This function is currently not implemented, but is required by the API.
		  Might only make sense for single material vectors, but not sure
		*/
		int get_data_size(Entity_kind on_what, std::string const& name) const {
			return 0; //fix if ever used
		}




		//////////////////////////////////////////////////
		// methods to register or access a complete state vector
		//////////////////////////////////////////////////
 
 		
		/*!
		  @brief Add a StateVectorBase to the state manager.
		  @param[in] sv		shared pointer to StateVectorBase
		  
		  Add a shared pointer to a state vector to the state manager. This specialization
		  works directly with a StateVectorBase
		*/
		void add(std::shared_ptr<StateVectorBase> sv){
		
			// create the key
			std::string  key{sv->get_name()};
			
			// if the map entry for this key already exists, throw an error
			if (state_vectors_.find(key)!=state_vectors_.end()) 
				throw std::runtime_error("Field " + sv->get_name() + 
					" already exists in the state manager");
					
			// copies the shared pointer
 			state_vectors_[key]=sv;
		}

		
		/*!
		  @brief Add a StateVectorUni<T> to the state manager.
		  @param[in] sv		shared pointer to StateVectorUni<T>
		  
		  Adds a shared pointer to a single material state vector to the state manager. 
		  This specialization works with a single material StateVectorUni templated
		  on type T. The code does some error checking, in particular making sure 
		  the the length of the data is the same as the number of mesh entities in the
		  underlying mesh.
		*/
		template <class T>
		void add(std::shared_ptr<StateVectorUni<T>> sv){
		
			// create the key
			std::string key{sv->get_name()};
			
			// if the map entry for this key already exists, throw an error
			if (state_vectors_.find(key)!=state_vectors_.end()) 
				throw std::runtime_error("Field " + sv->get_name() + 
					" already exists in the state manager");
					
			// check that the size is correct
			if (sv->get_data().size()!=mesh_.num_entities(sv->get_kind())) throw std::runtime_error(
				"The added data did not have the same number of elements (" + 
				std::to_string(sv->get_data().size()) +") as the number of cells ("+
				std::to_string(mesh_.num_owned_cells())+")");
				
			// copies the shared pointer
 			state_vectors_[key]=sv;
		}

		
		/*!
		  @brief Add a StateVectorMulti<T> to the state manager.
		  @param[in] sv		shared pointer to StateVectorMulti<T>
		  
		  Adds a shared pointer to a multi material state vector to the state manager. 
		  This specialization works with a multi material StateVectorMulti templated
		  on type T. The code does some error checking, in particular making sure 
		  the the shape of the data is the same as the shape of the material cells.
		*/
		template <class T>
		void add(std::shared_ptr<StateVectorMulti<T>> sv){
		
			// create the key
			std::string key{sv->get_name()};
			
			// if the map entry for this key already exists, throw an error
			// Note: we shouldn't use (state_vectors_[key]!=nullptr) as []creates the 
			// key even if the key doesn't exist
			if (state_vectors_.find(key)!=state_vectors_.end()) 
				throw std::runtime_error("Field " + sv->get_name() + 
					" already exists in the state manager");
					
			// check that the size is correct
			if (!shape_is_good<T>(sv)) throw std::runtime_error(
				"The shape of the data was not the same as the shape of the material cells");
			
			// copies the shared pointer
 			state_vectors_[key]=sv;
		}

		
		/*!
		  @brief Get a shared pointer to a StateVectorBase from the state manager.
		  @param[in] name		name of the registered state vector
		  @return						shared pointer to StateVectorBase
		  
		  Gets a shared pointer to a StateVectorBase from the state manager. The pointer
		  is to the base class and can be dynamically pointer cast to any type. The
		  name is checked with the state manager. If the name is not found, a nullptr
		  is returned so validity can be checked with a nullptr test.		  
		*/
 		std::shared_ptr<StateVectorBase> get(std::string name){
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return state_vectors_[key];
 			
		}	


		/*!
		  @brief Get a shared pointer to a templated type from the state manager.
		  @param[in] name		name of the registered state vector
		  @return						shared pointer to a type T
		  
		  Gets a shared pointer to a type T from the state manager. The pointer will
		  be dynamically cast to type T. If the name is not found, a nullptr
		  is returned so validity can be checked with a nullptr test.	T is the complete
		  type, which is likely itself templated, e.g. T=StateVectorMulti<double>.
		  With this signature we only need one specialization and it is used for
		  accessing both single material and multi material state vectors. 
		*/
 		template <typename T>
		std::shared_ptr<T> get(std::string name){
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return std::dynamic_pointer_cast<T>(state_vectors_[key]);
 			
		}	

		
		/*!
		  @brief Get a shared pointer to a templated type from the state manager.
		  @param[in] name			name of the registered state vector
		  @return							const shared pointer to a type T
		  
		  Gets a shared pointer to a type T from the state manager. The pointer will
		  be dynamically cast to type T. If the name is not found, a nullptr
		  is returned so validity can be checked with a nullptr test.	T is the complete
		  type, which is likely itself templated, e.g. T=StateVectorMulti<double>.
		  With this signature we only need one specialization and it is used for
		  accessing both single material and multi material state vectors. This is a
		  const version of the function, that is used for read only data. 
		  
		*/
 		template <typename T>
		std::shared_ptr<T> get(std::string name) const{
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return std::dynamic_pointer_cast<T>(state_vectors_.at(key));
 			
		}	

				
		
		//////////////////////////////////////////////////
		// material dominant accessor methods
		//////////////////////////////////////////////////

					
		
		/*!
		  @brief Return the number of cells for this material id.
		  @param[in] m		material id
		  @return					number of cells containing this material

		  Return the number of cells for this material id.
		*/
  	int num_material_cells(int m) const {return material_cells_.at(m).size();}


		/*!
		  @brief Return the cell ids for this material id.
		  @param[in] m		material id
		  @return					vector of cell indices containing this material

		  Return a read only reference to the cell ids for this material id. If no 
		  materials have been defined, return a reference to an empty vector.
		*/
		std::vector<int> const& get_material_cells(int m)  const {
			return material_cells_.at(m);
		}


		/*!
		  @brief Return the cell ids for this material id.
		  @param[in] m        		material id
		  @param[out] matcells		vector of cell indices containing this material

		  Return the cell ids for this material id. This is the in-place version of
		  the function. The cells are returned in a function argument rather than a
		  function return.
		*/
		void mat_get_cells(int m, std::vector<int> *matcells) const {
			// the "this" is for a gcc compiler error about assert and templates
		  assert(this->material_cells_.find(m)!=this->material_cells_.end());
		  matcells->clear();
		  *matcells = material_cells_.at(m);
		}

		
		/*!
		  @brief Return the cell data for this material id.
		  @param[in] name			name registered with the state manager
		  @param[in] m 				material id
		  @param[out] data		a read only pointer to the cell vector data 

		  Return the cell data for this material id. This is the in-place version of
		  the function. The data is returned in a location pointed to by the user
		  argument data. This is read-only version of the function. The vector data
		  cannot be modified.
		*/
		template <class T>
		void mat_get_celldata(std::string const& name, int m,
		                      T const **data) const {
		  *data= get<StateVectorMulti<T>>(name)->get_data(m).data();                
		}


		/*!
		  @brief Return the cell index in this material id.
		  @param[in] c		cell id
		  @param[in] m		material id
		  @return					index of cell in material 

		  Return the index (location) of this cell within the list of cells for this
		  material. 
		  
  		Get the cell index in material
  		since this uses material dominant structures, I wish the signature was
  		int cell_index_in_material(int m, int c)
		*/
  	int cell_index_in_material(int c, int m) const {
  	
  		// get the cell ids for a material
  		std::vector<int> cells{material_cells_.at(m)};
  		
  		// get an iterator to the element
  		auto it = std::find(cells.begin(),cells.end(),c);
  		
  		//return the distance
  		return std::distance(cells.begin(),it);
  		
  	}


		/*!
		  @brief Return the cell data for this material id.
		  @param[in] name			name registered with the state manager
		  @param[in] m 				material id
		  @param[out] data		a read only pointer to the cell vector data 

		  Return the cell data for this material id. This is the in-place version of
		  the function. The data is returned in a location pointed to by the user
		  argument data. This is the non-const version of the function. The data is
		  both writable and allocated if it doesn't exist.
		  
			Discussion points:
			This function has a side effect. When this function is called, the 
			resulting pointer needs to point to memory that is already allocated
			of sufficient size. It cannot be a nullptr. At some point in the code we
			need to do the allocation. In our case this is an std::vector.resize call. 
			This memory management could be done in add_material like Jali does. The
			downside is that this buries the allocation into a "hidden" location. I (DWS) 
			would personally like to see add_material go away. It is responsible for 
			creating a uid for the material in a distributed environment which is
			problematic and this memory management which is also problematic. The
			memory management is hidden here just the same,and my preferred solution
			is to make the allocate/resize explicit. We should add an allocate api that is
			called explicitly in mmdriver. That way we know exactly what we are doing
			when and where. It adds a step, but hopefully that will remove the segfault
			that developers such as myself got because we didn't know when and where
			to allocate the target state vectors. 
		  
		*/
		template <class T>
		void mat_get_celldata(std::string const& name, int m,
		                      T **data) {
		                      
			// get the correct row of the ragged right data structure
			// NOTE: (this bit me) the &, the local reference needs to refer to the 
			// actual data in the state manager, not a copy
			std::vector<T>& material_data = get<StateVectorMulti<T>>(name)->get_data(m);
			
			int n = material_cells_[m].size();
			
			// if the space in the row isn't already allocated, fix this
			if (material_data.size() != n) material_data.resize(n);	
		  
		  *data= material_data.data();                
		}



		//////////////////////////////////////////////////
		// cell dominant accessor methods
		//////////////////////////////////////////////////

					
		
		/*!
		  @brief Return the number of materials in this cell.
		  @param[in] c		cell id
		  @return					number of material in this cell

		  Return the number of materials in this cell.
		*/
  	int cell_get_num_mats(int c) const {return cell_materials_.at(c).size();}

  	
		/*!
		  @brief Return the unordered set of materials in this cell.
		  @param[in] c		cell id
		  @return					unordered set of material ids in this cell

		  Return the unordered set of materials in this cell. While cell id's need
		  to be kept in order within a material, the reverse is not true. The order
		  of materials in a cell is arbitrary, so use a data structure that reflects
		  this. Using a set makes comparing, interecting, adding etc. easy.
		*/
  	std::unordered_set<int> get_cell_materials(int c) const {
  		return cell_materials_.at(c);
  	}

		/*!
		  @brief Return the materials in this cell.
		  @param[in] c						cell id
		  @param[out] cellmats		pointer to vector of materials in this cell

		  Return the unordered set of materials in this cell. While cell id's need
		  to be kept in order within a material, the reverse is not true. Reflecting
		  the use of an unordered set for the underlying data structure, the order
		  of the materials is not guaranteed.
		*/
		void cell_get_mats(int c, std::vector<int> *cellmats) const {
		  cellmats->clear();
		  auto& mats=cell_materials_.at(c);
		  *cellmats = std::vector<int>{mats.begin(),mats.end()};
		} 



		//////////////////////////////////////////////////
		// single material accessor methods
		//////////////////////////////////////////////////


		
		/*!
		  @brief Return the const cell data for this material id.
		  @param[in] on_what		Portage::Entity_kind of data (unused, but required
		  	by the API)
		  @param[in] name				name registered with the state manager
		  @param[out] data			a read only pointer to the cell vector data 

		  Return the uni material data for this state vector. This is the in-place version of
		  the function. The data is returned in a location pointed to by the user
		  argument data. This is read-only version of the function. The vector data
		  cannot be modified.
		*/
		template <class T>
		void mesh_get_data(Entity_kind on_what, std::string const& name,
		                   T const **data) const {
		  *data = get<StateVectorUni<T>>(name)->get_data().data();		                 
		}

		
		/*!
		  @brief Return the cell data for this material id.
		  @param[in] on_what		Portage::Entity_kind of data (unused, but required
		  	by the API)
		  @param[in] name				name registered with the state manager
		  @param[out] data			a pointer to the cell vector data 

		  Return the uni material data for this state vector. This is the in-place version of
		  the function. The data is returned in a location pointed to by the user
		  argument data. This is the mutable form of the function. The data may be
		  modified by the caller.
		*/
		template <class T>
		void mesh_get_data(Entity_kind on_what, std::string const& name,
		                   T **data) {
		  *data = get<StateVectorUni<T>>(name)->get_data().data();		                 
		}


		
		//////////////////////////////////////////////////
		// material dominant mutator methods
		//////////////////////////////////////////////////


		/*!
			@brief Add cells for a new material to the material cells.
			@param[in] m					material id
			@param[in] newcells		vector of cell ids containing this material
			
			Add cells for a new material to the material cells. This function takes two
			arguments, the material id and the cells containing it. The function then
			set the cells for this material. Add is a bit of a misnomer. At the moment,
			this sets or replaces. No merge is done.
		*/
		void mat_add_cells(int m, std::vector<int> const& newcells) {
		
			// assign the cell ids to the material (copy constructor
		  material_cells_[m]=newcells;
		  
		  // need to update cell_materials_ inverse map
		  for (int c:newcells)
		  	cell_materials_[c].insert(m);
		  
		}


		/*!
			@brief Add cell data for a new material in a state vector to the state manager.
			@param[in] name 		name registered with the state manager
			@param[in] m				material id
			@param[in] values		vector of data for this material
			
			Add cell data for a new material in a state vector to the state manager. 
			This function takes three arguments, the name of the state vector, the 
			material id and the data itself. The function then copies the data into
			the state vector. Add is a bit of a misnomer. At the moment,
			this sets or replaces. No merge is done.
		*/
		template <class T>
		void mat_add_celldata(std::string const& name, int m,
		                      T const * values) {
		  // get the number of cells for this material
		  // should use api, but I'm going to change it
		  int ncells = material_cells_[m].size();
		  
		  // could use a std::copy, but I'll get to it later
		  // need a reference,because we are modifying the data in place
		  std::vector<T>& data=get<StateVectorMulti<T>>(name)->get_data(m);
		  data.clear();
		  for (int i=0;i<ncells;++i) data.push_back(*(values+i));
		}

	
		/*!
			@brief Add cell data for a new material in a state vector to the state manager.
			@param[in] matname 		name of the material
			@param[in] matcells		vector of material cells 
			
			NOT IMPLEMENTED AT PRESENT:
			Starting from scratch, add the material and cell id's
			I think this is problematic in a distributed sense, and I don't hit it
			in simple_mash app, so for now don't do anything
			the problem is that if a rank creates it's own material id, then those
			could conflict with other ranks. I believe the main caller should have a
			registry of all materials in the problem. It should not be up to a single
			node to create a name<->matid association.
		*/
		void add_material(
			std::string const& matname, std::vector<int> const& matcells) {
				std::cout<<"in add_material: "<<matname<<std::endl;
				throw std::runtime_error("StateManager::add_material is not implemented"); 
			}




		//////////////////////////////////////////////////
		// utility functions
		//////////////////////////////////////////////////
 		
		/*!
		  @brief Check that a StateVectorMulti<T> has a shape consistent with the 
		  	state manager.
		  @param[in] sv		shared pointer to StateVectorMulti<T>
		  
		  Check that a StateVectorMulti<T> has a shape consistent with the state 
		  manager. Both the number of materials and number of cells per material
		  must match in order to pass the test.
		*/
		template <class T=double>
		bool shape_is_good(const std::shared_ptr<StateVectorMulti<T>> sv){
		
			// get the data shape
			std::unordered_map<int,int> shape{get_material_shape()};
			
			// if we don't have cells yet, and don't know the shape, assume the
			// user knows what they are doing and the shape is good.
			if(shape.size()==0)return true; 
			
			// get the actual data out of the state vector
			const auto& data = sv->get_data();
			
			// check first that there are the correct number of materials
			if (shape.size()!= data.size()) return false;
			
			for (const auto& kv : data){			
				if (
					shape.find(kv.first)==shape.end() ||
					shape[kv.first]!= kv.second.size()
				) return false;
			}
			
			return true;
		}

		
	protected:
	
		// underlying mesh reference
		const MeshWrapper& mesh_;
		
		// state vectors
		std::unordered_map<std::string, std::shared_ptr<StateVectorBase>> state_vectors_;
		
		// material name to id
		std::unordered_map<std::string,int> material_ids_;
		
		// material id to name
		std::unordered_map<int, std::string> material_names_;
		
		// cell ids for a material
		std::unordered_map<int,std::vector<int>> material_cells_;
  
		// material id's for a cell
		std::unordered_map<int,std::unordered_set<int>> cell_materials_;
  
};

}

#endif //SRC_STATE_STATE_MANAGER_H_

