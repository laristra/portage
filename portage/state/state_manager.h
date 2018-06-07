/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_STATE_STATE_MANAGER_H_
#define SRC_STATE_STATE_MANAGER_H_

#include <string>
#include <iostream>
#include <map>
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
#include "portage/state/state_vector_uni_raw.h"
#include "portage/state/state_vector_multi_raw.h"

namespace Portage {  

template <class MeshWrapper>
class StateManager {

	public:
	 		
 		StateManager (const MeshWrapper& mesh, 
 			std::unordered_map<std::string,int> names={},
 			std::unordered_map<int,std::vector<int>> material_cells={}
 			):mesh_(mesh), num_cells_(mesh_.num_owned_cells()){
 				add_material_names(names);
 				add_material_cells(material_cells);
 		}
 		
		void add(std::shared_ptr<StateVectorBase> sv){
		
			// create the key
			std::string  key{sv->name()};
			
			// if the map entry for this key already exists, throw and erro
			if (state_vectors_.find(key)!=state_vectors_.end()) 
				throw std::runtime_error("Field " + sv->name() + 
					" already exists in the state manager");
					
			// copies the shared pointer
 			state_vectors_[key]=sv;
		}
		
		// specialization for adding StateVectorUni's
		// we need to do size checking and this is the easiest way
		template <class T>
		void add(std::shared_ptr<StateVectorUni<T>> sv){
		
			// create the key
			std::string key{sv->name()};
			
			// if the map entry for this key already exists, throw and erro
			if (state_vectors_.find(key)!=state_vectors_.end()) 
				throw std::runtime_error("Field " + sv->name() + 
					" already exists in the state manager");
					
			// check that the size is correct
			if (sv->get_data().size()!=num_cells_) throw std::runtime_error(
				"The added data did not have the same number of elements (" + 
				std::to_string(sv->get_data().size()) +") as the number of cells ("+
				std::to_string(num_cells_)+")");
				
			// copies the shared pointer
 			state_vectors_[key]=sv;
		}
		
		// specialization for adding StateVectorMulti's
		// we need to do size checking and this is the easiest way
		template <class T>
		void add(std::shared_ptr<StateVectorMulti<T>> sv){
		
			// create the key
			std::string key{sv->name()};
			
			// if the map entry for this key already exists, throw an error
			// Note: the follow is BAD, as it creates the key even if the data is bad
			//if (state_vectors_[key]!=nullptr) 
			if (state_vectors_.find(key)!=state_vectors_.end()) 
				throw std::runtime_error("Field " + sv->name() + 
					" already exists in the state manager");
					
			// check that the size is correct
			if (!shape_is_good<T>(sv)) throw std::runtime_error(
				"The shape of the data was not the same as the shape of the material cells");
			
			// copies the shared pointer
 			state_vectors_[key]=sv;
		}
		
		// is the data the correct shape
		template <class T=double>
		bool shape_is_good(const std::shared_ptr<StateVectorMulti<T>> sv){
		
			// if we don't have cells yet, and don't know the shape, assume the
			// user knows what they are doing and the shape is good.
			if(material_data_shape_.size()==0)return true; 
			
			// get the actual data out of the state vector
			const auto& data = sv->get_data();
			
			// check first that there are the correct number of materials
			if (material_data_shape_.size()!= data.size()) return false;
			
			for (const auto& kv : data){			
				if (
					material_data_shape_.find(kv.first)==material_data_shape_.end() ||
					material_data_shape_[kv.first]!= kv.second.size()
				) return false;
			}
			
			return true;
		}
					
 		// get the data and do the cast internal to the function
 		// note this call always works even if the key doesn't exist. We can check
 		// the return against a nullptr
 		template <class T, template<class> class StateVectorType>
		std::shared_ptr<StateVectorType<T>> get(std::string name){
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return std::dynamic_pointer_cast<StateVectorType<T>>(state_vectors_[key]);
 			
		}	
		
 		// get the data and do the cast internal to the function
 		// note this call always works even if the key doesn't exist. We can check
 		// the return against a nullptr
 		std::shared_ptr<StateVectorBase> get(std::string name){
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return state_vectors_[key];
 			
		}	
/**/		
 		// convenience function for getting StateVectorMulti
 		// note this call always works even if the key doesn't exist. We can check
 		// the return against a nullptr
/* 		template <class T=double>
		std::shared_ptr<StateVectorMulti<T>> get(std::string name){
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return std::dynamic_pointer_cast<StateVectorMulti<T>>(state_vectors_[key]);
 			
		}	
*/
		
		////////////////////
		// Const Versions
		////////////////
/* 		template <class T=double>
		std::shared_ptr<StateVectorMulti<T>> get(std::string name) const{
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return std::dynamic_pointer_cast<StateVectorMulti<T>>(state_vectors_.at(key));
 			
		}	
*/		
 		template <class T, template<class> class StateVectorType>
		std::shared_ptr<StateVectorType<T>> get(std::string name) const{
		
			// create the key
			std::string key{name};
			
			if (state_vectors_.find(key)==state_vectors_.end()){
				return nullptr;
			}			
			
			// do the return all in one step
			return std::dynamic_pointer_cast<StateVectorType<T>>(state_vectors_.at(key));
 			
		}	
				
		// add the names of the materials
		void add_material_names(std::unordered_map<std::string,int> names){
			material_ids_=names;
			for (auto& kv: names) material_names_[kv.second]=kv.first;				
		}
		
		// add the material cells
		void add_material_cells(std::unordered_map<int,std::vector<int>> cells){
		
			std::unordered_map<int,int> shape;
			
			// make sure all materials are known by index 
			for (auto &kv : cells){
				if ( material_names_.find(kv.first)==material_names_.end()){
					throw std::runtime_error("Tried to add cells for an unknown material: "
						+ std::to_string(kv.first));
				}
				
				// set the shape of the material cells
				shape[kv.first]=kv.second.size();
				
				// create the inverse map
				for (int cell : kv.second){
				
					// add the material to the cell
					cell_materials_[cell].insert(kv.first);
				}
			}
			
			// set the 
			material_cells_=cells;
			material_data_shape_=shape;
		}

		const std::unordered_map<int,int> get_material_shape() const {
				return material_data_shape_;
		}
			
		// get the number of materials
		int get_num_materials() const {return material_names_.size();}
		int num_materials() const {return material_names_.size();}

		// get a material name from its id
		std::string get_material_name(int id) {return material_names_[id];}
		
		// get the material id from its name
		int get_material_id(std::string name) {return material_ids_[name];}
		
		// Get the number of cells in the material
  	int num_material_cells(int m) {return material_data_shape_[m];}

		// Get the cells in the m'th material. If no materials have been defined
		// return reference to a dummy vector thats empty
		std::vector<int> const& get_material_cells(int m)  {
			return material_cells_[m];
		}
		
		// Get the number of materials in a cell
		int num_cell_materials(int c) {    	
    	return cell_materials_[c].size();
  	}
  	
  	// Get the materials in a cell
  	std::unordered_set<int> get_cell_materials(int c) {
  		return cell_materials_[c];
  	}
  	
  	// Get the cell index in material
  	int cell_index_in_material(int c, int m) const{
  	
  		// get the cell ids for a material
  		std::vector<int> cells{material_cells_.at(m)};
  		
  		// get an iterator to the element
  		auto it = std::find(cells.begin(),cells.end(),c);
  		
  		//return the distance
  		return std::distance(cells.begin(),it);
  		
  	}
  	
  	// Get the Entity_kind for a name
  	// NOTE: This may not be unique, a name can be reused, but for the moment...
  	Entity_kind get_entity(std::string name) const{
  	
  		// for right now, just loop and check that the name is correct
  		for (auto& kv: state_vectors_){
  			if (kv.first==name) return kv.second->kind();
  		}
  		
  		// default to CELL data, could also throw an error for not found...
  		return Entity_kind::CELL;
  	}
  	
  	Field_type field_type(Entity_kind kind, std::string name) const{
  		return state_vectors_.at(name)->type();
  	}

		// Get the states in the states in state manager
		std::vector<std::string> get_state_keys(){
			std::vector<std::string> keys;
			for (auto& kv: state_vectors_){
				keys.push_back(kv.first);
			} 
			return keys;
		}
		
		void mat_get_cells(int matid, std::vector<int> *matcells) const {
			// the "this" is for a gcc compiler error about assert and templates
		  assert(this->material_cells_.find(matid)!=this->material_cells_.end());
		  matcells->clear();
		  *matcells = material_cells_.at(matid);
		}
		
		template <class T>
		void mat_get_celldata(std::string const& var_name, int matid,
		                      T const **data) const {
		  *data= get<T,StateVectorMulti>(var_name)->get_data()[matid].data();                
		}

  	int cell_get_num_mats(int cellid) const {
    	return cell_materials_.at(cellid).size();
  	}
	
		void cell_get_mats(int cellid, std::vector<int> *cellmats) const {
		  cellmats->clear();
		  auto& mats=cell_materials_.at(cellid);
		  *cellmats = std::vector<int>{mats.begin(),mats.end()};
		} 
		
		template <class T>
		void mesh_get_data(Entity_kind on_what, std::string const& var_name,
		                   T const **data) const {
		  *data = get<T,StateVectorUni>(var_name)->get_data().data();		                 
		}
		
		template <class T>
		void mesh_get_data(Entity_kind on_what, std::string const& var_name,
		                   T **data) {
		  *data = get<T,StateVectorUni>(var_name)->get_data().data();		                 
		}
		
		std::string material_name(int matid) const {
     return material_names_.at(matid);
		} 
		
		void mat_add_cells(int matid, std::vector<int> const& newcells) {
		
			// assign the cell ids to the material (copy constructor
		  material_cells_[matid]=newcells;
		  
		  // need to update cell_materials_ inverse map
		  for (int cellid:newcells)
		  	cell_materials_[cellid].insert(matid);
		  
		}

		// starts from scratch, add the material and cell id's
		// I think this is problematic in a distributed sense, and I don't hit it
		// in simple_mash app, so for now don't do anythin
		void add_material(
			std::string const& matname, std::vector<int> const& matcells) {
				std::cout<<"in add_material: "<<matname<<std::endl;
			}

		template <class T>
		void mat_add_celldata(std::string const& var_name, int matid,
		                      T const * values) {
		  // get the number of cells for this material
		  // should use api, but I'm going to change it
		  int ncells = material_cells_[matid].size();
		  
		  // could use a std::copy, but I'll get to it later
		  // need a reference,because we are modifying the data in place
		  std::vector<T>& data=get<T,StateVectorMulti>(var_name)->get_data()[matid];
		  data.clear();
		  for (int i=0;i<ncells;++i) data.push_back(*(values+i));
		}
	
		// Discussion points:
		// This function has a side effect. When this function is called, the 
		// resulting pointer needs to point to memory that is already allocated
		// of sufficient size. It cannot be a nullptr. At some point in the code we
		// need to do the allocation. In our case is an std::vector.reserve call. 
		// This memory management could be done in add_material like Jali does. The
		// downside is that this buries the allocation into "hidden" location. I (DWS) 
		// would personally like to see add_material go away. It is responsible for 
		// creating a uid for the material in a distributed environment which is
		// problematic and this memory management which is also problematic. The
		// memory management is hidden here just the same,and my preferred solution
		// is to make the reserve explicit. We should add an allocate api that is
		// called explicitly in mmdriver. That way we know exactly what we are doing
		// when and where. It adds a step, but hopefully that will remove the segfault
		// that developers such as myself got because we didn't know when and where
		// to allocate the target state vectors. 
		template <class T>
		void mat_get_celldata(std::string const& var_name, int matid,
		                      T **data) {
		                      
			// get the correct row of the ragged right data structure
			// NOTE (this bit me) the &, the local reference needs to refer to the 
			// actual data in the state manager, not a copy
			std::vector<T>& material_data = get<T,StateVectorMulti>(var_name)->get_data()[matid];
			
			int n = material_cells_[matid].size();
			
			// if the space in the row isn't already allocated, fix this
			if (material_data.size() != n) material_data.resize(n);	
		  
		  *data= material_data.data();                
		}
		
		std::unordered_map<int,std::vector<int>> const get_material_cells(){
			return material_cells_;
		}

		int mat_get_num_cells(int matid) const {
			return material_cells_.at(matid).size();
		}
		
		int get_data_size(Entity_kind on_what, std::string const& var_name) const {
			return 0; //fix if ever used
		}

	private:
	
		const MeshWrapper& mesh_;
		
		// the number of cells
		int num_cells_;
		
		// I'm using a map instead of an unordered map, because there aren't
		// that many state vectors, and also, to use an unordered_map we need
		// to either provide a key function or a hashing function.
		std::unordered_map<std::string, std::shared_ptr<StateVectorBase>> state_vectors_;
		
		// material ids
		std::unordered_map<std::string,int> material_ids_;
		
		// material names
		std::unordered_map<int, std::string> material_names_;
		
		// cell ids for a material
		std::unordered_map<int,std::vector<int>> material_cells_;
  
		// material id's for a cell
		std::unordered_map<int,std::unordered_set<int>> cell_materials_;
  
  	// material data shape
  	std::unordered_map<int,int> material_data_shape_;
};

}

#endif //SRC_STATE_STATE_MANAGER_H_

