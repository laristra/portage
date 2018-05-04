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

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"

namespace Portage {  

template <class MeshWrapper>
class StateManager {

	public:
	
 		using pair_t = std::pair<Portage::Entity_kind, std::string>;
 		
 		StateManager (const MeshWrapper& mesh, 
 			std::unordered_map<std::string,int> names={},
 			std::unordered_map<int,std::vector<int>> material_cells={}
 			):mesh_(mesh), num_cells_(mesh_.num_owned_cells()){
 				add_material_names(names);
 				add_material_cells(material_cells);
 		}
 		
		void print_counts(){
			std::cout<< "happy days"<<std::endl;
			std::cout<< "mesh #cells: "<< mesh_.num_owned_cells()<<std::endl;
			std::cout<< "mesh #nodes: "<< mesh_.num_owned_nodes()<<std::endl;
		}
		
		void add(std::shared_ptr<StateVectorBase> sv){
		
			// create the key
			pair_t key{sv->kind(), sv->name()};
			
			// if the map entry for this key already exists, throw and erro
			if (state_vectors_[key]!=nullptr) 
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
			pair_t key{sv->kind(), sv->name()};
			
			// if the map entry for this key already exists, throw and erro
			if (state_vectors_[key]!=nullptr) 
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
			pair_t key{sv->kind(), sv->name()};
			
			// if the map entry for this key already exists, throw and erro
			if (state_vectors_[key]!=nullptr) 
				throw std::runtime_error("Field " + sv->name() + 
					" already exists in the state manager");
					
			// check that the size is correct
			if (!shape_is_good<T>(sv)) throw std::runtime_error(
				"The shape of the data was not the same as the shape of the material cells");
				
			// copies the shared pointer
 			state_vectors_[key]=sv;
		}
		
		// get with no template parameters returns the base vector, can be cast
		// by client
		bool get(Entity_kind kind, std::string name, 
			std::shared_ptr<StateVectorBase>& pv){
		
			// create the key
			pair_t key{kind, name};
			
			// lookup the vector
			pv = state_vectors_[key];
						
			return pv==nullptr;
 			
		}	
		
		// is the data the correct shape
		template <class T>
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
 		template <class T, template<class> class StateVectorType>
		bool get(Entity_kind kind, std::string name, 
			std::shared_ptr<StateVectorType<T>>& pv){
		
			// create the key
			pair_t key{kind, name};
			
			// lookup the vector
			auto bv = state_vectors_[key];
					
			// cast to the desired type
			pv = std::dynamic_pointer_cast<StateVectorType<T>>(bv);
			
			return pv==nullptr;
 			
		}	
		
 		// get the data and do the cast internal to the function
 		// note this call always works even if the key doesn't exist. We can check
 		// the return against a nullptr
 		template <class T, template<class> class StateVectorType>
		std::shared_ptr<StateVectorType<T>> get(Entity_kind kind, std::string name){
		
			// create the key
			pair_t key{kind, name};
			
			// do the return all in one step
			return std::dynamic_pointer_cast<StateVectorType<T>>(state_vectors_[key]);
 			
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


		
	private:
	
		const MeshWrapper& mesh_;
		
		// the number of cells
		int num_cells_;
		
		// I'm using a map instead of an unordered map, because there aren't
		// that many state vectors, and also, to use an unordered_map we need
		// to either provide a key function or a hashing function.
		std::map<pair_t, std::shared_ptr<StateVectorBase>> state_vectors_;
		
		// material names
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

