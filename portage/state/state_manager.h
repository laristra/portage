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
#include <memory>
#include <stdexcept>

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"

namespace Portage {  

template <class MeshWrapper>
class StateManager {

	public:
	
 		using pair_t = std::pair<Portage::Entity_kind,std::string>;
 		
 		StateManager (const MeshWrapper& mesh):mesh_(mesh){}
 		
 		StateManager (const MeshWrapper& mesh, std::shared_ptr<StateVectorBase> sv)
 			:mesh_(mesh) {
 				add(sv);
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
				throw std::runtime_error("Field already exists");
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

	private:
	
		const MeshWrapper& mesh_;
		
		// I'm using a map instead of an unordered map, because there aren't
		// that many state vectors, and also, to use an unordered_map we need
		// to either provide a key function or a hashing functions.
		std::map<pair_t, std::shared_ptr<StateVectorBase>> state_vectors_;
	
  
};

}

#endif //SRC_STATE_STATE_MANAGER_H_

