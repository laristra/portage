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
#include <memory>

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"

namespace Portage {  

template <class MeshWrapper>
class StateManager {

	public:
	
 		using pair_t = std::pair<std::string, Entity_kind>;
 		
 		StateManager (const MeshWrapper& mesh):mesh_(mesh){}
		
		void print_counts(){
			std::cout<< "happy days"<<std::endl;
			std::cout<< "mesh #cells: "<< mesh_.num_owned_cells()<<std::endl;
			std::cout<< "mesh #nodes: "<< mesh_.num_owned_nodes()<<std::endl;
		}
 
	private:
	
		const MeshWrapper& mesh_;
		
		//std::unordered_map<pair_t, std::shared_ptr<StateVectorBase>> state_vectors_;
	
  
};

}

#endif //SRC_STATE_STATE_MANAGER_H_

