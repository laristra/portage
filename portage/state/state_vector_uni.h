/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_STATE_STATE_VECTOR_UNI_H_
#define SRC_STATE_STATE_VECTOR_UNI_H_

#include <string>
#include <typeinfo>
#include <vector>

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"

namespace Portage {  

template <class T=double>
class StateVectorUni : public StateVectorBase {

 public:
  
		StateVectorUni(
			std::string name, 
			Entity_kind kind=Entity_kind::CELL,
			std::vector<T> data=std::vector<T>()
		) : StateVectorBase(name, Field_type::MESH_FIELD, kind),data_(data) {}
		

		//! Destructor
		~StateVectorUni() {}
		
		// print
		std::ostream & print(std::ostream & os) const {
		  os << "UniStateVector\n";
		  return os;
		}
		  
		// get the data type
		const std::type_info& data_type() {
			const std::type_info& ti = typeid(T);
			return ti;
		}
		
		/*!
			@brief Return a reference to the data in the state vector.
			@return a reference to the vector of data in the state vector

			Return a reference to the data in the state vector.
		*/
		std::vector<T>& get_data() { return data_; }
 	
 private:
 
 		std::vector<T> data_;
 
};

}

#endif //SRC_STATE_STATE_VECTOR_UNI_H_

