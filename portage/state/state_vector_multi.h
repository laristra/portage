/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_STATE_STATE_VECTOR_MULTI_H_
#define SRC_STATE_STATE_VECTOR_MULTI_H_

#include <string>
#include <typeinfo>
#include <vector>
#include <unordered_map>

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"

namespace Portage {  

template <class T=double>
class StateVectorMulti : public StateVectorBase {

	public:
  
		StateVectorMulti(
			std::string name, 
			std::unordered_map<int, std::vector<T>> data = std::unordered_map<int, std::vector<T>>()
		) : StateVectorBase(name, Field_type::MULTIMATERIAL_FIELD, Entity_kind::CELL), 
					data_(data) {}


		//! Destructor
		~StateVectorMulti() {}

		// print
		std::ostream & print(std::ostream & os) const {
			os << "StateVectorMulti\n";
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
		std::unordered_map<int, std::vector<T>>& get_data() { return data_; }
		
		
		/*!
			@brief Return a reference to the data in the state vector.
			@param[in] m	material id
			@return a reference to the vector of data in the state vector

			Return a reference to the data in the state vector. Note: this function
			has a potential side effect (which is a good thing, but users need to be
			aware). If the vector for key m does not exist, it will be created automatically.
			If you want to throw an error, use get_data().at(m) instead.
		*/
		std::vector<T>& get_data(int m) { return data_[m]; }

	private:
 
 	std::unordered_map<int, std::vector<T>> data_;
 
};

}

#endif //SRC_STATE_STATE_VECTOR_MULTI_H_

