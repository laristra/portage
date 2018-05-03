/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_STATE_STATE_VECTOR_MULTI_H_
#define SRC_STATE_STATE_VECTOR_MULTI_H_

#include <string>
#include <typeinfo>
#include <memory>

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"

namespace Portage {  

template <class T=double>
class StateVectorMultiRaw : public StateVectorBase {

 public:
  
  StateVectorMultiRaw(
  	std::string name, 
  	T** hdata=nullptr
  ) : StateVectorBase(name, Field_type::MULTIMATERIAL_FIELD, Entity_kind::CELL), 
  			hdata_(hdata) {}


  //! Destructor
  ~StateVectorMultiRaw() {}
  
  // print
  std::ostream & print(std::ostream & os) const {
    os << "StateVectorMultiRaw\n";
    return os;
  }
  
  // get the data type
  const std::type_info& data_type() {
  	const std::type_info& ti = typeid(T);
		return ti;
	}
		
	/// Get a shared pointer to the data
  T** get_data() { return hdata_; }
 	
 private:
 
 	T** hdata_;
 
};

}

#endif //SRC_STATE_STATE_VECTOR_MULTI_H_

