/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_STATE_STATE_VECTOR_UNI_H_
#define SRC_STATE_STATE_VECTOR_UNI_H_

#include <string>
#include <typeinfo>
#include <memory>

#include "portage/support/portage.h"
#include "portage/state/state_vector_base.h"

namespace Portage {  

template <class T=double>
class StateVectorUni : public StateVectorBase {

 public:
  
  StateVectorUni(
  	std::string name, 
  	Entity_kind kind=Entity_kind::CELL,
  	T* pdata=nullptr
  ) : StateVectorBase(name, Field_type::MESH_FIELD, kind),pdata_(pdata) {}
  

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
		
	/// Get a shared pointer to the data
  T* get_data() { return pdata_; }
 	
 private:
 
 	T* pdata_;
 
};

}

#endif //SRC_STATE_STATE_VECTOR_UNI_H_

