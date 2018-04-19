/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef SRC_STATE_STATE_VECTOR_BASE_H_
#define SRC_STATE_STATE_VECTOR_BASE_H_

#include <string>

#include "portage/support/portage.h"

namespace Portage {  

class StateVectorBase {

 public:
  
  /*!
    @brief Constructor with a name
    @param name   Name of the StateVector
    @param type   Type of the StateVector (Portage::Field_type)
  */  
  explicit StateVectorBase(std::string name, Portage::Field_type type) :
      name_(name), type_(type) {}


  //! Destructor

  virtual ~StateVectorBase() {}

  //! Virtual methods

  virtual std::ostream & print(std::ostream & os) const {
    os << "Print not implemented for data type of StateVectorBase\n";
    return os;
  }
  virtual const std::type_info& data_type() = 0;

  //! Query Metadata

	std::string name() const { return name_; }
	Portage::Field_type type() const { return type_; }
	
 protected:
 
  std::string name_;
  Portage::Field_type type_;

};

}

#endif //SRC_STATE_STATE_VECTOR_BASE_H_

