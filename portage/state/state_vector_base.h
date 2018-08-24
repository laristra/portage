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

/*!
	This class implements a base class state vector. A StateVectorBase instance
	holds metadata common for all field types, such as the name of the field and
	on what mesh entity the field is defined. The primary reason for this base class
	is so that we can implement a single container that holds different conceptual
	field types such as single or multi material fields with different internal
	representations. We need a common base class pointer so we can store all the
	different types together. We can recast the pointer to the correct type whenever
	needed.
*/class StateVectorBase {

 public:
  
		/*!
		  @brief Constructor with a name
		  @param name[in]   Name of the StateVector
		  @param type[in]   Type of the StateVector (Portage::Field_type) (Single/Multi)
		  @param kind[in]   Kind of StateVector (CELL, NODE)
		*/  
		explicit StateVectorBase(std::string name, Portage::Field_type type,
			Portage::Entity_kind kind=Portage::Entity_kind::CELL) :
		    name_(name), type_(type), kind_(kind) {}


		//! Destructor

		virtual ~StateVectorBase() {}

		//! Virtual methods

		virtual std::ostream & print(std::ostream & os) const {
		  os << "Print not implemented for data type of StateVectorBase\n";
		  return os;
		}
		
		
		virtual const std::type_info& data_type() = 0;

		//! Query Metadata

		/*!
			@brief Return the name of the state vector.
			@return the name of the state vector

			Return the name of the state vector.
		*/
		std::string get_name() const { return name_; }


		/*!
			@brief Return the field type [MESH_FIELD, MULTIMATERIAL_FIELD] of the state vector.
			@return the Portage::Field_type [MESH_FIELD, MULTIMATERIAL_FIELD] of the state vector

			Return the field type [MESH_FIELD, MULTIMATERIAL_FIELD] of the state vector.
		*/
		Portage::Field_type get_type() const { return type_; }
		
		
		/*!
			@brief Return the entity kind [CELL, NODE] of the state vector.
			@return the Portage::Entity_kind [CELL, NODE] of the state vector

			Return the entity kind [CELL, NODE] of the state vector.
		*/
		Portage::Entity_kind get_kind() const { return kind_; }
		
	
 protected:
 
		std::string name_;
		Portage::Field_type type_;
		Portage::Entity_kind kind_;

};

}

#endif //SRC_STATE_STATE_VECTOR_BASE_H_

