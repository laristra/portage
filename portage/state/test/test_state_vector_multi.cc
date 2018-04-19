/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/



#include "portage/support/portage.h"
#include "portage/state/state_vector_multi.h"

#include "gtest/gtest.h"

using namespace Portage;

TEST(State, BasicInt) {

	std::string name{"field"};
	Portage::Field_type field_type{Field_type::MESH_FIELD};
	
	StateVectorMulti<int> v(name, field_type);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), field_type);
  ASSERT_EQ(v.data_type(), typeid(int));
  ASSERT_NE(v.data_type(), typeid(double));
  
}

TEST(State, BasicDouble1) {

	std::string name{"field"};
	Portage::Field_type field_type{Field_type::MESH_FIELD};
	
	StateVectorMulti<> v(name, field_type);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), field_type);
  ASSERT_EQ(v.data_type(), typeid(double));
  ASSERT_EQ(v.get_data(), nullptr);
  
}



