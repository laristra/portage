/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "portage/support/portage.h"
#include "portage/state/state_vector_uni.h"

#include "gtest/gtest.h"

using namespace Portage;

TEST(StateUni, BasicInt) {

	std::string name{"field"};
	Portage::Field_type field_type{Field_type::MESH_FIELD};
	
	StateVectorUni<int> v(name, field_type);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), field_type);
  ASSERT_EQ(v.data_type(), typeid(int));
  ASSERT_NE(v.data_type(), typeid(double));
  
}

TEST(StateUni, BasicDouble1) {

	std::string name{"field"};
	Portage::Field_type field_type{Field_type::MESH_FIELD};
	
	StateVectorUni<> v(name, field_type);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), field_type);
  ASSERT_EQ(v.data_type(), typeid(double));
  ASSERT_EQ(v.get_data(), nullptr);
  
}



TEST(StateUni, DataAccess) {

	std::string name{"field"};
	Portage::Field_type field_type{Field_type::MESH_FIELD};	
	std::vector<double> data {1.,2.,3.};	
	
	StateVectorUni<double> sv(name, field_type, data.data());
	
	const double const * sv_data{sv.get_data()};

	for (int i=0; i<data.size(); i++) {
		ASSERT_EQ(sv_data[i],data[i]);
	}  
  
}


