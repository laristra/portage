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
	
	StateVectorUni<int> v(name);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), Field_type::MESH_FIELD);
  ASSERT_EQ(v.data_type(), typeid(int));
  ASSERT_NE(v.data_type(), typeid(double));
  
}

TEST(StateUni, BasicDouble1) {

	std::string name{"field"};
	
	StateVectorUni<> v(name);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), Field_type::MESH_FIELD);
  ASSERT_EQ(v.data_type(), typeid(double));
  ASSERT_EQ(v.get_data(), std::vector<double>());
  
}



TEST(StateUni, DataAccess) {

	std::string name{"field"};
	std::vector<double> data {1.,2.,3.};	
	
	// create the state vector
	StateVectorUni<double> sv(name, Entity_kind::CELL, data);
	
	// get the data
	std::vector<double>& out{sv.get_data()};

	for (int i=0; i<data.size(); i++) {
		ASSERT_EQ(out[i],data[i]);
	}  
  
}

TEST(StateUni, ModifyProtected) {

	std::string name{"field"};
	std::vector<double> data {1.,2.,3.};	
	
	// create the state vector
	StateVectorUni<double> sv(name, Entity_kind::CELL, data);
	
	// get the data by actual reference
	std::vector<double>& out{sv.get_data()};

	// test that the vector originally in the state vector is equal
	for (int i=0; i<data.size(); i++) {
		ASSERT_EQ(out[i],data[i]);
	}  
	
	// changes to out should change the state vector but not the original vector
	out[0]=-1.;
	ASSERT_EQ(out[0], sv.get_data()[0]);
	
	// check that we are isolated from the original data vector
	ASSERT_NE(sv.get_data()[0], data[0]);
  
}


