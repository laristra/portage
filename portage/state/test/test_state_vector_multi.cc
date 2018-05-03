/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/



#include "portage/support/portage.h"
#include "portage/state/state_vector_multi.h"

#include "gtest/gtest.h"

using namespace Portage;



TEST(StateMulti, BasicInt) {

	std::string name{"field"};
	
	StateVectorMulti<int> v(name);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), Portage::Field_type::MULTIMATERIAL_FIELD);
  ASSERT_EQ(v.data_type(), typeid(int));
  ASSERT_NE(v.data_type(), typeid(double));
  
}

TEST(StateMulti, BasicDouble1) {

	std::string name{"field"};
	
	StateVectorMulti<> v(name);
  
  ASSERT_EQ(v.name(), name);
  ASSERT_EQ(v.type(), Portage::Field_type::MULTIMATERIAL_FIELD);
  ASSERT_EQ(v.data_type(), typeid(double));
  
}


TEST(StateMulti, DataAccess) {

	std::string name{"field"};
	std::unordered_map<int, std::vector<double>> data {{1,{1.,2.,3.}},{2,{4.,5.}},{3,{6.}}};	
	
	StateVectorMulti<double> sv(name, data);
	
	std::unordered_map<int, std::vector<double>>& out{sv.get_data()};

	for (int i=0; i<data.size(); i++) {
		for (int j=0; j<data[i].size(); j++){
			ASSERT_EQ(out[i][j],data[i][j]);
		}
	}  
  
}

TEST(StateMulti, ModifyProtected) {

	std::string name{"field"};
	std::unordered_map<int, std::vector<double>> data {{1,{1.,2.,3.}},{2,{4.,5.}},{3,{6.}}};	
	
	StateVectorMulti<double> sv(name, data);
	
	std::unordered_map<int, std::vector<double>>& out{sv.get_data()};

	for (int i=0; i<data.size(); i++) {
		for (int j=0; j<data[i].size(); j++){
			ASSERT_EQ(out[i][j],data[i][j]);
		}
	}  

	// changes to out should change the state vector but not the original vector
	out[1][1]=-1.;
	ASSERT_EQ(out[1][1], sv.get_data()[1][1]);
	
	// check that we are isolated from the original data vector
	ASSERT_NE(sv.get_data()[1][1], data[1][1]);
  
}


