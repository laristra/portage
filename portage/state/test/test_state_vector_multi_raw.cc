/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/



#include "portage/support/portage.h"
#include "portage/state/state_vector_multi_raw.h"

#include "gtest/gtest.h"

using namespace Portage;



TEST(StateMulti, BasicInt) {

	std::string name{"field"};
	
	StateVectorMultiRaw<int> v(name);
  
  ASSERT_EQ(v.get_name(), name);
  ASSERT_EQ(v.get_type(), Portage::Field_type::MULTIMATERIAL_FIELD);
  ASSERT_EQ(v.data_type(), typeid(int));
  ASSERT_NE(v.data_type(), typeid(double));
  
}

TEST(StateMulti, BasicDouble1) {

	std::string name{"field"};
	
	StateVectorMultiRaw<> v(name);
  
  ASSERT_EQ(v.get_name(), name);
  ASSERT_EQ(v.get_type(), Portage::Field_type::MULTIMATERIAL_FIELD);
  ASSERT_EQ(v.data_type(), typeid(double));
  ASSERT_EQ(v.get_data(), nullptr);
  
}


TEST(StateMulti, DataAccess) {

	std::string name{"field"};
	std::vector<std::vector<double>> data {{1.,2.,3.},{4.,5.},{6.}};	
	double *temp[3];
	for (int i=0; i<3; i++)temp[i]=data[i].data();
	
	StateVectorMultiRaw<double> sv(name, temp);
	
	double const * const * const sv_data{sv.get_data()};

	for (int i=0; i<data.size(); i++) {
		for (int j=0; j<data[i].size(); j++){
			ASSERT_EQ(sv_data[i][j],data[i][j]);
		}
	}  
  
}


