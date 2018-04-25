/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <map>

#include "portage/support/portage.h"
#include "portage/state/state_vector_uni.h"
#include "portage/state/state_vector_multi.h"
#include "portage/state/state_manager.h"

#include "portage/simple_mesh/simple_mesh.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"

#include "gtest/gtest.h"

using namespace Portage;

TEST(StateManager, testPointerStateManager){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	StateVectorUni<> sv1(name1, data1.data());

	// check we can get the data out
	double *sv1_data = sv1.get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(sv1_data[i],data1[i]);
	}  

	// store the uni state vector in a map
	std::map<int, StateVectorBase*> state {{1,&sv1}};
	
	// get the data out of the map
	StateVectorUni<double>* out{static_cast<StateVectorUni<double>*>(state[1])};
	
	// get the name
	ASSERT_EQ("field1", out->name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}

}

TEST(StateManager, testSharedPointerStateManager){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	std::shared_ptr<StateVectorUni<>> pv = std::make_shared<StateVectorUni<>> (
		name1, data1.data());

	// check we can get the data out
	double *pv_data = pv->get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(pv_data[i],data1[i]);
	} 
	 
	// store the uni state vector in a map
	std::map<int, std::shared_ptr<StateVectorBase>> state {{1,pv}};
	
	std::shared_ptr<StateVectorUni<>> out=std::dynamic_pointer_cast<StateVectorUni<>>(state[1]);

	// get the name
	ASSERT_EQ("field1", out->name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}

}

TEST(StateManager, testSharedPointerStateManagerMap){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	std::shared_ptr<StateVectorUni<>> pv = std::make_shared<StateVectorUni<>> (
		name1, data1.data());

	// check we can get the data out
	double *pv_data = pv->get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(pv_data[i],data1[i]);
	} 

 	using pair_t = std::pair<Portage::Entity_kind,std::string>;
	 
	pair_t key{Entity_kind::CELL,name1};
	
	// store the uni state vector in a map
	std::map<pair_t, std::shared_ptr<StateVectorBase>> state {{key,pv}};
	
	std::shared_ptr<StateVectorUni<>> out=std::dynamic_pointer_cast<StateVectorUni<>>(state[key]);

	// get the name
	ASSERT_EQ("field1", out->name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}
	
}

TEST(StateManager, test1){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	StateVectorUni<> sv1(name1, data1.data());

	// check we can get the data out
	double *sv1_data = sv1.get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(sv1_data[i],data1[i]);
	}  

	// store the uni state vector in a map
	std::map<int, StateVectorBase*> state {{1,&sv1}};
	
	// get the data out of the map
	StateVectorUni<double>* out{static_cast<StateVectorUni<double>*>(state[1])};
	
	// get the name
	ASSERT_EQ("field1", out->name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}

}

TEST(StateManager, test2){

	std::string name1{"field1"};
	std::vector<std::vector<double>> data1 {{1.,2.,3.},{4.,5.},{6.}};	
	double *temp1[3];
	for (int i=0; i<3; i++)temp1[i]=data1[i].data();
	StateVectorMulti<double> sv1(name1, temp1);

	double **sv1_data = sv1.get_data();
	for (int i=0; i<data1.size(); i++) {
		for (int j=0; j<data1[i].size(); j++){
			ASSERT_EQ(sv1_data[i][j],data1[i][j]);
		}
	}  

	std::string name2{"field2"};
	std::vector<std::vector<int>> data2 {{-1,-2,-3},{-4,-5}};	
	int *temp2[2];
	for (int i=0; i<2; i++)temp2[i]=data2[i].data();
	StateVectorMulti<int> sv2(name2, temp2);

	std::vector<StateVectorBase*> state {&sv1, &sv2};

	int **sv2_data = sv2.get_data();
	for (int i=0; i<data2.size(); i++) {
		for (int j=0; j<data2[i].size(); j++){
			ASSERT_EQ(sv2_data[i][j],data2[i][j]);
		}
	}  
	
	// now check through state and dynamic casting
	double **sv1_data_state = (dynamic_cast<StateVectorMulti<double>*>(state[0]))->get_data();
	for (int i=0; i<data1.size(); i++) {
		for (int j=0; j<data1[i].size(); j++){
			ASSERT_EQ(sv1_data_state[i][j],data1[i][j]);
		}
	}  

	int **sv2_data_state = (dynamic_cast<StateVectorMulti<int>*>(state[1]))->get_data();
	for (int i=0; i<data2.size(); i++) {
		for (int j=0; j<data2[i].size(); j++){
			ASSERT_EQ(sv2_data_state[i][j],data2[i][j]);
		}
	}  

	// now check through state and static
	double **sv1_data_state2 = (static_cast<StateVectorMulti<double>*>(state[0]))->get_data();
	for (int i=0; i<data1.size(); i++) {
		for (int j=0; j<data1[i].size(); j++){
			ASSERT_EQ(sv1_data_state2[i][j],data1[i][j]);
		}
	}  

	int **sv2_data_state2 = (static_cast<StateVectorMulti<int>*>(state[1]))->get_data();
	for (int i=0; i<data2.size(); i++) {
		for (int j=0; j<data2[i].size(); j++){
			ASSERT_EQ(sv2_data_state2[i][j],data2[i][j]);
		}
	}  

}

TEST(StateManager, manageMeshField1){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper};
	

}

TEST(StateManager, manageMeshField2){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the data vector
	std::vector<double> data {1.};	

	// create a uni state vector
	std::shared_ptr<StateVectorUni<>> pv = std::make_shared<StateVectorUni<>> (
		"field", data.data());
	
	// make sure the data access is correct
	ASSERT_EQ(data[0], pv->get_data()[0]);
	
	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper, Entity_kind::CELL, pv};
	
	// make sure you can't clobber an existing field
	ASSERT_THROW(manager.add(Entity_kind::CELL, pv), std::runtime_error);

	// make sure the data access is correct
	std::shared_ptr<StateVectorBase> pv2;
	manager.get(Entity_kind::CELL, "field", pv2);
	
	// now do the pointer cast
	std::shared_ptr<StateVectorUni<>> out = std::dynamic_pointer_cast<StateVectorUni<>>(pv2);
	
	// test that the value is correct
	ASSERT_EQ(data[0], out->get_data()[0]);
	
	// now try with auto
	auto out2 = std::dynamic_pointer_cast<StateVectorUni<>>(pv2);
	
	// test that the value is correct
	ASSERT_EQ(data[0], out2->get_data()[0]);
	
	
}



TEST(StateManager, manageMeshFieldTemplate){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the data vector
	std::vector<double> data {1.};	

	// create a uni state vector
	std::shared_ptr<StateVectorUni<>> pv = std::make_shared<StateVectorUni<>> (
		"field", data.data());
	
	// make sure the data access is correct
	ASSERT_EQ(data[0], pv->get_data()[0]);
	
	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper, Entity_kind::CELL, pv};
	
	// make sure you can't clobber an existing field
	ASSERT_THROW(manager.add(Entity_kind::CELL, pv), std::runtime_error);

	// make sure the data access is correct
	std::shared_ptr<StateVectorUni<>> out;
	
	// get the data
	manager.get<double, StateVectorUni>(Entity_kind::CELL, "field", out);
	
	// test that the value is correct
	ASSERT_EQ(data[0], out->get_data()[0]);
	
	// try using the return reference
	auto out2 = manager.get<double, StateVectorUni>(Entity_kind::CELL, "field");
	
	// test that the value is correct
	ASSERT_EQ(data[0], out2->get_data()[0]);
	
	// try getting a non-existent field
	auto out3 = manager.get<double, StateVectorUni>(Entity_kind::CELL, "nonexistent");
	
	// test that the value is correct
	ASSERT_EQ(nullptr, out3);
	
	// add a second field 
	// create the data vector
	std::vector<double> data2 {10.};	

	// create a uni state vector
	std::shared_ptr<StateVectorUni<>> pv2 = std::make_shared<StateVectorUni<>> (
		"field2", data2.data());
		
	// add the field to the manager
	manager.add(Entity_kind::CELL, pv2);
	
	// check that the value is correct
	auto out4 = manager.get<double, StateVectorUni>(Entity_kind::CELL, "field2");
	
	// test that the value is correct
	ASSERT_EQ(data2[0], out4->get_data()[0]);
	
}






















