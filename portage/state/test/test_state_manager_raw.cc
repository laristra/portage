/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <map>

#include "portage/support/portage.h"
#include "portage/state/state_vector_uni_raw.h"
#include "portage/state/state_vector_multi_raw.h"
#include "portage/state/state_manager.h"

#include "portage/simple_mesh/simple_mesh.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"

#include "gtest/gtest.h"

using namespace Portage;

// test the pointer strategy to reference an object in a state manager
// but this does not in fact use the state manager
TEST(StateManager, testPointerStateManager){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	StateVectorUniRaw<> sv1(name1, Entity_kind::CELL, data1.data());

	// check we can get the data out
	double *sv1_data = sv1.get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(sv1_data[i],data1[i]);
	}  

	// store the uni state vector in a map
	std::map<int, StateVectorBase*> state {{1,&sv1}};
	
	// get the data out of the map
	StateVectorUniRaw<double>* out{static_cast<StateVectorUniRaw<double>*>(state[1])};
	
	// get the name
	ASSERT_EQ("field1", out->get_name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}

}

// test the pointer strategy to reference an object in a state manager
// but this does not in fact use the state manager
TEST(StateManager, testSharedPointerStateManager){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	std::shared_ptr<StateVectorUniRaw<>> pv = std::make_shared<StateVectorUniRaw<>> (
		name1, Entity_kind::CELL, data1.data());

	// check we can get the data out
	double *pv_data = pv->get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(pv_data[i],data1[i]);
	} 
	 
	// store the uni state vector in a map
	std::map<int, std::shared_ptr<StateVectorBase>> state {{1,pv}};
	
	std::shared_ptr<StateVectorUniRaw<>> out=std::dynamic_pointer_cast<StateVectorUniRaw<>>(state[1]);

	// get the name
	ASSERT_EQ("field1", out->get_name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}

}

// test the pointer strategy to reference an object in a state manager
// but this does not in fact use the state manager
TEST(StateManager, testSharedPointerStateManagerMap){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	std::shared_ptr<StateVectorUniRaw<>> pv = std::make_shared<StateVectorUniRaw<>> (
		name1, Entity_kind::CELL, data1.data());

	// check we can get the data out
	double *pv_data = pv->get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(pv_data[i],data1[i]);
	} 

 	using pair_t = std::pair<Portage::Entity_kind,std::string>;
	 
	pair_t key{Entity_kind::CELL,name1};
	
	// store the uni state vector in a map
	std::map<pair_t, std::shared_ptr<StateVectorBase>> state {{key,pv}};
	
	std::shared_ptr<StateVectorUniRaw<>> out=std::dynamic_pointer_cast<StateVectorUniRaw<>>(state[key]);

	// get the name
	ASSERT_EQ("field1", out->get_name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}
	
}

// test the pointer strategy to reference an object in a state manager
// but this does not in fact use the state manager
TEST(StateManager, test1){

	// create the uni state vector
	std::string name1{"field1"};
	std::vector<double> data1 {1.,2.,3.};	
	StateVectorUniRaw<> sv1(name1, Entity_kind::CELL, data1.data());

	// check we can get the data out
	double *sv1_data = sv1.get_data();
	for (int i=0; i<data1.size(); i++) {
		ASSERT_EQ(sv1_data[i],data1[i]);
	}  

	// store the uni state vector in a map
	std::map<int, StateVectorBase*> state {{1,&sv1}};
	
	// get the data out of the map
	StateVectorUniRaw<double>* out{static_cast<StateVectorUniRaw<double>*>(state[1])};
	
	// get the name
	ASSERT_EQ("field1", out->get_name());
	
	// make sure the data comes out correctly
	for (int i=0; i<data1.size(); ++i){
		ASSERT_EQ(data1[i], out->get_data()[i]);
	}

}

// test the pointer strategy to reference an object in a state manager
// but this does not in fact use the state manager
TEST(StateManager, test2){

	std::string name1{"field1"};
	std::vector<std::vector<double>> data1 {{1.,2.,3.},{4.,5.},{6.}};	
	double *temp1[3];
	for (int i=0; i<3; i++)temp1[i]=data1[i].data();
	StateVectorMultiRaw<double> sv1(name1, temp1);

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
	StateVectorMultiRaw<int> sv2(name2, temp2);

	std::vector<StateVectorBase*> state {&sv1, &sv2};

	int **sv2_data = sv2.get_data();
	for (int i=0; i<data2.size(); i++) {
		for (int j=0; j<data2[i].size(); j++){
			ASSERT_EQ(sv2_data[i][j],data2[i][j]);
		}
	}  
	
	// now check through state and dynamic casting
	double **sv1_data_state = (dynamic_cast<StateVectorMultiRaw<double>*>(state[0]))->get_data();
	for (int i=0; i<data1.size(); i++) {
		for (int j=0; j<data1[i].size(); j++){
			ASSERT_EQ(sv1_data_state[i][j],data1[i][j]);
		}
	}  

	int **sv2_data_state = (dynamic_cast<StateVectorMultiRaw<int>*>(state[1]))->get_data();
	for (int i=0; i<data2.size(); i++) {
		for (int j=0; j<data2[i].size(); j++){
			ASSERT_EQ(sv2_data_state[i][j],data2[i][j]);
		}
	}  

	// now check through state and static
	double **sv1_data_state2 = (static_cast<StateVectorMultiRaw<double>*>(state[0]))->get_data();
	for (int i=0; i<data1.size(); i++) {
		for (int j=0; j<data1[i].size(); j++){
			ASSERT_EQ(sv1_data_state2[i][j],data1[i][j]);
		}
	}  

	int **sv2_data_state2 = (static_cast<StateVectorMultiRaw<int>*>(state[1]))->get_data();
	for (int i=0; i<data2.size(); i++) {
		for (int j=0; j<data2[i].size(); j++){
			ASSERT_EQ(sv2_data_state2[i][j],data2[i][j]);
		}
	}  

}

TEST(StateManager, testConstructManager){

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
	std::shared_ptr<StateVectorUniRaw<>> pv = std::make_shared<StateVectorUniRaw<>> (
		"field", Entity_kind::CELL, data.data());
	
	// make sure the data access is correct
	ASSERT_EQ(data[0], pv->get_data()[0]);
	
	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper};
	
	// add the data
	manager.add(pv);
		
	// make sure you can't clobber an existing field
	ASSERT_THROW(manager.add(pv), std::runtime_error);

	// make sure the data access is correct
	auto pv2=manager.get("field");
	
	// now do the pointer cast
	std::shared_ptr<StateVectorUniRaw<>> out = std::dynamic_pointer_cast<StateVectorUniRaw<>>(pv2);
	
	// test that the value is correct
	ASSERT_EQ(data[0], out->get_data()[0]);
	
	// now try with auto
	auto out2 = std::dynamic_pointer_cast<StateVectorUniRaw<>>(pv2);
	
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
	std::shared_ptr<StateVectorUniRaw<>> pv = std::make_shared<StateVectorUniRaw<>> (
		"field", Entity_kind::CELL, data.data());
	
	// make sure the data access is correct
	ASSERT_EQ(data[0], pv->get_data()[0]);
	
	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper};
	
	// add the data
	manager.add(pv);
	
	// make sure you can't clobber an existing field
	ASSERT_THROW(manager.add(pv), std::runtime_error);

	// make sure the data access is correct
	auto out=manager.get<StateVectorUniRaw<double>>("field");
	
	// test that the value is correct
	ASSERT_EQ(data[0], out->get_data()[0]);
	
	// try using the return reference
	auto out2 = manager.get<StateVectorUniRaw<double>>("field");
	
	// test that the value is correct
	ASSERT_EQ(data[0], out2->get_data()[0]);
	
	// try getting a non-existent field
	auto out3 = manager.get<StateVectorUniRaw<double>>("nonexistent");
	
	// test that the value is correct
	ASSERT_EQ(nullptr, out3);
	
	// add a second field 
	// create the data vector
	std::vector<double> data2 {10.};	

	// create a uni state vector
	std::shared_ptr<StateVectorUniRaw<>> pv2 = std::make_shared<StateVectorUniRaw<>> (
		"field2", Entity_kind::CELL, data2.data());
		
	// add the field to the manager
	manager.add(pv2);
	
	// check that the value is correct
	auto out4 = manager.get<StateVectorUniRaw<double>>("field2");
	
	// test that the value is correct
	ASSERT_EQ(data2[0], out4->get_data()[0]);
	
}


TEST(StateManager, multiMatField){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	std::string name{"field"};
	std::vector<std::vector<double>> data {{10.},{10.1},{10.2}};	
	double *temp[3];
	for (int i=0; i<3; i++)temp[i]=data[i].data();
	
	// create a uni state vector
	auto pv= std::make_shared<StateVectorMultiRaw<>> ("field", temp);

	// make sure the data access is correct
	for (int i=0; i<data.size(); i++) {
		for (int j=0; j<data[i].size(); j++){
			ASSERT_EQ(data[i][j], pv->get_data()[i][j]);
		}
	}  

	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper};

	// add the data
	manager.add(pv);
		
	// try using the return reference
	auto out = manager.get<StateVectorMultiRaw<double>>("field");
	
	// test that the value is correct
	for (int i=0; i<data.size(); i++) {
		for (int j=0; j<data[i].size(); j++){
			ASSERT_EQ(data[i][j], out->get_data()[i][j]);
		}
	}  
	

}


TEST(StateManager, mixedFields){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	std::vector<std::vector<double>> data {{10.},{10.1},{10.2}};	
	double *temp[3];
	for (int i=0; i<3; i++)temp[i]=data[i].data();
	
	// create a uni state vector
	auto pv= std::make_shared<StateVectorMultiRaw<>> ("pressure", temp);

	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper};

	// add the data
	manager.add(pv);
		
	// try using the return reference
	auto out = manager.get<StateVectorMultiRaw<double>>("pressure");
	
	// test that the values are correct
	for (int i=0; i<data.size(); i++) {
		for (int j=0; j<data[i].size(); j++){
			ASSERT_EQ(data[i][j], out->get_data()[i][j]);
		}
	}  
	
	// create the data vector
	std::vector<double> sdata {1.};	

	// create a uni state vector
	auto puv = std::make_shared<StateVectorUniRaw<>> ("temperature", Entity_kind::CELL, sdata.data());
		
	// add the scalar field
	manager.add(puv);
	
	// get the data
	auto sout = manager.get<StateVectorUniRaw<double>>("temperature");
	
	// test that the values are correct
	for (int i=0; i<sdata.size();i++) {
		ASSERT_EQ(sdata[i], sout->get_data()[i]);
	}

}


TEST(StateManager, mixedFields2){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the state manager
	StateManager<Simple_Mesh_Wrapper> manager{wrapper};
	
	
	// scalar mm data
	std::vector<std::vector<double>> data1 {{10.},{10.1},{10.2}};	
	double *pdata1[3];
	for (int i=0; i<3; i++)pdata1[i]=data1[i].data();

	// create a uni state vector
	auto pv1= std::make_shared<StateVectorMultiRaw<>> ("pressure", pdata1);

	// add the mm field
	manager.add(pv1);
	
	// try using the return reference
	auto out1 = manager.get<StateVectorMultiRaw<double>>("pressure");
	
	// test that the values are correct
	for (int i=0; i<data1.size(); i++) {
		for (int j=0; j<data1[i].size(); j++){
			ASSERT_EQ(data1[i][j], out1->get_data()[i][j]);
		}
	} 
	
	 
	// vector mm data
	std::vector<std::vector<Point<2>>> data2 {
		{Point<2>{10.,-10.}},{Point<2>{10.,-10.}},{Point<2>{10.,-10.}}};	
	Point<2> *pdata2[3];
	for (int i=0; i<3; i++)pdata2[i]=data2[i].data();

	// create a uni state vector
	auto pv2= std::make_shared<StateVectorMultiRaw<Point<2>>> ("mm centroid", pdata2);

	// add the mm field
	manager.add(pv2);
	
	// try using the return reference
	auto out2 = manager.get<StateVectorMultiRaw<Point<2>>>("mm centroid");
	
	// test that the values are correct
	for (int i=0; i<data2.size(); i++) {
		for (int j=0; j<data2[i].size(); j++){
			ASSERT_EQ(data2[i][j][0], out2->get_data()[i][j][0]); //x coord
			ASSERT_EQ(data2[i][j][1], out2->get_data()[i][j][1]); //y coord
		}
	} 
	
	 
	
	// scalar uni data
	std::vector<double> data3 {1.};	

	// create a uni state vector
	auto pv3 = std::make_shared<StateVectorUniRaw<>> ("temperature", 
		Entity_kind::CELL, data3.data());
		
	// add the scalar field
	manager.add(pv3);
	
	// get the data
	auto out3 = manager.get<StateVectorUniRaw<double>>("temperature");
	
	// test that the values are correct
	for (int i=0; i<data3.size();i++) {
		ASSERT_EQ(data3[i], out3->get_data()[i]);
	}


	// vector uni data
	std::vector<Point<2>> data4 {Point<2>{1.,1.}};	

	// create a uni state vector
	auto pv4 = std::make_shared<StateVectorUniRaw<Point<2>>> ("cell centroid", 
		Entity_kind::CELL, data4.data());
		
	// add the scalar field
	manager.add(pv4);
	
	// get the data
	auto out4 = manager.get<StateVectorUniRaw<Point<2>>>("cell centroid");
	
	// test that the values are correct
	for (int i=0; i<data4.size();i++) {
		ASSERT_EQ(data4[i][0], out4->get_data()[i][0]); //x coord
		ASSERT_EQ(data4[i][1], out4->get_data()[i][1]); //y coord
	}

}
























