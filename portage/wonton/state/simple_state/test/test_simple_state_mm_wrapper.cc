/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <map>
#include <unordered_map>
#include <memory>

#include "portage/support/portage.h"
#include "portage/state/state_vector_uni.h"
#include "portage/state/state_vector_multi.h"
#include "portage/state/state_manager.h"

#include "portage/simple_mesh/simple_mesh.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/wonton/state/simple_state/simple_state_mm_wrapper.h"

#include "gtest/gtest.h"

using namespace Portage;




///////////////////////////////////////
// Begin the state manager tests proper
///////////////////////////////////////

TEST(Simple_State_Wrapper, manageMeshField2){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the data vector
	std::vector<double> data {1.};	

	// create a uni state vector
	std::shared_ptr<StateVectorUni<>> pv = std::make_shared<StateVectorUni<>> (
		"field", Entity_kind::CELL, data);
	
	// make sure the data access is correct
	ASSERT_EQ(data[0], pv->get_data()[0]);
	
	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};
	manager.add(pv);
	
	// make sure you can't clobber an existing field
	ASSERT_THROW(manager.add(pv), std::runtime_error);

	// make sure the data access is correct
	std::shared_ptr<StateVectorBase> pv2=manager.get("field");
	
	// now do the pointer cast
	std::shared_ptr<StateVectorUni<>> out = std::dynamic_pointer_cast<StateVectorUni<>>(pv2);
	
	// test that the value is correct
	ASSERT_EQ(data[0], out->get_data()[0]);
	
	// now try with auto
	auto out2 = std::dynamic_pointer_cast<StateVectorUni<>>(pv2);
	
	// test that the value is correct
	ASSERT_EQ(data[0], out2->get_data()[0]);
	
	
}



TEST(Simple_State_Wrapper, manageMeshFieldTemplate){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the data vector
	std::vector<double> data {1.};	

	// create a uni state vector
	std::shared_ptr<StateVectorUni<>> pv = std::make_shared<StateVectorUni<>> (
		"field", Entity_kind::CELL, data);
	
	// make sure the data access is correct
	ASSERT_EQ(data[0], pv->get_data()[0]);
	
	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};
	manager.add(pv);
	
	// make sure you can't clobber an existing field
	ASSERT_THROW(manager.add(pv), std::runtime_error);

	// make sure the data access is correct
	auto out=manager.get<StateVectorUni<double>>("field");
	
	// test that the value is correct
	ASSERT_EQ(data[0], out->get_data()[0]);
	
	// try using the return reference
	auto out2 = manager.get<StateVectorUni<double>>("field");
	
	// test that the value is correct
	ASSERT_EQ(data[0], out2->get_data()[0]);
	
	// try getting a non-existent field
	auto out3 = manager.get<StateVectorUni<double>>("nonexistent");
	
	// test that the value is correct
	ASSERT_EQ(nullptr, out3);
	
	// add a second field 
	// create the data vector
	std::vector<double> data2 {10.};	

	// create a uni state vector
	std::shared_ptr<StateVectorUni<>> pv2 = std::make_shared<StateVectorUni<>> (
		"field2", Entity_kind::CELL, data2);
		
	// add the field to the manager
	manager.add(pv2);
	
	// check that the value is correct
	auto out4 = manager.get<StateVectorUni<double>>("field2");
	
	// test that the value is correct
	ASSERT_EQ(data2[0], out4->get_data()[0]);
	
}


TEST(Simple_State_Wrapper, multiMatField){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	std::string name{"field"};
	std::unordered_map<int, std::vector<double>> data {{1,{10.}},{2,{10.1}},{5,{10.2}}};	
	
	auto pv= std::make_shared<StateVectorMulti<>> ("field", data);

	// make sure the data access is correct
	// now check data2 through state and dynamic casting
	for (auto& kv : pv->get_data()) {
		auto& out=kv.second;
		for (int j=0; j<out.size(); j++){
			ASSERT_EQ(data[kv.first][j], out[j]);
		}
	} 


	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};
	manager.add(pv);
	
	// try using the return reference
	auto out = manager.get<StateVectorMulti<double>>("field");

	// test that the value obtained through the state manager is correct
	for (auto& kv : out->get_data()) {
		auto& out=kv.second;
		for (int j=0; j<out.size(); j++){
			ASSERT_EQ(data[kv.first][j], out[j]);
		}
	} 
	
}

// this test mixes StateVectorMulti<double> and StateVectorUn<double> fields
TEST(Simple_State_Wrapper, mixedFields){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	std::unordered_map<int, std::vector<double>> data {{1,{10.}},{2,{10.1}},{5,{10.2}}};	
	
	// create a uni state vector
	auto spv= std::make_shared<StateVectorMulti<>> ("pressure", data);

	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};
	manager.add(spv);

	// try using the return reference
	auto out = manager.get<StateVectorMulti<double>>("pressure");
	
	// test that the value obtained through the state manager is correct
	for (auto& kv : out->get_data()) {
		auto& out=kv.second;
		for (int j=0; j<out.size(); j++){
			ASSERT_EQ(data[kv.first][j], out[j]);
		}
	} 
	
	// create the data vector
	std::vector<double> sdata {1.};	

	// create a uni state vector
	auto spsd = std::make_shared<StateVectorUni<>> ("temperature", Entity_kind::CELL, sdata);
		
	// add the scalar field
	manager.add(spsd);
	
	// get the data
	auto sout = manager.get<StateVectorUni<double>>("temperature");
	
	// test that the values are correct
	for (int i=0; i<sdata.size();i++) {
		ASSERT_EQ(sdata[i], sout->get_data()[i]);
	}

}


TEST(Simple_State_Wrapper, mixedFields2){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};
	
	
	// integer mm data
	std::unordered_map<int, std::vector<int>> data0 {{1,{10}},{2,{11}},{5,{12}}};	

	// create a uni state vector
	auto sp0= std::make_shared<StateVectorMulti<int>> ("indices", data0);

	// add the mm field
	manager.add(sp0);
	
	// try using the return reference
	auto out0 = manager.get<StateVectorMulti<int>>("indices");
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out0->get_data()) {
		auto& out=kv.second;
		for (int j=0; j<out.size(); j++){
			ASSERT_EQ(data0[kv.first][j], out[j]);
		}
	} 
	
	 
	// scalar mm data
	std::unordered_map<int, std::vector<double>> data1 {{1,{10.}},{2,{10.1}},{5,{10.2}}};	

	// create a uni state vector
	auto sp1= std::make_shared<StateVectorMulti<>> ("density", data1);

	// add the mm field
	manager.add(sp1);
	
	// try using the return reference
	auto out1 = manager.get<StateVectorMulti<double>>("density");
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out1->get_data()) {
		auto& out=kv.second;
		for (int j=0; j<out.size(); j++){
			ASSERT_EQ(data1[kv.first][j], out[j]);
		}
	} 
	
	 
	// vector mm data
	std::unordered_map<int, std::vector<Point<2>>> data2 {
		{1,{Point<2>{10.,-10.}}},{2,{Point<2>{11.,-11.}}},{5,{Point<2>{12.,-12.}}}};	

	// create a uni state vector
	auto sp2= std::make_shared<StateVectorMulti<Point<2>>> ("mm centroid", data2);

	// add the mm field
	manager.add(sp2);
	
	// try using the return reference
	auto out2 = manager.get<StateVectorMulti<Point<2>>>("mm centroid");
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out2->get_data()) {
		auto& out=kv.second;
		for (int j=0; j<out.size(); j++){
			ASSERT_EQ(data2[kv.first][j][0], out2->get_data()[kv.first][j][0]); //x coord
			ASSERT_EQ(data2[kv.first][j][1], out2->get_data()[kv.first][j][1]); //y coord
		}
	} 
	
	
	// scalar uni data
	std::vector<double> data3 {1.};	

	// create a uni state vector
	auto sp3 = std::make_shared<StateVectorUni<>> ("temperature", 
		Entity_kind::CELL, data3);
		
	// add the scalar field
	manager.add(sp3);
	
	// get the data
	auto out3 = manager.get<StateVectorUni<double>>("temperature");
	
	// test that the values are correct
	for (int i=0; i<data3.size();i++) {
		ASSERT_EQ(data3[i], out3->get_data()[i]);
	}

	 

	// vector uni data
	std::vector<Point<2>> data4 {Point<2>{1.,1.}};	

	// create a uni state vector
	auto sp4 = std::make_shared<StateVectorUni<Point<2>>> ("cell centroid", 
		Entity_kind::CELL, data4);
		
	// add the scalar field
	manager.add(sp4);
	
	// get the data
	auto out4 = manager.get<StateVectorUni<Point<2>>>("cell centroid");
	
	// test that the values are correct
	for (int i=0; i<data4.size();i++) {
		ASSERT_EQ(data4[i][0], out4->get_data()[i][0]); //x coord
		ASSERT_EQ(data4[i][1], out4->get_data()[i][1]); //y coord
	}

}

// Test the data isolation between the original data and the manager data
TEST(Simple_State_Wrapper, testUniIsolation){
	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};
	
	// scalar uni data
	std::vector<double> data {10.};	

	// add the scalar field
	manager.add(std::make_shared<StateVectorUni<>> ("temperature", 
		Entity_kind::CELL, data));
	
	// get the data
	auto& out = manager.get<StateVectorUni<double>>("temperature")->get_data();
	
	// test that the values are correct
	for (int i=0; i<data.size();i++) {
		ASSERT_EQ(data[i], out[i]);
	}
	
	// add data to the back of the managed data
	out.push_back(-50.);
	
	// make sure the sizes are correct and data is unaltered
	ASSERT_EQ(out.size(), 2);
	ASSERT_EQ(data.size(), 1);
	ASSERT_EQ(out[1], -50.);
}


// Test the data isolation between the original data and the manager data
TEST(Simple_State_Wrapper, testMultiIsolation){
	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};
	
	// scalar multi data
	std::unordered_map<int, std::vector<double>> data {{1,{10.}},{2,{10.1}},{5,{10.2}}};	

	// add the scalar field
	manager.add(std::make_shared<StateVectorMulti<>> ("pressure", data));
	
	// get the data
	auto& out = manager.get<StateVectorMulti<double>>("pressure")->get_data();
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out) {
		auto& d=kv.second;
		for (int j=0; j<d.size(); j++){
			ASSERT_EQ(data[kv.first][j], d[j]);
		}
	} 

	// make sure there are the correct number of materials in each map
	ASSERT_EQ(data.size(), out.size());
	
	// add data to the back of the managed data
	out[7]= std::vector<double>{11.};

	// make sure there are the correct number of materials in each map
	ASSERT_NE(data.size(), out.size()); // we have a shallow copy of the map vectors
	ASSERT_EQ(out.size(), 4);
	ASSERT_EQ(data.size(), 3);
	
	// now check for a deep copy of the vector data
	
	// unmodified
	ASSERT_EQ(data[1].size(), out[1].size());
	ASSERT_EQ(data[1].size(), 1);
	
	// add data to a single material
	out[1].push_back(-50.);
	
	// check sizes
	ASSERT_NE(data[1].size(), out[1].size());
	ASSERT_EQ(data[1].size(), 1);
	ASSERT_EQ(out[1].size(), 2);
}


// Test the data isolation between the original data and the manager data at 
// depth three (inside a point)
TEST(Simple_State_Wrapper, testMultiIsolationPoint){
	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 1., 1., 1, 1};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};

	// vector mm data
	std::unordered_map<int, std::vector<Point<2>>> data {
		{1,{Point<2>{10.,-10.}}},{2,{Point<2>{11.,-11.}}},{5,{Point<2>{12.,-12.}}}};	
	
	// add the scalar field
	manager.add(std::make_shared<StateVectorMulti<Point<2>>> ("centroid", data));
	
	// get the data
	auto& out = manager.get<StateVectorMulti<Point<2>>>("centroid")->get_data();
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out) {
		auto& d=kv.second;
		for (int j=0; j<d.size(); j++){
			ASSERT_EQ(data[kv.first][j][0], d[j][0]);
			ASSERT_EQ(data[kv.first][j][1], d[j][1]);
		}
	} 

	// make sure there are the correct number of materials in each map
	ASSERT_EQ(data.size(), out.size());
	
	// add data to the back of the managed data
	out[7]= std::vector<Point<2>>{Point<2>{15.,-15.}};

	// make sure there are the correct number of materials in each map
	ASSERT_NE(data.size(), out.size()); // we have a shallow copy of the map vectors
	ASSERT_EQ(out.size(), 4);
	ASSERT_EQ(data.size(), 3);

	
	// now check for a deep copy of the vector data
	
	// unmodified
	ASSERT_EQ(data[1].size(), out[1].size());
	ASSERT_EQ(data[1].size(), 1);
	
	// add data to a single material
	out[1].push_back(Point<2>{20.,-20.});
	
	// check sizes
	ASSERT_NE(data[1].size(), out[1].size());
	ASSERT_EQ(data[1].size(), 1);
	ASSERT_EQ(out[1].size(), 2);
	
	// finally check that the points themselves are copied
	
	// change coordinates
	out[5][0][0]=75.;
	out[5][0][1]=-75.;
	
	// make sure the Point components are different
	ASSERT_NE(out[5][0][0], data[5][0][0]);
	ASSERT_NE(out[5][0][1], data[5][0][1]);
}

// Test with a 4 cell mesh
TEST(Simple_State_Wrapper,test4Cell){

	using namespace Wonton;

	// create the mesh
	Simple_Mesh mesh{0., 0., 2., 2., 2, 2};
	
	// create the wrapper
	Simple_Mesh_Wrapper wrapper{mesh};
	
	// create the state manager
	Simple_State_Wrapper<Simple_Mesh_Wrapper> manager{wrapper};

	// map of names of materials
	std::unordered_map<std::string,int> matnames{{"mat0",1},{"mat2",3},{"mat3",5}};
	
	// add the names to the manager
	manager.add_material_names(matnames);
	
	// check get_nmats API function
	ASSERT_EQ(manager.num_materials(),3);
	for (const auto& kv: matnames){
		ASSERT_EQ(manager.get_material_id(kv.first),matnames[kv.first]);
		ASSERT_EQ(manager.material_name(kv.second),kv.first);
	}
	
	// create the material cells
	std::unordered_map<int,std::vector<int>> cells{{1,{0,1,2,3}},{3,{0,1}},{5,{0,2}}};
	
	// add the material cells
	manager.add_material_cells(cells);
	
	// create the material cells
	std::unordered_map<int,std::vector<int>> cells2{{10,{0,1,2,3}}};
	
	// add the cell ids for an unknown material
	ASSERT_THROW(manager.add_material_cells(cells2), std::runtime_error);
	
	// make sure the shape is correct
	for (const auto& kv : manager.get_material_shape()){
		ASSERT_EQ(cells[kv.first].size(), kv.second);
	}
	
	// create the multi material data
	std::unordered_map<int,std::vector<double>> density 
		{{1,{1.,1.1,1.2,1.3}},{3,{10.,10.1}},{5,{100.,100.1}}};
	
	// add the data, note the make_shared actually copies the data
	manager.add(std::make_shared<StateVectorMulti<>>("density",density));
	
	// try to add incorrectly sized data, missing a key
	std::unordered_map<int,std::vector<double>> density_bad_shape1
		{{1,{1.,1.1,1.2,1.3}},{3,{10.,10.1}}};		
	ASSERT_THROW(
		manager.add(std::make_shared<StateVectorMulti<>>("density_bad_shape1",density_bad_shape1)),
		std::runtime_error
	);
	
	// try to add incorrectly sized data, extra cell 
	std::unordered_map<int,std::vector<double>> density_bad_shape2
		{{1,{1.,1.1,1.2,1.3, -1.}},{3,{10.,10.1}},{5,{100.,100.1}}};
	ASSERT_THROW(
		manager.add(std::make_shared<StateVectorMulti<>>("density_bad_shape2",density_bad_shape2)),
		std::runtime_error
	);
	
	// try to add incorrectly sized data, extra material 
	std::unordered_map<int,std::vector<double>> density_bad_shape3
		{{1,{1.,1.1,1.2,1.3}},{3,{10.,10.1}},{5,{100.,100.1}},{7,{100.,100.1}}};
	ASSERT_THROW(
		manager.add(std::make_shared<StateVectorMulti<>>("density_bad_shape3",density_bad_shape3)),
		std::runtime_error
	);
	
	// check that the number of cells for each material is correct
	for (const auto& kv: cells){
			ASSERT_EQ(kv.second.size(), manager.num_material_cells(kv.first));
	}
	
	// check that the cells ids for each material are correct
	for (const auto& kv: cells){
		for (int i=0; i<kv.second.size(); i++){	
			ASSERT_EQ(kv.second[i],manager.get_material_cells(kv.first)[i]);
		}
	}

	// test that the number of materials in each cell is correct
	ASSERT_EQ(manager.cell_get_num_mats(0),3);
	ASSERT_EQ(manager.cell_get_num_mats(1),2);
	ASSERT_EQ(manager.cell_get_num_mats(2),2);
	ASSERT_EQ(manager.cell_get_num_mats(3),1);
	
	// check that the materials in the cell are correct
	std::unordered_set<int> dum{1,3,5};
	ASSERT_EQ(manager.get_cell_materials(0), dum);
	dum=std::unordered_set<int>{1,3};
	ASSERT_EQ(manager.get_cell_materials(1), dum);
	dum=std::unordered_set<int>{1,5};
	ASSERT_EQ(manager.get_cell_materials(2), dum);
	dum=std::unordered_set<int>{1};
	ASSERT_EQ(manager.get_cell_materials(3), dum);
	
	// check that the material cell indices work
	ASSERT_EQ(manager.cell_index_in_material(0,1),0);
	ASSERT_EQ(manager.cell_index_in_material(1,1),1);
	ASSERT_EQ(manager.cell_index_in_material(2,1),2);
	ASSERT_EQ(manager.cell_index_in_material(3,1),3);	
	ASSERT_EQ(manager.cell_index_in_material(0,3),0);
	ASSERT_EQ(manager.cell_index_in_material(1,3),1);
	ASSERT_EQ(manager.cell_index_in_material(0,5),0);
	ASSERT_EQ(manager.cell_index_in_material(2,5),1);

	// get the data
	auto out = manager.get<StateVectorMulti<double>>("density");
	
	ASSERT_EQ(out->get_type(),Portage::Field_type::MULTIMATERIAL_FIELD);

	auto& out_values= out->get_data();
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out_values) {
		auto& d=kv.second;
		for (int j=0; j<d.size(); j++){
			ASSERT_EQ(density[kv.first][j], d[j]);
		}
	} 
	
	// now add centroid data to the mesh
	// vector mm data
	std::unordered_map<int, std::vector<Point<2>>> centroid {
		{1,{Point<2>{.5,.5},Point<2>{1.5,.5},Point<2>{.5,1.5},Point<2>{1.5,1.5}}},
		{3,{Point<2>{.55,.55},Point<2>{1.55,.55}}},
		{5,{Point<2>{.45,.45},Point<2>{.45,1.45}}}};	
	
	// add the scalar field
	manager.add(std::make_shared<StateVectorMulti<Point<2>>> ("centroid", centroid));
	
	// Try different ways of accessing the data, some using template deduction
	
	// get the data
	auto& out2 = manager.get<StateVectorMulti<Point<2>>>("centroid")->get_data();
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out2) {
		auto& d=kv.second;
		for (int j=0; j<d.size(); j++){
			ASSERT_EQ(centroid[kv.first][j][0], d[j][0]);
			ASSERT_EQ(centroid[kv.first][j][1], d[j][1]);
		}
	} 

	// get the data using the simplified multi material get
	auto& out3 = manager.get<StateVectorMulti<Point<2>>>("centroid")->get_data();
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out3) {
		auto& d=kv.second;
		for (int j=0; j<d.size(); j++){
			ASSERT_EQ(centroid[kv.first][j][0], d[j][0]);
			ASSERT_EQ(centroid[kv.first][j][1], d[j][1]);
		}
	} 

	// get the data using the template parameter deduction form
	auto sv=manager.get<StateVectorMulti<Point<2>>>("centroid");
	auto& out4 =sv->get_data();
	
	// test that the values obtained through the state manager are correct
	for (auto& kv : out4) {
		auto& d=kv.second;
		for (int j=0; j<d.size(); j++){
			ASSERT_EQ(centroid[kv.first][j][0], d[j][0]);
			ASSERT_EQ(centroid[kv.first][j][1], d[j][1]);
		}
	} 

	// what is in the state manager?
	//for (auto& kv: manager.get_state_keys()){
	//	std::cout<<"type: "<<kv.first<<" name: " << kv.second<<std::endl;
	//}
}





















