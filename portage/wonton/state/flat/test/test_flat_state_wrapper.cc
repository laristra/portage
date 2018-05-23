/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#include "portage/wonton/state/jali/jali_state_wrapper.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"

#include "portage/support/portage.h"

#include <iostream>
#include <exception>
#include <string>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wonton/state/flat/flat_state_wrapper.h"

TEST(Flat_State_Wrapper, VectorInit) {
  std::vector<double> vertx={0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> verty={0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
  std::vector<double>
        &nfield1=*new std::vector<double>(16,0.),
	&nfield2=*new std::vector<double>(16,0.),
	&cfield1=*new std::vector<double>(9,0.),
	&cfield2=*new std::vector<double>(9,0.);
  std::vector<std::string> names = {"nf1", "cf1", "cf2", "nf2"};
  for (size_t i=0; i<16; i++) {
      nfield1[i] = vertx[i]*vertx[i] + verty[i]*verty[i];
      nfield2[i] = vertx[i]+verty[i];
  }
  for (size_t i=0; i<3; i++) for (size_t j=0; j<3; j++) {
      double xc = .25*(vertx[i+j*4] + vertx[i+1+j*4] + vertx[i+j*4] + vertx[i+1+j*4]);
      double yc = .5*(verty[i+(j+1)*4] + verty[i+j*4]);
      cfield1[i+j*3] = xc*xc + yc*yc;
      cfield2[i+j*3] = xc + yc;
  }

  std::shared_ptr<std::vector<double>> nf1(&nfield1),nf2(&nfield2),cf1(&cfield1),cf2(&cfield2);
  std::vector<std::shared_ptr<std::vector<double>>> data ={nf1,cf1,cf2,nf2};

  using namespace Portage;
  using namespace Wonton;

  std::vector<Entity_kind> entities = {NODE, CELL, CELL, NODE};

  Flat_State_Wrapper<double> flatwrap;
  flatwrap.initialize(names, entities, data);

  // check clearing function for reinitialize
  flatwrap.initialize(names, entities, data);
}

TEST(Flat_State_Wrapper, VectorInitFail1) {
  std::vector<double> vertx={0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> verty={0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
  std::vector<double>
        &nfield1=*new std::vector<double>(16,0.),
	&nfield2=*new std::vector<double>(16,0.),
	&cfield1=*new std::vector<double>(9,0.),
	&cfield2=*new std::vector<double>(9,0.),
	&badfield=*new std::vector<double>(23,0.);
  std::vector<std::string> names = {"nf1", "cf1", "cf2", "nf2", "bad"};
  for (size_t i=0; i<16; i++) {
      nfield1[i] = vertx[i]*vertx[i] + verty[i]*verty[i];
      nfield2[i] = vertx[i]+verty[i];
  }
  for (size_t i=0; i<3; i++) for (size_t j=0; j<3; j++) {
      double xc = .25*(vertx[i+j*4] + vertx[i+1+j*4] + vertx[i+j*4] + vertx[i+1+j*4]);
      double yc = .5*(verty[i+(j+1)*4] + verty[i+j*4]);
      cfield1[i+j*3] = xc*xc + yc*yc;
      cfield2[i+j*3] = xc + yc;
  }

  std::shared_ptr<std::vector<double>> nf1(&nfield1),nf2(&nfield2),cf1(&cfield1),cf2(&cfield2),bf(&badfield);
  std::vector<std::shared_ptr<std::vector<double>>> data ={nf1,cf1,cf2,nf2,bf};

  using namespace Portage;
  using namespace Wonton;

  std::vector<Entity_kind> entities = {NODE, CELL, CELL, NODE, NODE};

  Flat_State_Wrapper<double> flatwrap;
  try {
      flatwrap.initialize(names, entities, data);
  } catch (std::exception& err) {
      std::string var_name = "bad";
      std::string msg=std::string("variable ")+var_name+" has incompatible size on add";
      ASSERT_STREQ(err.what(), msg.c_str());
  }

  std::vector<std::string> enames;
  flatwrap.get_names(NODE,enames);
  ASSERT_EQ(enames.size(), 2);
  ASSERT_EQ(enames[0], "nf1");
  ASSERT_EQ(enames[1], "nf2");

  flatwrap.get_names(CELL,enames);
  ASSERT_EQ(enames.size(), 2);
  ASSERT_EQ(enames[0], "cf1");
  ASSERT_EQ(enames[1], "cf2");
}

TEST(Flat_State_Wrapper, VectorInitFail2) {
  std::vector<double> vertx={0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> verty={0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
  std::vector<double>
        &nfield1=*new std::vector<double>(16,0.),
	&nfield2=*new std::vector<double>(16,0.),
	&cfield1=*new std::vector<double>(9,0.),
	&cfield2=*new std::vector<double>(9,0.);
  std::vector<std::string> names = {"nf1", "cf1", "cf2", "nf2"},
                           badnames = {"nf1", "cf1", "cf2"};
  for (size_t i=0; i<16; i++) {
      nfield1[i] = vertx[i]*vertx[i] + verty[i]*verty[i];
      nfield2[i] = vertx[i]+verty[i];
  }
  for (size_t i=0; i<3; i++) for (size_t j=0; j<3; j++) {
      double xc = .25*(vertx[i+j*4] + vertx[i+1+j*4] + vertx[i+j*4] + vertx[i+1+j*4]);
      double yc = .5*(verty[i+(j+1)*4] + verty[i+j*4]);
      cfield1[i+j*3] = xc*xc + yc*yc;
      cfield2[i+j*3] = xc + yc;
  }

  std::shared_ptr<std::vector<double>> nf1(&nfield1),nf2(&nfield2),cf1(&cfield1),cf2(&cfield2);
  std::vector<std::shared_ptr<std::vector<double>>> data ={nf1,cf1,cf2,nf2},
      baddata= {nf1,cf1,cf2};

  using namespace Portage;
  using namespace Wonton;

  std::vector<Entity_kind> entities = {NODE, CELL, CELL, NODE}, badentities={NODE, CELL, CELL};

  try {
      Flat_State_Wrapper<double> flatwrap;
      flatwrap.initialize(badnames, entities, data);
  } catch (std::exception& err) {
      std::string msg="argument sizes do not agree";
      ASSERT_STREQ(err.what(), msg.c_str());
  }

  try {
      Flat_State_Wrapper<double> flatwrap;
      flatwrap.initialize(names, badentities, data);
  } catch (std::exception& err) {
      std::string msg="argument sizes do not agree";
      ASSERT_STREQ(err.what(), msg.c_str());
  }

  try {
      Flat_State_Wrapper<double> flatwrap;
      flatwrap.initialize(names, entities, baddata);
  } catch (std::exception& err) {
      std::string msg="argument sizes do not agree";
      ASSERT_STREQ(err.what(), msg.c_str());
  }
}

TEST(Flat_State_Wrapper, DataTypes2D) {

  // Add multiple state vector types
  int const n_cells = 4;
  double dtest1[] = {1.1, 2.2, 3.3, 4.4};
  double dtest2[] = {1.2, 2.2, 2.3, 2.4};
  double dtest3[] = {1.3, 3.2, 3.3, 3.4};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
  Wonton::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  std::shared_ptr<Jali::State> jali_state(Jali::State::create(inputMesh));
  Wonton::Jali_State_Wrapper jali_state_wrapper(*jali_state);

  jali_state->add("d1", inputMesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, dtest1);
  jali_state->add("d2", inputMesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, dtest2);
  jali_state->add("d3", inputMesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, dtest3);

  Wonton::Flat_State_Wrapper<> flat_state;
  flat_state.initialize(jali_state_wrapper, {"d1", "d2", "d3"});

  // check clear function
  flat_state.initialize(jali_state_wrapper, {"d1", "d2", "d3"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using native wrapper
  double* ddata = nullptr;
  jali_state_wrapper.mesh_get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d1)
  ddata = nullptr;
  flat_state.mesh_get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d2)
  ddata = nullptr;
  flat_state.mesh_get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest2[i]);

  // Get raw float data using the flat mesh wrapper (d3)
  ddata = nullptr;
  flat_state.mesh_get_data(Portage::CELL, "d3", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest3[i]);

  // Check entity types
  Portage::Entity_kind entity;
  entity = flat_state.get_entity((int)0);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity((int)1);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity((int)2);
  ASSERT_EQ(Portage::CELL, entity);
}

TEST(Flat_State_Wrapper, DataTypes3D) {

  // Add multiple state vector types
  int const n_cells = 8;
  double dtest1[] = {1.1, 2.2, 3.3, 4.4, 5, 6, 7, 8};
  double dtest2[] = {1.2, 2.2, 2.3, 2.4, 5.1, 6.1, 7.1, 8.1};
  double dtest3[] = {1.3, 3.2, 3.3, 3.4, 5.2, 6.2, 7.2, 8.2};

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  std::shared_ptr<Jali::Mesh> inputMesh = mf(0.0, 0.0, 0.0,
                                             1.0, 1.0, 1.0,
                                             2, 2, 2);
  Wonton::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  std::shared_ptr<Jali::State> jali_state(Jali::State::create(inputMesh));
  Wonton::Jali_State_Wrapper jali_state_wrapper(*jali_state);

  jali_state->add("d1", inputMesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, dtest1);
  jali_state->add("d2", inputMesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, dtest2);
  jali_state->add("d3", inputMesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, dtest3);

  Wonton::Flat_State_Wrapper<> flat_state;
  flat_state.initialize(jali_state_wrapper, {"d1", "d2", "d3"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using wrapper
  double* ddata = nullptr;
  jali_state_wrapper.mesh_get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d1)
  ddata = nullptr;
  flat_state.mesh_get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d2)
  ddata = nullptr;
  flat_state.mesh_get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest2[i]);

  // Get raw float data using the flat mesh wrapper (d3)
  ddata = nullptr;
  flat_state.mesh_get_data(Portage::CELL, "d3", &ddata);
  for (unsigned int i = 0; i < n_cells; i++) ASSERT_EQ(ddata[i], dtest3[i]);

  // Check entity types
  Portage::Entity_kind entity;
  entity = flat_state.get_entity((int)0);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity((int)1);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity((int)2);
  ASSERT_EQ(Portage::CELL, entity);
}

// test a mesh with one cell containing 3 materials
TEST(Flat_State_Wrapper, mm1) {
	
	// allocate the Jali mesh and state
	std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;
  
  // create the mesh and state
  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 1, 1);
	sourceState = Jali::State::create(sourceMesh);
	
	// create the mesh and state wrappers
	Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  
  // number of materials in the problem
  constexpr int nmats = 3;
  
  // define the material names
  std::string matnames[nmats] = {"mat0", "mat1", "mat2"};
  
  // add the fixed length vector (nmats) of variable length cells containing
  // that material
  // this initializer says all three materials are in our one cell
  // the shape of matcells_src defines the shape of all multimaterial fields
	std::vector<int> matcells_src[nmats]{{0},{0},{0}};  
	
  // add the materials to the state wrapper
  // note that this step implicitly defines the order (and indices) of the materials
 	for (int m = 0; m < nmats; m++)
		sourceStateWrapper.add_material(matnames[m], matcells_src[m]);
		
	// the raw data
	// note that the shape of the data needs to be the same shape as matcells_src
	std::vector<std::vector<double>> density_data {{1.},{10.},{100.}};
	std::vector<std::vector<double>> volume_fraction_data {{.5},{.25},{.25}};
	std::vector<std::vector<Portage::Point<2>>> centroid_data {{{.25,.25}},{{.5,.5}},{{.75,.75}}};
		
	// add the data 
	for (int m = 0; m < nmats; m++){
		sourceStateWrapper.mat_add_celldata("density", m, &(density_data[m][0]));
		sourceStateWrapper.mat_add_celldata("volume fraction", m, &(volume_fraction_data[m][0]));
		sourceStateWrapper.mat_add_celldata("centroid", m, &(centroid_data[m][0]));
	}
		
	// create the flat state wrapper
	Wonton::Flat_State_Wrapper<> flat_state;
	
	// initialize the flat state wrapper
	std::vector<std::string> field_names = {"density", "volume fraction"};
  flat_state.initialize(sourceStateWrapper, field_names);// , "centroid"


	//////////////////////////////////////////////		
  // diagnostics
  //////////////////////////////////////////////
  
	// check number of materials
	ASSERT_EQ(nmats, flat_state.num_materials());
	ASSERT_EQ(sourceStateWrapper.num_materials(), flat_state.num_materials());
	
	// check material names
	for (int m=0; m<nmats; ++m){
		ASSERT_EQ(matnames[m], flat_state.material_name(m));
		ASSERT_EQ(sourceStateWrapper.material_name(m), flat_state.material_name(m));
	}
		
	// check number of cells for each material
	for (int m=0; m<nmats; ++m){
		ASSERT_EQ(matcells_src[m].size(), flat_state.mat_get_num_cells(m));
	}
	
	// check that the cell id's are correct each material
	for (int m=0; m<nmats; ++m){
		
		// temp storage for the cell id's
		std::vector<int> cellids;
		
		// get the cell id's for this material
		flat_state.mat_get_cells(m, &cellids);
		
	  for (int i=0; i<flat_state.mat_get_num_cells(m); ++i){
			ASSERT_EQ(matcells_src[m][i], cellids[i]);
	  }
	}
	
	// check that there are the correct number of materials in the cell
	ASSERT_EQ(flat_state.cell_get_num_mats(0),3);
	
	// check that the materials in the cell are correct
	// these are ugly tests, but work since we know the answer
	std::vector<int> cell_mats;
	flat_state.cell_get_mats(0, &cell_mats);
	ASSERT_EQ(cell_mats[0],0);
	ASSERT_EQ(cell_mats[1],1);
	ASSERT_EQ(cell_mats[2],2);
	
	// check that the material cell indices work
	ASSERT_EQ(flat_state.cell_index_in_material(0,0),0);
	ASSERT_EQ(flat_state.cell_index_in_material(0,1),0);
	ASSERT_EQ(flat_state.cell_index_in_material(0,2),0);
	
	// test the field types
	ASSERT_EQ(
		flat_state.field_type(Portage::Entity_kind::CELL, field_names[0]),
		Portage::Field_type::MULTIMATERIAL_FIELD
	);
	ASSERT_EQ(
		flat_state.field_type(Portage::Entity_kind::CELL, field_names[1]),
		Portage::Field_type::MULTIMATERIAL_FIELD
	);
	
	// test getting the data
	double *data;
	flat_state.mat_get_celldata("density", 0, & data);
	ASSERT_EQ(data[0],1.);
	flat_state.mat_get_celldata("density", 1, &data);
	ASSERT_EQ(data[0],10.);
	flat_state.mat_get_celldata("density", 2, &data);
	ASSERT_EQ(data[0],100.);
	flat_state.mat_get_celldata("volume fraction", 0, &data);
	ASSERT_EQ(data[0],.5);
	flat_state.mat_get_celldata("volume fraction", 1, &data);
	ASSERT_EQ(data[0],.25);
	flat_state.mat_get_celldata("volume fraction", 2, &data);
	ASSERT_EQ(data[0],.25);

}

// test a mesh with 4 cells containing 3 materials
TEST(Flat_State_Wrapper, mm2) {
	
	// allocate the Jali mesh and state
	std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::State> sourceState;
  
  // create the mesh and state
  sourceMesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, 1.0, 1.0, 2, 2);
	sourceState = Jali::State::create(sourceMesh);
	
	// create the mesh and state wrappers
	Wonton::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Wonton::Jali_State_Wrapper sourceStateWrapper(*sourceState);
  
  // number of materials in the problem
  constexpr int nmats = 3;
  
  // define the material names
  std::string matnames[nmats] = {"mat0", "mat1", "mat2"};
  
  // add the fixed length vector (nmats) of variable length cells containing
  // that material
  // this initializer says all three materials are in our one cell
  // the shape of matcells_src defines the shape of all multimaterial fields
	std::vector<int> matcells_src[nmats]{{0,1,2,3},{0,1},{0,2}};  
	
  // add the materials to the state wrapper
  // note that this step implicitly defines the order (and indices) of the materials
 	for (int m = 0; m < nmats; m++)
		sourceStateWrapper.add_material(matnames[m], matcells_src[m]);
		
	// the raw data
	// note that the shape of the data needs to be the same shape as matcells_src
	std::vector<std::vector<double>> density_data {{1.,1.1,1.2,1.3},{10.,10.1},{100.,100.1}};
		
	// add the data 
	for (int m = 0; m < nmats; m++){
		sourceStateWrapper.mat_add_celldata("density", m, &(density_data[m][0]));
	}
	
	// create the flat state wrapper
	Wonton::Flat_State_Wrapper<> flat_state;
	
	// initialize the flat state wrapper
	std::vector<std::string> field_names = {"density"};
  flat_state.initialize(sourceStateWrapper, field_names);


	//////////////////////////////////////////////		
  // diagnostics
  //////////////////////////////////////////////
  
	// check number of materials
	ASSERT_EQ(nmats, flat_state.num_materials());
	ASSERT_EQ(sourceStateWrapper.num_materials(), flat_state.num_materials());
	
	// check material names
	for (int m=0; m<nmats; ++m){
		ASSERT_EQ(matnames[m], flat_state.material_name(m));
		ASSERT_EQ(sourceStateWrapper.material_name(m), flat_state.material_name(m));
	}
			
	// check number of cells for each material
	for (int m=0; m<nmats; ++m){
		ASSERT_EQ(matcells_src[m].size(), flat_state.mat_get_num_cells(m));
	}
	
	// check that the cell id's are correct each material
	for (int m=0; m<nmats; ++m){
		
		// temp storage for the cell id's
		std::vector<int> cellids;
		
		// get the cell id's for this material
		flat_state.mat_get_cells(m, &cellids);
		
	  for (int i=0; i<flat_state.mat_get_num_cells(m); ++i){
			ASSERT_EQ(matcells_src[m][i], cellids[i]);
	  }
	}
	
	// check that there are the correct number of materials in the cell
	ASSERT_EQ(flat_state.cell_get_num_mats(0),3);
	ASSERT_EQ(flat_state.cell_get_num_mats(1),2);
	ASSERT_EQ(flat_state.cell_get_num_mats(2),2);
	ASSERT_EQ(flat_state.cell_get_num_mats(3),1);
	
	// check that the materials in the cell are correct
	// these are ugly tests, but work since we know the answer
	std::vector<int> cell_mats;
	flat_state.cell_get_mats(0, &cell_mats);
	ASSERT_EQ(cell_mats[0],0);
	ASSERT_EQ(cell_mats[1],1);
	ASSERT_EQ(cell_mats[2],2);
	flat_state.cell_get_mats(1, &cell_mats);
	ASSERT_EQ(cell_mats[0],0);
	ASSERT_EQ(cell_mats[1],1);
	flat_state.cell_get_mats(2, &cell_mats);
	ASSERT_EQ(cell_mats[0],0);
	ASSERT_EQ(cell_mats[1],2);
	flat_state.cell_get_mats(3, &cell_mats);
	ASSERT_EQ(cell_mats[0],0);
	
	
	// check that the material cell indices work
	ASSERT_EQ(flat_state.cell_index_in_material(0,0),0);
	ASSERT_EQ(flat_state.cell_index_in_material(1,0),1);
	ASSERT_EQ(flat_state.cell_index_in_material(2,0),2);
	ASSERT_EQ(flat_state.cell_index_in_material(3,0),3);	
	ASSERT_EQ(flat_state.cell_index_in_material(0,1),0);
	ASSERT_EQ(flat_state.cell_index_in_material(1,1),1);
	ASSERT_EQ(flat_state.cell_index_in_material(0,2),0);
	ASSERT_EQ(flat_state.cell_index_in_material(2,2),1);
	
	// test the field types
	ASSERT_EQ(
		flat_state.field_type(Portage::Entity_kind::CELL, field_names[0]),
		Portage::Field_type::MULTIMATERIAL_FIELD
	);

	// test getting the data
	double *data;
	flat_state.mat_get_celldata("density", 0, &data);
	ASSERT_EQ(data[0],density_data[0][0]);
	ASSERT_EQ(data[1],density_data[0][1]);
	ASSERT_EQ(data[2],density_data[0][2]);
	ASSERT_EQ(data[3],density_data[0][3]);
	flat_state.mat_get_celldata("density", 1, &data);
	ASSERT_EQ(data[0],density_data[1][0]);
	ASSERT_EQ(data[1],density_data[1][1]);
	flat_state.mat_get_celldata("density", 2, &data);
	ASSERT_EQ(data[0],density_data[2][0]);
	ASSERT_EQ(data[1],density_data[2][1]);

}
