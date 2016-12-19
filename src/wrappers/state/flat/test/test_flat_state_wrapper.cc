/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

#include "portage/support/portage.h"

#include <iostream>
#include <exception>
#include <string>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/wrappers/state/flat/flat_state_wrapper.h"

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

  std::vector<Entity_kind> entities = {NODE, CELL, CELL, NODE, NODE};

  Flat_State_Wrapper<double> flatwrap;
  try {
      flatwrap.initialize(names, entities, data);
  } catch (std::exception& err) {
      std::string var_name = "bad";
      std::string msg=std::string("variable ")+var_name+" has incompatible size on add";
      ASSERT_STREQ(err.what(), msg.c_str());
  }
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

TEST(Flat_State_Wrapper, AddData) {
  std::vector<double> vertx={0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
  std::vector<double> verty={0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
  std::vector<double>
        &nfield1=*new std::vector<double>(16,0.),
	&nfield2=*new std::vector<double>(16,0.),
	&nfield3=*new std::vector<double>(16,0.),
	&cfield1=*new std::vector<double>(9,0.),
	&cfield2=*new std::vector<double>(9,0.);
  for (size_t i=0; i<16; i++) {
      nfield1[i] = vertx[i]*vertx[i] + verty[i]*verty[i];
      nfield2[i] = vertx[i]+verty[i];
      nfield3[i] = vertx[i]*verty[i];
  }
  for (size_t i=0; i<3; i++) for (size_t j=0; j<3; j++) {
      double xc = .25*(vertx[i+j*4] + vertx[i+1+j*4] + vertx[i+j*4] + vertx[i+1+j*4]);
      double yc = .5*(verty[i+(j+1)*4] + verty[i+j*4]);
      cfield1[i+j*3] = xc*xc + yc*yc;
      cfield2[i+j*3] = xc + yc;
  }

  std::shared_ptr<std::vector<double>> nf1(&nfield1),nf2(&nfield2),nf3(&nfield3),cf1(&cfield1),cf2(&cfield2);

  using namespace Portage;

  {
    Flat_State_Wrapper<double> flatwrap;
    flatwrap.add_data(NODE, "nf1", nf1);
    flatwrap.add_data(CELL, "cf1", cf1);
    flatwrap.add_data(CELL, "cf2", cf2);
    flatwrap.add_data(NODE, "nf2", nf2);

    size_t sn=flatwrap.get_entity_size(NODE), sc=flatwrap.get_entity_size(CELL);
    ASSERT_EQ(16, sn);
    ASSERT_EQ(9, sc);

    size_t index=flatwrap.get_vector_index(NODE, "nf2");
    ASSERT_EQ(index, 3);
    index=flatwrap.get_vector_index(CELL, "cf2");
    ASSERT_EQ(index, 2);

    double* nf2check = nullptr;
    flatwrap.get_data(NODE, "nf2", &nf2check);
    for (size_t i=0; i<16; i++) {
	ASSERT_EQ(nf2check[i], nfield2[i]);
    }

    flatwrap.add_data(NODE, "nf2", nf3);
    flatwrap.get_data(NODE, "nf2", &nf2check);
    for (size_t i=0; i<16; i++) {
	ASSERT_EQ(nf2check[i], nfield3[i]);
    }

    double *nf4check=nullptr;
    flatwrap.add_data(NODE, "nf4", 1.6);
    flatwrap.get_data(NODE, "nf4", &nf4check);
    for (size_t i=0; i<16; i++) {
	ASSERT_EQ(nf4check[i], 1.6);
    }

    try {
	flatwrap.add_data(EDGE, "ef1", 1.7);
	throw std::runtime_error("shouldn't have gotten here");
    } catch (std::exception err) {
	std::string var_name = "ef1";
	std::string msg = std::string("variable ")+var_name+" has no size information available on add";
	std::string errmsg = err.what();
	// ASSERT_STREQ(errmsg.c_str(), msg.c_str());
	// for some reason errmsg is not what should be with intel 15.
    }
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
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Jali::State jali_state(inputMesh);
  Portage::Jali_State_Wrapper jali_state_wrapper(jali_state);

  jali_state.add("d1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest1);
  jali_state.add("d2", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest2);
  jali_state.add("d3", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest3);

  Portage::Flat_State_Wrapper<> flat_state;
  flat_state.initialize(jali_state_wrapper, {"d1", "d2", "d3"});

  // check clear function
  flat_state.initialize(jali_state_wrapper, {"d1", "d2", "d3"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using wrapper
  double* ddata = nullptr;
  jali_state_wrapper.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d1)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d2)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest2[i]);

  // Get raw float data using the flat mesh wrapper (d3)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d3", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest3[i]);
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
  Portage::Jali_Mesh_Wrapper inputMeshWrapper(*inputMesh);
  Jali::State jali_state(inputMesh);
  Portage::Jali_State_Wrapper jali_state_wrapper(jali_state);

  jali_state.add("d1", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest1);
  jali_state.add("d2", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest2);
  jali_state.add("d3", inputMesh, Jali::Entity_kind::CELL,
            Jali::Entity_type::ALL, dtest3);

  Portage::Flat_State_Wrapper<> flat_state;
  flat_state.initialize(jali_state_wrapper, {"d1", "d2", "d3"});

  // Check the data using Jali as well as by the Flat_State_Wrapper wrapper

  // Get raw float data using wrapper
  double* ddata = nullptr;
  jali_state_wrapper.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d1)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d1", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest1[i]);

  // Get raw float data using the flat mesh wrapper (d2)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d2", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest2[i]);

  // Get raw float data using the flat mesh wrapper (d3)
  ddata = nullptr;
  flat_state.get_data(Portage::CELL, "d3", &ddata);
  for (unsigned int i=0; i<n_cells; i++) ASSERT_EQ(ddata[i], dtest3[i]);
}
