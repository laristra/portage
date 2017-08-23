/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
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

  // Check entity types
  Portage::Entity_kind entity;
  entity = flat_state.get_entity(0);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity(1);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity(2);
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

  // Check entity types
  Portage::Entity_kind entity;
  entity = flat_state.get_entity(0);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity(1);
  ASSERT_EQ(Portage::CELL, entity);
  entity = flat_state.get_entity(2);
  ASSERT_EQ(Portage::CELL, entity);
}
