/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#include "gtest/gtest.h"

#include "portage/wonton/state/jali/jali_state_wrapper.h"
#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"
#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/wonton/state/flat/flat_state_wrapper.h"

TEST(SwarmState, basic) {
  using std::make_shared;

  const size_t npoints=10;
  std::vector<Portage::Point<3>> points(npoints), extents(npoints);

    // set up swarm
    double h = 0.01;
    srand(time(NULL));
    for (int i = 0; i < npoints; i++) {
      points[i] = Portage::Point<3>(((double)rand()/RAND_MAX),
                                    ((double)rand()/RAND_MAX),
                                    ((double)rand()/RAND_MAX));
    }
    auto p_ptr = std::make_shared<std::vector<Portage::Point<3>>>(points);
    auto swarm = Portage::Meshfree::Swarm<3>(p_ptr);

    // create state
    using namespace Portage::Meshfree;
    SwarmState<3> state(swarm);
    ASSERT_EQ(state.get_size(), npoints);

    // create state fields
    SwarmState<3>::DblVecPtr dbl_field1 =
        make_shared<SwarmState<3>::DblVec>(npoints,0.);
    SwarmState<3>::DblVecPtr dbl_field2 =
        make_shared<SwarmState<3>::DblVec>(npoints,0.);
    SwarmState<3>::DblVecPtr bad_dbl_field =
        make_shared<SwarmState<3>::DblVec>(npoints+5,0.);
    SwarmState<3>::IntVecPtr int_field1 =
        make_shared<SwarmState<3>::IntVec>(npoints,0.);
    SwarmState<3>::IntVecPtr int_field2 =
        make_shared<SwarmState<3>::IntVec>(npoints,0.);
    SwarmState<3>::IntVecPtr bad_int_field =
        make_shared<SwarmState<3>::IntVec>(npoints+5,0.);

    // fill in fields
    for (size_t i=0; i<npoints; i++) {
      (*dbl_field1)[i] = i+.1;
      (*dbl_field2)[i] = i+.01;
      (*int_field1)[i] = i+10;
      (*int_field2)[i] = i+100;
    }

    // add the fields to the state
    state.add_field("d1", dbl_field1);
    state.add_field("d2", dbl_field2);
    state.add_field("i1", int_field1);
    state.add_field("i2", int_field2);

    // check that fields are correct
    SwarmState<3>::DblVecPtr d1p, d2p;
    SwarmState<3>::IntVecPtr i1p, i2p;
    state.get_field("d1",d1p);
    state.get_field("d2",d2p);
    state.get_field("i1",i1p);
    state.get_field("i2",i2p);
    for (size_t i=0; i<npoints; i++) {
      ASSERT_EQ((*d1p)[i], (*dbl_field1)[i]);
      ASSERT_EQ((*d2p)[i], (*dbl_field2)[i]);
      ASSERT_EQ((*i1p)[i], (*int_field1)[i]);
      ASSERT_EQ((*i2p)[i], (*int_field2)[i]);
    }

    // check names lists are correct
    std::vector<std::string>
      dnames = state.field_names_double();
    std::vector<std::string>
      inames = state.field_names_int();
    ASSERT_EQ(dnames.size(), 2);
    ASSERT_EQ(dnames[0], "d1");
    ASSERT_EQ(dnames[1], "d2");
    ASSERT_EQ(inames.size(), 2);
    ASSERT_EQ(inames[0], "i1");
    ASSERT_EQ(inames[1], "i2");

    // check failure on adding field twice
    try {
      state.add_field("d1", dbl_field1);
      ASSERT_FALSE(true);
    } catch (std::exception err) {
      ASSERT_TRUE(true);
    }
    try {
      state.add_field("i1", int_field1);
      ASSERT_FALSE(true);
    } catch (std::exception err) {
      ASSERT_TRUE(true);
    }

    // check failure on adding bad fields
    try {
      state.add_field("bad", bad_dbl_field);
      ASSERT_FALSE(true);
    } catch (...) {
      ASSERT_TRUE(true);
    }
    try {
      state.add_field("bad", bad_int_field);
      ASSERT_FALSE(true);
    } catch (...) {
      ASSERT_TRUE(true);
    }
}


/*!
  @brief Unit test for constructor with Flat_Mesh_Wrapper in 3D using cells
 */
TEST(Swarm, Flat_Mesh_Wrapper_Cell) {
  using namespace Portage;
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  ASSERT_TRUE(mesh != NULL);
  Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);
  Portage::Flat_Mesh_Wrapper<double> mesh_flat;
  mesh_flat.initialize(mesh_wrapper);

  int ncells = mesh_flat.num_owned_cells();
  int nnodes = mesh_flat.num_owned_nodes();
  std::vector<double> &nfield1 = *new std::vector<double>(nnodes);
  std::vector<double> &cfield1 = *new std::vector<double>(ncells);
  for (int i=0; i<ncells; i++) {
    cfield1[i] = 1.; 
  }
  for (int i=0; i<nnodes; i++) {
    nfield1[i] = 2.; 
  }
  std::shared_ptr<std::vector<double>> nf1(&nfield1),cf1(&cfield1);
  std::vector<std::shared_ptr<std::vector<double>>> data ={nf1,cf1};
  std::vector<Entity_kind> entities = {NODE, CELL};
  std::vector<std::string> names = {"nf1", "cf1"};
  Flat_State_Wrapper<double> state_flat;
  state_flat.initialize(names, entities, data);

  // create swarm state from mesh and state wrappers for cells
  Portage::Meshfree::SwarmState<3> state(mesh_flat, Portage::CELL, state_flat);

  // test size
  ASSERT_EQ(state.get_size(), ncells);

  // check data fields
  std::vector<std::string> intnames = state.field_names_int();
  ASSERT_EQ(intnames.size(), 0);
  std::vector<std::string> dblnames = state.field_names_double();
  ASSERT_EQ(dblnames.size(), 1);
  ASSERT_EQ(dblnames[0], "cf1");
  Meshfree::SwarmState<3>::DblVecPtr field;
  state.get_field("cf1", field);
  for (int i=0; i<ncells; i++) {
    ASSERT_EQ((*field)[i], 1.0);
  }
}


/*!
  @brief Unit test for constructor with Flat_Mesh_Wrapper in 3D using cells
 */
TEST(Swarm, Flat_Mesh_Wrapper_Node) {
  using namespace Portage;
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
                        Jali::Entity_kind::FACE,
                        Jali::Entity_kind::WEDGE,
                        Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  ASSERT_TRUE(mesh != NULL);
  Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);
  Portage::Flat_Mesh_Wrapper<double> mesh_flat;
  mesh_flat.initialize(mesh_wrapper);

  int ncells = mesh_flat.num_owned_cells();
  int nnodes = mesh_flat.num_owned_nodes();
  std::vector<double> &nfield1 = *new std::vector<double>(nnodes);
  std::vector<double> &cfield1 = *new std::vector<double>(ncells);
  for (int i=0; i<ncells; i++) {
    cfield1[i] = 1.; 
  }
  for (int i=0; i<nnodes; i++) {
    nfield1[i] = 2.; 
  }
  std::shared_ptr<std::vector<double>> nf1(&nfield1),cf1(&cfield1);
  std::vector<std::shared_ptr<std::vector<double>>> data ={nf1,cf1};
  std::vector<Entity_kind> entities = {NODE, CELL};
  std::vector<std::string> names = {"nf1", "cf1"};
  Flat_State_Wrapper<double> state_flat;
  state_flat.initialize(names, entities, data);

  // create swarm state from mesh and state wrappers for nodes
  Portage::Meshfree::SwarmState<3> state(mesh_flat, Portage::NODE, state_flat);

  // test size
  ASSERT_EQ(state.get_size(), nnodes);

  // check data fields
  std::vector<std::string> intnames = state.field_names_int();
  ASSERT_EQ(intnames.size(), 0);
  auto dblnames = state.field_names_double();
  ASSERT_EQ(dblnames.size(), 1);
  ASSERT_EQ(dblnames[0], "nf1");
  Meshfree::SwarmState<3>::DblVecPtr field;
  state.get_field("nf1", field);
  for (int i=0; i<nnodes; i++) {
    ASSERT_EQ((*field)[i], 2.0);
  }
}
