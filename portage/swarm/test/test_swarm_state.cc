/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/
#include <vector>
#include <memory>
#include <cstdlib>
#include <string>
#include <exception>

#include "gtest/gtest.h"

// portage includes
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_mm_wrapper.h"
#include "wonton/state/state_vector_uni.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/support/Point.h"

TEST(SwarmState, basic) {
//  using std::make_shared;
  using Portage::Meshfree::SwarmState;
  using namespace Portage::Meshfree;

  // set up a random swarm
  std::random_device device;
  std::mt19937 engine { device() };
  std::uniform_real_distribution<double> generator(0.0, 1.0);

  int const num_points = 10;
  Portage::vector<Wonton::Point<3>> points(num_points);

  for (int i = 0; i < num_points; i++) {
    // point coordinates are not always initialized in order
    // so enforce random number picking sequence
    double const noise[] = { generator(engine),
                             generator(engine),
                             generator(engine) };

    points[i] = Wonton::Point<3>(noise[0], noise[1], noise[2]);
  }

  // create state
  Swarm<3> swarm(points);
  SwarmState<3> state(swarm);
  ASSERT_EQ(state.get_size(), num_points);

  // create state fields
  Portage::vector<double> dbl_field1(num_points, 0.);
  Portage::vector<double> dbl_field2(num_points, 0.);
  Portage::vector<int>    int_field1(num_points, 0.);
  Portage::vector<int>    int_field2(num_points, 0.);
  std::vector<double>     dbl_field3(num_points, 0.);
  std::vector<int>        int_field3(num_points, 0.);

  // fill in fields
  for (int i = 0; i < num_points; i++) {
    dbl_field1[i] = i + 0.10;
    dbl_field2[i] = i + 0.01;
    int_field1[i] = i + 10;
    int_field2[i] = i + 100;
  }

  std::copy(dbl_field1.begin(), dbl_field1.end(), dbl_field3.begin());
  std::copy(int_field1.begin(), int_field1.end(), int_field3.begin());

  // add the fields to the state
  state.add_field("d1", dbl_field1);
  state.add_field("d2", dbl_field2);
  state.add_field("i1", int_field1);
  state.add_field("i2", int_field2);
  state.add_field("d3", dbl_field3);
  state.add_field("i3", int_field3);

  // check that fields are correct
  int i3[num_points];
  double d3[num_points];

  auto d1 = state.get_field_dbl("d1");
  auto d2 = state.get_field_dbl("d2");
  auto i1 = state.get_field_int("i1");
  auto i2 = state.get_field_int("i2");
  state.copy_field("d3", d3);
  state.copy_field("i3", i3);

  for (int i = 0; i < num_points; i++) {
    ASSERT_EQ(i1[i], int_field1[i]);
    ASSERT_EQ(i2[i], int_field2[i]);
    ASSERT_EQ(i3[i], i1[i]);
    ASSERT_DOUBLE_EQ(d1[i], dbl_field1[i]);
    ASSERT_DOUBLE_EQ(d2[i], dbl_field2[i]);
    ASSERT_DOUBLE_EQ(d3[i], d1[i]);
  }

  // check names lists are correct
  auto dnames = state.get_field_names<double>();
  auto inames = state.get_field_names<int>();
  ASSERT_EQ(dnames.size(), 3);
  ASSERT_EQ(dnames[0], "d1");
  ASSERT_EQ(dnames[1], "d2");
  ASSERT_EQ(dnames[2], "d3");
  ASSERT_EQ(inames.size(), 3);
  ASSERT_EQ(inames[0], "i1");
  ASSERT_EQ(inames[1], "i2");
  ASSERT_EQ(inames[2], "i3");

  // check creation by size alone
  SwarmState<3> state2(num_points);
  state2.add_field("d1", dbl_field1);
  auto d1p2 = state2.get_field_dbl("d1");
  ASSERT_EQ(d1p2.size(), num_points);
}


/*!
  @brief Unit test for constructor with Simple_State_Wrapper in 3D using cells
*/
TEST(SwarmState, Simple_State_Wrapper) {

  using namespace Portage::Meshfree;
  using Mesh = Wonton::Simple_Mesh;
  using Wrapper = Wonton::Simple_Mesh_Wrapper;
  using State = Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>;

  Mesh mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  Wrapper mesh_wrapper(mesh);
  State mesh_state(mesh_wrapper);

  int ncells = mesh_wrapper.num_owned_cells();
  int nnodes = mesh_wrapper.num_owned_nodes();
  std::vector<double> cfield1(ncells, 1.);
  std::vector<double> nfield1(nnodes, 2.);

  mesh_state.add(std::make_shared<Wonton::StateVectorUni<>>("cf1", Wonton::CELL, cfield1));
  mesh_state.add(std::make_shared<Wonton::StateVectorUni<>>("nf1", Wonton::NODE, nfield1));

  {
    SwarmState<3> state(mesh_state, Wonton::CELL);

    auto intnames = state.get_field_names<int>();
    auto dblnames = state.get_field_names<double>();

    ASSERT_EQ(state.get_size(), ncells);
    ASSERT_EQ(intnames.size(), 0);
    ASSERT_EQ(dblnames.size(), 1);
    ASSERT_EQ(dblnames[0], "cf1");

    auto field = state.get_field_dbl("cf1");
    for (int i = 0; i < ncells; i++)
      ASSERT_EQ(field[i], 1.0);
  }

  // create swarm state from mesh state wrapper for nodes
  {
    SwarmState<3> state(mesh_state, Wonton::NODE);

    auto intnames = state.get_field_names<int>();
    auto dblnames = state.get_field_names<double>();

    ASSERT_EQ(state.get_size(), nnodes);
    ASSERT_EQ(intnames.size(), 0);
    ASSERT_EQ(dblnames.size(), 1);
    ASSERT_EQ(dblnames[0], "nf1");

    auto field = state.get_field_dbl("nf1");
    for (int i=0; i < nnodes; i++)
      ASSERT_EQ(field[i], 2.0);
  }
}

