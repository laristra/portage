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
#include "wonton/state/simple/simple_state_wrapper.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/support/Point.h"

TEST(SwarmState, Multiple2D) {

  using namespace Portage::Meshfree;

  Wonton::Simple_Mesh mesh0(0.0, 0.0, 1.0, 1.0, 4, 4);
  Wonton::Simple_Mesh mesh1(1.0, 0.0, 2.0, 1.0, 4, 5);
  Wonton::Simple_Mesh mesh2(1.0, 1.0, 2.0, 2.0, 4, 6);

  Wonton::Simple_State state0(std::make_shared<Wonton::Simple_Mesh>(mesh0));
  Wonton::Simple_State state1(std::make_shared<Wonton::Simple_Mesh>(mesh1));
  Wonton::Simple_State state2(std::make_shared<Wonton::Simple_Mesh>(mesh2));

  double f00[25], f01[25], f10[30], f11[30], f20[35], f21[35];
  for (int i=0; i<90; i++) {
    if (i<25) {
      f00[i] = i;
      f01[i] = 1000.+i;
    } else if (i<55) {
      f10[i-25] = i;
      f11[i-25] = 1000.+i;
    } else if (i<90) {
      f20[i-55] = i;
      f21[i-55] = 1000.+i;
    }
  }

  state0.add("f0", Wonton::NODE, f00);
  state0.add("f1", Wonton::NODE, f01);
  state1.add("f0", Wonton::NODE, f10);
  state1.add("f1", Wonton::NODE, f11);
  state2.add("f0", Wonton::NODE, f20);
  state2.add("f1", Wonton::NODE, f21);

  Wonton::Simple_Mesh_Wrapper  mesh_wrapper0(mesh0);
  Wonton::Simple_Mesh_Wrapper  mesh_wrapper1(mesh1);
  Wonton::Simple_Mesh_Wrapper  mesh_wrapper2(mesh2);
  Wonton::Simple_State_Wrapper state_wrapper0(state0);
  Wonton::Simple_State_Wrapper state_wrapper1(state1);
  Wonton::Simple_State_Wrapper state_wrapper2(state2);

  std::vector<Wonton::Simple_Mesh_Wrapper*>  mesh_wrappers = { &mesh_wrapper0,
                                                               &mesh_wrapper1,
                                                               &mesh_wrapper2 };

  std::vector<Wonton::Simple_State_Wrapper*> state_wrappers = { &state_wrapper0,
                                                                &state_wrapper1,
                                                                &state_wrapper2 };

  Swarm<2> swarm(mesh_wrappers, Wonton::NODE);
  SwarmState<2> swarm_state(state_wrappers, Wonton::NODE);

  auto g0 = swarm_state.get_field_double("f0");
  auto g1 = swarm_state.get_field_double("f1");

  for (int i = 0; i < 90; i++) {
    ASSERT_EQ(g0[i], i);
    ASSERT_EQ(g1[i], 1000.+i);
  }
}
