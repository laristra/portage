/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/



#include "gtest/gtest.h"

#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"

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

    // check failure on adding field twice
    /* deferred
    std::string msg1 = "tried to add double field "+string("d1")+
    " when it already existed";
    try {
      state.add_field("d1", dbl_field1);
      ASSERT_FALSE(true);
    } catch (std::exception err) {
      ASSERT_EQ(err.what(), msg1);
    }
    std::string msg2 = "tried to add int field "+string("i1")+
    " when it already existed";
    try {
      state.add_field("i1", int_field1);
      ASSERT_FALSE(true);
    } catch (std::exception err) {
      ASSERT_EQ(err.what(), msg2);
    }
    */

    // check failure on adding bad fields
    /* deferred
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
    */
}
