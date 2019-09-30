/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <iostream>
#include <memory>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <cassert>
#include <string>
#include <limits>

#include "gtest/gtest.h"
#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

// portage includes
#include "portage/driver/driver_mesh_swarm_mesh.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/operator.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_wrapper.h"

namespace {

  double TOL = 1e-12;

  TEST(PartByParticle, 2D) {
    const size_t NCELLS = 4;
    std::shared_ptr<Wonton::Simple_Mesh> smesh_ptr = 
      std::make_shared<Wonton::Simple_Mesh>(-1., -1., 1., 1., NCELLS, NCELLS);
    Wonton::Simple_Mesh &smesh = *smesh_ptr;
    Wonton::Simple_Mesh_Wrapper smwrapper(smesh);
    Portage::vector<std::vector<std::vector<double>>> smoothing;
    Portage::vector<Wonton::Point<2>> extents;

    size_t nsource = smesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
    assert(nsource == NCELLS*NCELLS);
    std::vector<double> values(nsource, 1.0);
    for (int i=0; i<nsource; i++) {
      Wonton::Point<2> pnt;
      smwrapper.cell_centroid(i, &pnt);
      if      (pnt[0]<0. and pnt[1]<0.) values[i] = 2.;
      else if (pnt[0]>0. and pnt[1]>0.) values[i] = 2.;
    }
    Wonton::Simple_State sstate(smesh_ptr); 
    double *valptr = values.data();
    std::vector<double> &sadded = sstate.add("indicate", Wonton::CELL, valptr);
    Wonton::Simple_State_Wrapper sswrapper(sstate);

    std::shared_ptr<Wonton::Simple_Mesh> tmesh_ptr = 
      std::make_shared<Wonton::Simple_Mesh>(-1., -1., 1., 1., 2*NCELLS+2, 2*NCELLS+2);
    Wonton::Simple_Mesh &tmesh = *tmesh_ptr;
    Wonton::Simple_Mesh_Wrapper tmwrapper(tmesh);

    Wonton::Simple_State tstate(tmesh_ptr); 
    Wonton::Simple_State_Wrapper tswrapper(tstate);

    size_t ntarget = tmesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
    std::vector<double> tvalues(ntarget);
    double *tvalues_ptr=tvalues.data();
    const double** tvpp=const_cast<const double**>(&tvalues_ptr);
    tswrapper.mesh_add_data(Wonton::CELL, "indicate", tvpp);

    using MSM_Driver_Type =
    Portage::MSM_Driver<
      Portage::SearchPointsByCells,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      2,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper
      >;

    MSM_Driver_Type *msmdriver_ptr;

    msmdriver_ptr = new MSM_Driver_Type
      (smwrapper, sswrapper,
       tmwrapper, tswrapper,
       1.0, 
       1.0, 
       Portage::Meshfree::Weight::FACETED, 
       Portage::Meshfree::Weight::POLYRAMP, 
       Portage::Meshfree::Scatter, 
       std::string("indicate"), 
       0.25);

    MSM_Driver_Type &msmdriver(*msmdriver_ptr);

    msmdriver.set_remap_var_names({"indicate"}, {"indicate"});

    msmdriver.run();

    double *trvalues;
    tswrapper.mesh_get_data(Wonton::CELL, "indicate", &trvalues);

    for (int i=0; i<ntarget; i++) {
      Wonton::Point<2> pnt;
      tmwrapper.cell_centroid(i, &pnt);
      double value=1.0;
      if      (pnt[0]<0. and pnt[1]<0.) value = 2.;
      else if (pnt[0]>0. and pnt[1]>0.) value = 2.;
      ASSERT_NEAR(value, trvalues[i], TOL);
    }
  }

  TEST(PartByParticle, 3D) {
    const size_t NCELLS = 4;
    std::shared_ptr<Wonton::Simple_Mesh> smesh_ptr = 
      std::make_shared<Wonton::Simple_Mesh>(-1., -1., -1., 1., 1., 1., NCELLS, NCELLS, NCELLS);
    Wonton::Simple_Mesh &smesh = *smesh_ptr;
    Wonton::Simple_Mesh_Wrapper smwrapper(smesh);
    Portage::vector<std::vector<std::vector<double>>> smoothing;
    Portage::vector<Wonton::Point<3>> extents;

    size_t nsource = smesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
    assert(nsource == NCELLS*NCELLS*NCELLS);
    std::vector<double> values(nsource, 1.0);
    for (int i=0; i<nsource; i++) {
      Wonton::Point<3> pnt;
      smwrapper.cell_centroid(i, &pnt);
      if      (pnt[0]<0. and pnt[1]<0. and pnt[2]<0.) values[i] = 2.;
      else if (pnt[0]>0. and pnt[1]>0. and pnt[2]>0.) values[i] = 2.;
    }
    Wonton::Simple_State sstate(smesh_ptr); 
    double *valptr = values.data();
    std::vector<double> &sadded = sstate.add("indicate", Wonton::CELL, valptr);
    Wonton::Simple_State_Wrapper sswrapper(sstate);

    std::shared_ptr<Wonton::Simple_Mesh> tmesh_ptr = 
      std::make_shared<Wonton::Simple_Mesh>(-1., -1., -1., 1., 1., 1., 2*NCELLS+2, 2*NCELLS+2, 2*NCELLS+2);
    Wonton::Simple_Mesh &tmesh = *tmesh_ptr;
    Wonton::Simple_Mesh_Wrapper tmwrapper(tmesh);

    Wonton::Simple_State tstate(tmesh_ptr); 
    Wonton::Simple_State_Wrapper tswrapper(tstate);

    size_t ntarget = tmesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
    std::vector<double> tvalues(ntarget);
    double *tvalues_ptr=tvalues.data();
    const double** tvpp=const_cast<const double**>(&tvalues_ptr);
    tswrapper.mesh_add_data(Wonton::CELL, "indicate", tvpp);

    using MSM_Driver_Type =
    Portage::MSM_Driver<
      Portage::SearchPointsByCells,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      3,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper
      >;

    MSM_Driver_Type *msmdriver_ptr;

    msmdriver_ptr = new MSM_Driver_Type
      (smwrapper, sswrapper,
       tmwrapper, tswrapper,
       1.0, 
       1.0, 
       Portage::Meshfree::Weight::FACETED, 
       Portage::Meshfree::Weight::POLYRAMP, 
       Portage::Meshfree::Scatter, 
       std::string("indicate"), 
       0.25);

    MSM_Driver_Type &msmdriver(*msmdriver_ptr);

    msmdriver.set_remap_var_names({"indicate"}, {"indicate"});

    msmdriver.run();

    double *trvalues;
    tswrapper.mesh_get_data(Wonton::CELL, "indicate", &trvalues);

    for (int i=0; i<ntarget; i++) {
      Wonton::Point<3> pnt;
      tmwrapper.cell_centroid(i, &pnt);
      double value=1.0;
      if      (pnt[0]<0. and pnt[1]<0. and pnt[2]<0.) value = 2.;
      else if (pnt[0]>0. and pnt[1]>0. and pnt[2]>0.) value = 2.;
      ASSERT_NEAR(value, trvalues[i], TOL);
    }
  }


}  // namespace
