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



#include <iostream>
#include <memory>

#include "gtest/gtest.h"
#include "mpi.h"

#include "portage/driver/driver.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"


double TOL = 1e-12;

//This is a set of integration tests based off of main.cc.  There will be at least one test corresponding to each case found in main.cc.  

//amh: FIXME!! much of this code can be shared using the setup/teardown functions--refactor to remove the duplication of this information between tests.
TEST(Driver, 2D_1stOrderLinearCellCntrConformal1proc){

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;

  mf.included_entities({Jali::Entity_kind::FACE});
  // 2d quad source mesh from (0,0) to (1,1) with nxn zones
  sourceMesh = mf(0.0, 0.0, 1.0, 1.0, 11, 11);
  // 2d quad output mesh from (0,0) to (1,1) with (n+1)x(n+1) zones
  targetMesh = mf(0.0, 0.0, 1.0, 1.0, 3, 3);
  // Wrappers for interfacing with the underlying mesh data structures
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper(*sourceMesh);
  Portage::Jali_Mesh_Wrapper targetMeshWrapper(*targetMesh);

  const int nsrccells = sourceMeshWrapper.num_owned_cells() + sourceMeshWrapper.num_ghost_cells();
  const int ntarcells = targetMeshWrapper.num_owned_cells();
    
  // Fill the source state data with the specified profile
  Jali::State sourceState(sourceMesh);
  std::vector<double> sourceData(nsrccells);

  for (unsigned int c = 0; c < nsrccells; ++c) {
    JaliGeometry::Point cen = sourceMesh->cell_centroid(c);
    sourceData[c] = cen[0]+cen[1];
  }
  sourceState.add("celldata", sourceMesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, &(sourceData[0]));
  Portage::Jali_State_Wrapper sourceStateWrapper(sourceState);
    
  // Build the target state storage
  Jali::State targetState(targetMesh);
  std::vector<double> targetData(ntarcells, 0.0);
  auto& cellvecout = targetState.add("celldata", targetMesh,
                                     Jali::Entity_kind::CELL,
                                     Jali::Entity_type::ALL,
                                     &(targetData[0]));
  Portage::Jali_State_Wrapper targetStateWrapper(targetState);

  // Build the main driver data for this mesh type

  // Register the variable name and interpolation order with the driver
  std::vector<std::string> remap_fields;
  remap_fields.push_back("celldata");
    
  Portage::Driver<Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_1stOrder, 2,
      Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
      Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>  
  d(sourceMeshWrapper, sourceStateWrapper,targetMeshWrapper,
    targetStateWrapper);
  d.set_remap_var_names(remap_fields);    
  //run on one processor
  d.run(false);

  Portage::Point<2> nodexy;
  const int ntarnodes = targetMeshWrapper.num_owned_nodes();
  double stdval, err;
  double toterr=0.;
  for (int c = 0; c < ntarcells; ++c) {
    JaliGeometry::Point ccen = targetMesh->cell_centroid(c);

    double error;
    error = ccen[0]+ccen[1] - cellvecout[c];
    // dump diagnostics for each cell    
    std::printf("Cell=% 4d Centroid = (% 5.3lf,% 5.3lf)", c,
                ccen[0], ccen[1]); 
    std::printf("  Value = % 10.6lf  Err = % lf\n",
              cellvecout[c], error);
    toterr += error*error;
  }

  //amh: FIXME!!  Compare individual, per-node  values/ error norms here
  //amh: FIXME !! Right now this assertion fails--not enough digits

  std::printf("\n\nL2 NORM OF ERROR = %lf\n\n", sqrt(toterr));
  ASSERT_NEAR(0.0095429796560267122, sqrt(toterr), TOL);
}
