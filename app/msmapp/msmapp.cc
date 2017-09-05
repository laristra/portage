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
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <fstream>

#include "mpi.h"

#include "portage/driver/driver.h"
#include "portage/driver/driver_mesh_swarm_mesh.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/simple_mesh/simple_state.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"
#include "portage/support/Point.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/wonton/state/simple_state/simple_state_wrapper.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"

template<size_t dim>
struct Controls {
  double smin[dim], smax[dim], tmin[dim], tmax[dim];
  int scells[dim], tcells[dim];
  int order;
  int example;
  double smoothing_factor;
  int print_detail;
};

template<template<class, class, class, Portage::Entity_kind, long> class T>
class consistent_order{
public:
  static bool check(Portage::Meshfree::Basis::Type type) {return false;}
};

template<>
class consistent_order<Portage::Interpolate_1stOrder>{
public:
  static bool check(Portage::Meshfree::Basis::Type type){
    if (type == Portage::Meshfree::Basis::Linear) return true;
    else return false;
  }
};

template<>
class consistent_order<Portage::Interpolate_2ndOrder>{
public:
  static bool check(Portage::Meshfree::Basis::Type type){
    if (type == Portage::Meshfree::Basis::Quadratic) return true;
    else return false;
  }
};


//methods for computing initial field values
template<unsigned int D>
double field_func(int example, Portage::Point<D> coord) {
  double value = 0.0;
  switch (example) {
  case -1: {
    double rsqr = 0.0;
    for (int i = 0; i < D; i++) 
      rsqr += (coord[i]-.5)*(coord[i]-.5);
    value = 1.0;
    for (int i = 0; i < D; i++)
      value *= sin(0.9*2.*3.1415926535898*(coord[i]-.5));
    value *= exp(-1.5*sqrt(rsqr));
    break;
  }
  case 0:
    value = 25.3;
    break;
  case 1:
    for (int i = 0; i < D; i++)
      value += coord[i];
    break;
  case 2:
    for (int i = 0; i < D; i++) 
      for (int j = i; j < D; j++) 
        value += coord[i]*coord[j];
    break;
  case 3:
    for (int i = 0; i < D; i++)
      for (int j = i; j < D; j++) 
        for (int k = j; k < D; k++) 
          value += coord[i]*coord[j]*coord[k];
    break;
  default:
    throw std::runtime_error("Unknown example!");
  }

  return value;
}


// class for running examples
class runMSM {
protected:
  
  //Source and target meshes
  std::shared_ptr<Portage::Simple_Mesh> sourceMesh;
  std::shared_ptr<Portage::Simple_Mesh> targetMesh;
  //Source and target mesh state
  Portage::Simple_State sourceState;
  Portage::Simple_State targetState;
  Portage::Simple_State targetState2;
  // Wrappers for interfacing with the underlying mesh data structures
  Portage::Simple_Mesh_Wrapper sourceMeshWrapper;
  Portage::Simple_Mesh_Wrapper targetMeshWrapper;
  Portage::Simple_State_Wrapper sourceStateWrapper;
  Portage::Simple_State_Wrapper targetStateWrapper;
  Portage::Simple_State_Wrapper targetStateWrapper2;
  // run parameters
  Controls<3> controls_;

public:
  
  // This is the basic test method to be called for each unit test. It will work
  // for 2-D and 3-D, coincident and non-coincident cell-centered remaps.
  template <
    template<class, class> class Intersect,
    template<class, class, class, Portage::Entity_kind, long> class Interpolate,
    template <int, class, class> class SwarmSearch,
    int Dimension=3
    >
  void runit()
  {
    if (Dimension != 3) {
      throw std::runtime_error("2D not available yet");
    }

    double smoothing_factor = controls_.smoothing_factor;
    Portage::Meshfree::Basis::Type basis;
    if (controls_.order == 0) basis = Portage::Meshfree::Basis::Unitary;
    if (controls_.order == 1) basis = Portage::Meshfree::Basis::Linear;
    if (controls_.order == 2) basis = Portage::Meshfree::Basis::Quadratic;
    assert(consistent_order<Interpolate>::check(basis));      

    // Fill the source state data with the specified profile
    const int nsrccells = sourceMeshWrapper.num_owned_cells();
    std::vector<double> sourceData(nsrccells);
    const int nsrcnodes = sourceMeshWrapper.num_owned_nodes();
    std::vector<double> sourceDataNode(nsrcnodes);

    //Create the source data for given function
    Portage::Flat_Mesh_Wrapper<double> sourceFlatMesh;
    sourceFlatMesh.initialize(sourceMeshWrapper);
    for (unsigned int c = 0; c < nsrccells; ++c) {
      Portage::Point<Dimension> cen;
      sourceFlatMesh.cell_centroid(c, &cen);
      sourceData[c] = field_func<3>(controls_.example, cen);
    }
    Portage::Simple_State::vec &sourceVec(sourceState.add("celldata",
                                                          Portage::Entity_kind::CELL, &(sourceData[0])));
    
    for (unsigned int c = 0; c < nsrcnodes; ++c) {
      Portage::Point<Dimension> cen;
      sourceFlatMesh.node_get_coordinates(c, &cen);
      sourceDataNode[c] = field_func<3>(controls_.example, cen);
    }
    Portage::Simple_State::vec &sourceVecNode(sourceState.add("nodedata",
                                                              Portage::Entity_kind::NODE, &(sourceDataNode[0])));

    //Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();
    std::vector<double> targetData(ntarcells), targetData2(ntarcells);
    std::vector<double> targetDataNode(ntarnodes), targetData2Node(ntarnodes);
    Portage::Simple_State::vec &targetVec(targetState.add("celldata",
                                                          Portage::Entity_kind::CELL, &(targetData[0])));
    Portage::Simple_State::vec &targetVec2(targetState2.add("celldata",
                                                            Portage::Entity_kind::CELL, &(targetData2[0])));
    Portage::Simple_State::vec &targetVecNode(targetState.add("nodedata",
                                                              Portage::Entity_kind::NODE, &(targetDataNode[0])));
    Portage::Simple_State::vec &targetVec2Node(targetState2.add("nodedata",
                                                                Portage::Entity_kind::NODE, &(targetData2Node[0])));

    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    remap_fields.push_back("nodedata");

    // Build the mesh-mesh driver data for this mesh type
    Portage::Driver<Portage::SearchKDTree,
                    Intersect,
                    Interpolate,
                    Dimension,
                    Portage::Simple_Mesh_Wrapper, Portage::Simple_State_Wrapper,
                    Portage::Simple_Mesh_Wrapper, Portage::Simple_State_Wrapper>
      mmdriver(sourceMeshWrapper, sourceStateWrapper,
               targetMeshWrapper, targetStateWrapper);
    mmdriver.set_remap_var_names(remap_fields);
    //run on one processor
    mmdriver.run(false);

    // Build the mesh-swarm-mesh driver data for this mesh type
    Portage::MSM_Driver<
      SwarmSearch,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      Dimension,
      Portage::Simple_Mesh_Wrapper, Portage::Simple_State_Wrapper,
      Portage::Simple_Mesh_Wrapper, Portage::Simple_State_Wrapper
      >
      msmdriver(sourceMeshWrapper, sourceStateWrapper,
                targetMeshWrapper, targetStateWrapper2,
                smoothing_factor, basis);
    msmdriver.set_remap_var_names(remap_fields);
    //run on one processor
    msmdriver.run(false);

    //Check the answer
    double stdval, err;
    double totmerr=0., totserr=0.;

    Portage::Simple_State::vec &cellvecout(targetState.get("celldata", Portage::CELL));
    Portage::Simple_State::vec &cellvecout2(targetState2.get("celldata", Portage::CELL));
    Portage::Simple_State::vec &nodevecout(targetState.get("nodedata", Portage::NODE));
    Portage::Simple_State::vec &nodevecout2(targetState2.get("nodedata", Portage::NODE));

    Portage::Flat_Mesh_Wrapper<double> targetFlatMesh;
    targetFlatMesh.initialize(targetMeshWrapper);
    if (controls_.print_detail == 1) {
      std::printf("Cell Centroid-coord-1-2-3 Exact Mesh-Mesh Error Mesh-Swarm-Mesh Error\n");
    }
    for (int c = 0; c < ntarcells; ++c) {
      Portage::Point<Dimension> ccen;
      targetFlatMesh.cell_centroid(c, &ccen);
      double value = field_func<3>(controls_.example, ccen);
      double merror, serror;
      merror = cellvecout[c] - value;
      serror = cellvecout2[c] - value;
      // dump diagnostics for each cell
      if (controls_.print_detail == 1){
	std::printf("cell-data %8d %19.13le %19.13le %19.13le ", c, ccen[0], ccen[1], ccen[2]);
	std::printf("%19.13le %19.13le %19.13le %19.13le %19.13le\n",
		    value, cellvecout[c], merror, cellvecout2[c], serror);
      }
      totmerr = std::max(totmerr, std::fabs(merror));
      totserr = std::max(totserr, std::fabs(serror));
    }

    std::printf("\n\nLinf NORM OF MM CELL ERROR: %le\n\n", totmerr);
    std::printf("\n\nLinf NORM OF MSM CELL ERROR: %le\n\n", totserr);

    if (controls_.print_detail == 1) {
      std::printf("Node Node-coord-1-2-3 Exact Mesh-Mesh Error Mesh-Swarm-Mesh Error\n");
    }
    totmerr = totserr = 0.;
    for (int n = 0; n < ntarnodes; ++n) {
      Portage::Point<Dimension> node;
      targetFlatMesh.node_get_coordinates(n, &node);
      double value = field_func<3>(controls_.example, node);
      double merror, serror;
      merror = nodevecout[n] - value;
      serror = nodevecout2[n] - value;
      // dump diagnostics for each node
      if (controls_.print_detail == 1){
	std::printf("node-data %8d %19.13le %19.13le %19.13le ", n, node[0], node[1], node[2]);
	std::printf("%19.13le %19.13le %19.13le %19.13le %19.13le\n",
		    value, cellvecout[n], merror, cellvecout2[n], serror);
      }
      totmerr = std::max(totmerr, std::fabs(merror));
      totserr = std::max(totserr, std::fabs(serror));
    }

    std::printf("\n\nLinf NORM OF MM NODE ERROR: %lf\n\n", totmerr);
    std::printf("\n\nLinf NORM OF MSM NODE ERROR: %lf\n\n", totserr);
  }

  //Constructor
  runMSM(Controls<3> controls, std::shared_ptr<Portage::Simple_Mesh> s,std::shared_ptr<Portage::Simple_Mesh> t) :
    sourceMesh(s), targetMesh(t), 
    sourceState(sourceMesh), targetState(targetMesh), targetState2(targetMesh),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(sourceState), targetStateWrapper(targetState), 
    targetStateWrapper2(targetState2), 
    controls_(controls)
  {}
};

void usage() {
  std::cout << "Usage: msmapp file" << std::endl;
  std::cout << "\
    example file format: \n\
    note: do not include this line or text before and including \":\" in lines below\n\
    source mesh min, double:          0. 0. 0.\n\
    source mesh max, double:          1. 1. 1.\n\
    target mesh min, double:         .2 .2 .2\n\
    target mesh max, double:         .8 .8 .8\n\
    number source cells, int:        28 32 20\n\
    number target cells, int:        26 30 23\n\
    order, int, choice of 1 or 2:           2\n\
    example, int, choice of -1...3:         3\n\
    double, smoothing factor:             1.5\n\
    print detail, choice of 0 or 1:         1\n\
  ";
}


int main(int argc, char** argv) {
  if (argc < 2) {
    usage();
    return 0;
  }

  Controls<3> ctl;
  std::string filename(argv[1]);
  std::fstream file;
  file.open(filename);
  file >> ctl.smin[0] >> ctl.smin[1] >> ctl.smin[2];
  file >> ctl.smax[0] >> ctl.smax[1] >> ctl.smax[2];
  file >> ctl.tmin[0] >> ctl.tmin[1] >> ctl.tmin[2];
  file >> ctl.tmax[0] >> ctl.tmax[1] >> ctl.tmax[2];
  file >> ctl.scells[0] >> ctl.scells[1] >> ctl.scells[2];
  file >> ctl.tcells[0] >> ctl.tcells[1] >> ctl.tcells[2];
  file >> ctl.order;
  file >> ctl.example;
  file >> ctl.smoothing_factor;
  file >> ctl.print_detail;

  bool error = false;
  for (int i; i<3; i++) {
    if (ctl.smin[i]>=ctl.smax[i]) error = true;
    if (ctl.tmin[i]>=ctl.tmax[i]) error = true;
    if (ctl.scells[i]<0) error = true;
    if (ctl.tcells[i]<0) error = true;
  }
  if (ctl.order<1 or ctl.order>2) error = true;
  if (ctl.example<-1 or ctl.example>3) error = true;
  if (ctl.smoothing_factor<0.) error = true;
  if (ctl.print_detail<0 or ctl.print_detail>1) error = true;
  if (error) {
    throw std::runtime_error("error in input file");
  }

#ifdef ENABLE_MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
    MPI_Init(&argc, &argv);
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  if (numpe > 1) {
    std::cerr << "Particle Remap is only designed for a single process!"
              << std::endl;
    return 1;
  }
#endif

  std::shared_ptr<Portage::Simple_Mesh> src_mesh = std::make_shared<Portage::Simple_Mesh>
    (ctl.smin[0], ctl.smin[1], ctl.smin[2], ctl.smax[0], ctl.smax[1], ctl.smax[2],
     ctl.scells[0], ctl.scells[1], ctl.scells[2]);
  std::shared_ptr<Portage::Simple_Mesh> tgt_mesh = std::make_shared<Portage::Simple_Mesh>
    (ctl.tmin[0], ctl.tmin[1], ctl.tmin[2], ctl.tmax[0], ctl.tmax[1], ctl.tmax[2],
     ctl.tcells[0], ctl.tcells[1], ctl.tcells[2]);
  
  runMSM msmguy(ctl, src_mesh, tgt_mesh);

  if (ctl.order == 1) {
  
    msmguy.runit<Portage::IntersectR3D, 
                 Portage::Interpolate_1stOrder, 
                 Portage::SearchPointsByCells, 3>();
  
  } else if (ctl.order == 2) {
 
    msmguy.runit<Portage::IntersectR3D, 
		 Portage::Interpolate_2ndOrder, 
		 Portage::SearchPointsByCells, 3>();

  }

  std::cout << "starting msmapp..." << std::endl;
  std::cout << "running with input file " << filename << std::endl;
}
