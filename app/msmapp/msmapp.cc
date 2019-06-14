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
#include <fstream>
#include <vector>

#ifdef PORTAGE_ENABLE_MPI
#include <mpi.h>
#else
#define PORTAGE_SERIAL_ONLY
#endif

// portage includes
#include "portage/support/portage.h"
#include "portage/driver/mmdriver.h"
#include "portage/driver/driver_mesh_swarm_mesh.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/search/search_points_by_cells.h"
#include "portage/accumulate/accumulate.h"
#include "portage/estimate/estimate.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_mm_wrapper.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliState.h"

template<size_t dim>
struct Controls {
  // See the usage() function below for the meaning of these quantities.
  double smin[dim], smax[dim], tmin[dim], tmax[dim];
  int scells[dim], tcells[dim];
  int order;
  int example;
  double smoothing_factor;
  int print_detail;
  std::string source_file="none", target_file="none";
  std::string oper8tor="none";
  bool domeshmesh=true;
};

template<template<int, Portage::Entity_kind, class, class, class,
                  template<class, int, class, class> class,
                  class, class, class> class T>
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


// Class for performing mesh-mesh and mesh-swarm-mesh remapping comparisons using simple mesh.
class runMSM {
protected:
  //Source and target meshes
  std::shared_ptr<Wonton::Simple_Mesh> sourceMesh;
  std::shared_ptr<Wonton::Simple_Mesh> targetMesh;
  // Wrappers for interfacing with the underlying mesh data structures
  Wonton::Simple_Mesh_Wrapper sourceMeshWrapper;
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper;
  //Source and target mesh state
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> sourceStateWrapper;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> targetStateWrapper;
  Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper> targetStateWrapper2;
  // run parameters
  Controls<3> controls_;

public:

  // Perform comparison. First remap using traditional mesh techniques. Second remap using particles as intermediary.
  // Report timing and accuracy.
  template <
  template<Portage::Entity_kind, class, class, class,
  template <class, int, class, class> class,
  class, class> class Intersect,
  template<int, Portage::Entity_kind, class, class, class,
  template<class, int, class, class> class,
  class, class, class> class Interpolate,
  template <int, class, class> class SwarmSearch,
  int Dimension=3
  >
  void runit()
  {
    if (Dimension != 3) {
      throw std::runtime_error("2D not available yet");
    }

    // process controls
    double smoothing_factor = controls_.smoothing_factor;
    Portage::Meshfree::Basis::Type basis;
    if (controls_.order == 0) basis = Portage::Meshfree::Basis::Unitary;
    if (controls_.order == 1) basis = Portage::Meshfree::Basis::Linear;
    if (controls_.order == 2) basis = Portage::Meshfree::Basis::Quadratic;
    assert(consistent_order<Interpolate>::check(basis));
    Portage::Meshfree::Operator::Type oper8tor;
    if (controls_.oper8tor == "VolumeIntegral") {
      oper8tor = Portage::Meshfree::Operator::VolumeIntegral;
    } else if (controls_.oper8tor == "none") {
      oper8tor = Portage::Meshfree::Operator::LastOperator;
    } else {
      throw std::runtime_error("illegal operator specified");
    }

    // Fill the source state data with the specified profile
    const int nsrccells = sourceMeshWrapper.num_owned_cells();
    std::vector<double> sourceData(nsrccells);
    const int nsrcnodes = sourceMeshWrapper.num_owned_nodes();
    std::vector<double> sourceDataNode(nsrcnodes);

    //Create the source data for given function
    for (unsigned int c = 0; c < nsrccells; ++c) {
      Portage::Point<Dimension> cen;
      sourceMeshWrapper.cell_centroid(c, &cen);
      sourceData[c] = field_func<3>(controls_.example, cen);
    }
    sourceStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"celldata", Portage::Entity_kind::CELL, sourceData)
    );
    
    for (unsigned int c = 0; c < nsrcnodes; ++c) {
      Portage::Point<Dimension> cen;
      sourceMesh->node_get_coordinates(c, &cen);
      sourceDataNode[c] = field_func<3>(controls_.example, cen);
    }
    sourceStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"nodedata", Portage::Entity_kind::NODE, sourceDataNode));

    //Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();
    std::vector<double> targetData(ntarcells), targetData2(ntarcells);
    std::vector<double> targetDataNode(ntarnodes), targetData2Node(ntarnodes);
    targetStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"celldata", Portage::Entity_kind::CELL, targetData)
    );
    targetStateWrapper.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"nodedata", Portage::Entity_kind::NODE, targetDataNode)
    );

    targetStateWrapper2.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"celldata", Portage::Entity_kind::CELL, targetData2)
    );
    targetStateWrapper2.add(std::make_shared<Wonton::StateVectorUni<>>(
    	"nodedata", Portage::Entity_kind::NODE, targetData2Node)
    );

    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    remap_fields.push_back("nodedata");

    // If an operator is requested, collect the information required.
    Portage::vector<std::vector<Portage::Point<Dimension>>> data;
    Portage::vector<Portage::Meshfree::Operator::Domain> domains;
    if (controls_.oper8tor == "VolumeIntegral") {
      int numcells = targetMesh->num_entities(Portage::Entity_kind::CELL,
                                              Portage::Entity_type::ALL);
      domains.resize(numcells);
      data.resize(numcells);

      for (int c=0; c<numcells; c++) {
	// get integration domains
        std::vector<Portage::Point<Dimension>> cellverts;
        targetMesh->cell_get_coordinates(c, &cellverts);
        data[c] = cellverts;
        domains[c] = Portage::Meshfree::Operator::domain_from_points<Dimension>(cellverts);
      }
    }

    // Build the mesh-mesh driver data for this mesh type
    Portage::MMDriver<Portage::SearchKDTree,
                    Intersect,
                    Interpolate,
                    Dimension,
                    Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>,
                    Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>>
      driver(sourceMeshWrapper, sourceStateWrapper,
             targetMeshWrapper, targetStateWrapper);
    driver.set_remap_var_names(remap_fields);
    driver.run();  // executor argument defaults to nullptr -> run in serial


    // Build the mesh-swarm-mesh driver data for this mesh type
    Portage::MSM_Driver<
      SwarmSearch,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      Dimension,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>,
      Wonton::Simple_Mesh_Wrapper, Wonton::Simple_State_Wrapper<Wonton::Simple_Mesh_Wrapper>
      >
      msmdriver(sourceMeshWrapper, sourceStateWrapper,
                targetMeshWrapper, targetStateWrapper2,
                smoothing_factor);
    Portage::Meshfree::EstimateType estimate=Portage::Meshfree::LocalRegression;
    if (oper8tor == Portage::Meshfree::Operator::VolumeIntegral)
      estimate=Portage::Meshfree::OperatorRegression;
    msmdriver.set_remap_var_names(remap_fields, remap_fields,
                                  estimate, basis, oper8tor, domains, data);
    //run on one processor
    msmdriver.run();

    //Check the answer
    double stdval, err;
    double totmerr=0., totserr=0., totint=0.;

    auto& cellvecout = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper.get("celldata"))->get_data();
    auto& cellvecout2 = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper2.get("celldata"))->get_data();
    auto& nodevecout = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper.get("nodedata"))->get_data();
    auto& nodevecout2 = std::static_pointer_cast<Wonton::StateVectorUni<>>(targetStateWrapper2.get("nodedata"))->get_data();

    if (controls_.domeshmesh) {
      if (controls_.print_detail == 1) {
        std::printf("Cell Centroid-coord-1-2-3 Exact Mesh-Mesh Error Mesh-Swarm-Mesh Error\n");
      }
      for (int c = 0; c < ntarcells; ++c) {
        Portage::Point<Dimension> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);
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
      totmerr = totserr = totint = 0.;
      for (int n = 0; n < ntarnodes; ++n) {
        Portage::Point<Dimension> node;
        targetMesh->node_get_coordinates(n, &node);
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
    } else {
      if (controls_.print_detail == 1) {
        std::printf("Cell Centroid-coord-1-2-3 Exact Mesh-Swarm-Mesh Error\n");
      }
      totmerr = totserr = 0.;
      for (int c = 0; c < ntarcells; ++c) {
        Portage::Point<Dimension> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);
        double value = field_func<3>(controls_.example, ccen);
        double serror;
        serror = cellvecout2[c] - value;
        // dump diagnostics for each cell
        if (controls_.print_detail == 1){
          std::printf("cell-data %8d %19.13le %19.13le %19.13le ", c, ccen[0], ccen[1], ccen[2]);
          if (oper8tor != Portage::Meshfree::Operator::VolumeIntegral) {
            std::printf("%19.13le %19.13le %19.13le\n", value, cellvecout2[c], serror);
          } else {
            std::printf("%19.13le %19.13le\n", value, cellvecout2[c]);
          }
        }
        totserr = std::max(totserr, std::fabs(serror));
      }
      // accumulate integral
      totint = 0.;
      if (oper8tor == Portage::Meshfree::Operator::VolumeIntegral) {
        for (int c = 0; c < ntarcells; ++c) {
          totint += cellvecout2[c];
        }
      }

      std::printf("\n\nLinf NORM OF MSM CELL ERROR: %le\n\n", totserr);
      if (oper8tor == Portage::Meshfree::Operator::VolumeIntegral)
        std::printf("\n\nTOTAL INTEGRAL: %le\n\n", totint);

      if (controls_.print_detail == 1) {
        std::printf("Node Node-coord-1-2-3 Exact Mesh-Swarm-Mesh Error\n");
      }
      totserr = 0.;
      for (int n = 0; n < ntarnodes; ++n) {
        Portage::Point<Dimension> node;
        targetMesh->node_get_coordinates(n, &node);
        double value = field_func<3>(controls_.example, node);
        double serror;
        serror = nodevecout2[n] - value;
        // dump diagnostics for each node
        if (controls_.print_detail == 1){
          std::printf("node-data %8d %19.13le %19.13le %19.13le ", n, node[0], node[1], node[2]);
          std::printf("%19.13le %19.13le %19.13le\n", value, cellvecout2[n], serror);
        }
        totserr = std::max(totserr, std::fabs(serror));
      }

      std::printf("\n\nLinf NORM OF MSM NODE ERROR: %lf\n\n", totserr);
    }
  }

  //Constructor
  runMSM(Controls<3> controls, std::shared_ptr<Wonton::Simple_Mesh> s,std::shared_ptr<Wonton::Simple_Mesh> t) :
    sourceMesh(s), targetMesh(t), 
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(sourceMeshWrapper), targetStateWrapper(targetMeshWrapper),
    targetStateWrapper2(targetMeshWrapper),
    controls_(controls)
  {}
};

// Class for performing mesh-mesh and mesh-swarm-mesh remapping comparisons using jali mesh.
class runMSMJali {
protected:
  //Source and target meshes
  std::shared_ptr<Jali::Mesh> sourceMesh;
  std::shared_ptr<Jali::Mesh> targetMesh;
  //Source and target mesh state
  std::shared_ptr<Jali::State> sourceState;
  std::shared_ptr<Jali::State> targetState;
  std::shared_ptr<Jali::State> targetState2;
protected:
  //Source and target mesh and state wrappers
  Portage::Jali_Mesh_Wrapper sourceMeshWrapper;
  Portage::Jali_Mesh_Wrapper targetMeshWrapper;
  Portage::Jali_State_Wrapper sourceStateWrapper;
  Portage::Jali_State_Wrapper targetStateWrapper;
  Portage::Jali_State_Wrapper targetStateWrapper2;
  // run parameters
  Controls<3> controls_;

 public:

  // Perform comparison. First remap using traditional mesh techniques. Second remap using particles as intermediary.
  // Report timing and accuracy.
  template <
    template<Portage::Entity_kind, class, class, class,
    template<class, int, class, class> class,
    class, class> class Intersect,
    template<int, Portage::Entity_kind, class, class, class,
    template<class, int, class, class> class,
    class, class, class> class Interpolate,
    template <int, class, class> class SwarmSearch,
    int Dimension=3
    >
  void runit()
  {
    if (Dimension != 3) {
      throw std::runtime_error("2D not available yet");
    }

    // process controls
    double smoothing_factor = controls_.smoothing_factor;
    Portage::Meshfree::Basis::Type basis;
    if (controls_.order == 0) basis = Portage::Meshfree::Basis::Unitary;
    if (controls_.order == 1) basis = Portage::Meshfree::Basis::Linear;
    if (controls_.order == 2) basis = Portage::Meshfree::Basis::Quadratic;
    assert(consistent_order<Interpolate>::check(basis));
    Portage::Meshfree::Operator::Type oper8tor;
    if (controls_.oper8tor == "VolumeIntegral") {
      oper8tor = Portage::Meshfree::Operator::VolumeIntegral;
    } else if (controls_.oper8tor != "none") {
      throw std::runtime_error("illegal operator specified");
    }

    // Fill the source state data with the specified profile
    const int nsrccells = sourceMeshWrapper.num_owned_cells();
    std::vector<double> sourceData(nsrccells);
    const int nsrcnodes = sourceMeshWrapper.num_owned_nodes();
    std::vector<double> sourceDataNode(nsrcnodes);

    //Create the source data for given function
    for (unsigned int c = 0; c < nsrccells; ++c) {
      Portage::Point<Dimension> cen;
      sourceMeshWrapper.cell_centroid(c, &cen);
      sourceData[c] = field_func<3>(controls_.example, cen);
    }
    Jali::UniStateVector<double> &sourceVec(sourceState->add("celldata",
      sourceMesh, Jali::Entity_kind::CELL, Jali::Entity_type::ALL, &(sourceData[0])));

    for (unsigned int c = 0; c < nsrcnodes; ++c) {
      Portage::Point<Dimension> cen;
      sourceMeshWrapper.node_get_coordinates(c, &cen);
      sourceDataNode[c] = field_func<3>(controls_.example, cen);
    }
    Jali::UniStateVector<double> &sourceVecNode(sourceState->add("nodedata",
      sourceMesh, Jali::Entity_kind::NODE, Jali::Entity_type::ALL, &(sourceDataNode[0])));

    //Build the target state storage
    const int ntarcells = targetMeshWrapper.num_owned_cells();
    const int ntarnodes = targetMeshWrapper.num_owned_nodes();
    std::vector<double> targetData(ntarcells), targetData2(ntarcells);
    std::vector<double> targetDataNode(ntarnodes), targetData2Node(ntarnodes);
    targetStateWrapper.mesh_add_data<double>(Portage::Entity_kind::CELL, "celldata", &(targetData[0]));
    targetStateWrapper2.mesh_add_data<double>(Portage::Entity_kind::CELL, "celldata", &(targetData2[0]));
    targetStateWrapper.mesh_add_data<double>(Portage::Entity_kind::NODE, "nodedata", &(targetDataNode[0]));
    targetStateWrapper2.mesh_add_data<double>(Portage::Entity_kind::NODE, "nodedata", &(targetData2Node[0]));

    // Register the variable name and interpolation order with the driver
    std::vector<std::string> remap_fields;
    remap_fields.push_back("celldata");
    remap_fields.push_back("nodedata");

    // If an operator is requested, collect the information required.
    Portage::vector<std::vector<Portage::Point<Dimension>>> data;
    Portage::vector<Portage::Meshfree::Operator::Domain> domains;
    if (controls_.oper8tor == "VolumeIntegral") {
      int numcells = targetMeshWrapper.num_owned_cells();
      domains.resize(numcells);
      data.resize(numcells);

      for (int c=0; c<numcells; c++) {
	// get integration domains
        std::vector<int> cellnodes;
        targetMesh->cell_get_nodes(c, &cellnodes);
        std::vector<Portage::Point<Dimension>> cellverts(cellnodes.size());
        for (int i=0; i<cellnodes.size(); i) {
          targetMeshWrapper.node_get_coordinates(cellnodes[i], &cellverts[i]);
        }
        data[c] = cellverts;
        domains[c] = Portage::Meshfree::Operator::domain_from_points<Dimension>(cellverts);
      }
    }

    // Build the mesh-mesh driver data for this mesh type
    Portage::MMDriver<Portage::SearchKDTree,
                    Intersect,
                    Interpolate,
                    Dimension,
                    Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
                    Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper>
    driver(sourceMeshWrapper, sourceStateWrapper,
           targetMeshWrapper, targetStateWrapper);
    driver.set_remap_var_names(remap_fields);
    driver.run();  // executor argument to driver defaults to nullptr -> serial

    // Build the mesh-swarm-mesh driver data for this mesh type
    Portage::MSM_Driver<
      SwarmSearch,
      Portage::Meshfree::Accumulate,
      Portage::Meshfree::Estimate,
      Dimension,
      Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper,
      Portage::Jali_Mesh_Wrapper, Portage::Jali_State_Wrapper
      >
      msmdriver(sourceMeshWrapper, sourceStateWrapper,
                targetMeshWrapper, targetStateWrapper2,
                smoothing_factor);
    msmdriver.set_remap_var_names(remap_fields, remap_fields,
                                  Portage::Meshfree::LocalRegression, basis,
                                  oper8tor, domains, data);
    // run on one processor
    msmdriver.run();

    //Check the answer
    double stdval, err;
    double totmerr=0., totserr=0., totint=0.;

    std::vector<double> cellvecout(ntarcells), cellvecout2(ntarcells);
    std::vector<double> nodevecout(ntarnodes), nodevecout2(ntarnodes);
    Jali::UniStateVector<double, Jali::Mesh> cvp, cv2p, nvp, nv2p;
    targetState->get("celldata", targetMesh, Jali::Entity_kind::CELL, Jali::Entity_type::ALL, &cvp);
    targetState2->get("celldata", targetMesh, Jali::Entity_kind::CELL, Jali::Entity_type::ALL, &cv2p);
    targetState->get("nodedata", targetMesh, Jali::Entity_kind::NODE, Jali::Entity_type::ALL, &nvp);
    targetState2->get("nodedata", targetMesh, Jali::Entity_kind::NODE, Jali::Entity_type::ALL, &nv2p);
    for (int i=0; i<ntarcells; i++) {cellvecout[i]=cvp[i]; cellvecout2[i]=cv2p[i];}
    for (int i=0; i<ntarnodes; i++) {nodevecout[i]=nvp[i]; nodevecout2[i]=nv2p[i];}

    if (controls_.domeshmesh) {
      if (controls_.print_detail == 1) {
        std::printf("Cell Centroid-coord-1-2-3 Exact Mesh-Mesh Error Mesh-Swarm-Mesh Error\n");
      }
      for (int c = 0; c < ntarcells; ++c) {
        Portage::Point<Dimension> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);
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
        targetMeshWrapper.node_get_coordinates(n, &node);
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
    } else {
      if (controls_.print_detail == 1) {
        std::printf("Cell Centroid-coord-1-2-3 Exact Mesh-Swarm-Mesh Error\n");
      }
      totserr = totint = 0.;
      for (int c = 0; c < ntarcells; ++c) {
        Portage::Point<Dimension> ccen;
        targetMeshWrapper.cell_centroid(c, &ccen);
        double value = field_func<3>(controls_.example, ccen);
        double serror;
        serror = cellvecout2[c] - value;
        // accumulate integral
        if (oper8tor == Portage::Meshfree::Operator::VolumeIntegral) {
          totint += cellvecout2[c];
        }
        // dump diagnostics for each cell
        if (controls_.print_detail == 1){
          std::printf("cell-data %8d %19.13le %19.13le %19.13le ", c, ccen[0], ccen[1], ccen[2]);
          if (oper8tor != Portage::Meshfree::Operator::VolumeIntegral) {
            std::printf("%19.13le %19.13le %19.13le\n", value, cellvecout2[c], serror);
          } else {
            std::printf("%19.13le %19.13le\n", value, cellvecout2[c]);
          }
        }
        totserr = std::max(totserr, std::fabs(serror));
      }
      std::printf("\n\nLinf NORM OF MSM CELL ERROR: %le\n\n", totserr);
      if (oper8tor == Portage::Meshfree::Operator::VolumeIntegral)
        std::printf("\n\nTOTAL INTEGRAL: %le\n\n", totint);

      if (controls_.print_detail == 1) {
        std::printf("Node Node-coord-1-2-3 Exact Mesh-Swarm-Mesh Error\n");
      }
      totserr = 0.;
      for (int n = 0; n < ntarnodes; ++n) {
        Portage::Point<Dimension> node;
        targetMeshWrapper.node_get_coordinates(n, &node);
        double value = field_func<3>(controls_.example, node);
        double serror;
        serror = nodevecout2[n] - value;
        // dump diagnostics for each node
        if (controls_.print_detail == 1){
          std::printf("node-data %8d %19.13le %19.13le %19.13le ", n, node[0], node[1], node[2]);
          std::printf("%19.13le %19.13le %19.13le\n",
                      value, cellvecout2[n], serror);
        }
        totserr = std::max(totserr, std::fabs(serror));
      }

      std::printf("\n\nLinf NORM OF MSM NODE ERROR: %lf\n\n", totserr);
    }
  }

  //Constructor
  runMSMJali(Controls<3> controls, std::shared_ptr<Jali::Mesh> s,std::shared_ptr<Jali::Mesh> t) :
    sourceMesh(s), targetMesh(t),
    sourceState(Jali::State::create(sourceMesh)),
    targetState(Jali::State::create(targetMesh)),
    targetState2(Jali::State::create(targetMesh)),
    sourceMeshWrapper(*sourceMesh), targetMeshWrapper(*targetMesh),
    sourceStateWrapper(*sourceState), targetStateWrapper(*targetState),
    targetStateWrapper2(*targetState2),
    controls_(controls)
  {}
};

void usage() {
  std::cout << "Usage: msmapp file [-domm | -nomm]" << std::endl;
  std::cout << "\
    Uses specifications in \"file\" to perform a mesh-mesh and also a mesh-swarm-mesh remap and compare.\n\
    Options:\n\
      -domm do the direct mesh-mesh remap\n\
      -nomm do not do the direct mesh-mesh remap\n\
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
    (optional) source mesh EXODUS file:  src_file_name \n\
    (optional) target mesh EXODUS file:  tgt_file_name \n\
    (optional) operator specification: VolumeIntegral \n\
    \n\
    must specify both source and target mesh files if so using \n\
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
  try {
    file >> ctl.source_file;
    if (ctl.source_file != "VolumeIntegral") {
      file >> ctl.target_file;
      try {
        file >> ctl.oper8tor;
      } catch (...) {
      }
    } else {
      ctl.oper8tor = "VolumeIntegral";
      ctl.source_file = "none";
    }
  } catch (...) {
  }

  if (argc >= 3) {
    if (argv[2] == std::string("-domm")) ctl.domeshmesh=true;
    else if (argv[2] == std::string("-nomm")) ctl.domeshmesh=false;
  }

  bool error = false;
  for (size_t i=0; i<3; i++) {
    if (ctl.smin[i]>=ctl.smax[i]) error = true;
    if (ctl.tmin[i]>=ctl.tmax[i]) error = true;
    if (ctl.scells[i]<0) error = true;
    if (ctl.tcells[i]<0) error = true;
  }
  if (ctl.order<1 or ctl.order>2) error = true;
  if (ctl.example<-1 or ctl.example>3) error = true;
  if (ctl.smoothing_factor<0.) error = true;
  if (ctl.print_detail<0 or ctl.print_detail>1) error = true;
  if (ctl.source_file == "none" and ctl.target_file != "none") error = true;
  if (ctl.target_file == "none" and ctl.source_file != "none") error = true;
  if (error) {
    throw std::runtime_error("error in input file");
  }

#ifdef PORTAGE_ENABLE_MPI
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

  if (ctl.source_file == "none") {
    std::shared_ptr<Wonton::Simple_Mesh> src_mesh = std::make_shared<Wonton::Simple_Mesh>
      (ctl.smin[0], ctl.smin[1], ctl.smin[2], ctl.smax[0], ctl.smax[1], ctl.smax[2],
       ctl.scells[0], ctl.scells[1], ctl.scells[2]);
    std::shared_ptr<Wonton::Simple_Mesh> tgt_mesh = std::make_shared<Wonton::Simple_Mesh>
      (ctl.tmin[0], ctl.tmin[1], ctl.tmin[2], ctl.tmax[0], ctl.tmax[1], ctl.tmax[2],
       ctl.tcells[0], ctl.tcells[1], ctl.tcells[2]);

    runMSM msmguy(ctl, src_mesh, tgt_mesh);

    std::cout << "starting msmapp..." << std::endl;
    std::cout << "running with input file " << filename << std::endl;

    if (ctl.order == 1) {

      msmguy.runit<Portage::IntersectR3D,
		   Portage::Interpolate_1stOrder,
		   Portage::SearchPointsByCells, 3>();

    } else if (ctl.order == 2) {
      msmguy.runit<Portage::IntersectR3D,
      Portage::Interpolate_2ndOrder,
      Portage::SearchPointsByCells, 3>();
    }
  } else {
    Jali::MeshFactory jmf(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> src_mesh = jmf(ctl.source_file);
    std::shared_ptr<Jali::Mesh> tgt_mesh = jmf(ctl.target_file);

    runMSMJali msmguy(ctl, src_mesh, tgt_mesh);

    std::cout << "starting msmapp..." << std::endl;
    std::cout << "running with input file " << filename << std::endl;

    if (ctl.order == 1) {

      msmguy.runit<Portage::IntersectR3D,
		   Portage::Interpolate_1stOrder,
		   Portage::SearchPointsByCells, 3>();

    } else if (ctl.order == 2) {

      //      msmguy.runit<Portage::IntersectR3D,
      //		   Portage::Interpolate_2ndOrder,
      //		   Portage::SearchPointsByCells, 3>();

    }
  }
}
