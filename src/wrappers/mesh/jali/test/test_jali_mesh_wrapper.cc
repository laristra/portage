/*---------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *---------------------------------------------------------------------------~*/

#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

using std::abs;

/*!
  @file test_jali_mesh_wrapper.cc
  @brief Unit tests for the Jali mesh wrapper class
 */

/*!
  @brief Returns true if a == b within the accuracy 'eps'
  @param[in] a First input vector of pairs of doubles
  @param[in] b Second input vector of pairs of doubles
  @param[in] eps Tolerance
  @return Whether a == b within the accuracy 'eps'
 */
bool vdd_eq(const std::vector<std::pair<double,double>> &a,
    const std::vector<std::pair<double,double>> &b, const double eps=1e-12)
{
    // Can't be equal if # of entries differ:
    if (a.size() != b.size()) return false;
    // Loop over elements in "a" and "b":
    for (size_t i = 0; i < a.size(); i++) {
        // values not equal
        if (abs(std::get<0>(a[i]) - std::get<0>(b[i])) > eps or
            abs(std::get<1>(a[i]) - std::get<1>(b[i])) > eps) return false;
    }
    return true;
}

/*!
  @brief Rotate the 'xylist' vector into a canonical (unique) form. The first point
  will be the one with the lowest angle between it, the nodeid and the x-axis.
  @param[in] center_node The center node to test
  @param[in,out] xylist      The xylist vector
 */
void coordinates_canonical_rotation(
      const std::pair<double, double> center_node,
      std::vector<std::pair<double, double>> * const xylist)
{
    int i = 0;
    auto angle = [&]() {
        return std::atan2(std::get<1>((*xylist)[i])-std::get<1>(center_node),
                std::get<0>((*xylist)[i])-std::get<0>(center_node));
    };
    double a = angle();
    while (a >= 0) { i++; i = i % xylist->size(); a = angle(); }
    while (a < 0) { i++; i = i % xylist->size(); a = angle(); }
    std::rotate(xylist->begin(), xylist->begin()+i, xylist->end());
}

/*!
  @brief Rotate the 'xylist' vector into a canonical (unique) form. The first point
  will be the one with the lowest angle between it, the nodeid and the x-axis.
  @param[in] mesh_wrapper The mesh wrapper
  @param[in] nodeid The node id
  @param[in,out] xylist The xylist
 */
void dual_cell_coordinates_canonical_rotation(
    const Portage::Jali_Mesh_Wrapper &mesh_wrapper,
    int const nodeid,
    std::vector<std::pair<double,double> > * const xylist) {
  std::pair<double, double> center_node;
  mesh_wrapper.node_get_coordinates(nodeid, &center_node);
  coordinates_canonical_rotation(center_node, xylist);
}

/*!
  @brief Unit test for equality comparisons
 */
TEST(Jali_Mesh, vdd_eq) {
    std::vector<std::pair<double,double>> a, b, c;
    a = {{0.25, 0}, {0.25, 0.25}};
    b = {{0.25, 0}, {0.25, 0.25}};
    c = {{0.25, 0.25}, {0.25, 0}};
    ASSERT_TRUE(vdd_eq(a, b));
    ASSERT_TRUE(vdd_eq(a, {{0.25, 0}, {0.25, 0.25}}));
    ASSERT_TRUE(not vdd_eq(a, {{0.25, 0}}));
    ASSERT_TRUE(not vdd_eq(a, c));
    ASSERT_TRUE(not vdd_eq(a, {{0.25, 0.25}, {0.25, 0}}));
    ASSERT_TRUE(vdd_eq(a, c, 0.5));
}

/*!
  @brief Unit test for canonical rotations of coordinates
 */
TEST(Jali_Mesh, coordinates_canonical_rotation) {
    std::vector<std::pair<double,double>> xylist, xylist_canonical;
    xylist_canonical = {
                {0.75, 0.5},
                {0.75, 0.75},
                {0.5, 0.75},
                {0.25, 0.75},
                {0.25, 0.5},
                {0.25, 0.25},
                {0.5, 0.25},
                {0.75, 0.25},
                };

    xylist = { {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));
}

/*!
  @brief Unit test for counter-clockwise winding of coordinates
 */
TEST(Jali_Mesh, ccw) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);
    ASSERT_TRUE(mesh != NULL);
    Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

    ASSERT_TRUE(mesh_wrapper.ccw({-1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(not mesh_wrapper.ccw({1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(not mesh_wrapper.ccw({-1, 0}, {0, 0}, {1, 0}));
    ASSERT_TRUE(mesh_wrapper.ccw({0.25, 0}, {0.5, 0}, {0.25, 0.25}));
    ASSERT_TRUE(not mesh_wrapper.ccw({0.25, 0}, {0, 0}, {0.25, 0.25}));
    ASSERT_TRUE(mesh_wrapper.ccw({0.25, 0.25}, {0, 0}, {0.25, 0}));
}

/*!
  @brief Unit test for getting dual cell coordinates
 */
TEST(Jali_Mesh, dual_cell_get_coordinates) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2, NULL, true, true, true, true);
    ASSERT_TRUE(mesh != NULL);
    Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);
    double eps = 1e-12;

    std::vector<std::pair<double,double>> xylist;
    mesh_wrapper.dual_cell_get_coordinates(0, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {0, 0},
                {0.25, 0},
                {0.25, 0.25},
                {0, 0.25},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(1, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {0.5, 0},
                {0.75, 0},
                {0.75, 0.25},
                {0.5, 0.25},
                {0.25, 0.25},
                {0.25, 0},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(2, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {1, 0},
                {1, 0.25},
                {0.75, 0.25},
                {0.75, 0},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(3, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {0, 0.5},
                {0, 0.25},
                {0.25, 0.25},
                {0.25, 0.5},
                {0.25, 0.75},
                {0, 0.75},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(4, &xylist);
    dual_cell_coordinates_canonical_rotation(mesh_wrapper, 4, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {0.75, 0.5},
                {0.75, 0.75},
                {0.5, 0.75},
                {0.25, 0.75},
                {0.25, 0.5},
                {0.25, 0.25},
                {0.5, 0.25},
                {0.75, 0.25},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(5, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {1, 0.5},
                {1, 0.75},
                {0.75, 0.75},
                {0.75, 0.5},
                {0.75, 0.25},
                {1, 0.25},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(6, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {0, 1},
                {0, 0.75},
                {0.25, 0.75},
                {0.25, 1},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(7, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {0.5, 1},
                {0.25, 1},
                {0.25, 0.75},
                {0.5, 0.75},
                {0.75, 0.75},
                {0.75, 1},
                }));
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(8, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, {
                {1, 1},
                {0.75, 1},
                {0.75, 0.75},
                {1, 0.75},
                }));
    xylist.clear();

    /*
    // Uncomment this code to print the xylist
    std::cout << "xylist:" << std::endl;
    std::cout << xylist.size() << std::endl;
    for (auto &v: xylist) {
        std::cout << "{" << v.first << ", " << v.second << "}," << std::endl;
    }
    */
}

/*!
  @brief Unit test for getting neighbor cells
 */
TEST(Jali_Mesh, Get_Neighbor_Cells) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  Jali::Mesh *mesh = mf(0.0,0.0,0.0,1.0,1.0,1.0,2,2,2);
  ASSERT_TRUE(mesh != NULL);
  Portage::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  // This is a regular mesh with 2 cells (3 nodes) in each direction 

  // If we ask for the node adjacent neighbors of any cell, we should
  // get back all the other cells in the mesh

  int anycell = 3;

  // 1 (2nd argument) means retrieve only owned cells

  std::vector<int> adjcellids;
  mesh_wrapper.cell_get_node_adj_cells(anycell, Portage::ALL, &adjcellids);

  int ncells = mesh_wrapper.num_owned_cells();
  EXPECT_EQ(ncells-1,adjcellids.size()) <<
      "Not the right number of adjacent cells" << std::endl;

  for (int i = 0; i < ncells; i++) {
    if (i == anycell) continue;

    ASSERT_TRUE(std::find(adjcellids.begin(),adjcellids.end(),i) != 
                adjcellids.end()) << "Cell " << i << 
        " not found in adjacent cell list" << std::endl;
  }

  // The center node of this 3 node x 3 node x 3 node mesh corresponds
  // to a dual cell that is surrounded completely by the other dual
  // cells of the mesh. If we ask for the dual cell neighbors of the central
  // node, we should get back all the other dual cells in the mesh
  
  int center_node = 13;
  std::tuple<double,double,double> cxyz;
  mesh_wrapper.node_get_coordinates(center_node,&cxyz);

  ASSERT_EQ(0.5,std::get<0>(cxyz));
  ASSERT_EQ(0.5,std::get<1>(cxyz));
  ASSERT_EQ(0.5,std::get<2>(cxyz));

  std::vector<int> adjdualcellids; 

  // 1 (2nd argument) means owned dual cells
  mesh_wrapper.dual_cell_get_node_adj_cells(center_node, Portage::ALL, &adjdualcellids); 

  // List of adjacent dual cell ids should contain all node ids except
  // the center one

  int nnodes = mesh_wrapper.num_owned_nodes();

  EXPECT_EQ(nnodes-1,adjdualcellids.size()) << 
      "Not the right number of adjacent dual cells" << std::endl;

  for (int i = 0; i < nnodes; i++) {
    if (i == center_node) continue;

    ASSERT_TRUE(std::find(adjdualcellids.begin(),adjdualcellids.end(),i) != 
                adjdualcellids.end()) << "Dual cell " << i << 
        " not found in adjacent dual cell list" << std::endl;
  }
}

/*!
  @brief Unit test for creating a mesh from the shotshell exo file
 */
TEST(Jali_Mesh, mesh_shotshell) {
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);

  // Make sure we request faces, edges, wedges and corners
  Jali::Mesh *mesh = mesh_factory("test_data/shotshell.exo",NULL,true,true,true,true);
  ASSERT_TRUE(mesh != NULL);

  // Make sure we request faces, edges, wedges and corners
  mesh = mesh_factory("test_data/shotshell-v.exo",NULL,true,true,true,true);
  ASSERT_TRUE(mesh != NULL);
}
