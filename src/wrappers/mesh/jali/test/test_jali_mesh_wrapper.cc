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
#include "FrameworkTraits.hh"

using std::abs;

// Returns true if a == b within the accuracy 'eps'
bool vdd_eq(const std::vector<std::pair<double,double>> &a,
    const std::vector<std::pair<double,double>> &b, double eps=1e-12)
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

// Rotate the 'xylist' vector into a canonical (unique) form. The first point
// will be the one with the lowest angle between it, the nodeid and the
// x-axis.
void coordinates_canonical_rotation(
      const std::pair<double, double> center_node,
      std::vector<std::pair<double, double>> *xylist)
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

// Rotate the 'xylist' vector into a canonical (unique) form. The first point
// will be the one with the lowest angle between it, the nodeid and the
// x-axis.
void dual_cell_coordinates_canonical_rotation(
        const Jali_Mesh_Wrapper &mesh_wrapper,
        int const nodeid,
        std::vector<std::pair<double,double> > *xylist) {
    std::pair<double, double> center_node;
    mesh_wrapper.node_get_coordinates(nodeid, &center_node);
    coordinates_canonical_rotation(center_node, xylist);
}

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

TEST(Jali_Mesh, ccw) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);
    ASSERT_TRUE(mesh != NULL);
    Jali_Mesh_Wrapper mesh_wrapper(*mesh);

    ASSERT_TRUE(mesh_wrapper.ccw({-1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(not mesh_wrapper.ccw({1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(not mesh_wrapper.ccw({-1, 0}, {0, 0}, {1, 0}));
    ASSERT_TRUE(mesh_wrapper.ccw({0.25, 0}, {0.5, 0}, {0.25, 0.25}));
    ASSERT_TRUE(not mesh_wrapper.ccw({0.25, 0}, {0, 0}, {0.25, 0.25}));
    ASSERT_TRUE(mesh_wrapper.ccw({0.25, 0.25}, {0, 0}, {0.25, 0}));
}

TEST(Jali_Mesh, dual_cell_get_coordinates) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2, NULL, true, true, true, true);
    ASSERT_TRUE(mesh != NULL);
    Jali_Mesh_Wrapper mesh_wrapper(*mesh);
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

TEST(Jali_Mesh, mesh_shotshell) {
  int nproc, myrank;
  MPI_Comm_size(MPI_COMM_WORLD,&nproc);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

  ASSERT(Jali::framework_available(Jali::MSTK));

  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);
  //mesh_factory.preferences(Jali::MSTK);

  // Make sure we request faces, edges, wedges and corners
  Jali::Mesh *mesh = mesh_factory("../test_data/shotshell.exo",NULL,true,true,true,true);
  ASSERT_TRUE(mesh != NULL);
}
