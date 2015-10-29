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
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, };
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, };
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, };
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, };
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, };
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, };
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, };
    Jali_Mesh_Wrapper::coordinates_canonical_rotation({0.5, 0.5}, &xylist);
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
    mesh_wrapper.dual_cell_coordinates_canonical_rotation(4, &xylist);
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
