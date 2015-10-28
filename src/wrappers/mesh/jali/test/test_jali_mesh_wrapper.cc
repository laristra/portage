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

TEST(Jali_Mesh, ccw) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2);
    ASSERT_TRUE(mesh != NULL);
    Jali_Mesh_Wrapper mesh_wrapper(*mesh);

    ASSERT_TRUE(mesh_wrapper.ccw({-1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(not mesh_wrapper.ccw({1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(not mesh_wrapper.ccw({-1, 0}, {0, 0}, {1, 0}));
}

TEST(Jali_Mesh, dual_cell_get_coordinates) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    Jali::Mesh *mesh = mf(0.0,0.0,1.0,1.0,2,2, NULL, true, true, true, true);
    ASSERT_TRUE(mesh != NULL);
    Jali_Mesh_Wrapper mesh_wrapper(*mesh);
    double eps = 1e-12;

    std::vector<std::pair<double,double>> xylist;
    mesh_wrapper.dual_cell_get_coordinates(0, &xylist);
    ASSERT_TRUE(xylist.size() == 2);
    ASSERT_TRUE(abs(std::get<0>(xylist[0]) - 0.25) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[0]) - 0.25) < eps);
    ASSERT_TRUE(abs(std::get<0>(xylist[1]) - 0.25) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[1]) - 0   ) < eps);
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(1, &xylist);
    ASSERT_TRUE(xylist.size() == 4);
    ASSERT_TRUE(abs(std::get<0>(xylist[0]) - 0.25) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[0]) - 0   ) < eps);
    ASSERT_TRUE(abs(std::get<0>(xylist[1]) - 0.25) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[1]) - 0.25) < eps);
    ASSERT_TRUE(abs(std::get<0>(xylist[2]) - 0.75) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[2]) - 0.25) < eps);
    ASSERT_TRUE(abs(std::get<0>(xylist[3]) - 0.75) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[3]) - 0   ) < eps);
    xylist.clear();

    mesh_wrapper.dual_cell_get_coordinates(2, &xylist);
    ASSERT_TRUE(xylist.size() == 2);
    ASSERT_TRUE(abs(std::get<0>(xylist[0]) - 0.75) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[0]) - 0   ) < eps);
    ASSERT_TRUE(abs(std::get<0>(xylist[1]) - 0.75) < eps);
    ASSERT_TRUE(abs(std::get<1>(xylist[1]) - 0.25) < eps);
    xylist.clear();

    /*
    // Uncomment this code to print the xylist
    std::cout << "xylist:" << std::endl;
    std::cout << xylist.size() << std::endl;
    for (auto &v: xylist) {
        std::cout << v.first << " " << v.second << std::endl;
    }
    */
}
