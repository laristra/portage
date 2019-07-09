/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <vector>
#include <array>

#include "gtest/gtest.h"

#include "tangram/support/MatPoly.h"

#include "portage/intersect/intersect_r3d.h"

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

double eps = 1.e-8;

TEST(TANGRAM_3D, test_matpoly_succeeds) {
  // test that we can create a Tangram MatPoly

  Tangram::MatPoly<3> matpoly;
  SUCCEED();
}

TEST(TANGRAM_3D, test_matpoly_create) {
  // test that we can construct a real matpoly (lifted from
  // tangram/src/support/test/test_MatPoly_3D.cc)
  int mat_id = 1;

  // Test for a right triangular prism
  std::vector<Wonton::Point<3>> prism_points = {
      Wonton::Point<3>(1.0, 0.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(0.0, 1.0, 1.0),
      Wonton::Point<3>(1.0, 0.0, 1.0), Wonton::Point<3>(0.0, 0.0, 1.0)};
  std::vector<std::vector<int>> prism_faces = {
      {0, 2, 1}, {2, 0, 4, 5}, {3, 4, 0, 1}, {5, 2, 1, 3}, {3, 4, 5}};
  std::vector<Wonton::Point<3>> face_centroids = {
      Wonton::Point<3>(1.0 / 3.0, 1.0 / 3.0, 0.0),
      Wonton::Point<3>(0.5, 0.0, 0.5), Wonton::Point<3>(0.5, 0.5, 0.5),
      Wonton::Point<3>(0.0, 0.5, 0.5),
      Wonton::Point<3>(1.0 / 3.0, 1.0 / 3.0, 1.0)};

  // Check material ID correctness
  Tangram::MatPoly<3> prism_matpoly(mat_id);
  ASSERT_EQ(mat_id, prism_matpoly.mat_id());

  // Initialization
  prism_matpoly.initialize(prism_points, prism_faces);

  // Verify coordinates
  const std::vector<Wonton::Point<3>>& matpoly_points = prism_matpoly.points();
  ASSERT_EQ(prism_points.size(), prism_matpoly.num_vertices());
  for (int ivrt = 0; ivrt < prism_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(prism_points[ivrt], matpoly_points[ivrt], 1.0e-15));

  // Verify faces
  ASSERT_EQ(prism_faces.size(), prism_matpoly.num_faces());
  for (int iface = 0; iface < prism_faces.size(); iface++) {
    const std::vector<int>& face_vertices = prism_matpoly.face_vertices(iface);
    ASSERT_EQ(prism_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < prism_faces[iface].size(); ivrt++)
      ASSERT_EQ(prism_faces[iface][ivrt], face_vertices[ivrt]);
  }

  // Verify centroids
  for (int iface = 0; iface < prism_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         prism_matpoly.face_centroid(iface), 1.0e-15));

  std::vector<Wonton::Point<3>> faceted_prism_points = {
      Wonton::Point<3>(1.0, 0.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(0.0, 1.0, 1.0),
      Wonton::Point<3>(1.0, 0.0, 1.0), Wonton::Point<3>(0.0, 0.0, 1.0),
      Wonton::Point<3>(0.5, 0.0, 0.5), Wonton::Point<3>(0.5, 0.5, 0.5),
      Wonton::Point<3>(0.0, 0.5, 0.5),
  };
  std::vector<std::vector<int>> faceted_prism_faces = {
      {0, 2, 1}, {6, 2, 0}, {6, 0, 4}, {6, 4, 5}, {6, 5, 2},
      {7, 3, 4}, {7, 4, 0}, {7, 0, 1}, {7, 1, 3}, {8, 5, 2},
      {8, 2, 1}, {8, 1, 3}, {8, 3, 5}, {3, 4, 5}};

  // Create faceted poly
  Tangram::MatPoly<3> faceted_prism_matpoly;
  prism_matpoly.faceted_matpoly(&faceted_prism_matpoly);

  // Check material ID correctness
  ASSERT_EQ(mat_id, faceted_prism_matpoly.mat_id());

  // Verify facetization
  // Verify node coordinates
  const std::vector<Wonton::Point<3>>& faceted_matpoly_points =
      faceted_prism_matpoly.points();
  ASSERT_EQ(faceted_prism_points.size(), faceted_prism_matpoly.num_vertices());
  for (int ivrt = 0; ivrt < faceted_prism_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(faceted_prism_points[ivrt],
                         faceted_matpoly_points[ivrt], 1.0e-15));

  // Verify facets
  ASSERT_EQ(faceted_prism_faces.size(), faceted_prism_matpoly.num_faces());
  for (int iface = 0; iface < faceted_prism_faces.size(); iface++) {
    const std::vector<int>& face_vertices =
        faceted_prism_matpoly.face_vertices(iface);
    ASSERT_EQ(faceted_prism_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < faceted_prism_faces[iface].size(); ivrt++)
      ASSERT_EQ(faceted_prism_faces[iface][ivrt], face_vertices[ivrt]);
  }

  // Non-convex distorted prism
  std::vector<Wonton::Point<3>> ncv_prism_points = {
      Wonton::Point<3>(1.0, 0.0, 0.0), Wonton::Point<3>(0.4, 0.8, 0.2),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.5, 0.1, 0.1), Wonton::Point<3>(0.0, 0.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 1.0, 1.0),
      Wonton::Point<3>(0.5, 0.9, 1.1), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(0.0, 1.0, 1.0), Wonton::Point<3>(0.6, 0.2, 0.8)};
  std::vector<std::vector<int>> ncv_prism_faces = {
      {3, 1, 2, 0, 4, 5}, {5, 4, 11, 10}, {0, 9, 11, 4}, {7, 9, 0, 2},
      {2, 1, 8, 7},       {8, 1, 3, 10},  {5, 6, 10, 3}, {7, 11, 10, 6, 8, 9}};
  std::vector<Wonton::Point<3>> ncv_prism_face_centroids = {
      Wonton::Point<3>(0.49982029799691324, 0.48793283780407681,
                        0.038845671672909921),
      Wonton::Point<3>(0.2519019047393049, 0.36565382681732062,
                        0.51522129026172814),
      Wonton::Point<3>(0.78729807432173626, 0.070061777159007799,
                        0.46574083274162237),
      Wonton::Point<3>(1.0, 0.5, 0.5),
      Wonton::Point<3>(0.72714164844017881, 0.92490808257951729,
                        0.55650182505587276),
      Wonton::Point<3>(0.22120243710721832, 0.92700420108481663,
                        0.58407100072352491),
      Wonton::Point<3>(0.0, 0.5, 0.5),
      Wonton::Point<3>(0.50915092483338054, 0.51261128330217054,
                        0.98738871669782957)};

  Tangram::MatPoly<3> ncv_prism_matpoly(mat_id);
  // Initialization
  ncv_prism_matpoly.initialize(ncv_prism_points, ncv_prism_faces);

  // Verify centroids
  for (int iface = 0; iface < ncv_prism_faces.size(); iface++)
    ASSERT_TRUE(approxEq(ncv_prism_face_centroids[iface],
                         ncv_prism_matpoly.face_centroid(iface), 1.0e-15));
}

TEST(TANGRAM_3D, test_matpoly_cube) {
  // test that we can construct a cube matpoly
  int mat_id = 1;

  // Test for a cube
  std::vector<Wonton::Point<3>> cube_points = {
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(1.0, 0.0, 0.0),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(1.0, 1.0, 1.0), Wonton::Point<3>(0.0, 1.0, 1.0)};
  std::vector<std::vector<int>> cube_faces = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                              {2, 3, 7, 6}, {3, 0, 4, 7},
                                              {0, 3, 2, 1}, {4, 5, 6, 7}};
  std::vector<Wonton::Point<3>> face_centroids = {
      Wonton::Point<3>(0.5, 0.0, 0.5), Wonton::Point<3>(1.0, 0.5, 0.5),
      Wonton::Point<3>(0.5, 1.0, 0.5), Wonton::Point<3>(0.0, 0.5, 0.5),
      Wonton::Point<3>(0.5, 0.5, 0.0), Wonton::Point<3>(0.5, 0.5, 1.0)};

  // Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  ASSERT_EQ(mat_id, matpoly.mat_id());

  // Initialization
  matpoly.initialize(cube_points, cube_faces);

  // Verify coordinates
  const std::vector<Wonton::Point<3>>& matpoly_points = matpoly.points();
  ASSERT_EQ(cube_points.size(), matpoly.num_vertices());
  for (int ivrt = 0; ivrt < cube_points.size(); ivrt++)
    for (int j = 0; j < 3; j++)
      ASSERT_NEAR(cube_points[ivrt][j], matpoly_points[ivrt][j], 1.0e-15);

  // Verify faces
  ASSERT_EQ(cube_faces.size(), matpoly.num_faces());
  for (int iface = 0; iface < cube_faces.size(); iface++) {
    const std::vector<int>& face_vertices = matpoly.face_vertices(iface);
    ASSERT_EQ(cube_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < cube_faces[iface].size(); ivrt++)
      ASSERT_EQ(cube_faces[iface][ivrt], face_vertices[ivrt]);
  }

  // Verify centroids
  for (int iface = 0; iface < cube_faces.size(); iface++)
    ASSERT_TRUE(
        approxEq(face_centroids[iface], matpoly.face_centroid(iface), 1.0e-15));
}

TEST(TANGRAM_3D, test_matpoly_faceted_cube_by_hand) {
  // test that we can construct a faceted cube matpoly by hand
  int mat_id = 1;

  // Test for a cube
  std::vector<Wonton::Point<3>> faceted_cube_points = {
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(1.0, 0.0, 0.0),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(1.0, 1.0, 1.0), Wonton::Point<3>(0.0, 1.0, 1.0),
      Wonton::Point<3>(0.5, 0.0, 0.5), Wonton::Point<3>(1.0, 0.5, 0.5),
      Wonton::Point<3>(0.5, 1.0, 0.5), Wonton::Point<3>(0.0, 0.5, 0.5),
      Wonton::Point<3>(0.5, 0.5, 0.0), Wonton::Point<3>(0.5, 0.5, 1.0)};
  std::vector<std::vector<int>> faceted_cube_faces = {{0, 1, 8}, {1, 5, 8},
      {5, 4, 8}, {4, 0, 8}, {1, 2, 9}, {2, 6, 9}, {6, 5, 9}, {5, 1, 9},
      {2, 3, 10}, {3, 7, 10}, {7, 6, 10}, {6, 2, 10}, {3, 0, 11}, {0, 4, 11},
      {4, 7, 11}, {7, 3, 11}, {0, 3, 12}, {3, 2, 12}, {2, 1, 12}, {1, 0, 12},
      {4, 5, 13}, {5, 6, 13}, {6, 7, 13}, {7, 4, 12}};

  // Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  ASSERT_EQ(mat_id, matpoly.mat_id());

  // Initialization
  matpoly.initialize(faceted_cube_points, faceted_cube_faces);

  // Verify coordinates
  const std::vector<Wonton::Point<3>>& matpoly_points = matpoly.points();
  ASSERT_EQ(faceted_cube_points.size(), matpoly.num_vertices());
  for (int ivrt = 0; ivrt < faceted_cube_points.size(); ivrt++)
    for (int j = 0; j < 3; j++)
      ASSERT_NEAR(faceted_cube_points[ivrt][j], matpoly_points[ivrt][j],
                  1.0e-15);

  // Verify faces
  ASSERT_EQ(faceted_cube_faces.size(), matpoly.num_faces());
  for (int iface = 0; iface < faceted_cube_faces.size(); iface++) {
    const std::vector<int>& face_vertices = matpoly.face_vertices(iface);
    ASSERT_EQ(faceted_cube_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < faceted_cube_faces[iface].size(); ivrt++)
      ASSERT_EQ(faceted_cube_faces[iface][ivrt], face_vertices[ivrt]);
  }
}

TEST(TANGRAM_3D, test_matpoly_intersect) {
  // test that we can construct a cube matpoly using .faceted_matpoly
  // and intersect with a portag
  int mat_id = 1;

  // Test for a cube
  std::vector<Wonton::Point<3>> cube_points = {
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(1.0, 0.0, 0.0),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(1.0, 1.0, 1.0), Wonton::Point<3>(0.0, 1.0, 1.0)};
  std::vector<std::vector<int>> cube_faces = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                              {2, 3, 7, 6}, {3, 0, 4, 7},
                                              {0, 3, 2, 1}, {4, 5, 6, 7}};
  std::vector<Wonton::Point<3>> face_centroids = {
      Wonton::Point<3>(0.5, 0.0, 0.5), Wonton::Point<3>(1.0, 0.5, 0.5),
      Wonton::Point<3>(0.5, 1.0, 0.5), Wonton::Point<3>(0.0, 0.5, 0.5),
      Wonton::Point<3>(0.5, 0.5, 0.0), Wonton::Point<3>(0.5, 0.5, 1.0)};

  // Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  ASSERT_EQ(mat_id, matpoly.mat_id());

  // Initialization
  matpoly.initialize(cube_points, cube_faces);

  // Verify coordinates
  const std::vector<Wonton::Point<3>>& matpoly_points = matpoly.points();
  ASSERT_EQ(cube_points.size(), matpoly.num_vertices());
  for (int ivrt = 0; ivrt < cube_points.size(); ivrt++)
    for (int j = 0; j < 3; j++)
      ASSERT_NEAR(cube_points[ivrt][j], matpoly_points[ivrt][j], 1.0e-15);

  // Verify faces
  ASSERT_EQ(cube_faces.size(), matpoly.num_faces());
  for (int iface = 0; iface < cube_faces.size(); iface++) {
    const std::vector<int>& face_vertices = matpoly.face_vertices(iface);
    ASSERT_EQ(cube_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < cube_faces[iface].size(); ivrt++)
      ASSERT_EQ(cube_faces[iface][ivrt], face_vertices[ivrt]);
  }

  // Verify centroids
  for (int iface = 0; iface < cube_faces.size(); iface++)
    ASSERT_TRUE(
        approxEq(face_centroids[iface], matpoly.face_centroid(iface), 1.0e-15));

  // Facet the cube
  Tangram::MatPoly<3> faceted_matpoly;
  matpoly.faceted_matpoly(&faceted_matpoly);

  // Check that the faceted matpoly has 24 faces and 14 vertices
  ASSERT_EQ(faceted_matpoly.num_faces(), 24);
  ASSERT_EQ(faceted_matpoly.num_vertices(), 14);

  // Convert the points
  std::vector<Wonton::Point<3>> _source_points = faceted_matpoly.points();
  std::vector<Wonton::Point<3>> source_points;
  for (auto p : _source_points) source_points.push_back(Wonton::Point<3>(p));

  // Convert the face indices from the matpoly to a Portage style
  std::vector<std::vector<int>> facetpoints;
  for (int iface = 0; iface < faceted_matpoly.num_faces(); ++iface) {
    const std::vector<int>& face_vertices =
        faceted_matpoly.face_vertices(iface);
    // for (int f: face_vertices) std::cout << f << ", ";
    // std::cout <<std::endl;
    facetpoints.push_back(face_vertices);
  }

  // Create the facetedpoly_t structure
  Portage::facetedpoly_t srcpoly{facetpoints, source_points};

  // Create a length 1 vector of tets by hand
  std::vector<std::array<Wonton::Point<3>, 4>> target_tet_coords{
      {Wonton::Point<3>{0., 0., 0.}, Wonton::Point<3>{1., 0., 0.},
       Wonton::Point<3>{0., 1., 0.}, Wonton::Point<3>{0., 0., 1.}}};

  // use default tolerances
  Portage::NumericTolerances_t num_tols;

  // Hope for a miracle with intersection
  std::vector<double> moments(
      Portage::intersect_polys_r3d(srcpoly, target_tet_coords, num_tols));

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 1. / 6., eps);
  ASSERT_NEAR(moments[1] / moments[0], .25, eps);
  ASSERT_NEAR(moments[2] / moments[0], .25, eps);
  ASSERT_NEAR(moments[3] / moments[0], .25, eps);
}

TEST(TANGRAM_3D, test_intersect_matpoly_gold) {
  // test that we can construct a cube matpoly and intersect with a constructed
  // tet
  int mat_id = 1;

  // Test for a cube
  std::vector<Wonton::Point<3>> cube_points = {
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(1.0, 0.0, 0.0),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(1.0, 1.0, 1.0), Wonton::Point<3>(0.0, 1.0, 1.0)};
  std::vector<std::vector<int>> cube_faces = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                              {2, 3, 7, 6}, {3, 0, 4, 7},
                                              {0, 3, 2, 1}, {4, 5, 6, 7}};

  // Create the MatPoly
  Tangram::MatPoly<3> matpoly(mat_id);

  // Initialization
  matpoly.initialize(cube_points, cube_faces);

  // Create a length 1 vector of tets by hand
  std::vector<std::array<Wonton::Point<3>, 4>> target_tet_coords{
      {Wonton::Point<3>{0., 0., 0.}, Wonton::Point<3>{1., 0., 0.},
       Wonton::Point<3>{0., 1., 0.}, Wonton::Point<3>{0., 0., 1.}}};

  // use default tolerances
  Portage::NumericTolerances_t num_tols;

  // Intersect
  std::vector<double> moments(Portage::intersect_polys_r3d(
      Portage::get_faceted_matpoly(matpoly), target_tet_coords, num_tols));

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 1. / 6., eps);
  ASSERT_NEAR(moments[1] / moments[0], .25, eps);
  ASSERT_NEAR(moments[2] / moments[0], .25, eps);
  ASSERT_NEAR(moments[3] / moments[0], .25, eps);
}

TEST(TANGRAM_3D, test_intersect_matpoly_gold2) {
  // test that we can construct a cube matpoly and intersect with a cell
  int mat_id = 1;

  // Test for a cube
  std::vector<Wonton::Point<3>> cube_points = {
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(1.0, 0.0, 0.0),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(1.0, 1.0, 1.0), Wonton::Point<3>(0.0, 1.0, 1.0)};
  std::vector<std::vector<int>> cube_faces = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                              {2, 3, 7, 6}, {3, 0, 4, 7},
                                              {0, 3, 2, 1}, {4, 5, 6, 7}};

  // Create the MatPoly
  Tangram::MatPoly<3> matpoly(mat_id);

  // Initialization
  matpoly.initialize(cube_points, cube_faces);

  // create a simple mesh with a single cell
  Wonton::Simple_Mesh target_mesh(0., 0., 0., 1., 1., 1., 1, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(target_mesh);

  // decompose the cell into tets
  std::vector<std::array<Wonton::Point<3>, 4>> tcoords;
  bool planar_hex = true;
  target_mesh_wrapper.decompose_cell_into_tets(0, &tcoords, planar_hex);

  // use default tolerances
  Portage::NumericTolerances_t num_tols;

  // Intersect
  std::vector<double> moments(Portage::intersect_polys_r3d(
      Portage::get_faceted_matpoly(matpoly), tcoords, num_tols));

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 1., eps);
  ASSERT_NEAR(moments[1] / moments[0], .5, eps);
  ASSERT_NEAR(moments[2] / moments[0], .5, eps);
  ASSERT_NEAR(moments[3] / moments[0], .5, eps);
}

TEST(TANGRAM_3D, test_intersect_matpoly_gold3) {
  // test that we can construct a cube matpoly and intersect with a cell using
  // the non-planar flag
  int mat_id = 1;

  // Test for a cube
  std::vector<Wonton::Point<3>> cube_points = {
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(1.0, 0.0, 0.0),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(1.0, 1.0, 1.0), Wonton::Point<3>(0.0, 1.0, 1.0)};
  std::vector<std::vector<int>> cube_faces = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                              {2, 3, 7, 6}, {3, 0, 4, 7},
                                              {0, 3, 2, 1}, {4, 5, 6, 7}};

  // Create the MatPoly
  Tangram::MatPoly<3> matpoly(mat_id);

  // Initialization
  matpoly.initialize(cube_points, cube_faces);

  // create a simple mesh with a single cell
  Wonton::Simple_Mesh target_mesh(0., 0., 0., 1., 1., 1., 1, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(target_mesh);

  // decompose the cell into tets
  std::vector<std::array<Wonton::Point<3>, 4>> tcoords;
  bool planar_hex = false;
  target_mesh_wrapper.decompose_cell_into_tets(0, &tcoords, planar_hex);

  // use default tolerances
  Portage::NumericTolerances_t num_tols;

  // Intersect
  std::vector<double> moments(Portage::intersect_polys_r3d(
      Portage::get_faceted_matpoly(matpoly), tcoords, num_tols));

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 1., eps);
  ASSERT_NEAR(moments[1] / moments[0], .5, eps);
  ASSERT_NEAR(moments[2] / moments[0], .5, eps);
  ASSERT_NEAR(moments[3] / moments[0], .5, eps);
}

TEST(TANGRAM_3D, test_intersect_matpoly_gold4) {
  // test that we can intersect a cube matpoly and intersect with a
  // non-coincident cell
  int mat_id = 1;

  // Test for a cube
  std::vector<Wonton::Point<3>> cube_points = {
      Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(1.0, 0.0, 0.0),
      Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
      Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 0.0, 1.0),
      Wonton::Point<3>(1.0, 1.0, 1.0), Wonton::Point<3>(0.0, 1.0, 1.0)};
  std::vector<std::vector<int>> cube_faces = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                              {2, 3, 7, 6}, {3, 0, 4, 7},
                                              {0, 3, 2, 1}, {4, 5, 6, 7}};

  // Create the MatPoly
  Tangram::MatPoly<3> matpoly(mat_id);

  // Initialization
  matpoly.initialize(cube_points, cube_faces);

  // create a simple mesh with a single cell
  Wonton::Simple_Mesh target_mesh(0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(target_mesh);

  // decompose the cell into tets
  std::vector<std::array<Wonton::Point<3>, 4>> tcoords;
  bool planar_hex = false;
  target_mesh_wrapper.decompose_cell_into_tets(0, &tcoords, planar_hex);

  // use default tolerances
  Portage::NumericTolerances_t num_tols;

  // Intersect
  std::vector<double> moments(Portage::intersect_polys_r3d(
      Portage::get_faceted_matpoly(matpoly), tcoords, num_tols));

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 1. / 8., eps);
  ASSERT_NEAR(moments[1] / moments[0], .75, eps);
  ASSERT_NEAR(moments[2] / moments[0], .75, eps);
  ASSERT_NEAR(moments[3] / moments[0], .75, eps);
}
