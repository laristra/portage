/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <vector>

#include "gtest/gtest.h"

#include "tangram/support/MatPoly.h"

#include "portage/intersect/intersect_r2d.h"
#include "portage/support/portage.h"

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

double eps = 1.e-8;

TEST(TANGRAM_2D, test_matpoly_succeeds) {
  // test that we can create a Tangram MatPoly

  Tangram::MatPoly<2> matpoly;
  SUCCEED();
}

TEST(TANGRAM_2D, test_matpoly_create) {
  // test that we can construct a real matpoly (lifted from
  // tangram/src/support/test/test_MatPoly_2D.cc)
  int mat_id = 1;

  // create data  for a unit square
  std::vector<Wonton::Point2> square_points = {
      Wonton::Point2(0.0, 0.0), Wonton::Point2(1.0, 0.0),
      Wonton::Point2(1.0, 1.0), Wonton::Point2(0.0, 1.0)};
  std::vector<std::vector<int>> square_faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  std::vector<Wonton::Point2> face_centroids = {
      Wonton::Point2(0.5, 0.0), Wonton::Point2(1.0, 0.5),
      Wonton::Point2(0.5, 1.0), Wonton::Point2(0.0, 0.5)};

  // Check material ID correctness
  Tangram::MatPoly<2> square_matpoly;
  ASSERT_EQ(-1, square_matpoly.mat_id());
  square_matpoly.set_mat_id(mat_id);
  ASSERT_EQ(mat_id, square_matpoly.mat_id());
  square_matpoly.reset_mat_id();
  ASSERT_EQ(-1, square_matpoly.mat_id());

  // Initialization from ccw ordered vertices
  square_matpoly.initialize(square_points);

  // Verify coordinates
  const std::vector<Wonton::Point2>& matpoly_points = square_matpoly.points();
  ASSERT_EQ(square_points.size(), square_matpoly.num_vertices());
  for (int ivrt = 0; ivrt < square_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(square_points[ivrt], matpoly_points[ivrt], 1.0e-15));

  // Verify faces
  ASSERT_EQ(square_faces.size(), square_matpoly.num_faces());
  for (int iface = 0; iface < square_faces.size(); iface++) {
    const std::vector<int>& face_vertices = square_matpoly.face_vertices(iface);
    ASSERT_EQ(2, face_vertices.size());
    ASSERT_EQ(square_faces[iface][0], face_vertices[0]);
    ASSERT_EQ(square_faces[iface][1], face_vertices[1]);
  }

  // Verify centroids
  for (int iface = 0; iface < square_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         square_matpoly.face_centroid(iface), 1.0e-15));
}

TEST(TANGRAM_2D, test_matpoly_intersect_unit_cells) {
  // test that we can construct a matpoly and intersect with a cell
  // the source and target geometries are both unit cells
  int mat_id = 1;

  double xl = 0., xh = 1., yl = 0., yh = 1.;

  // create data  for a unit square
  std::vector<Wonton::Point<2>> square_points = {
      Wonton::Point<2>(xl, yl), Wonton::Point<2>(xh, yl),
      Wonton::Point<2>(xh, yh), Wonton::Point<2>(xl, yh)};
  std::vector<std::vector<int>> square_faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  // create the matpoly
  Tangram::MatPoly<2> square_matpoly;
  square_matpoly.set_mat_id(mat_id);
  square_matpoly.initialize(square_points);

  // extract the matpoly points
  std::vector<Wonton::Point<2>> _source_points = square_matpoly.points();

  // unfortunately we seem to need to use portage points only in intersection
  // so we need to convert from Tangram points to Portage points
  std::vector<Wonton::Point<2>> source_points;
  for (auto p : _source_points) source_points.push_back(Wonton::Point<2>(p));

  // create a simple mesh with a single cell
  Wonton::Simple_Mesh mesh(xl, yl, xh, yh, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper meshWrapper(mesh);

  // get the coordinates of the single cell
  std::vector<Wonton::Point<2>> target_points;
  meshWrapper.cell_get_coordinates(0, &target_points);

  // use default tolerances
  Portage::NumericTolerances_t num_tols;
  num_tols.use_default();

  // actually intersect
  std::vector<double> moments =
      Portage::intersect_polys_r2d(source_points, target_points, num_tols);

  // test that the moments are correct
  ASSERT_NEAR(moments[0], (xh - xl) * (yh - yl), eps);
  ASSERT_NEAR(moments[1], moments[0] * (xh - xl) / 2., eps);
  ASSERT_NEAR(moments[2], moments[0] * (yh - yl) / 2., eps);
}

TEST(TANGRAM_2D, test_matpoly_intersect_non_coincident) {
  // test that we can construct a matpoly and intersect with a cell
  // the source and target geometries are side 4 squares that intersect
  // at a corner
  int mat_id = 1;

  double xl = 0., xh = 4., yl = 0., yh = 4., xoffset = 2, yoffset = 2;

  // create data  for a unit square
  std::vector<Wonton::Point<2>> square_points = {
      Wonton::Point<2>(xl, yl), Wonton::Point<2>(xh, yl),
      Wonton::Point<2>(xh, yh), Wonton::Point<2>(xl, yh)};
  std::vector<std::vector<int>> square_faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  // create the matpoly
  Tangram::MatPoly<2> square_matpoly;
  square_matpoly.set_mat_id(mat_id);
  square_matpoly.initialize(square_points);

  // extract the matpoly points
  std::vector<Wonton::Point<2>> _source_points = square_matpoly.points();

  // unfortunately we seem to need to use portage points only in intersection
  // so we need to convert from Tangram points to Portage points
  std::vector<Wonton::Point<2>> source_points;
  for (auto p : _source_points) source_points.push_back(Wonton::Point<2>(p));

  // create a simple mesh with a single cell
  Wonton::Simple_Mesh mesh(xl + xoffset, yl + yoffset, xh + xoffset,
                            yh + yoffset, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper meshWrapper(mesh);

  // get the coordinates of the single cell
  std::vector<Wonton::Point<2>> target_points;
  meshWrapper.cell_get_coordinates(0, &target_points);

  // use default tolerances
  Portage::NumericTolerances_t num_tols;
  num_tols.use_default();

  // actually intersect
  std::vector<double> moments =
      Portage::intersect_polys_r2d(source_points, target_points, num_tols);

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 4., eps);
  ASSERT_NEAR(moments[1], 12., eps);
  ASSERT_NEAR(moments[2], 12., eps);
}
