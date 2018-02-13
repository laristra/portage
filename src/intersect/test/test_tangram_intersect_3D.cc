/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifdef HAVE_TANGRAM

#include <iostream>
#include "gtest/gtest.h"

#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"
//#include "tangram/simple_mesh/simple_mesh.h"

#include "portage/intersect/intersect_r3d.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "portage/support/Point.h"
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"

double eps = 1.e-8;
TEST(TANGRAM, test_tangram_exists) {
  // test that we can instantiate a Tangram Point even though we don't include
  // the Tangram header

  float x{1.}, y{2}, z{3};
  Tangram::Point<3> p{x, y, z};
  ASSERT_FLOAT_EQ(x, p[0]);
  ASSERT_FLOAT_EQ(y, p[1]);
  ASSERT_FLOAT_EQ(z, p[2]);
}

TEST(TANGRAM, test_portage_exists) {
  // test that we can instantiate a Portage Point

  float x{1.}, y{2}, z{3};
  Portage::Point<3> p{x, y, z};
  ASSERT_FLOAT_EQ(x, p[0]);
  ASSERT_FLOAT_EQ(y, p[1]);
  ASSERT_FLOAT_EQ(z, p[2]);
}

TEST(TANGRAM, test_tangram_to_portage) {
  // test that we can create a Portage Point from a Tangram Point

  float x{1.}, y{2}, z{3};
  Tangram::Point<3> pt{x, y,z};
  Portage::Point<3> pp{pt};
  ASSERT_FLOAT_EQ(x, pp[0]);
  ASSERT_FLOAT_EQ(y, pp[1]);
  ASSERT_FLOAT_EQ(z, pp[2]);
}

TEST(TANGRAM, test_portage_to_tangram) {
  // test that we can create a Tangram Point from a Portage Point

  float x{1.}, y{2}, z{3};
  Portage::Point<3> pp{x, y, z};
  Tangram::Point<3> pt{pp};
  ASSERT_FLOAT_EQ(x, pt[0]);
  ASSERT_FLOAT_EQ(y, pt[1]);
  ASSERT_FLOAT_EQ(z, pt[2]);
}

TEST(TANGRAM, test_matpoly_succeeds) {
  // test that we can create a Tangram MatPoly

  Tangram::MatPoly<3> matpoly;
  SUCCEED();
}

/*
TEST(TANGRAM, test_matpoly_create) {
  // test that we can construct a real matpoly (lifted from
  // tangram/src/support/test/test_MatPoly_2D.cc)
  int mat_id = 1;

  // create data  for a unit square
  std::vector<Tangram::Point2> square_points = {
      Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
      Tangram::Point2(1.0, 1.0), Tangram::Point2(0.0, 1.0)};
  std::vector<std::vector<int>> square_faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
  std::vector<Tangram::Point2> face_centroids = {
      Tangram::Point2(0.5, 0.0), Tangram::Point2(1.0, 0.5),
      Tangram::Point2(0.5, 1.0), Tangram::Point2(0.0, 0.5)};

  // Check material ID correctness
  Tangram::MatPoly<3> square_matpoly;
  ASSERT_EQ(-1, square_matpoly.mat_id());
  square_matpoly.set_mat_id(mat_id);
  ASSERT_EQ(mat_id, square_matpoly.mat_id());
  square_matpoly.reset_mat_id();
  ASSERT_EQ(-1, square_matpoly.mat_id());

  // Initialization from ccw ordered vertices
  square_matpoly.initialize(square_points);

  // Verify coordinates
  const std::vector<Tangram::Point2>& matpoly_points = square_matpoly.points();
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

TEST(TANGRAM, test_matpoly_extract_points) {
  // test that we can construct a matpoly and extract its points
  int mat_id = 1;

  // create data  for a unit square
  std::vector<Tangram::Point2> square_points = {
      Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
      Tangram::Point2(1.0, 1.0), Tangram::Point2(0.0, 1.0)};
  std::vector<std::vector<int>> square_faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  // create the matpoly
  Tangram::MatPoly<3> square_matpoly;
  square_matpoly.set_mat_id(mat_id);
  square_matpoly.initialize(square_points);

  // extract the matpoly points
  std::vector<Tangram::Point<3>> points = square_matpoly.points();

  // test they have the correct values
  ASSERT_FLOAT_EQ(0., points[0][0]);
  ASSERT_FLOAT_EQ(0., points[0][1]);

  ASSERT_FLOAT_EQ(1., points[1][0]);
  ASSERT_FLOAT_EQ(0., points[1][1]);

  ASSERT_FLOAT_EQ(1., points[2][0]);
  ASSERT_FLOAT_EQ(1., points[2][1]);

  ASSERT_FLOAT_EQ(0., points[3][0]);
  ASSERT_FLOAT_EQ(1., points[3][1]);
}

TEST(TANGRAM, test_matpoly_intersect_unit_cells) {
  // test that we can construct a matpoly and intersect with a cell
  // the source and target geometries are both unit cells
  int mat_id = 1;

  double xl = 0., xh = 1., yl = 0., yh = 1.;

  // create data  for a unit square
  std::vector<Tangram::Point<3>> square_points = {
      Tangram::Point<3>(xl, yl), Tangram::Point<3>(xh, yl),
      Tangram::Point<3>(xh, yh), Tangram::Point<3>(xl, yh)};
  std::vector<std::vector<int>> square_faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  // create the matpoly
  Tangram::MatPoly<3> square_matpoly;
  square_matpoly.set_mat_id(mat_id);
  square_matpoly.initialize(square_points);

  // extract the matpoly points
  std::vector<Tangram::Point<3>> _source_points = square_matpoly.points();

  // unfortunately we seem to need to use portage points only in intersection
  // so we need to convert from Tangram points to Portage points
  std::vector<Portage::Point<3>> source_points;
  for (auto p : _source_points) source_points.push_back(Portage::Point<3>(p));

  // create a simple mesh with a single cell
  Portage::Simple_Mesh mesh(xl, yl, xh, yh, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper meshWrapper(mesh);

  // get the coordinates of the single cell
  std::vector<Portage::Point<3>> target_points;
  meshWrapper.cell_get_coordinates(0, &target_points);

  // actually intersect
  std::vector<double> moments =
      Portage::intersect_2Dpolys(source_points, target_points);

  // test that the moments are correct
  ASSERT_NEAR(moments[0], (xh - xl) * (yh - yl), eps);
  ASSERT_NEAR(moments[1], moments[0] * (xh - xl) / 2., eps);
  ASSERT_NEAR(moments[2], moments[0] * (yh - yl) / 2., eps);
}

TEST(TANGRAM, test_matpoly_intersect_non_coincident) {
  // test that we can construct a matpoly and intersect with a cell
  // the source and target geometries are side 4 squares that intersect
  // at a corner
  int mat_id = 1;

  double xl = 0., xh = 4., yl = 0., yh = 4., xoffset = 2, yoffset = 2;

  // create data  for a unit square
  std::vector<Tangram::Point<3>> square_points = {
      Tangram::Point<3>(xl, yl), Tangram::Point<3>(xh, yl),
      Tangram::Point<3>(xh, yh), Tangram::Point<3>(xl, yh)};
  std::vector<std::vector<int>> square_faces = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};

  // create the matpoly
  Tangram::MatPoly<3> square_matpoly;
  square_matpoly.set_mat_id(mat_id);
  square_matpoly.initialize(square_points);

  // extract the matpoly points
  std::vector<Tangram::Point<3>> _source_points = square_matpoly.points();

  // unfortunately we seem to need to use portage points only in intersection
  // so we need to convert from Tangram points to Portage points
  std::vector<Portage::Point<3>> source_points;
  for (auto p : _source_points) source_points.push_back(Portage::Point<3>(p));

  // create a simple mesh with a single cell
  Portage::Simple_Mesh mesh(xl + xoffset, yl + yoffset, xh + xoffset,
                            yh + yoffset, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper meshWrapper(mesh);

  // get the coordinates of the single cell
  std::vector<Portage::Point<3>> target_points;
  meshWrapper.cell_get_coordinates(0, &target_points);

  // actually intersect
  std::vector<double> moments =
      Portage::intersect_2Dpolys(source_points, target_points);

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 4., eps);
  ASSERT_NEAR(moments[1], 12., eps);
  ASSERT_NEAR(moments[2], 12., eps);
}

TEST(TANGRAM, test_cellmatpoly_succeeds) {
  // test that we can create a Tangram CellMatPoly

  Tangram::CellMatPoly<3> cellmatpoly;
  SUCCEED();
}

TEST(TANGRAM, test_cellmatpoly_create) {
  // create a CellMatPoly by hand
  // taken from tangram/src/driver/test/test_CellMatPoly_2D

  // dimensions
  double xl = 0., xh = 1., yl = 0., yh = 1.;

  // create a simple mesh with a single cell
  Portage::Simple_Mesh mesh(xl, yl, xh, yh, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper meshWrapper(mesh);

  // Make a 2-material CellMatPoly object for the cell

  Tangram::CellMatPoly<3> cellmatpoly;

  // Cell info
  int cellid = 0;
  cellmatpoly.set_cell(cellid);

  std::vector<int> cnodes(4);
  meshWrapper.cell_get_nodes(cellid, &cnodes);

  // Two material polygons
  //
  //   o----------o
  //   |         /|
  //   |        / |
  //   |       /  |
  //   |      /   |
  //   o-----o----o

  // Specify the first polygon with parent info for points and faces
  std::vector<Tangram::Point<3>> points0(4);
  std::vector<Tangram::Entity_kind> vparentkind0(4);
  std::vector<int> vparentid0(4);

  points0[0] = Tangram::Point<3>(0.0, 0.0);
  vparentkind0[0] = Tangram::Entity_kind::NODE;
  vparentid0[0] = 0;  // parent is node 0
  points0[1] = Tangram::Point<3>(0.5, 0.0);
  vparentkind0[1] = Tangram::Entity_kind::FACE;
  vparentid0[1] = 0;  // parent is face 0
  points0[2] = Tangram::Point<3>(1.0, 1.0);
  vparentkind0[2] = Tangram::Entity_kind::NODE;
  vparentid0[2] = 2;  // parent is node 2
  points0[3] = Tangram::Point<3>(0.0, 1.0);
  vparentkind0[3] = Tangram::Entity_kind::NODE;
  vparentid0[3] = 3;  // parent is node 3

  std::vector<Tangram::Entity_kind> fparentkind0(4);
  std::vector<int> fparentid0(4);

  fparentkind0[0] = Tangram::Entity_kind::FACE;
  fparentid0[0] = 0;  // parent is face 0
  fparentkind0[1] = Tangram::Entity_kind::CELL;
  fparentid0[1] = 0;  // parent is cell 0
  fparentkind0[2] = Tangram::Entity_kind::FACE;
  fparentid0[2] = 2;  // parent is face 2
  fparentkind0[3] = Tangram::Entity_kind::FACE;
  fparentid0[3] = 3;  // parent is face 3

  // expected vertices of faces for matpoly 0

  std::vector<std::vector<int>> fverts0;
  fverts0.push_back({0, 1});
  fverts0.push_back({1, 2});
  fverts0.push_back({2, 3});
  fverts0.push_back({3, 0});

  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys0;
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, 1});  // 2nd face is internal; cncted to matpolys 0,1
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, -1});

  cellmatpoly.add_matpoly(0, 4, &(points0[0]), &(vparentkind0[0]),
                          &(vparentid0[0]), &(fparentkind0[0]),
                          &(fparentid0[0]));

  // Specify the second polygon with parent info only for points.
  // Trigger a particular check (and possibility of failure) by
  // specifying the internal face points first

  std::vector<Tangram::Point<3>> points1(3);
  std::vector<Tangram::Entity_kind> vparentkind1(3);
  std::vector<int> vparentid1(3);

  points1[0] = Tangram::Point<3>(1.0, 1.0);
  vparentkind1[0] = Tangram::Entity_kind::NODE;
  vparentid1[0] = 2;  // parent is node 2
  points1[1] = Tangram::Point<3>(0.5, 0.0);
  vparentkind1[1] = Tangram::Entity_kind::FACE;
  vparentid1[1] = 0;  // parent is face 0
  points1[2] = Tangram::Point<3>(1.0, 0.0);
  vparentkind1[2] = Tangram::Entity_kind::NODE;
  vparentid1[2] = 1;  // parent is node 1

  // expected vertices of faces for matpoly 1

  std::vector<std::vector<int>> fverts1;
  fverts1.push_back({1, 2});  // this face got defined earlier as 1, 2
  fverts1.push_back({1, 4});
  fverts1.push_back({4, 2});

  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys1;
  fmatpolys1.push_back({0, 1});  // 1st face is internal; cncted to matpolys 0,1
  fmatpolys1.push_back({1, -1});
  fmatpolys1.push_back({1, -1});

  cellmatpoly.add_matpoly(1, 3, &(points1[0]), &(vparentkind1[0]),
                          &(vparentid1[0]), nullptr, nullptr);

  // Verify the matpoly info in cell 0

  ASSERT_EQ(cellid, cellmatpoly.cell());

  // There should be two material polygons, with material IDs 0 and 1

  ASSERT_EQ(2, cellmatpoly.num_matpolys());
  ASSERT_EQ(0, cellmatpoly.matpoly_matid(0));
  ASSERT_EQ(1, cellmatpoly.matpoly_matid(1));

  // Verify info for matpoly 0
  {
    std::vector<int> mverts_out = cellmatpoly.matpoly_vertices(0);
    ASSERT_EQ(4, mverts_out.size());
    ASSERT_EQ(0, mverts_out[0]);
    ASSERT_EQ(1, mverts_out[1]);
    ASSERT_EQ(2, mverts_out[2]);
    ASSERT_EQ(3, mverts_out[3]);

    for (int i = 0; i < 4; i++) {
      int v = mverts_out[i];
      ASSERT_EQ(vparentkind0[i], cellmatpoly.matvertex_parent_kind(v));
      ASSERT_EQ(vparentid0[i], cellmatpoly.matvertex_parent_id(v));
    }

    std::vector<Tangram::Point<3>> mpoints = cellmatpoly.matpoly_points(0);
    ASSERT_EQ(4, mpoints.size());
    for (int i = 0; i < 4; i++) {
      ASSERT_TRUE(approxEq(points0[i], mpoints[i], 1.0e-8));
    }

    std::vector<int> const& mfaces_out = cellmatpoly.matpoly_faces(0);
    ASSERT_EQ(4, mfaces_out.size());
    ASSERT_EQ(0, mfaces_out[0]);
    ASSERT_EQ(1, mfaces_out[1]);
    ASSERT_EQ(2, mfaces_out[2]);
    ASSERT_EQ(3, mfaces_out[3]);

    for (int i = 0; i < 4; i++) {
      int f = mfaces_out[i];
      std::vector<int> const& mfverts_out = cellmatpoly.matface_vertices(f);
      for (int j = 0; j < mfverts_out.size(); j++)
        ASSERT_EQ(fverts0[i][j], mfverts_out[j]);

      int fmatpolys_out[2];
      cellmatpoly.matface_matpolys(f, &(fmatpolys_out[0]), &(fmatpolys_out[1]));
      ASSERT_EQ(fmatpolys0[i][0], fmatpolys_out[0]);
      ASSERT_EQ(fmatpolys0[i][1], fmatpolys_out[1]);

      if (fmatpolys_out[0] != -1 && fmatpolys_out[1] != -1)
        ASSERT_TRUE(cellmatpoly.matface_is_interface(mfaces_out[i]));
      else
        ASSERT_TRUE(!cellmatpoly.matface_is_interface(i));

      ASSERT_EQ(fparentkind0[i], cellmatpoly.matface_parent_kind(f));
      ASSERT_EQ(fparentid0[i], cellmatpoly.matface_parent_id(f));
    }

    // Verify volume

    ASSERT_NEAR(0.75 * meshWrapper.cell_volume(0),
                cellmatpoly.matpoly_volume(0), 1.0e-08);

    // Verify centroid

    Tangram::Point<3> expcen0(0.388888889, 0.55555556);
    Tangram::Point<3> mcen0 = cellmatpoly.matpoly_centroid(0);
    ASSERT_TRUE(approxEq(expcen0, mcen0, 1.0e-08));
  }

  // Verify info for matpoly 1

  {
    std::vector<int> mverts_out = cellmatpoly.matpoly_vertices(1);
    ASSERT_EQ(3, mverts_out.size());
    ASSERT_EQ(2, mverts_out[0]);
    ASSERT_EQ(1, mverts_out[1]);
    ASSERT_EQ(4, mverts_out[2]);

    for (int i = 0; i < 3; i++) {
      int v = mverts_out[i];
      ASSERT_EQ(vparentkind1[i], cellmatpoly.matvertex_parent_kind(v));
      ASSERT_EQ(vparentid1[i], cellmatpoly.matvertex_parent_id(v));
    }

    std::vector<Tangram::Point<3>> mpoints = cellmatpoly.matpoly_points(1);
    ASSERT_EQ(3, mpoints.size());
    for (int i = 0; i < 3; i++) {
      ASSERT_TRUE(approxEq(points1[i], mpoints[i], 1.0e-08));
    }

    std::vector<int> const& mfaces_out = cellmatpoly.matpoly_faces(1);
    ASSERT_EQ(3, mfaces_out.size());
    ASSERT_EQ(1, mfaces_out[0]);
    ASSERT_EQ(4, mfaces_out[1]);
    ASSERT_EQ(5, mfaces_out[2]);

    for (int i = 0; i < 3; i++) {
      int f = mfaces_out[i];
      std::vector<int> const& mfverts_out = cellmatpoly.matface_vertices(f);
      for (int j = 0; j < mfverts_out.size(); j++)
        ASSERT_EQ(fverts1[i][j], mfverts_out[j]);

      int fmatpolys_out[2];
      cellmatpoly.matface_matpolys(f, &(fmatpolys_out[0]), &(fmatpolys_out[1]));
      ASSERT_EQ(fmatpolys1[i][0], fmatpolys_out[0]);
      ASSERT_EQ(fmatpolys1[i][1], fmatpolys_out[1]);

      if (fmatpolys_out[0] != -1 && fmatpolys_out[1] != -1)
        ASSERT_TRUE(cellmatpoly.matface_is_interface(f));
      else
        ASSERT_TRUE(!cellmatpoly.matface_is_interface(f));

      if (f != 1)
        ASSERT_EQ(Tangram::Entity_kind::UNKNOWN_KIND,
                  cellmatpoly.matface_parent_kind(f));
      else
        ASSERT_EQ(Tangram::Entity_kind::CELL,
                  cellmatpoly.matface_parent_kind(f));

      if (f != 1)
        ASSERT_EQ(-1, cellmatpoly.matface_parent_id(f));
      else
        ASSERT_EQ(0, cellmatpoly.matface_parent_id(f));
    }

    ASSERT_NEAR(0.25 * meshWrapper.cell_volume(0),
                cellmatpoly.matpoly_volume(1), 1.0e-08);

    Tangram::Point<3> expcen1(2.5 / 3.0, 1.0 / 3.0);
    Tangram::Point<3> mcen1 = cellmatpoly.matpoly_centroid(1);
    ASSERT_TRUE(approxEq(expcen1, mcen1, 1.0e-08));
  }

  // Extract matpoly 1 as a MatPoly object
  {
    Tangram::MatPoly<3> MatPoly1 = cellmatpoly.get_ith_matpoly(1);

    // Verify material ID
    ASSERT_EQ(1, MatPoly1.mat_id());

    // Verify vertices
    const std::vector<Tangram::Point2>& matpoly_points = MatPoly1.points();
    ASSERT_EQ(3, MatPoly1.num_vertices());
    for (int ivrt = 0; ivrt < 3; ivrt++)
      ASSERT_TRUE(approxEq(points1[ivrt], matpoly_points[ivrt], 1.0e-15));

    std::vector<std::vector<int>> expected_mp1_faces = {{0, 1}, {1, 2}, {2, 0}};

    // Verify faces
    ASSERT_EQ(3, MatPoly1.num_faces());
    for (int iface = 0; iface < 3; iface++) {
      const std::vector<int>& face_vertices = MatPoly1.face_vertices(iface);
      ASSERT_EQ(2, face_vertices.size());
      ASSERT_EQ(expected_mp1_faces[iface][0], face_vertices[0]);
      ASSERT_EQ(expected_mp1_faces[iface][1], face_vertices[1]);
    }
  }
}

TEST(TANGRAM, test_cellmatpoly_intersect) {
  // create a CellMatPoly by hand
  // taken from tangram/src/driver/test/test_CellMatPoly_2D

  // dimensions
  double xl = 0., xh = 1., yl = 0., yh = 1., xoffset = 0., yoffset = 0.;

  // create a simple mesh with a single cell
  Portage::Simple_Mesh mesh(xl, yl, xh, yh, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper meshWrapper(mesh);

  // Make a 2-material CellMatPoly object for the cell

  Tangram::CellMatPoly<3> cellmatpoly;

  // Cell info
  int cellid = 0;
  cellmatpoly.set_cell(cellid);

  std::vector<int> cnodes(4);
  meshWrapper.cell_get_nodes(cellid, &cnodes);

  // Two material polygons
  //
  //   o----------o
  //   |         /|
  //   |        / |
  //   |       /  |
  //   |      /   |
  //   o-----o----o

  // Specify the first polygon with parent info for points and faces
  std::vector<Tangram::Point<3>> points0(4);
  std::vector<Tangram::Entity_kind> vparentkind0(4);
  std::vector<int> vparentid0(4);

  points0[0] = Tangram::Point<3>(0.0, 0.0);
  vparentkind0[0] = Tangram::Entity_kind::NODE;
  vparentid0[0] = 0;  // parent is node 0
  points0[1] = Tangram::Point<3>(0.5, 0.0);
  vparentkind0[1] = Tangram::Entity_kind::FACE;
  vparentid0[1] = 0;  // parent is face 0
  points0[2] = Tangram::Point<3>(1.0, 1.0);
  vparentkind0[2] = Tangram::Entity_kind::NODE;
  vparentid0[2] = 2;  // parent is node 2
  points0[3] = Tangram::Point<3>(0.0, 1.0);
  vparentkind0[3] = Tangram::Entity_kind::NODE;
  vparentid0[3] = 3;  // parent is node 3

  std::vector<Tangram::Entity_kind> fparentkind0(4);
  std::vector<int> fparentid0(4);

  fparentkind0[0] = Tangram::Entity_kind::FACE;
  fparentid0[0] = 0;  // parent is face 0
  fparentkind0[1] = Tangram::Entity_kind::CELL;
  fparentid0[1] = 0;  // parent is cell 0
  fparentkind0[2] = Tangram::Entity_kind::FACE;
  fparentid0[2] = 2;  // parent is face 2
  fparentkind0[3] = Tangram::Entity_kind::FACE;
  fparentid0[3] = 3;  // parent is face 3

  // expected vertices of faces for matpoly 0

  std::vector<std::vector<int>> fverts0;
  fverts0.push_back({0, 1});
  fverts0.push_back({1, 2});
  fverts0.push_back({2, 3});
  fverts0.push_back({3, 0});

  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys0;
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, 1});  // 2nd face is internal; cncted to matpolys 0,1
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, -1});

  cellmatpoly.add_matpoly(0, 4, &(points0[0]), &(vparentkind0[0]),
                          &(vparentid0[0]), &(fparentkind0[0]),
                          &(fparentid0[0]));

  // Specify the second polygon with parent info only for points.
  // Trigger a particular check (and possibility of failure) by
  // specifying the internal face points first

  std::vector<Tangram::Point<3>> points1(3);
  std::vector<Tangram::Entity_kind> vparentkind1(3);
  std::vector<int> vparentid1(3);

  points1[0] = Tangram::Point<3>(1.0, 1.0);
  vparentkind1[0] = Tangram::Entity_kind::NODE;
  vparentid1[0] = 2;  // parent is node 2
  points1[1] = Tangram::Point<3>(0.5, 0.0);
  vparentkind1[1] = Tangram::Entity_kind::FACE;
  vparentid1[1] = 0;  // parent is face 0
  points1[2] = Tangram::Point<3>(1.0, 0.0);
  vparentkind1[2] = Tangram::Entity_kind::NODE;
  vparentid1[2] = 1;  // parent is node 1

  // expected vertices of faces for matpoly 1

  std::vector<std::vector<int>> fverts1;
  fverts1.push_back({1, 2});  // this face got defined earlier as 1, 2
  fverts1.push_back({1, 4});
  fverts1.push_back({4, 2});

  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys1;
  fmatpolys1.push_back({0, 1});  // 1st face is internal; cncted to matpolys 0,1
  fmatpolys1.push_back({1, -1});
  fmatpolys1.push_back({1, -1});

  cellmatpoly.add_matpoly(1, 3, &(points1[0]), &(vparentkind1[0]),
                          &(vparentid1[0]), nullptr, nullptr);

  // create a simple targetMesh with a single cell
  Portage::Simple_Mesh targetMesh(xl + xoffset, yl + yoffset, xh + xoffset,
                                  yh + yoffset, 1, 1);

  // Create mesh wrappers
  Wonton::Simple_Mesh_Wrapper targetMeshWrapper(targetMesh);

  // get the coordinates of the single cell
  std::vector<Portage::Point<3>> target_points;
  targetMeshWrapper.cell_get_coordinates(0, &target_points);

  // extract the MatPoly's from the CellMatPoly
  int nMatPoly = cellmatpoly.num_matpolys();

  for (int i = 0; i < nMatPoly; ++i) {
    // extract the i'th matpoly
    Tangram::MatPoly<3> matPoly = cellmatpoly.get_ith_matpoly(i);

    // extract the matpoly points
    std::vector<Tangram::Point<3>> _matPoly_points = matPoly.points();

    // unfortunately we seem to need to use portage points only in intersection
    // so we need to convert from Tangram points to Portage points
    std::vector<Portage::Point<3>> matPoly_points;
    for (auto p : _matPoly_points)
      matPoly_points.push_back(Portage::Point<3>(p));

    // actually intersect
    std::vector<double> moments =
        Portage::intersect_2Dpolys(matPoly_points, target_points);

    if (i == 0) {
      // Verify volume
      ASSERT_NEAR(0.75 * meshWrapper.cell_volume(0), moments[0], eps);

      // calculated centroid
      ASSERT_NEAR(.38888889, moments[1] / moments[0], eps);
      ASSERT_NEAR(.55555556, moments[2] / moments[0], eps);

    } else {
      // Verify volume
      ASSERT_NEAR(0.25 * meshWrapper.cell_volume(0), moments[0], eps);

      // calculated centroid
      ASSERT_NEAR(.83333333, moments[1] / moments[0], eps);
      ASSERT_NEAR(.33333333, moments[2] / moments[0], eps);
    }
  }
}
*/
#endif
