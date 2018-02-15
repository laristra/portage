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
TEST(TANGRAM_3D, test_tangram_exists) {
  // test that we can instantiate a Tangram Point even though we don't include
  // the Tangram header

  float x{1.}, y{2}, z{3};
  Tangram::Point<3> p{x, y, z};
  ASSERT_FLOAT_EQ(x, p[0]);
  ASSERT_FLOAT_EQ(y, p[1]);
  ASSERT_FLOAT_EQ(z, p[2]);
}

TEST(TANGRAM_3D, test_portage_exists) {
  // test that we can instantiate a Portage Point

  float x{1.}, y{2}, z{3};
  Portage::Point<3> p{x, y, z};
  ASSERT_FLOAT_EQ(x, p[0]);
  ASSERT_FLOAT_EQ(y, p[1]);
  ASSERT_FLOAT_EQ(z, p[2]);
}

TEST(TANGRAM_3D, test_tangram_to_portage) {
  // test that we can create a Portage Point from a Tangram Point

  float x{1.}, y{2}, z{3};
  Tangram::Point<3> pt{x, y,z};
  Portage::Point<3> pp{pt};
  ASSERT_FLOAT_EQ(x, pp[0]);
  ASSERT_FLOAT_EQ(y, pp[1]);
  ASSERT_FLOAT_EQ(z, pp[2]);
}

TEST(TANGRAM_3D, test_portage_to_tangram) {
  // test that we can create a Tangram Point from a Portage Point

  float x{1.}, y{2}, z{3};
  Portage::Point<3> pp{x, y, z};
  Tangram::Point<3> pt{pp};
  ASSERT_FLOAT_EQ(x, pt[0]);
  ASSERT_FLOAT_EQ(y, pt[1]);
  ASSERT_FLOAT_EQ(z, pt[2]);
}

TEST(TANGRAM_3D, test_matpoly_succeeds) {
  // test that we can create a Tangram MatPoly

  Tangram::MatPoly<3> matpoly;
  SUCCEED();
}

TEST(TANGRAM_3D, test_matpoly_create) {
  // test that we can construct a real matpoly (lifted from
  // tangram/src/support/test/test_MatPoly_3D.cc)
  int mat_id = 1;
  
  //Test for a right triangular prism
  std::vector<Tangram::Point<3>> prism_points = {
    Tangram::Point<3>(1.0, 0.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 0.0), Tangram::Point<3>(0.0, 1.0, 1.0),
    Tangram::Point<3>(1.0, 0.0, 1.0), Tangram::Point<3>(0.0, 0.0, 1.0) };
  std::vector< std::vector<int> > prism_faces = {
    {0, 2, 1}, {2, 0, 4, 5}, {3, 4, 0, 1}, {5, 2, 1, 3}, {3, 4, 5} };
  std::vector<Tangram::Point<3>> face_centroids = {
    Tangram::Point<3>(1.0/3.0, 1.0/3.0, 0.0), Tangram::Point<3>(0.5, 0.0, 0.5),
    Tangram::Point<3>(0.5, 0.5, 0.5), Tangram::Point<3>(0.0, 0.5, 0.5),
    Tangram::Point<3>(1.0/3.0, 1.0/3.0, 1.0) };
  
  //Check material ID correctness
  Tangram::MatPoly<3> prism_matpoly(mat_id);
  ASSERT_EQ(mat_id, prism_matpoly.mat_id());
  
  //Initialization
  prism_matpoly.initialize(prism_points, prism_faces);
  
  //Verify coordinates
  const std::vector<Tangram::Point<3>>& matpoly_points = prism_matpoly.points();
  ASSERT_EQ(prism_points.size(), prism_matpoly.num_vertices());
  for (int ivrt = 0; ivrt < prism_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(prism_points[ivrt], matpoly_points[ivrt], 1.0e-15));
  
  //Verify faces
  ASSERT_EQ(prism_faces.size(), prism_matpoly.num_faces());
  for (int iface = 0; iface < prism_faces.size(); iface++) {
    const std::vector<int>& face_vertices = prism_matpoly.face_vertices(iface);
    ASSERT_EQ(prism_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < prism_faces[iface].size(); ivrt++)
      ASSERT_EQ(prism_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Verify centroids
  for (int iface = 0; iface < prism_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         prism_matpoly.face_centroid(iface), 1.0e-15));
  
  std::vector<Tangram::Point<3>> faceted_prism_points = {
    Tangram::Point<3>(1.0, 0.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 0.0), Tangram::Point<3>(0.0, 1.0, 1.0),
    Tangram::Point<3>(1.0, 0.0, 1.0), Tangram::Point<3>(0.0, 0.0, 1.0),
    Tangram::Point<3>(0.5, 0.0, 0.5), Tangram::Point<3>(0.5, 0.5, 0.5),
    Tangram::Point<3>(0.0, 0.5, 0.5),
  };
  std::vector< std::vector<int> > faceted_prism_faces = {
    {0, 2, 1},
    {6, 2, 0}, {6, 0, 4}, {6, 4, 5}, {6, 5, 2},
    {7, 3, 4}, {7, 4, 0}, {7, 0, 1}, {7, 1, 3},
    {8, 5, 2}, {8, 2, 1}, {8, 1, 3}, {8, 3, 5},
    {3, 4, 5} };
  
  //Create faceted poly
  Tangram::MatPoly<3> faceted_prism_matpoly;
  prism_matpoly.faceted_matpoly(&faceted_prism_matpoly);
  
  //Check material ID correctness
  ASSERT_EQ(mat_id, faceted_prism_matpoly.mat_id());
  
  //Verify facetization
  //Verify node coordinates
  const std::vector<Tangram::Point<3>>& faceted_matpoly_points =
    faceted_prism_matpoly.points();
  ASSERT_EQ(faceted_prism_points.size(), faceted_prism_matpoly.num_vertices());
  for (int ivrt = 0; ivrt < faceted_prism_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(faceted_prism_points[ivrt],
                         faceted_matpoly_points[ivrt], 1.0e-15));
  
  //Verify facets
  ASSERT_EQ(faceted_prism_faces.size(), faceted_prism_matpoly.num_faces());
  for (int iface = 0; iface < faceted_prism_faces.size(); iface++) {
    const std::vector<int>& face_vertices = faceted_prism_matpoly.face_vertices(iface);
    ASSERT_EQ(faceted_prism_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < faceted_prism_faces[iface].size(); ivrt++)
      ASSERT_EQ(faceted_prism_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Non-convex distorted prism
  std::vector<Tangram::Point<3>> ncv_prism_points = {
    Tangram::Point<3>(1.0, 0.0, 0.0), Tangram::Point<3>(0.4, 0.8, 0.2),
    Tangram::Point<3>(1.0, 1.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.5, 0.1, 0.1), Tangram::Point<3>(0.0, 0.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 1.0), Tangram::Point<3>(1.0, 1.0, 1.0),
    Tangram::Point<3>(0.5, 0.9, 1.1), Tangram::Point<3>(1.0, 0.0, 1.0),
    Tangram::Point<3>(0.0, 1.0, 1.0), Tangram::Point<3>(0.6, 0.2, 0.8) };
  std::vector< std::vector<int> > ncv_prism_faces = {
    {3, 1, 2, 0, 4, 5},
    {5, 4, 11, 10}, {0, 9, 11, 4}, {7, 9, 0, 2},
    {2, 1, 8, 7}, {8, 1, 3, 10}, {5, 6, 10, 3},
    {7, 11, 10, 6, 8, 9} };
  std::vector<Tangram::Point<3>> ncv_prism_face_centroids = {
    Tangram::Point<3>(0.49982029799691324, 0.48793283780407681, 0.038845671672909921),
    Tangram::Point<3>(0.2519019047393049, 0.36565382681732062, 0.51522129026172814),
    Tangram::Point<3>(0.78729807432173626, 0.070061777159007799, 0.46574083274162237),
    Tangram::Point<3>(1.0, 0.5, 0.5),
    Tangram::Point<3>(0.72714164844017881, 0.92490808257951729, 0.55650182505587276),
    Tangram::Point<3>(0.22120243710721832, 0.92700420108481663, 0.58407100072352491),
    Tangram::Point<3>(0.0, 0.5, 0.5),
    Tangram::Point<3>(0.50915092483338054, 0.51261128330217054, 0.98738871669782957) };
  
  Tangram::MatPoly<3> ncv_prism_matpoly(mat_id);
  //Initialization
  ncv_prism_matpoly.initialize(ncv_prism_points, ncv_prism_faces);
  
  //Verify centroids
  for (int iface = 0; iface < ncv_prism_faces.size(); iface++)
    ASSERT_TRUE(approxEq(ncv_prism_face_centroids[iface],
                         ncv_prism_matpoly.face_centroid(iface), 1.0e-15));
}


TEST(TANGRAM_3D, test_matpoly_cube) {
  // test that we can construct a cube matpoly and extract its points
  int mat_id = 1;
  
  //Test for a cube
  std::vector<Tangram::Point<3>> cube_points = {
    Tangram::Point<3>(0.0, 0.0, 0.0), Tangram::Point<3>(1.0, 0.0, 0.0),
    Tangram::Point<3>(1.0, 1.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 1.0), Tangram::Point<3>(1.0, 0.0, 1.0),
    Tangram::Point<3>(1.0, 1.0, 1.0), Tangram::Point<3>(0.0, 1.0, 1.0)};
  std::vector< std::vector<int> >cube_faces = {
    {0, 1, 5, 4},{1, 2, 6, 5},{2, 3, 7, 6},{3, 0, 4, 7},{0, 3, 2, 1},
    {4, 5, 6, 7} };
  std::vector<Tangram::Point<3>> face_centroids = {
    Tangram::Point<3>(0.5, 0.0, 0.5), Tangram::Point<3>(1.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 1.0, 0.5), Tangram::Point<3>(0.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 0.5, 0.0), Tangram::Point<3>(0.5, 0.5, 1.0) };
  
  //Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  ASSERT_EQ(mat_id, matpoly.mat_id());
  
  //Initialization
  matpoly.initialize(cube_points, cube_faces);
  
  //Verify coordinates
  const std::vector<Tangram::Point<3>>& matpoly_points = matpoly.points();
  ASSERT_EQ(cube_points.size(), matpoly.num_vertices());
  for (int ivrt = 0; ivrt < cube_points.size(); ivrt++)
    for (int j=0; j<3; j++)
    	ASSERT_NEAR(cube_points[ivrt][j], matpoly_points[ivrt][j], 1.0e-15);
    	
  //Verify faces
  ASSERT_EQ(cube_faces.size(), matpoly.num_faces());
  for (int iface = 0; iface < cube_faces.size(); iface++) {
    const std::vector<int>& face_vertices = matpoly.face_vertices(iface);
    ASSERT_EQ(cube_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < cube_faces[iface].size(); ivrt++)
      ASSERT_EQ(cube_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Verify centroids
  for (int iface = 0; iface < cube_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         matpoly.face_centroid(iface), 1.0e-15));
  
    	
  std::cout <<"finished"<<std::endl;
}


TEST(TANGRAM_3D, test_matpoly_faceted_cube_by_hand) {
  // test that we can construct a faceted cube matpoly by hand and extract its points
  int mat_id = 1;
  
  //Test for a cube
  std::vector<Tangram::Point<3>> faceted_cube_points = {
    Tangram::Point<3>(0.0, 0.0, 0.0), Tangram::Point<3>(1.0, 0.0, 0.0),
    Tangram::Point<3>(1.0, 1.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 1.0), Tangram::Point<3>(1.0, 0.0, 1.0),
    Tangram::Point<3>(1.0, 1.0, 1.0), Tangram::Point<3>(0.0, 1.0, 1.0),
    Tangram::Point<3>(0.5, 0.0, 0.5), Tangram::Point<3>(1.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 1.0, 0.5), Tangram::Point<3>(0.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 0.5, 0.0), Tangram::Point<3>(0.5, 0.5, 1.0)     
    };
  std::vector< std::vector<int> >faceted_cube_faces = {
    {0, 1, 8}, {1,5,8}, {5,4,8}, {4,0,8}, 
    {1,2,9}, {2,6,9}, {6,5,9}, {5,1,9}, 
    {2,3,10}, {3,7,10}, {7,6,10}, {6,2,10}, 
    {3,0,11}, {0,4,11}, {4,7,11}, {7,3,11}, 
    {0,3,12}, {3,2,12}, {2,1,12,}, {1,0,12}, 
    {4,5,13}, {5,6,13}, {6,7,13}, {7,4,12}};

  
  //Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  ASSERT_EQ(mat_id, matpoly.mat_id());
  
  //Initialization
  matpoly.initialize(faceted_cube_points, faceted_cube_faces);
  
  //Verify coordinates
  const std::vector<Tangram::Point<3>>& matpoly_points = matpoly.points();
  ASSERT_EQ(faceted_cube_points.size(), matpoly.num_vertices());
  for (int ivrt = 0; ivrt < faceted_cube_points.size(); ivrt++)
    for (int j=0; j<3; j++)
    	ASSERT_NEAR(faceted_cube_points[ivrt][j], matpoly_points[ivrt][j], 1.0e-15);
    	
  //Verify faces
  ASSERT_EQ(faceted_cube_faces.size(), matpoly.num_faces());
  for (int iface = 0; iface < faceted_cube_faces.size(); iface++) {
    const std::vector<int>& face_vertices = matpoly.face_vertices(iface);
    ASSERT_EQ(faceted_cube_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < faceted_cube_faces[iface].size(); ivrt++)
      ASSERT_EQ(faceted_cube_faces[iface][ivrt], face_vertices[ivrt]);
  }
      	
  std::cout <<"finished 2"<<std::endl;
}

TEST(TANGRAM_3D, test_matpoly_faceted_cube) {
  // test that we can construct a cube matpoly and extract its points
  int mat_id = 1;
  
  //Test for a cube
  std::vector<Tangram::Point<3>> cube_points = {
    Tangram::Point<3>(0.0, 0.0, 0.0), Tangram::Point<3>(1.0, 0.0, 0.0),
    Tangram::Point<3>(1.0, 1.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 1.0), Tangram::Point<3>(1.0, 0.0, 1.0),
    Tangram::Point<3>(1.0, 1.0, 1.0), Tangram::Point<3>(0.0, 1.0, 1.0)};
  std::vector< std::vector<int> >cube_faces = {
    {0, 1, 5, 4},{1, 2, 6, 5},{2, 3, 7, 6},{3, 0, 4, 7},{0, 3, 2, 1},
    {4, 5, 6, 7} };
  std::vector<Tangram::Point<3>> face_centroids = {
    Tangram::Point<3>(0.5, 0.0, 0.5), Tangram::Point<3>(1.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 1.0, 0.5), Tangram::Point<3>(0.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 0.5, 0.0), Tangram::Point<3>(0.5, 0.5, 1.0) };
  
  //Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  ASSERT_EQ(mat_id, matpoly.mat_id());
  
  //Initialization
  matpoly.initialize(cube_points, cube_faces);
  
  //Verify coordinates
  const std::vector<Tangram::Point<3>>& matpoly_points = matpoly.points();
  ASSERT_EQ(cube_points.size(), matpoly.num_vertices());
  for (int ivrt = 0; ivrt < cube_points.size(); ivrt++)
    for (int j=0; j<3; j++)
    	ASSERT_NEAR(cube_points[ivrt][j], matpoly_points[ivrt][j], 1.0e-15);
    	
  //Verify faces
  ASSERT_EQ(cube_faces.size(), matpoly.num_faces());
  for (int iface = 0; iface < cube_faces.size(); iface++) {
    const std::vector<int>& face_vertices = matpoly.face_vertices(iface);
    ASSERT_EQ(cube_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < cube_faces[iface].size(); ivrt++)
      ASSERT_EQ(cube_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Verify centroids
  for (int iface = 0; iface < cube_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         matpoly.face_centroid(iface), 1.0e-15));
  
  //Facet the cube
  Tangram::MatPoly<3> faceted_matpoly;
  matpoly.faceted_matpoly(&faceted_matpoly);
  
  //Check that the faceted matpoly has 24 faces and 14 vertices
  ASSERT_EQ(faceted_matpoly.num_faces(), 24);
  ASSERT_EQ(faceted_matpoly.num_vertices(), 14);
  
  //Convert the points
  std::vector<Tangram::Point<3>> _source_points = faceted_matpoly.points();
  std::vector<Portage::Point<3>> source_points;
  for (auto p : _source_points) source_points.push_back(Portage::Point<3>(p));
  
  //Convert the face indices from the matpoly to a Portage style
  std::vector<std::vector<int>> facetpoints;
  for (int iface = 0; iface < faceted_matpoly.num_faces(); ++iface){
    const std::vector<int>& face_vertices = faceted_matpoly.face_vertices(iface);
    //for (int f: face_vertices) std::cout << f << ", ";
  	//std::cout <<std::endl;
    facetpoints.push_back(face_vertices);
  }
  
  //Create the facetedpoly_t structure
  Portage::facetedpoly_t srcpoly{facetpoints, source_points};
  
  //Create a length 1 vector of tets by hand
  std::vector<std::array<Portage::Point<3>, 4>> target_tet_coords{{
    Portage::Point<3>{0.,0.,0.},
    Portage::Point<3>{1.,0.,0.},
    Portage::Point<3>{0.,1.,0.},
    Portage::Point<3>{0.,0.,1.}}};
    
  //Hope for a miracle with intersection
  std::vector<double> moments(Portage::intersect_3Dpolys(srcpoly, target_tet_coords));

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 1./6., eps); 
  ASSERT_NEAR(moments[1]/moments[0], .25, eps);
  ASSERT_NEAR(moments[2]/moments[0], .25, eps);
  ASSERT_NEAR(moments[3]/moments[0], .25, eps);
  
  for (auto m: moments) std::cout << m << ", ";
  std::cout <<std::endl;

    
  std::cout <<"finished 3"<<std::endl;
  
}

TEST(TANGRAM_3D, test_matpoly_faceted_cube2) {
  // test that we can construct a cube matpoly and extract its points
  int mat_id = 1;
  
  //Test for a cube
  std::vector<Tangram::Point<3>> cube_points = {
    Tangram::Point<3>(0.0, 0.0, 0.0), Tangram::Point<3>(1.0, 0.0, 0.0),
    Tangram::Point<3>(1.0, 1.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 1.0), Tangram::Point<3>(1.0, 0.0, 1.0),
    Tangram::Point<3>(1.0, 1.0, 1.0), Tangram::Point<3>(0.0, 1.0, 1.0)};
  std::vector< std::vector<int> >cube_faces = {
    {0, 1, 5, 4},{1, 2, 6, 5},{2, 3, 7, 6},{3, 0, 4, 7},{0, 3, 2, 1},
    {4, 5, 6, 7} };
  std::vector<Tangram::Point<3>> face_centroids = {
    Tangram::Point<3>(0.5, 0.0, 0.5), Tangram::Point<3>(1.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 1.0, 0.5), Tangram::Point<3>(0.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 0.5, 0.0), Tangram::Point<3>(0.5, 0.5, 1.0) };
  
  //Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  ASSERT_EQ(mat_id, matpoly.mat_id());
  
  //Initialization
  matpoly.initialize(cube_points, cube_faces);
  
  //Verify coordinates
  const std::vector<Tangram::Point<3>>& matpoly_points = matpoly.points();
  ASSERT_EQ(cube_points.size(), matpoly.num_vertices());
  for (int ivrt = 0; ivrt < cube_points.size(); ivrt++)
    for (int j=0; j<3; j++)
    	ASSERT_NEAR(cube_points[ivrt][j], matpoly_points[ivrt][j], 1.0e-15);
    	
  //Verify faces
  ASSERT_EQ(cube_faces.size(), matpoly.num_faces());
  for (int iface = 0; iface < cube_faces.size(); iface++) {
    const std::vector<int>& face_vertices = matpoly.face_vertices(iface);
    ASSERT_EQ(cube_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < cube_faces[iface].size(); ivrt++)
      ASSERT_EQ(cube_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Verify centroids
  for (int iface = 0; iface < cube_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         matpoly.face_centroid(iface), 1.0e-15));
  
  //Facet the cube
  Tangram::MatPoly<3> faceted_matpoly;
  matpoly.faceted_matpoly(&faceted_matpoly);
  
  //Check that the faceted matpoly has 24 faces and 14 vertices
  ASSERT_EQ(faceted_matpoly.num_faces(), 24);
  ASSERT_EQ(faceted_matpoly.num_vertices(), 14);
  
  //Convert the points
  std::vector<Tangram::Point<3>> _source_points = faceted_matpoly.points();
  std::vector<Portage::Point<3>> source_points;
  for (auto p : _source_points) source_points.push_back(Portage::Point<3>(p));
  
  //Convert the face indices from the matpoly to a Portage style
  /* WITH THE NEW API FUNCTION THIS IS NOW UNNECESSARY
  std::vector<std::vector<int>> facetpoints;
  for (int iface = 0; iface < faceted_matpoly.num_faces(); ++iface){
    const std::vector<int>& face_vertices = faceted_matpoly.face_vertices(iface);
    //for (int f: face_vertices) std::cout << f << ", ";
  	//std::cout <<std::endl;
    facetpoints.push_back(face_vertices);
  }
  */
  
  //Create the facetedpoly_t structure
  Portage::facetedpoly_t srcpoly{faceted_matpoly.face_vertices(), source_points};
  //NOW DIFFERENT
  //Portage::facetedpoly_t srcpoly{facetpoints, source_points};
  
  //Create a length 1 vector of tets by hand
  std::vector<std::array<Portage::Point<3>, 4>> target_tet_coords{{
    Portage::Point<3>{0.,0.,0.},
    Portage::Point<3>{1.,0.,0.},
    Portage::Point<3>{0.,1.,0.},
    Portage::Point<3>{0.,0.,1.}}};
    
  //Hope for a miracle with intersection
  std::vector<double> moments(Portage::intersect_3Dpolys(srcpoly, target_tet_coords));

  // test that the moments are correct
  ASSERT_NEAR(moments[0], 1./6., eps); 
  ASSERT_NEAR(moments[1]/moments[0], .25, eps);
  ASSERT_NEAR(moments[2]/moments[0], .25, eps);
  ASSERT_NEAR(moments[3]/moments[0], .25, eps);
  
  for (auto m: moments) std::cout << m << ", ";
  std::cout <<std::endl;

    
  std::cout <<"finished 3"<<std::endl;
  
}

/*
TEST(TANGRAM_3D, test_matpoly_faceted_cube_intersection) {
  // test that we can construct a cube matpoly and extract its points
  int mat_id = 1;
  
  //Test for a cube
  std::vector<Tangram::Point<3>> faceted_cube_points = {
    Tangram::Point<3>(0.0, 0.0, 0.0), Tangram::Point<3>(1.0, 0.0, 0.0),
    Tangram::Point<3>(1.0, 1.0, 0.0), Tangram::Point<3>(0.0, 1.0, 0.0),
    Tangram::Point<3>(0.0, 0.0, 1.0), Tangram::Point<3>(1.0, 0.0, 1.0),
    Tangram::Point<3>(1.0, 1.0, 1.0), Tangram::Point<3>(0.0, 1.0, 1.0),
    Tangram::Point<3>(0.5, 0.0, 0.5), Tangram::Point<3>(1.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 1.0, 0.5), Tangram::Point<3>(0.0, 0.5, 0.5),
    Tangram::Point<3>(0.5, 0.5, 0.0), Tangram::Point<3>(0.5, 0.5, 1.0)     
    };
  std::vector< std::vector<int> >faceted_cube_faces = {
    {0, 1, 8}, {1,5,8}, {5,4,8}, {4,0,8}, 
    {1,2,9}, {2,6,9}, {6,5,9}, {5,1,9}, 
    {2,3,10}, {3,7,10}, {7,6,10}, {6,2,10}, 
    {3,0,11}, {0,4,11}, {4,7,11}, {7,3,11}, 
    {0,3,12}, {3,2,12}, {2,1,12,}, {1,0,12}, 
    {4,5,13}, {5,6,13}, {6,7,13}, {7,4,12}};

  
  //Check material ID correctness
  Tangram::MatPoly<3> matpoly(mat_id);
  
  //Initialization
  matpoly.initialize(faceted_cube_points, faceted_cube_faces);
  
  //Convert the matpoly into a facetedpoly_t
  std::vector<std::vector<int>> facetpoints;
  
  std::cout <<"finished 3"<<std::endl;
}
*/
  
/*
TEST(TANGRAM_3D, test_matpoly_intersect_unit_cells) {
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

TEST(TANGRAM_3D, test_matpoly_intersect_non_coincident) {
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

TEST(TANGRAM_3D, test_cellmatpoly_succeeds) {
  // test that we can create a Tangram CellMatPoly

  Tangram::CellMatPoly<3> cellmatpoly;
  SUCCEED();
}

TEST(TANGRAM_3D, test_cellmatpoly_create) {
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

TEST(TANGRAM_3D, test_cellmatpoly_intersect) {
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
