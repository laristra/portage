/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/





#include "portage/wonton/mesh/jali/jali_mesh_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "portage/support/Point.h"

using std::abs;

/*!
  file test_jali_mesh_wrapper.cc
  @brief Unit tests for the Jali mesh wrapper class
 */

/*!
  @brief Returns true if a == b within the accuracy 'eps'
  @param[in] a First input vector of pairs of doubles
  @param[in] b Second input vector of pairs of doubles
  @param[in] eps Tolerance
  @return Whether a == b within the accuracy @c eps
 */
bool vdd_eq(const std::vector<Portage::Point<2>> &a,
    const std::vector<Portage::Point<2>> &b, const double eps = 1e-12)
{
    // Can't be equal if # of entries differ:
    if (a.size() != b.size()) return false;
    // Loop over elements in "a" and "b":
    for (size_t i = 0; i < a.size(); i++) {
        // values not equal
        if (abs(a[i][0] - b[i][0]) > eps ||
            abs(a[i][1] - b[i][1]) > eps) return false;
    }
    return true;
}

/*!
  @brief Rotate the @c xylist vector into a canonical (unique) form. The first point
  will be the one with the lowest angle between it, the nodeid and the x-axis.
  @param[in] center_node The center node to test
  @param[in,out] xylist      The xylist vector
 */
void coordinates_canonical_rotation(
      const Portage::Point<2>& center_node,
      std::vector<Portage::Point<2>> * const xylist)
{
    int i = 0;
    auto angle = [&]() {
        return std::atan2((*xylist)[i][1] - center_node[1],
                          (*xylist)[i][0] - center_node[0]);
    };
    double a = angle();
    while (a >= 0) { i++; i = i % xylist->size(); a = angle(); }
    while (a < 0) { i++; i = i % xylist->size(); a = angle(); }
    std::rotate(xylist->begin(), xylist->begin()+i, xylist->end());
}

/*!
  @brief Rotate the @c xylist vector into a canonical (unique) form. The first point
  will be the one with the lowest angle between it, the nodeid and the x-axis.
  @param[in] mesh_wrapper The mesh wrapper
  @param[in] nodeid The node id
  @param[in,out] xylist The xylist
 */
void dual_cell_coordinates_canonical_rotation(
    const Wonton::Jali_Mesh_Wrapper &mesh_wrapper,
    int const nodeid,
    std::vector<Portage::Point<2>> * const xylist) {
  Portage::Point<2> center_node;
  mesh_wrapper.node_get_coordinates(nodeid, &center_node);
  coordinates_canonical_rotation(center_node, xylist);
}

/*!
  @brief Check validity of faceted polyhedron
  @param fctpoints Vector of unique points of facetization
  @param facets  Vector of vector of facet points (indexing fctpoints)
  @param centroid  Centroid of polyhedron
  @param expected_area  Expected area of polyhedron surface (-1 means don't check)
  @param expected_volume  Expected volume of polyhedron (-1 means don't check)
  @param expected_centroid Expected centroid of polyhedron ([-99, -99, -99] means don't check)
*/

bool faceted_poly_ok(std::vector<Portage::Point<3>> const & fctpoints,
                     std::vector<std::vector<int>> const & facets,
                     double expected_area, double expected_volume,
                     Portage::Point<3> expected_centroid,
                     double tol) {

  int nfacets = facets.size();

  Portage::Point<3> gcen(0.0, 0.0, 0.0);  // Geometric center
  int npoints = fctpoints.size();
  for (int i = 0; i < npoints; i++)
    gcen += fctpoints[i];
  gcen /= npoints;

  double ctetvolume = 0.0, cfctvolume = 0.0, cfctarea = 0.0;
  Portage::Point<3> cfctcen, ctetcen;
  for (int i = 0; i < nfacets; i++) {
    int nfpoints = facets[i].size();
    for (int j = 1; j < nfpoints-1; j++) {
      Portage::Point<3> p0 = fctpoints[facets[i][0]];
      Portage::Point<3> p1 = fctpoints[facets[i][j]];
      Portage::Point<3> p2 = fctpoints[facets[i][j+1]];

      Portage::Vector<3> v0 = p1 - p0;
      Portage::Vector<3> v1 = p2 - p0;
      Portage::Vector<3> outnormal = cross(v0, v1);
      double fctarea = outnormal.norm()/2.0;

      cfctarea += fctarea;

      // Compute contribution to the volume by divergence theorem
      Portage::Vector<3> vec;
      for (int k = 0; k < 3; k++)
        vec[k] = p0[k];
      cfctvolume += dot(vec, outnormal)/6.0;

      // Compute contribution to the centroid by divergence theorem
      for (int k = 0; k < 3; k++)
        cfctcen[k] += ((p0[k]+p1[k])*(p0[k]+p1[k]) +
                       (p1[k]+p2[k])*(p1[k]+p2[k]) +
                       (p2[k]+p0[k])*(p2[k]+p0[k]))*outnormal[k]/24.0;

      // Also compute volume of a tet formed by the facet and the centroid
      Portage::Vector<3> v2 = gcen - p0;
      double tetvolume = -dot(v2, outnormal)/6.0;
      if (tetvolume < 1.0e-12) {
        std::cerr << "Tet formed by facet and centroid has volume " <<
            tetvolume << std::endl;
        return false;
      }
      ctetvolume += tetvolume;

      // Also compute centroid by volume weighted sum of tet centroids
      ctetcen += tetvolume*(p0 + p1 + p2 + gcen)/4.0;
    }
  }
  cfctcen /= (2*ctetvolume);
  ctetcen /= ctetvolume;

  if (expected_area > 0.0 && std::abs(expected_area-cfctarea) > 1.0e-12) {
    std::cerr << "Expected area: " << expected_area << std::endl;
    std::cerr << "Computed area: " << cfctarea << std::endl;
    return false;
  }
  if (expected_volume > 0.0 && std::abs(expected_volume-cfctvolume) > 1.0e-12) {
    std::cerr << "Expected volume: " << expected_volume << std::endl;
    std::cerr << "Computed volume by divergence theorem: " << cfctvolume <<
        std::endl;
    return false;
  }
  if (expected_volume > 0.0 && std::abs(expected_volume-ctetvolume) > 1.0e-12) {
    std::cerr << "Expected volume: " << expected_volume << std::endl;
    std::cerr << "Computed volume by summing tet volumes: " << ctetvolume <<
        std::endl;
    return false;
  }
  Portage::Point<3> dummy_point(-99, -99, -99);
  if (!Portage::approxEq(expected_centroid, dummy_point) &&
      !Portage::approxEq(expected_centroid, cfctcen, 1.0e-08)) {
    std::cerr << "Expected centroid " << expected_centroid[0] << " " <<
        expected_centroid[1] << " " << expected_centroid[2] << "\n";
    std::cerr << "Computed centroid by divergence theorem " << cfctcen[0] <<
        " " << cfctcen[1] << " " << cfctcen[2] << "\n";
    return false;
  }
  if (!Portage::approxEq(expected_centroid, ctetcen, 1.0e-08)) {
    std::cerr << "Expected centroid " << expected_centroid[0] << " " <<
        expected_centroid[1] << " " << expected_centroid[2] << "\n";
    std::cerr << "Computed centroid by vol weighted sum of tet centroids " <<
        ctetcen[0] << " " << ctetcen[1] << " " << ctetcen[2] << "\n";
    return false;
  }

  return true;
}

/*!
  @brief Unit test for equality comparisons
 */
TEST(Jali_Mesh, vdd_eq) {
    std::vector<Portage::Point<2>> a, b, c;
    a = {{0.25, 0}, {0.25, 0.25}};
    b = {{0.25, 0}, {0.25, 0.25}};
    c = {{0.25, 0.25}, {0.25, 0}};
    ASSERT_TRUE(vdd_eq(a, b));
    ASSERT_TRUE(vdd_eq(a, {{0.25, 0}, {0.25, 0.25}}));
    ASSERT_TRUE(!vdd_eq(a, {{0.25, 0}}));
    ASSERT_TRUE(!vdd_eq(a, c));
    ASSERT_TRUE(!vdd_eq(a, {{0.25, 0.25}, {0.25, 0}}));
    ASSERT_TRUE(vdd_eq(a, c, 0.5));
}

/*!
  @brief Unit test for canonical rotations of coordinates
 */
TEST(Jali_Mesh, coordinates_canonical_rotation) {
    std::vector<Portage::Point<2>> xylist, xylist_canonical;
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

    xylist = { {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75},
               {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5},
               {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25},
               {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25},
               {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25},
               {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5},
               {0.75, 0.75}, {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.5, 0.25}, {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75},
               {0.5, 0.75}, {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));

    xylist = { {0.75, 0.25}, {0.75, 0.5}, {0.75, 0.75}, {0.5, 0.75},
               {0.25, 0.75}, {0.25, 0.5}, {0.25, 0.25}, {0.5, 0.25} };
    coordinates_canonical_rotation({0.5, 0.5}, &xylist);
    ASSERT_TRUE(vdd_eq(xylist, xylist_canonical));
}

/*!
  @brief Unit test for counter-clockwise winding of coordinates
 */
TEST(Jali_Mesh_Wrapper, ccw) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);
    std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
    ASSERT_TRUE(mesh != nullptr);
    Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

    ASSERT_TRUE(mesh_wrapper.ccw({-1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(!mesh_wrapper.ccw({1, 0}, {0, 0}, {0, 1}));
    ASSERT_TRUE(!mesh_wrapper.ccw({-1, 0}, {0, 0}, {1, 0}));
    ASSERT_TRUE(mesh_wrapper.ccw({0.25, 0}, {0.5, 0}, {0.25, 0.25}));
    ASSERT_TRUE(!mesh_wrapper.ccw({0.25, 0}, {0, 0}, {0.25, 0.25}));
    ASSERT_TRUE(mesh_wrapper.ccw({0.25, 0.25}, {0, 0}, {0.25, 0}));
}

/*!
  @brief Unit test for getting dual cell coordinates
 */
TEST(Jali_Mesh_Wrapper, dual_cell_get_coordinates) {
    Jali::MeshFactory mf(MPI_COMM_WORLD);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size > 1) return;

    mf.included_entities({Jali::Entity_kind::EDGE,
                          Jali::Entity_kind::FACE,
                          Jali::Entity_kind::WEDGE,
                          Jali::Entity_kind::CORNER});
    std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);
    ASSERT_TRUE(mesh != nullptr);
    Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);
    double eps = 1e-12;

    std::vector<Portage::Point<2>> xylist;
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
TEST(Jali_Mesh_Wrapper, Get_Neighbor_Cells) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);

  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size > 1) return;

  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  ASSERT_TRUE(mesh != NULL);
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  // This is a regular mesh with 2 cells (3 nodes) in each direction

  // If we ask for the node adjacent neighbors of any cell, we should
  // get back all the other cells in the mesh

  int anycell = 3;

  // 1 (2nd argument) means retrieve only owned cells

  std::vector<int> adjcellids;
  mesh_wrapper.cell_get_node_adj_cells(anycell, Portage::ALL, &adjcellids);

  int ncells = mesh_wrapper.num_owned_cells();
  EXPECT_EQ(ncells-1, adjcellids.size()) <<
      "Not the right number of adjacent cells" << std::endl;

  for (int i = 0; i < ncells; i++) {
    if (i == anycell) continue;

    ASSERT_NE(std::find(adjcellids.begin(), adjcellids.end(), i),
              adjcellids.end()) << "Cell " << i <<
        " not found in adjacent cell list" << std::endl;
  }

  // The center node of this 3 node x 3 node x 3 node mesh corresponds
  // to a dual cell that is surrounded completely by the other dual
  // cells of the mesh. If we ask for the dual cell neighbors of the central
  // node, we should get back all the other dual cells in the mesh

  int center_node = 13;
  Portage::Point<3> cxyz;
  mesh_wrapper.node_get_coordinates(center_node, &cxyz);

  ASSERT_EQ(0.5, cxyz[0]);
  ASSERT_EQ(0.5, cxyz[1]);
  ASSERT_EQ(0.5, cxyz[2]);

  std::vector<int> adjdualcellids;

  // 1 (2nd argument) means owned dual cells
  mesh_wrapper.dual_cell_get_node_adj_cells(center_node, Portage::ALL,
                                            &adjdualcellids);

  // List of adjacent dual cell ids should contain all node ids except
  // the center one

  int nnodes = mesh_wrapper.num_owned_nodes();

  EXPECT_EQ(nnodes-1, adjdualcellids.size()) <<
      "Not the right number of adjacent dual cells" << std::endl;

  for (int i = 0; i < nnodes; i++) {
    if (i == center_node) continue;

    ASSERT_NE(std::find(adjdualcellids.begin(), adjdualcellids.end(), i),
              adjdualcellids.end()) << "Dual cell " << i <<
        " not found in adjacent dual cell list" << std::endl;
  }
}

/*!
  @brief Unit test for getting entities on the exterior boundary
 */
TEST(Jali_Mesh_Wrapper, Get_Exterior_Flag) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  ASSERT_TRUE(mesh != NULL);
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  int nfaces = mesh_wrapper.num_entities(Portage::Entity_kind::FACE,
                                         Portage::Entity_type::ALL);
  for (int f = 0; f < nfaces; f++) {
    std::vector<int> fcells;
    mesh_wrapper.face_get_cells(f, Portage::Entity_type::ALL, &fcells);
    if (mesh_wrapper.on_exterior_boundary(Portage::Entity_kind::FACE, f))
      ASSERT_EQ(1, fcells.size());
    else
      ASSERT_EQ(2, fcells.size());
  }

  int ncells = mesh_wrapper.num_entities(Portage::Entity_kind::CELL,
                                         Portage::Entity_type::ALL);
  for (int c = 0; c < ncells; c++) {
    // if this is an exterior cell, it must have at least one exterior face
    std::vector<int> cfaces;
    std::vector<int> cfdirs;
    mesh_wrapper.cell_get_faces_and_dirs(c, &cfaces, &cfdirs);
    bool exterior_face_found = false;
    for (auto const &f : cfaces) {
      if (mesh_wrapper.on_exterior_boundary(Portage::Entity_kind::FACE, f)) {
        exterior_face_found = true;
        break;
      }
    }
    if (mesh_wrapper.on_exterior_boundary(Portage::Entity_kind::CELL, c))
      ASSERT_TRUE(exterior_face_found);
    else
      ASSERT_TRUE(!exterior_face_found);
  }

  int nnodes = mesh_wrapper.num_entities(Portage::Entity_kind::NODE,
                                         Portage::Entity_type::ALL);
  for (int n = 0; n < nnodes; n++) {
    // If this is an exterior node, it must have at least one exterior cell
    // If its an interior node, it may still have an exterior cell
    // coonnected to it - so we can't check for that. We would have to check
    // for faces connected to the node but thats not done here
    std::vector<int> nodecorners;
    mesh_wrapper.node_get_corners(n, Portage::Entity_type::ALL, &nodecorners);
    bool exterior_cell_found = false;
    for (auto const &cn : nodecorners) {
      int c = mesh_wrapper.corner_get_cell(cn);
      if (mesh_wrapper.on_exterior_boundary(Portage::Entity_kind::CELL, c))
        exterior_cell_found = true;
    }
    if (mesh_wrapper.on_exterior_boundary(Portage::Entity_kind::NODE, n))
      ASSERT_TRUE(exterior_cell_found);
  }
}

/*!
  @brief Unit test for 5-tet hex decomposition
 */
TEST(Jali_Mesh_Wrapper, Decompose_Cell_Into_Tets) {
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::EDGE,
          Jali::Entity_kind::FACE,
          Jali::Entity_kind::WEDGE,
          Jali::Entity_kind::CORNER});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
  ASSERT_NE(mesh, nullptr);
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  std::vector<std::array<Portage::Point<3>, 4>> tcoords;

  // The standard decomposition has 24 tets:
  mesh_wrapper.decompose_cell_into_tets(0, &tcoords, false);
  EXPECT_EQ(tcoords.size(), 24);

  // The special decomposition has 5 tets:
  tcoords.clear();
  mesh_wrapper.decompose_cell_into_tets(0, &tcoords, true);
  EXPECT_EQ(tcoords.size(), 5);

}


/*!  @brief Unit test for 2D sides construction and queries in the
  wrapper class (and not natively in Jali).

  For this we only need cells (always present), faces and nodes
  (always present) from Jali
*/
TEST(Jali_Mesh_Wrapper, MESH_SIDES_2D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  int nsides_owned = mesh_wrapper.num_entities(Portage::SIDE,
                                               Portage::PARALLEL_OWNED);
  int nsides_ghost = mesh_wrapper.num_entities(Portage::SIDE,
                                               Portage::PARALLEL_GHOST);
  ASSERT_GT(nsides_owned, 0);
  if (nproc > 1)
    ASSERT_TRUE(nsides_ghost);
  else
    ASSERT_TRUE(!nsides_ghost);

  double dp;

  int ncells = mesh_wrapper.num_entities(Portage::CELL, Portage::ALL);
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> csides;
    mesh_wrapper.cell_get_sides(c, &csides);

    // Quad elements have 4 sides

    ASSERT_EQ(csides.size(), 4);

    double cellvol = mesh->cell_volume(c);

    for (auto const & s : csides) {

      // Since the side came from the cell 'c' the cell of the
      // side must be 'c'

      int wc = mesh_wrapper.side_get_cell(s);
      ASSERT_EQ(c, wc);

      Portage::Point<2> ccen;
      mesh_wrapper.cell_centroid(c, &ccen);

      // Make sure the side knows which face its associated with

      int f = mesh_wrapper.side_get_face(s);
      ASSERT_GE(f, 0);
      Portage::Point<2> fcen;
      mesh_wrapper.face_centroid(f, &fcen);

      // Make sure the side knows which nodes its associated with

      int n0 = mesh_wrapper.side_get_node(s, 0);
      ASSERT_GE(n0, 0);
      Portage::Point<2> npnt0;
      mesh_wrapper.node_get_coordinates(n0, &npnt0);

      // Make sure the side knows which nodes its associated with

      int n1 = mesh_wrapper.side_get_node(s, 1);
      ASSERT_GE(n1, 0);
      Portage::Point<2> npnt1;
      mesh_wrapper.node_get_coordinates(n1, &npnt1);

      // Make sure the link between the side and the opposite side
      // are correct

      int s2 = mesh_wrapper.side_get_opposite_side(s);
      ASSERT_GE(s2, -1);
      if (s2 != -1) {
        ASSERT_EQ(s, mesh_wrapper.side_get_opposite_side(s2));
        ASSERT_EQ(f, mesh_wrapper.side_get_face(s2));
        ASSERT_EQ(n0, mesh_wrapper.side_get_node(s2, 1));
        ASSERT_EQ(n1, mesh_wrapper.side_get_node(s2, 0));
        ASSERT_NE(c, mesh_wrapper.side_get_cell(s2));

        // Get side coordinates and make sure they match up with the
        // expected coordinates (node, edge center, face center, cell
        // center)

        std::array<Portage::Point<2>, 3> scoords;
        mesh_wrapper.side_get_coordinates(s, &scoords);

        for (int i = 0; i < 2; ++i) {
          ASSERT_EQ(scoords[0][i], npnt0[i]);
          ASSERT_EQ(scoords[1][i], npnt1[i]);
          ASSERT_EQ(scoords[2][i], ccen[i]);
        }

        // Since there are 4 sides in a quad element, its volume should
        // be 1/4th that of the cell

        double volume = mesh_wrapper.side_volume(s);
        ASSERT_NEAR(volume, cellvol/csides.size(), 1.0e-06);

        // Now get the side coordinate in the positive volume order
        std::array<Portage::Point<2>, 3> scoords2;
        mesh_wrapper.side_get_coordinates(s, &scoords2, true);

        if (scoords[0][0] == scoords2[1][0] &&
            scoords[0][1] == scoords2[1][1] &&
            scoords[1][0] == scoords2[0][0] &&
            scoords[1][1] == scoords2[0][1]) {

          // coordinates are flipped - verify that the coordinate flipping
          // is warranted by computing the side volume using the natural
          // ordering and checking that its the opposite sign of the volume
          // returned by side_volume operator

          Portage::Vector<2> vec0 = scoords[1]-scoords[0];
          Portage::Vector<2> vec1 = scoords[2]-scoords[0];
          double cpvec = cross(vec0, vec1);
          double altvolume = cpvec;
          ASSERT_TRUE(altvolume*volume < 0.0);

          vec0 = scoords2[1]-scoords2[0];
          vec1 = scoords2[2]-scoords2[0];
          cpvec = cross(vec0, vec1);
          altvolume = cpvec;
          ASSERT_TRUE(altvolume*volume > 0.0);
        }
      }  // for (s : csides)
    }
  }

}

/*!  @brief Unit test for 3D sides construction and queries in the
  wrapper class (and not natively in Jali).

  For this we only need cells (always present), faces and nodes
  (always present) from Jali
*/
TEST(Jali_Mesh_Wrapper, MESH_SIDES_3D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE});

  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  int ncells = mesh_wrapper.num_entities(Portage::CELL, Portage::ALL);
  double dp;
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> csides;
    mesh_wrapper.cell_get_sides(c, &csides);

    // Hex elements have 24 sides

    ASSERT_EQ(24, csides.size());

    double cellvol = mesh_wrapper.cell_volume(c);

    for (auto const & s : csides) {

      // Since the side came from the cell 'c' the cell of of the
      // side must be 'c'

      int sc = mesh_wrapper.side_get_cell(s);
      ASSERT_EQ(c, sc);

      Portage::Point<3> ccen;
      mesh_wrapper.cell_centroid(c, &ccen);

      // Make sure the side knows which face its associated with

      int f = mesh_wrapper.side_get_face(s);
      ASSERT_GE(f, 0);
      Portage::Point<3> fcen;
      mesh_wrapper.face_centroid(f, &fcen);

      // Make sure the side knows which nodes its associated with

      int n0 = mesh_wrapper.side_get_node(s, 0);
      ASSERT_GE(n0, 0);
      Portage::Point<3> npnt0;
      mesh_wrapper.node_get_coordinates(n0, &npnt0);

      int n1 = mesh_wrapper.side_get_node(s, 1);
      ASSERT_GE(n1, 0);
      Portage::Point<3> npnt1;
      mesh_wrapper.node_get_coordinates(n1, &npnt1);

      // Make sure the link between the side and the opposite side
      // is correct

      int s2 = mesh_wrapper.side_get_opposite_side(s);
      ASSERT_GE(s2, -1);
      if (s2 != -1) {
        ASSERT_EQ(s, mesh_wrapper.side_get_opposite_side(s2));
        ASSERT_EQ(f, mesh_wrapper.side_get_face(s2));
        ASSERT_EQ(n0, mesh_wrapper.side_get_node(s2, 1));
        ASSERT_EQ(n1, mesh_wrapper.side_get_node(s2, 0));
        ASSERT_TRUE(c != mesh_wrapper.side_get_cell(s2));
      }

      // Get side coordinates in the natural order and make sure
      // they match up with the expected coordinates (node 0, node
      // 1, face center, cell center)

      std::array<Portage::Point<3>, 4> scoords;
      mesh_wrapper.side_get_coordinates(s, &scoords);

      for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(scoords[0][i], npnt0[i]);
        ASSERT_EQ(scoords[1][i], npnt1[i]);
        ASSERT_EQ(scoords[2][i], fcen[i]);
        ASSERT_EQ(scoords[3][i], ccen[i]);
      }

      // Since there are 24 sides in a hex element, and the hex is
      // cubical, each side volume should be 1/24th that of the cell

      double volume = mesh_wrapper.side_volume(s);
      ASSERT_NEAR(volume, cellvol/csides.size(), 1.0e-06);

      // Now get the side coordinate in the positive volume order
      std::array<Portage::Point<3>, 4> scoords2;
      mesh_wrapper.side_get_coordinates(s, &scoords2, true);

      if (scoords[0][0] == scoords2[1][0] &&
          scoords[0][1] == scoords2[1][1] &&
          scoords[0][2] == scoords2[1][2] &&
          scoords[1][0] == scoords2[0][0] &&
          scoords[1][1] == scoords2[0][1] &&
          scoords[1][2] == scoords2[0][2]) {
        // coordinates are flipped - verify that the coordinate flipping
        // is warranted by computing the side volume using the natural
        // ordering and checking that its the opposite sign of the volume
        // returned by side_volume operator

        Portage::Vector<3> vec0 = scoords[1]-scoords[0];
        Portage::Vector<3> vec1 = scoords[2]-scoords[0];
        Portage::Vector<3> vec2 = scoords[3]-scoords[0];
        Portage::Vector<3> cpvec = cross(vec0, vec1);
        double altvolume = dot(cpvec, vec2);
        ASSERT_TRUE(altvolume < 0.0);

        // Check that we get a +ve volume using the corrected ordering

        vec0 = scoords2[1]-scoords2[0];
        vec1 = scoords2[2]-scoords2[0];
        vec2 = scoords2[3]-scoords2[0];
        cpvec = cross(vec0, vec1);
        altvolume = dot(cpvec, vec2);
        ASSERT_TRUE(altvolume > 0.0);
      }
    }
  }

}



/*!  @brief Unit test for 2D wedges construction and queries in the
  wrapper class (and not natively in Jali).

  For this we only need cells (always present), faces and nodes
  (always present) from Jali
*/
TEST(Jali_Mesh_Wrapper, MESH_WEDGES_2D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  // Create the mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  int nwedges_owned = mesh_wrapper.num_entities(Portage::WEDGE,
                                                Portage::PARALLEL_OWNED);
  int nwedges_ghost = mesh_wrapper.num_entities(Portage::WEDGE,
                                                Portage::PARALLEL_GHOST);
  ASSERT_GT(nwedges_owned, 0);
  if (nproc > 1)
    ASSERT_TRUE(nwedges_ghost);
  else
    ASSERT_TRUE(!nwedges_ghost);


  double dp;
  int ncells = mesh_wrapper.num_entities(Portage::CELL, Portage::ALL);
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cwedges;
    mesh_wrapper.cell_get_wedges(c, &cwedges);

    // Quad elements have 8 wedges

    ASSERT_EQ(cwedges.size(), 8);

    double cellvol = mesh_wrapper.cell_volume(c);

    for (auto const & w : cwedges) {

      // Since the wedge came from the cell 'c' the cell of of the
      // wedge must be 'c'

      int wc = mesh_wrapper.wedge_get_cell(w);
      ASSERT_EQ(c, wc);

      Portage::Point<2> ccen;
      mesh_wrapper.cell_centroid(c, &ccen);

      // Make sure the wedge knows which face its associated with

      int f = mesh_wrapper.wedge_get_face(w);
      ASSERT_GE(f, 0);
      Portage::Point<2> fcen;
      mesh_wrapper.face_centroid(f, &fcen);

      // Make sure the wedge knows which node its associated with

      int n = mesh_wrapper.wedge_get_node(w);
      ASSERT_GE(n, 0);
      Portage::Point<2> npnt;
      mesh_wrapper.node_get_coordinates(n, &npnt);

      // Make sure the link between the wedge and the opposite wedge
      // are correct

      int w2 = mesh_wrapper.wedge_get_opposite_wedge(w);
      ASSERT_GE(w2, -1);
      if (w2 != -1) {
        ASSERT_EQ(w, mesh_wrapper.wedge_get_opposite_wedge(w2));
        ASSERT_EQ(f, mesh_wrapper.wedge_get_face(w2));
        //        ASSERT_EQ(e, mesh_wrapper.wedge_get_edge(w2));
        ASSERT_EQ(n, mesh_wrapper.wedge_get_node(w2));
        ASSERT_TRUE(c != mesh_wrapper.wedge_get_cell(w2));
      }

      // Make sure the link between the edge and the adjacent wedge
      // are correct

      int w3 = mesh_wrapper.wedge_get_adjacent_wedge(w);
      ASSERT_EQ(w, mesh_wrapper.wedge_get_adjacent_wedge(w3));
      ASSERT_EQ(f, mesh_wrapper.wedge_get_face(w3));
      //      ASSERT_EQ(e, mesh_wrapper.wedge_get_edge(w3));
      ASSERT_EQ(c, mesh_wrapper.wedge_get_cell(w3));

      // Get the opposite wedge (w4) of the adjacent wedge (w3) of w
      // and make sure that it is adjacent to the opposite wedge
      // (w2) of w.

      int w4 = mesh_wrapper.wedge_get_opposite_wedge(w3);
      ASSERT_GE(w4, -1);
      if (w4 != -1 && w2 != -1)
        ASSERT_EQ(w4, mesh_wrapper.wedge_get_adjacent_wedge(w2));

      // Get wedge coordinates and make sure they match up with the
      // expected coordinates (node, edge center, face center, cell
      // center)

      std::array<Portage::Point<2>, 3> wcoords;
      mesh_wrapper.wedge_get_coordinates(w, &wcoords);

      for (int i = 0; i < 2; ++i) {
        ASSERT_EQ(wcoords[0][i], npnt[i]);
        ASSERT_EQ(wcoords[1][i], fcen[i]);
        ASSERT_EQ(wcoords[2][i], ccen[i]);
      }

      // Since there are 8 wedges in a quad element, its volume should
      // be 1/8th that of the cell

      double volume = mesh_wrapper.wedge_volume(w);
      ASSERT_NEAR(volume, cellvol/cwedges.size(), 1.0e-06);

      // Now get the wedge coordinate in the positive volume order
      std::array<Portage::Point<2>, 3> wcoords2;
      mesh_wrapper.wedge_get_coordinates(w, &wcoords2, true);

      bool flipped = true;

      if (wcoords[1][0] == wcoords2[2][0] &&
          wcoords[1][1] == wcoords2[2][1] &&
          wcoords[2][0] == wcoords2[1][0] &&
          wcoords[2][1] == wcoords2[1][1]) {

        // coordinates are flipped - verify that the coordinate flipping
        // is warranted by computing the wedge volume using the natural
        // ordering and checking that its the opposite sign of the volume
        // returned by wedge_volume operator

        Portage::Vector<2> vec0 = wcoords[1]-wcoords[0];
        Portage::Vector<2> vec1 = wcoords[2]-wcoords[0];
        double cpvec = cross(vec0, vec1);
        double altvolume = cpvec;
        ASSERT_TRUE(altvolume*volume < 0.0);

        vec0 = wcoords2[1]-wcoords2[0];
        vec1 = wcoords2[2]-wcoords2[0];
        cpvec = cross(vec0, vec1);
        altvolume = cpvec;
        ASSERT_TRUE(altvolume*volume > 0.0);
      }
    }  // for (w : cwedges)
  }
}

/*!  @brief Unit test for 3D wedges construction and queries in the
  wrapper class (and not natively in Jali).

  For this we only need cells (always present), faces and nodes
  (always present) from Jali
*/
TEST(Jali_Mesh_Wrapper, MESH_WEDGES_3D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  // Create the mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  double dp;
  int ncells = mesh_wrapper.num_entities(Portage::CELL, Portage::ALL);
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> cwedges;
    mesh_wrapper.cell_get_wedges(c, &cwedges);

    // Hex elements have 48 wedges

    ASSERT_EQ(cwedges.size(), 48);

    double cellvol = mesh_wrapper.cell_volume(c);

    for (auto const & w : cwedges) {

      // Since the wedge came from the cell 'c' the cell of of the
      // wedge must be 'c'

      int wc = mesh_wrapper.wedge_get_cell(w);
      ASSERT_EQ(c, wc);

      Portage::Point<3> ccen;
      mesh_wrapper.cell_centroid(c, &ccen);

      // Make sure the wedge knows which face its associated with

      int f = mesh_wrapper.wedge_get_face(w);
      ASSERT_GE(f, 0);
      Portage::Point<3> fcen;
      mesh_wrapper.face_centroid(f, &fcen);

      // Make sure the wedge knows which node its associated with

      int n = mesh_wrapper.wedge_get_node(w);
      ASSERT_GE(n, 0);
      Portage::Point<3> npnt;
      mesh_wrapper.node_get_coordinates(n, &npnt);

      // Make sure the link between the wedge and the opposite wedge
      // are correct

      int w2 = mesh_wrapper.wedge_get_opposite_wedge(w);
      ASSERT_GE(w2, -1);
      if (w2 != -1) {
        ASSERT_EQ(w, mesh_wrapper.wedge_get_opposite_wedge(w2));
        ASSERT_EQ(f, mesh_wrapper.wedge_get_face(w2));
        //        ASSERT_EQ(e, mesh_wrapper.wedge_get_edge(w2));
        ASSERT_EQ(n, mesh_wrapper.wedge_get_node(w2));
        ASSERT_TRUE(c != mesh_wrapper.wedge_get_cell(w2));
      }

      // Make sure the link between the edge and the adjacent wedge
        // are correct

      int w3 = mesh_wrapper.wedge_get_adjacent_wedge(w);
      ASSERT_EQ(w, mesh_wrapper.wedge_get_adjacent_wedge(w3));
      ASSERT_EQ(f, mesh_wrapper.wedge_get_face(w3));
      //      ASSERT_EQ(e, mesh_wrapper.wedge_get_edge(w3));
      ASSERT_EQ(c, mesh_wrapper.wedge_get_cell(w3));

      // Get the opposite wedge (w4) of the adjacent wedge (w3) of w
      // and make sure that it is adjacent to the opposite wedge
      // (w2) of w.

      int w4 = mesh_wrapper.wedge_get_opposite_wedge(w3);
      ASSERT_GE(w4, -1);
      if (w4 != -1 && w2 != -1)
        ASSERT_EQ(w4, mesh_wrapper.wedge_get_adjacent_wedge(w2));

      // Get wedge coordinates in the natural order and make sure
      // they match up with the expected coordinates (node, edge
      // center, face center, cell center)

      std::array<Portage::Point<3>, 4> wcoords;
      mesh_wrapper.wedge_get_coordinates(w, &wcoords);

      for (int i = 0; i < 3; ++i) {
        ASSERT_EQ(wcoords[0][i], npnt[i]);
        //        ASSERT_EQ(wcoords[1][i], ecen[i]);
        ASSERT_EQ(wcoords[2][i], fcen[i]);
        ASSERT_EQ(wcoords[3][i], ccen[i]);
      }

      // Since there are 48 wedges in a hex element, its volume should
      // be 1/48th that of the cell

      double volume = mesh_wrapper.wedge_volume(w);
      ASSERT_NEAR(volume, cellvol/cwedges.size(), 1.0e-06);

      // Now get the wedge coordinate in the positive volume order
      std::array<Portage::Point<3>, 4> wcoords2;
      mesh_wrapper.wedge_get_coordinates(w, &wcoords2, true);

      if (wcoords[1][0] == wcoords2[2][0] &&
          wcoords[1][1] == wcoords2[2][1] &&
          wcoords[1][2] == wcoords2[2][2] &&
          wcoords[2][0] == wcoords2[1][0] &&
          wcoords[2][1] == wcoords2[1][1] &&
          wcoords[2][2] == wcoords2[1][2]) {
        // coordinates are flipped - verify that the coordinate flipping
        // is warranted by computing the wedge volume using the natural
        // ordering and checking that its the opposite sign of the volume
        // returned by wedge_volume operator

        Portage::Vector<3> vec0 = wcoords[1]-wcoords[0];
        Portage::Vector<3> vec1 = wcoords[2]-wcoords[0];
        Portage::Vector<3> vec2 = wcoords[3]-wcoords[0];
        Portage::Vector<3> cpvec = cross(vec0, vec1);
        double altvolume = dot(cpvec, vec2);
        ASSERT_TRUE(altvolume < 0.0);

        // Check that we get a +ve volume using the corrected ordering

        vec0 = wcoords2[1]-wcoords2[0];
        vec1 = wcoords2[2]-wcoords2[0];
        vec2 = wcoords2[3]-wcoords2[0];
        cpvec = cross(vec0, vec1);
        altvolume = dot(cpvec, vec2);
        ASSERT_TRUE(altvolume > 0.0);
      }
    }
  }

}


/*!  @brief Unit test for 2D corners construction and queries in the
  wrapper class (and not natively in Jali).

  For this we only need cells (always present), faces and nodes
  (always present) from Jali
*/
TEST(Jali_Mesh_Wrapper, MESH_CORNERS_2D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  // Create the mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);

  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  int ncorners_owned = mesh_wrapper.num_entities(Portage::CORNER,
                                                 Portage::PARALLEL_OWNED);
  int ncorners_ghost = mesh_wrapper.num_entities(Portage::CORNER,
                                                 Portage::PARALLEL_GHOST);
  ASSERT_GT(ncorners_owned, 0);
  if (nproc > 1)
    ASSERT_TRUE(ncorners_ghost);
  else
    ASSERT_TRUE(!ncorners_ghost);

  double totalvol = 0.0;  // total volume of domain
  int ncells = mesh_wrapper.num_entities(Portage::CELL, Portage::ALL);
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> ccorners;
    mesh_wrapper.cell_get_corners(c, &ccorners);

    // Quad elements have 4 corners

    ASSERT_EQ(4, ccorners.size());

    double cellvol = mesh_wrapper.cell_volume(c);
    totalvol += cellvol;

    for (auto const & cn : ccorners) {

      // Since the corner came from the cell 'c' the cell of of the
      // corner must be 'c'

      int wc = mesh_wrapper.corner_get_cell(cn);
      ASSERT_EQ(c, wc);

      Portage::Point<2> ccen;
      mesh_wrapper.cell_centroid(c, &ccen);

      // Make sure the corner knows which node its associated with

      int n = mesh_wrapper.corner_get_node(cn);
      ASSERT_NE(n, -1);
      Portage::Point<2> npnt;
      mesh_wrapper.node_get_coordinates(n, &npnt);

      // Get the wedges of the corner

      std::vector<int> cnwedges;

      mesh_wrapper.corner_get_wedges(cn, &cnwedges);
      ASSERT_EQ(2, cnwedges.size());   // corner in a 2D cell has 2 wedges

      double volwedges = 0;

      for (auto const & w : cnwedges) {
        // Make sure that the node of the wedge is the same as the
        // node of the corner

        ASSERT_EQ(n, mesh_wrapper.wedge_get_node(w));

        // Add up the volume of the wedges

        volwedges += mesh_wrapper.wedge_volume(w);
      }


      double cnvolume = mesh_wrapper.corner_volume(cn);

      ASSERT_NEAR(cellvol/4.0, cnvolume, 1.0e-06);
      ASSERT_NEAR(volwedges, cnvolume, 1.0e-06);
    }  // for (cn : ccorners)


    // Cross check in a different way. Get corner of cell at each
    // node of the cell and add up the volumes of the corners
    // obtained this way. Compare to cell volume

    std::vector<int> cnodes;
    mesh_wrapper.cell_get_nodes(c, &cnodes);

    double cellvol2 = 0.0;
    for (auto const & n : cnodes) {
      int corner = mesh_wrapper.cell_get_corner_at_node(c, n);
      cellvol2 += mesh_wrapper.corner_volume(corner);
    }

    ASSERT_NEAR(cellvol, cellvol2, 1.0e-06);

  }  // for (c : mesh_wrapper.cells)

  // Now get corners of nodes, add up their volumes and make sure
  // it compares accurately to the total volume of the domain

  double totalvol2 = 0.0;
  int nnodes = mesh_wrapper.num_entities(Portage::NODE, Portage::ALL);
  for (int n = 0; n < nnodes; ++n) {
    std::vector<int> corners;
    mesh_wrapper.node_get_corners(n, Portage::ALL, &corners);

    for (auto const & cn : corners)
      totalvol2 += mesh_wrapper.corner_volume(cn);
  }

  ASSERT_NEAR(totalvol, totalvol2, 1.0e-06);
}


/*!  @brief Unit test for 3D sides construction and queries in the
  wrapper class (and not natively in Jali).

  For this we only need cells (always present), faces and nodes
  (always present) from Jali
*/
TEST(Jali_Mesh_Wrapper, MESH_CORNERS_3D) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  // Create the mesh

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);

  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  double totalvol = 0.0;  // total volume of domain
  int ncells = mesh_wrapper.num_entities(Portage::CELL, Portage::ALL);
  for (int c = 0; c < ncells; ++c) {
    std::vector<int> ccorners;
    mesh_wrapper.cell_get_corners(c, &ccorners);

    // Hex elements have 8 corners

    ASSERT_EQ(8, ccorners.size());

    double cellvol = mesh_wrapper.cell_volume(c);
    totalvol += cellvol;

    for (auto const & cn : ccorners) {

      // Since the corner came from the cell 'c' the cell of of the
      // corner must be 'c'

      int wc = mesh_wrapper.corner_get_cell(cn);
      ASSERT_EQ(c, wc);

      Portage::Point<3> ccen;
      mesh_wrapper.cell_centroid(c, &ccen);

      // Make sure the corner knows which node its associated with

      int n = mesh_wrapper.corner_get_node(cn);
      ASSERT_NE(n, -1);
      Portage::Point<3> npnt;
      mesh_wrapper.node_get_coordinates(n, &npnt);

      // Get the wedges of the corner

      std::vector<int> cnwedges;

      mesh_wrapper.corner_get_wedges(cn, &cnwedges);
      ASSERT_EQ(6, cnwedges.size());   // 6 wedges in corner at trivalent
      //                                 // node of 3D cell

      double volwedges = 0;

      for (auto const & w : cnwedges) {
        // Make sure that the node of the wedge is the same as the
        // node of the corner

        ASSERT_EQ(n, mesh_wrapper.wedge_get_node(w));

        // Add up the volume of the wedges

        volwedges += mesh_wrapper.wedge_volume(w);
      }


      double cnvolume = mesh_wrapper.corner_volume(cn);

      ASSERT_NEAR(cellvol/8.0, cnvolume, 1.0e-06);
      ASSERT_NEAR(volwedges, cnvolume, 1.0e-06);
    }  // for (cn : ccorners)

    // Cross check in a different way. Get corner of cell at each
    // node of the cell and add up the volumes of the corners
    // obtained this way. Compare to cell volume

    std::vector<int> cnodes;
    mesh_wrapper.cell_get_nodes(c, &cnodes);

    double cellvol2 = 0.0;
    for (auto const & n : cnodes) {
      int corner = mesh_wrapper.cell_get_corner_at_node(c, n);
      cellvol2 += mesh_wrapper.corner_volume(corner);
    }

    ASSERT_NEAR(cellvol, cellvol2, 1.0e-06);

  }  // for c = 0, ncells


  // Now get corners of nodes, add up their volumes and make sure
  // it compares accurately to the total volume of the domain

  double totalvol2 = 0.0;

  int nnodes = mesh_wrapper.num_entities(Portage::NODE, Portage::ALL);
  for (int n = 0; n < nnodes; ++n) {
    std::vector<int> corners;
    mesh_wrapper.node_get_corners(n, Portage::ALL, &corners);

    for (auto const & cn : corners)
      totalvol2 += mesh_wrapper.corner_volume(cn);
  }

  ASSERT_NEAR(totalvol, totalvol2, 1.0e-06);

}  // MESH_CORNERS_3D


/*!  @brief Test instantiation of Jali mesh wrapper with non-default options

  For this we only need cells (always present), faces and nodes
  (always present) from Jali
*/
TEST(Jali_Mesh_Wrapper, MESH_NON_DEFAULT_OPTS) {

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 2, 2);


  // Should get positive number of sides and 0 for wedges and corners

  bool request_sides = true;
  bool request_wedges = false;
  bool request_corners = false;

  Wonton::Jali_Mesh_Wrapper mesh_wrapper1(*mesh, request_sides, request_wedges,
                                           request_corners);

  ASSERT_GT(mesh_wrapper1.num_entities(Portage::SIDE, Portage::ALL), 0);
  ASSERT_EQ(mesh_wrapper1.num_entities(Portage::WEDGE, Portage::ALL), 0);
  ASSERT_EQ(mesh_wrapper1.num_entities(Portage::CORNER, Portage::ALL), 0);



  // Should get positive number of wedges but also sides (because wedges are
  // implicitly derived from sides) but 0 for corners

  request_wedges = true;
  request_corners = false;

  Wonton::Jali_Mesh_Wrapper mesh_wrapper2(*mesh, request_sides, request_wedges,
                                           request_corners);

  ASSERT_GT(mesh_wrapper2.num_entities(Portage::WEDGE, Portage::ALL), 0);
  ASSERT_EQ(mesh_wrapper2.num_entities(Portage::CORNER, Portage::ALL), 0);



  // Should get positive number of corners but also sides and wedges (because
  // some corner data is implicitly derived from wedges)

  request_wedges = false;
  request_corners = true;

  Wonton::Jali_Mesh_Wrapper mesh_wrapper3(*mesh, request_sides, request_wedges,
                                           request_corners);

  ASSERT_GT(mesh_wrapper3.num_entities(Portage::WEDGE, Portage::ALL), 0);
  ASSERT_GT(mesh_wrapper3.num_entities(Portage::CORNER, Portage::ALL), 0);


  // Bare bones mesh wrapper - should get zero for auxiliary entity counts

  request_wedges = false;
  request_corners = false;

  Wonton::Jali_Mesh_Wrapper mesh_wrapper4(*mesh, request_sides, request_wedges,
                                           request_corners);

  ASSERT_EQ(mesh_wrapper4.num_entities(Portage::WEDGE, Portage::ALL), 0);
  ASSERT_EQ(mesh_wrapper4.num_entities(Portage::CORNER, Portage::ALL), 0);

}



// Check that we can facetize cells and dual cells in 3D correctly

TEST(Jali_Mesh_Wrapper, MultiCell_Facetization) {
  // Create a 8x8x8 mesh
  double xmin(0.0), ymin(0.0), zmin(0.0);
  double xmax(5.0), ymax(5.0), zmax(5.0);
  int nx(5), ny(5), nz(5);

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  std::shared_ptr<Jali::Mesh> mesh = mf(xmin, ymin, zmin,
                                        xmax, ymax, zmax,
                                        nx, ny, nz);
  mesh->write_to_gmv_file("test.gmv");

  bool request_sides = true;
  bool request_wedges = true;
  bool request_corners = true;
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh, request_sides, request_wedges,
                                         request_corners);

  int nnodes = mesh_wrapper.num_owned_nodes() + mesh_wrapper.num_ghost_nodes();

  // The volume of any cell in the mesh is known to be 1
  double cvolume = 1.0;

  // Sum of areas of the cell faces is known to be 6
  double cfarea = 6.0;

  // Get the facetization of any one cell and check surface area and
  // volume (by divergence theorem). By definition, each facet is
  // planar and has only one normal

  std::vector<Portage::Point<3>> fctpoints;
  std::vector<std::vector<int>> facets;
  mesh_wrapper.cell_get_facetization(1, &facets, &fctpoints);

  Portage::Point<3> cen;
  mesh_wrapper.cell_centroid(1, &cen);

  int nfacets = facets.size();
  ASSERT_EQ(24, nfacets);

  ASSERT_TRUE(faceted_poly_ok(fctpoints, facets, cfarea, cvolume, cen, 1e-12));

  // Get the facetization of the control volume around an interior
  // node and check surface area and volume. The volume and surface
  // area should be the same as a cell

  // Find index of node in interior of domain but on partition
  // boundary without assuming any numbering scheme
  int inode = -1;
  for (int n = 0; n < nnodes; n++) {
    std::vector<int> ncells;
    mesh_wrapper.node_get_cells(n, Portage::Entity_type::ALL, &ncells);
    if (ncells.size() == 8) {
      if (nproc == 1) {
        inode = n;
        break;
      } else {
        bool owned_found = false, ghost_found = false;
        for (int i = 0; i < ncells.size(); i++) {
          Portage::Entity_type ctype = mesh_wrapper.cell_get_type(ncells[i]);
          if (ctype == Portage::Entity_type::PARALLEL_OWNED)
            owned_found = true;
          else
            ghost_found = true;
        }
        if (owned_found && ghost_found) {
          inode = n;
          break;
        }
      }
    }
  }
  ASSERT_NE(-1, inode);
  mesh_wrapper.dual_cell_get_facetization(inode, &facets, &fctpoints);

  nfacets = facets.size();
  ASSERT_EQ(48, nfacets);

  mesh_wrapper.dual_cell_centroid(inode, &cen);

  ASSERT_TRUE(faceted_poly_ok(fctpoints, facets, cfarea, cvolume, cen, 1e-12));


  // Get the facetization of the control volume around a corner node
  // and check surface area and volume. The volume should be 1/8th
  // the volume of a cell and the surface area should be 1/4th the
  // surface area of a cell.

  // Find index of corner node without assuming any numbering scheme
  inode = -1;
  for (int n = 0; n < nnodes; n++) {
    std::vector<int> ncells;
    mesh_wrapper.node_get_cells(n, Portage::Entity_type::ALL, &ncells);
    if (ncells.size() == 1) {
      inode = n;
      break;
    }
  }
  ASSERT_NE(-1, inode);
  mesh_wrapper.dual_cell_get_facetization(inode, &facets, &fctpoints);

  nfacets = facets.size();
  ASSERT_EQ(12, nfacets);

  mesh_wrapper.dual_cell_centroid(inode, &cen);

  ASSERT_TRUE(faceted_poly_ok(fctpoints, facets, cfarea/4.0, cvolume/8.0, cen,
                              1e-12));

  // Get the facetization of the control volume around a corner node
  // and check surface area and volume. The volume should be 1/8th
  // the volume of a cell and the surface area should be 1/4th the
  // surface area of a cell.

  // Find index of a node on a **domain edge** without assuming any
  // numbering scheme
  inode = -1;
  for (int n = 0; n < nnodes; n++) {
    std::vector<int> ncells;
    mesh_wrapper.node_get_cells(n, Portage::Entity_type::ALL, &ncells);
    if (ncells.size() == 2) {
      if (nproc == 1) {
        inode = n;
        break;
      } else {
        bool owned_found = false, ghost_found = false;
        for (int i = 0; i < ncells.size(); i++) {
          Portage::Entity_type ctype = mesh_wrapper.cell_get_type(ncells[i]);
          if (ctype == Portage::Entity_type::PARALLEL_OWNED)
            owned_found = true;
          else
            ghost_found = true;
        }
        if (owned_found && ghost_found) {
          inode = n;
          break;
        }
      }
    }
  }
  ASSERT_NE(-1, inode);
  mesh_wrapper.dual_cell_get_facetization(inode, &facets, &fctpoints);

  nfacets = facets.size();
  ASSERT_EQ(20, nfacets);

  mesh_wrapper.dual_cell_centroid(inode, &cen);

  ASSERT_TRUE(faceted_poly_ok(fctpoints, facets, 5*cfarea/12, cvolume/4, cen,
                              1e-12));

  // Get the facetization of the control volume around a corner node
  // and check surface area and volume. The volume should be 1/8th
  // the volume of a cell and the surface area should be 1/4th the
  // surface area of a cell.

  // Find index of node on a **domain face** without assuming any
  // numbering scheme
  inode = -1;
  for (int n = 0; n < nnodes; n++) {
    std::vector<int> ncells;
    mesh_wrapper.node_get_cells(n, Portage::Entity_type::ALL, &ncells);
    if (ncells.size() == 4) {
      if (nproc == 1) {
        inode = n;
        break;
      } else {
        bool owned_found = false, ghost_found = false;
        for (int i = 0; i < ncells.size(); i++) {
          Portage::Entity_type ctype = mesh_wrapper.cell_get_type(ncells[i]);
          if (ctype == Portage::Entity_type::PARALLEL_OWNED)
            owned_found = true;
          else
            ghost_found = true;
        }
        if (owned_found && ghost_found) {
          inode = n;
          break;
        }
      }
    }
  }
  ASSERT_NE(-1, inode);
  mesh_wrapper.dual_cell_get_facetization(inode, &facets, &fctpoints);

  nfacets = facets.size();
  ASSERT_EQ(32, nfacets);

  mesh_wrapper.dual_cell_centroid(inode, &cen);

  ASSERT_TRUE(faceted_poly_ok(fctpoints, facets, 2*cfarea/3, cvolume/2, cen,
                              1e-12));

}



// Check that we can compute volumes and centroids of skewed cells correctly

TEST(Jali_Mesh_Wrapper, Skewed_2DCell_Geometry) {
  // Create a 1 cell mesh
  double xmin(0.0), ymin(0.0);
  double xmax(1.0), ymax(1.0);
  int nx(1), ny(1);

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  std::shared_ptr<Jali::Mesh> mesh = mf(xmin, ymin, xmax, ymax, nx, ny);

  // Move the node that's at 1.0, 1.0 to 1.25, 1.25
  int nnodes = mesh->num_nodes();
  for (int i = 0; i < nnodes; i++) {
    JaliGeometry::Point ncoord;
    mesh->node_get_coordinates(i, &ncoord);
    if (fabs(ncoord[0]-1.0) < 1.0e-12 && fabs(ncoord[1]-1.0) < 1.0e-12) {
      ncoord[0] = 1.25; ncoord[1] = 1.25;
      mesh->node_set_coordinates(i, ncoord);
      break;
    }
  }

  // We can hand-calculate the area and centroid of this skewed 2D element
  double carea = 1.25;
  Portage::Point<2> cen(0.58333333333333333, 0.58333333333333333);

  bool request_sides = true;
  bool request_wedges = true;
  bool request_corners = true;
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh, request_sides, request_wedges,
                                         request_corners);

  ASSERT_NEAR(carea, mesh_wrapper.cell_volume(0), 1.0e-12);

  Portage::Point<2> cen_computed;
  mesh_wrapper.cell_centroid(0, &cen_computed);
  ASSERT_TRUE(Portage::approxEq(cen, cen_computed, 1.0e-8));
}


// Check that we can compute volumes and centroids of skewed cells correctly

TEST(Jali_Mesh_Wrapper, Skewed_3DCell_Geometry) {
  // Create a 1 cell mesh
  double xmin(0.0), ymin(0.0), zmin(0.0);
  double xmax(1.0), ymax(1.0), zmax(1.0);
  int nx(1), ny(1), nz(1);

  int nproc, me;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);

  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});
  mf.partitioner(Jali::Partitioner_type::BLOCK);
  std::shared_ptr<Jali::Mesh> mesh = mf(xmin, ymin, zmin,
                                        xmax, ymax, zmax,
                                        nx, ny, nz);

  // Move the node that's at 1.0, 1.0, 1.0 to 1.25, 1.25, 1.25
  int nnodes = mesh->num_nodes();
  for (int i = 0; i < nnodes; i++) {
    JaliGeometry::Point ncoord;
    mesh->node_get_coordinates(i, &ncoord);
    if (fabs(ncoord[0]-1.0) < 1.0e-12 && fabs(ncoord[1] - 1.0) < 1.0e-12 &&
        fabs(ncoord[2]-1.0) < 1.0e-12) {
      ncoord[0] = 1.25; ncoord[1] = 1.25; ncoord[2] = 1.25;
      mesh->node_set_coordinates(i, ncoord);
      break;
    }
  }

  bool request_sides = true;
  bool request_wedges = true;
  bool request_corners = true;
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh, request_sides, request_wedges,
                                         request_corners);

  // Get the facetization of any one cell and check surface area and
  // volume (by divergence theorem). By definition, each facet is
  // planar and has only one normal

  std::vector<Portage::Point<3>> fctpoints;
  std::vector<std::vector<int>> facets;
  mesh_wrapper.cell_get_facetization(0, &facets, &fctpoints);

  Portage::Point<3> cen;
  mesh_wrapper.cell_centroid(0, &cen);

  double cvolume = mesh_wrapper.cell_volume(0);
  
  // The mesh wrappers don't have face area computation and the Jali
  // facets use a different face point than Portage mesh wrapper. So
  // don't compare this
  double cfarea = -1.0;

  ASSERT_TRUE(faceted_poly_ok(fctpoints, facets, cfarea, cvolume, cen, 1e-12));
}
