/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
// \file
// \brief Tests general features of the portage mesh wrapper.
////////////////////////////////////////////////////////////////////////////////

// user includes
#include "flecsi_2d_test.h"

// system includes
#include <vector>

namespace math = flecsi::sp::math;
namespace mesh = flecsi::sp::burton;

////////////////////////////////////////////////////////////////////////////////
//! \brief dump the mesh to std out
////////////////////////////////////////////////////////////////////////////////
TEST_F(flecsi_2d, first_order) 
{

  // number of cells wide
  constexpr size_t num_x = 10;
  // number of cells high
  constexpr size_t num_y = 10;
  // width of mesh
  constexpr size_t len_x = 1.;
  // length of mesh
  constexpr size_t len_y = 1.;

  auto mesh_a = mesh::box<mesh_t>( num_x, num_y, 0, 0, len_x, len_y );
  auto mesh_b = mesh::box<mesh_t>( num_x*2, num_y*2, 0, 0, len_x, len_y );

  flecsi_mesh_t mesh_wrapper_a( mesh_a );
  flecsi_mesh_t mesh_wrapper_b( mesh_b );

  flecsi_state_t state_wrapper_a( mesh_a );
  flecsi_state_t state_wrapper_b( mesh_b );

  // create the remapper object
  portage_1st_order_driver_t remapper( 
    mesh_wrapper_a, 
    state_wrapper_a, 
    mesh_wrapper_b, 
    state_wrapper_b
  );

  // register data on the source mesh
  flecsi_register_data(mesh_a, hydro, cell_data, real_t, dense, 1, cells);
  flecsi_register_data(mesh_b, hydro, cell_data, real_t, dense, 1, cells);

  auto a = flecsi_get_accessor(mesh_a, hydro, cell_data, real_t, dense, 0);
  auto b = flecsi_get_accessor(mesh_b, hydro, cell_data, real_t, dense, 0);

  a.attributes().set(persistent);
  b.attributes().set(persistent);


  // set some initial distribution in mesh A
  for ( auto c : mesh_a.cells() ) a[c] = 1.0;
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Declare which variables are remapped
  std::vector<std::string> var_names(1, "cell_data");
  remapper.set_remap_var_names(var_names);

  // Do the remap
  remapper.run(false);

  // the result should the same
  for ( auto c : mesh_b.cells() ) EXPECT_NEAR( b[c], 1.0, test_tolerance );

  // now try something more complicated

  // set some initial distribution in mesh A
  real_t total_a = 0;
  for ( auto c : mesh_a.cells() ) {
    const auto & xc = c->centroid();
    auto val = math::sqr(xc[0]) + math::sqr(xc[1]);
    a[c] = val;
    total_a += c->volume() * val;
  }
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Do the remap
  remapper.run(false);

  // the result should the same
  EXPECT_FALSE( mesh::write_mesh(output_prefix()+"_a.vtk", mesh_a, false) );
  EXPECT_FALSE( mesh::write_mesh(output_prefix()+"_b.vtk", mesh_b, false) );

  // Check the result
  real_t total_b = 0;
  for ( auto c : mesh_b.cells() ) 
    total_b += c->volume() * b[c];
 
  EXPECT_NEAR( total_a, total_b, test_tolerance );

}

#ifdef HAVE_EXODUS

////////////////////////////////////////////////////////////////////////////////
//! \brief Remap between a triangular mesh and a quad mesh
////////////////////////////////////////////////////////////////////////////////
TEST_F(flecsi_2d, unstruct) 
{

  mesh_t mesh_a, mesh_b;
  mesh::read_mesh( "box-quad.g", mesh_a); 
  mesh::read_mesh( "box-tri.g",  mesh_b); 

  flecsi_mesh_t mesh_wrapper_a( mesh_a );
  flecsi_mesh_t mesh_wrapper_b( mesh_b );

  flecsi_state_t state_wrapper_a( mesh_a );
  flecsi_state_t state_wrapper_b( mesh_b );

  // create the remapper object
  portage_1st_order_driver_t remapper( 
    mesh_wrapper_a, 
    state_wrapper_a, 
    mesh_wrapper_b, 
    state_wrapper_b
  );

  // register data on the source mesh
  flecsi_register_data(mesh_a, hydro, cell_data, real_t, dense, 1, cells);
  flecsi_register_data(mesh_b, hydro, cell_data, real_t, dense, 1, cells);

  auto a = flecsi_get_accessor(mesh_a, hydro, cell_data, real_t, dense, 0);
  auto b = flecsi_get_accessor(mesh_b, hydro, cell_data, real_t, dense, 0);

  a.attributes().set(persistent);
  b.attributes().set(persistent);


  // set some initial distribution in mesh A
  for ( auto c : mesh_a.cells() ) a[c] = 1.0;
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Declare which variables are remapped
  std::vector<std::string> var_names(1, "cell_data");
  remapper.set_remap_var_names(var_names);

  // Do the remap
  remapper.run(false);

  // the result should the same
  for ( auto c : mesh_b.cells() ) EXPECT_NEAR( b[c], 1.0, test_tolerance );

  // now try something more complicated

  // set some initial distribution in mesh A
  real_t total_a = 0;
  for ( auto c : mesh_a.cells() ) {
    const auto & xc = c->centroid();
    auto val = math::sqr(xc[0]) + math::sqr(xc[1]);
    a[c] = val;
    total_a += c->volume() * val;
  }
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Do the remap
  remapper.run(false);

  // the result should the same
  EXPECT_FALSE( mesh::write_mesh(output_prefix()+"_a.vtk", mesh_a, false) );
  EXPECT_FALSE( mesh::write_mesh(output_prefix()+"_b.vtk", mesh_b, false) );

  // Check the result
  real_t total_b = 0;
  for ( auto c : mesh_b.cells() ) 
    total_b += c->volume() * b[c];
 
  EXPECT_NEAR( total_a, total_b, test_tolerance );
}

#endif // HAVE_EXODUS

