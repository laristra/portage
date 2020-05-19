/*
This file is part of the Ristra Portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

////////////////////////////////////////////////////////////////////////////////
// \file
// \brief Demo app using  wonton flecsi mesh/state wrappers.
////////////////////////////////////////////////////////////////////////////////

// flecsi-sp includes
#include "flecsi-sp.h"
#include "flecsi-sp/burton/burton.h"
#include "flecsi-sp/burton/factory.h"

//wonton includes
#include "wonton/support/wonton.h"
#include "wonton/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "wonton/state/flecsi/flecsi_state_wrapper.h"

//portage includes
#include "portage/support/portage.h"
#include "portage/driver/mmdriver.h"

// system includes
#include <cassert>
#include <cmath>
#include <string>
#include <vector>

namespace math = flecsi::sp::math;
namespace mesh = flecsi::sp::burton;

class flecsi_2d {
public:

  //---------------------------------------------------------------------------
  // Types From Mesh
  //---------------------------------------------------------------------------

  //! \brief the mesh type
  using mesh_t = mesh::burton_mesh_2d_t;
  //! \brief the wonton mesh wrapper type
  using flecsi_mesh_t = Wonton::flecsi_mesh_t<mesh_t>;
  //! \brief the wonton state wrapper typename
  using flecsi_state_t = Wonton::flecsi_state_t<mesh_t>;
  //! \brief The portage driver type
  using portage_1st_order_driver_t =
      Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR2D,
      Portage::Interpolate_1stOrder,
      mesh_t::num_dimensions,
      flecsi_mesh_t,
      flecsi_state_t
    >;

  //! \brief the size type
  using size_t= typename mesh_t::size_t;
  //! \brief the counter type
  using counter_t= typename mesh_t::counter_t;
  //! \brief the mesh float type
  using real_t   = typename mesh_t::real_t;
  //! \brief the mesh int type
  using integer_t= typename mesh_t::integer_t;
  //! \brief the mesh dimensions
  static constexpr size_t num_dimensions = mesh_t::num_dimensions;

  //! \brief the point
  using point_t  = typename mesh_t::point_t;
  //! \brief the vector type
  using vector_t = typename mesh_t::vector_t;
  //! \brief the vertex type
  using vertex_t = typename mesh_t::vertex_t;
  //! \brief the vertex type
  using edge_t   = typename mesh_t::edge_t;
  //! \brief the cell type
  using cell_t   = typename mesh_t::cell_t;


  //---------------------------------------------------------------------------
  // Types
  //---------------------------------------------------------------------------

  //! \brief some test tolerance
  static constexpr real_t test_tolerance = 10*flecsi::sp::common::test_tolerance;
};

class flecsi_3d {
public:

  //---------------------------------------------------------------------------
  // Types From Mesh
  //---------------------------------------------------------------------------

  //! \brief the mesh type
  using mesh_t = mesh::burton_mesh_3d_t;
  //! \brief the wonton mesh wrapper type
  using flecsi_mesh_t = Wonton::flecsi_mesh_t<mesh_t>;
  //! \brief the wonton state wrapper typename
  using flecsi_state_t = Wonton::flecsi_state_t<mesh_t>;
  //! \brief The portage driver type
  using portage_1st_order_driver_t =
      Portage::MMDriver<
      Portage::SearchKDTree,
      Portage::IntersectR3D,
      Portage::Interpolate_1stOrder,
      mesh_t::num_dimensions,
      flecsi_mesh_t,
      flecsi_state_t
    >;

  //! \brief the size type
  using size_t= typename mesh_t::size_t;
  //! \brief the counter type
  using counter_t= typename mesh_t::counter_t;
  //! \brief the mesh float type
  using real_t   = typename mesh_t::real_t;
  //! \brief the mesh int type
  using integer_t= typename mesh_t::integer_t;
  //! \brief the mesh dimensions
  static constexpr size_t num_dimensions = mesh_t::num_dimensions;

  //! \brief the point
  using point_t  = typename mesh_t::point_t;
  //! \brief the vector type
  using vector_t = typename mesh_t::vector_t;
  //! \brief the vertex type
  using vertex_t = typename mesh_t::vertex_t;
  //! \brief the vertex type
  using edge_t   = typename mesh_t::edge_t;
  //! \brief the cell type
  using cell_t   = typename mesh_t::cell_t;


  //---------------------------------------------------------------------------
  // Types
  //---------------------------------------------------------------------------

  //! \brief some test tolerance
  static constexpr real_t test_tolerance = 10*flecsi::sp::common::test_tolerance;
};

void run_2d(size_t num_x, size_t num_y, std::string& output_prefix)
{
  // width of mesh
  constexpr size_t len_x = 1.;
  // length of mesh
  constexpr size_t len_y = 1.;

  auto mesh_a = mesh::box<flecsi_2d::mesh_t>( num_x, num_y, 0, 0, len_x, len_y );
  auto mesh_b = mesh::box<flecsi_2d::mesh_t>( num_x*2, num_y*2, 0, 0, len_x, len_y );

  flecsi_2d::flecsi_mesh_t mesh_wrapper_a( mesh_a );
  flecsi_2d::flecsi_mesh_t mesh_wrapper_b( mesh_b );

  flecsi_2d::flecsi_state_t state_wrapper_a( mesh_a );
  flecsi_2d::flecsi_state_t state_wrapper_b( mesh_b );

  // create the remapper object
  flecsi_2d::portage_1st_order_driver_t remapper( 
    mesh_wrapper_a, 
    state_wrapper_a, 
    mesh_wrapper_b, 
    state_wrapper_b
  );

  // register data on the source mesh
  flecsi_register_data(mesh_a, hydro, cell_data, flecsi_2d::real_t, dense, 1, cells);
  flecsi_register_data(mesh_b, hydro, cell_data, flecsi_2d::real_t, dense, 1, cells);

  auto a = flecsi_get_accessor(mesh_a, hydro, cell_data, flecsi_2d::real_t, dense, 0);
  auto b = flecsi_get_accessor(mesh_b, hydro, cell_data, flecsi_2d::real_t, dense, 0);

  a.attributes().set(persistent);
  b.attributes().set(persistent);

  // set some initial distribution in mesh A
  for ( auto c : mesh_a.cells() ) a[c] = 1.0;
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Declare which variables are remapped
  std::vector<std::string> var_names(1, "cell_data");
  remapper.set_remap_var_names(var_names);

  // Do the remap
  remapper.run();  // executor argument defaults to nullptr -> serial run

  // the result should the same
  for ( auto c : mesh_b.cells() )
    assert( std::fabs(b[c]- 1.0) <= flecsi_2d::test_tolerance );
  
  // now try something more complicated

  // set some initial distribution in mesh A
  flecsi_2d::real_t total_a = 0;
  for ( auto c : mesh_a.cells() ) {
    const auto & xc = c->centroid();
    auto val = math::sqr(xc[0]) + math::sqr(xc[1]);
    a[c] = val;
    total_a += c->volume() * val;
  }
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Do the remap
  remapper.run();

  // Check the result
  flecsi_2d::real_t total_b = 0;
  for ( auto c : mesh_b.cells() ) 
    total_b += c->volume() * b[c];
 
  assert(std::fabs(total_a-total_b) <= flecsi_2d::test_tolerance);

  // write out results  
  mesh::write_mesh(output_prefix+"_a.vtk", mesh_a, false);
  mesh::write_mesh(output_prefix+"_b.vtk", mesh_b, false);
};


void run_3d(size_t num_x, size_t num_y, size_t num_z, std::string& output_prefix)
{
 // width of mesh
  constexpr size_t len_x = 1.;
  // length of mesh
  constexpr size_t len_y = 1.;
  // depth of mesh
  constexpr size_t len_z = 1.;

  auto mesh_a = mesh::box<flecsi_3d::mesh_t>( num_x, num_y, num_z,  0, 0, 0, len_x, len_y, len_z );
  auto mesh_b = mesh::box<flecsi_3d::mesh_t>( num_x*2, num_y*2, num_z*2, 0, 0, 0, len_x, len_y, len_z );

  flecsi_3d::flecsi_mesh_t mesh_wrapper_a( mesh_a );
  flecsi_3d::flecsi_mesh_t mesh_wrapper_b( mesh_b );

  flecsi_3d::flecsi_state_t state_wrapper_a( mesh_a );
  flecsi_3d::flecsi_state_t state_wrapper_b( mesh_b );

  // create the remapper object
  flecsi_3d::portage_1st_order_driver_t remapper(
    mesh_wrapper_a,
    state_wrapper_a,
    mesh_wrapper_b,
    state_wrapper_b
  );

  // register data on the source mesh
  flecsi_register_data(mesh_a, hydro, cell_data, flecsi_3d::real_t, dense, 1, cells);
  flecsi_register_data(mesh_b, hydro, cell_data, flecsi_3d::real_t, dense, 1, cells);

  auto a = flecsi_get_accessor(mesh_a, hydro, cell_data, flecsi_3d::real_t, dense, 0);
  auto b = flecsi_get_accessor(mesh_b, hydro, cell_data, flecsi_3d::real_t, dense, 0);

  a.attributes().set(persistent);
  b.attributes().set(persistent);

  // set some initial distribution in mesh A
  for ( auto c : mesh_a.cells() ) a[c] = 1.0;
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Declare which variables are remapped
  std::vector<std::string> var_names(1, "cell_data");
  remapper.set_remap_var_names(var_names);

  // Do the remap
  remapper.run();

 // the result should the same
  for ( auto c : mesh_b.cells() ) {
    auto val = b[c];
    assert( std::fabs(val - 1.0) <= flecsi_3d::test_tolerance );
  }

  // now try something more complicated

  // set some initial distribution in mesh A
  flecsi_3d::real_t total_a = 0;
  for ( auto c : mesh_a.cells() ) {
    const auto & xc = c->centroid();
    auto val = math::sqr(xc[0]) + math::sqr(xc[1]) + math::sqr(xc[2]);
    a[c] = val;
    total_a += c->volume() * val;
  }
  for ( auto c : mesh_b.cells() ) b[c] = 0.0;

  // Do the remap
  remapper.run();

  // Check the result
  flecsi_3d::real_t total_b = 0;
  for ( auto c : mesh_b.cells() )
    total_b += c->volume() * b[c];

  assert(std::fabs(total_a-total_b) <= flecsi_3d::test_tolerance);

  // write out results
  mesh::write_mesh(output_prefix+"_a.vtk", mesh_a, false);
  mesh::write_mesh(output_prefix+"_b.vtk", mesh_b, false);

};


////////////////////////////////////////////////////////////////////////////////
//! \brief dump the mesh to std out
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
#ifdef WONTON_ENABLE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");
#endif
 
  if ((argc < 5) || (argc > 6)) {
      std::ostringstream os;
      os << std::endl <<
      "Correct usage: flecsiapp <dim> <nx> <ny> <output_filename> or" 
      " flecsiapp <dim> <nx> <ny> <nz> <output_filename>"  << std::endl;
      throw std::runtime_error(os.str());
  }
  
  int dim = atoi(argv[1]); 

  if (dim == 2)
  { 
    // number of cells wide
    size_t num_x = atoi(argv[2]);
    // number of cells high
    size_t num_y = atoi(argv[3]);
    std::string output_prefix = argv[4];

    run_2d(num_x, num_y, output_prefix);
  }
  else if (dim == 3)
  {
    // number of cells wide
    size_t nx = atoi(argv[2]);
    // number of cells high
    size_t ny = atoi(argv[3]);
    // number of cells deep
    size_t nz = atoi(argv[4]);
    std::string output_prefix = argv[5];

    run_3d(nx, ny, nz, output_prefix);
  }
  else
    std::cerr<<"Dimensions supported are 2 and 3 !"<<std::endl; 
}

/*#ifdef HAVE_EXODUS

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
  wonton_1st_order_driver_t remapper( 
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
*/
