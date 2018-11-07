/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/
////////////////////////////////////////////////////////////////////////////////
///
/// \file
/// 
/// \brief Defines a test fixture.
///
////////////////////////////////////////////////////////////////////////////////
#pragma once

// test include
#include "flecsi_test_base.h"

// user includes
#include "flecsi-sp/burton/factory.h"

// some general using statements
using std::vector;


////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh
////////////////////////////////////////////////////////////////////////////////
class flecsi_3d : public flecsi_test_base {
public: 

  //---------------------------------------------------------------------------
  // Types From Mesh
  //---------------------------------------------------------------------------

  //! \brief the mesh type
  using mesh_t = mesh_3d_t;

  //! \brief the wonton mesh wrapper type
  using flecsi_mesh_t = Wonton::flecsi_mesh_t<mesh_t>;
  //! \brief the wonton state wrapper typename
  using flecsi_state_t = Wonton::flecsi_state_t<mesh_t>;
  //! \brief The wonton driver type
  using wonton_1st_order_driver_t =
    Wonton::Driver<
      Wonton::SearchKDTree,
      Wonton::IntersectR3D,
      Wonton::Interpolate_1stOrder,
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
  static constexpr real_t test_tolerance = 1000*flecsi::sp::common::test_tolerance;



protected:
  
  //---------------------------------------------------------------------------
  //! \brief the test setup function
  //! \remark this function is called before each test
  //---------------------------------------------------------------------------
  virtual void SetUp() { } // SetUp

  //---------------------------------------------------------------------------
  //! \brief the test teardown function
  //! \remark this function is called after each test
  //---------------------------------------------------------------------------
  virtual void TearDown() { }



};
