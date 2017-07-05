/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
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
#include "flecsi-sp/specializations/burton/factory.h"

// some general using statements
using std::vector;


////////////////////////////////////////////////////////////////////////////////
//! \brief test fixture for creating the mesh
////////////////////////////////////////////////////////////////////////////////
class flecsi_2d : public flecsi_test_base {
public: 

  //---------------------------------------------------------------------------
  // Types From Mesh
  //---------------------------------------------------------------------------

  //! \brief the mesh type
  using mesh_t = mesh_2d_t;

  //! \brief the portage mesh wrapper type
  using flecsi_mesh_t = wonton::flecsi_mesh_t<mesh_t>;
  //! \brief the portage state wrapper typename
  using flecsi_state_t = wonton::flecsi_state_t<mesh_t>;
  //! \brief The portage driver type
  using portage_1st_order_driver_t =
    Portage::Driver<
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
  static constexpr real_t test_tolerance = 10*flecsale::common::test_tolerance;



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
