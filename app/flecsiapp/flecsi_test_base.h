/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

////////////////////////////////////////////////////////////////////////////////
/// \file burton_test_base.h
/// \brief Defines a base test fixture.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "flecsi-sp.h"
#include "flecsi-sp/burton/burton.h"

#include "wonton/mesh/flecsi/flecsi_mesh_wrapper.h"
#include "wonton/state/flecsi/flecsi_state_wrapper.h"

// system includes
#include <cinchtest.h>

//! \brief the mesh type
using mesh_2d_t   = flecsi::sp::burton::burton_mesh_2d_t;
//! \brief the mesh type
using mesh_3d_t   = flecsi::sp::burton::burton_mesh_3d_t;

// some general using statements
using std::string;
using flecsi::sp::burton::write_mesh;
using flecsi::sp::burton::read_mesh;

////////////////////////////////////////////////////////////////////////////////
//! \brief base test fixture for wonton
////////////////////////////////////////////////////////////////////////////////
class flecsi_test_base : public ::testing::Test {
protected:

  //---------------------------------------------------------------------------
  //! \brief get the output prefix for files
  //---------------------------------------------------------------------------
  auto output_prefix() 
  {
    auto test_info =
      ::testing::UnitTest::GetInstance()->current_test_info();
    string test_name = test_info->name();
    string case_name = test_info->test_case_name();
    return case_name + "." + test_name;
  }

};
