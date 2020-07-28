/* This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#include "gtest/gtest.h"

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// portage includes
#include "portage/support/portage.h"
#include "portage/distributed/mpi_update_ghosts.h"

// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"

TEST(UpdateGhosts, CommunicationMatrix) {

}