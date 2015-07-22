/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "driver.h"

#include <cstdio>

#include "search.h"
#include "intersect.h"
#include "remap.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

Driver::Driver()
{
  // Create some default meshes to test with
  // 2d quad mesh from (0,0) to (1,1) with 61x61 zones
  inputMesh_ = Jali::MeshFactory::create(0.0, 0.0,
					 1.0, 1.0,
					 61, 61);
  // 2d quad mesh from (0.11, 0.17) to (0.71, 0.67) with 23x41 zones
  targetMesh_ = Jali::MeshFactory::create(0.11, 0.17,
					  0.71, 0.67,
					  23, 41);
} // Driver::Driver

void Driver::run()
{
    std::printf("in Driver::run()...\n");

    Search s;
    s.search(inputMesh_, targetMesh_);

    Intersect isect(s.candidates);
    isect.intersect();

    Remap r(isect.moments);
    r.remap();

} // Driver::run

/*-------------------------------------------------------------------------~--*
 * Formatting options for Emacs and vim.
 *
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *-------------------------------------------------------------------------~--*/
