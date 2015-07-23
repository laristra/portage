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


void Driver::run()
{
    std::printf("in Driver::run()...\n");

    // Search for possible intersections.
    Search s;
    s.search(inputMesh_, targetMesh_);

    // // Calculate the overlap of actual intersections.
    Intersect isect(&s);
    isect.intersect();

    // // Remap from inputMesh_ to targetMesh_
     Remap r(&isect);
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
