/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "driver.h"

#include <cstdio>

#include "portage/search/search.h"
#include "portage/intersect/intersect.h"
#include "portage/remap/remap.h"

#include "Mesh.hh"
#include "MeshFactory.hh"


namespace Portage {

void Driver::run()
{
    std::printf("in Driver::run()...\n");

    // Search for possible intersections.
//    Search s;
//    s.search(0.0, 0.0);

    // Calculate the overlap of actual intersections.
    Intersect isect(sourceMesh_, targetMesh_);
    isect.intersect();

    // Remap from sourceMesh_ to targetMesh_
    Remap r(sourceMesh_, sourceState_, targetMesh_, targetState_);
    r.remap(remap_var_names_);

} // Driver::run

} // namespace Portage

/*-------------------------------------------------------------------------~--*
 * Formatting options for Emacs and vim.
 *
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *-------------------------------------------------------------------------~--*/
