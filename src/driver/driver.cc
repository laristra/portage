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

void Driver::run()
{
    std::printf("in Driver::run()...\n");

    Search s;
    s.search(0.5, 0.866);

    Intersect isect;
    isect.intersect();

    Remap r;
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
