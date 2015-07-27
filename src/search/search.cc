/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "search.h"

#include <cstdio>

namespace Portage {

void Search::search(double arg0, double arg1)
{
    std::printf("in Search::search(%g, %g)...\n", arg0, arg1);
} // Search::search

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
