/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "intersect.h"

#include <cstdio>

#include "clipper.hpp"

namespace Portage {

void Intersect::intersect(Jali::Entity_ID cellId, Jali::Entity_ID_List* candidates, std::vector<float>* moments)
{
    std::printf("in Intersect::intersect()...\n");
    for (unsigned int i=0; i<candidates->size(); i++)
    {
        moments->push_back(1.0f);
    }
} // Intersect::intersect

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
