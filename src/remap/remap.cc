/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "remap.h"

#include <cstdio>

namespace Portage {

double Remap::remap(std::string const & remap_var_name,
                    Jali::Entity_ID cellId, Jali::Entity_ID_List candidates, 
					std::vector<std::vector<std::vector<double> > > moments) const
{
    std::printf("in Remap::remap()...\n");

    std::vector<StateVector>::const_iterator field = sourceState_.find(remap_var_name, Jali::CELL);
    Portage::StateVector stateVector = *field;
    double value = 0.0;
    for (unsigned int j=0; j<candidates.size(); j++)
    {
       double x = *(stateVector.begin() + candidates[j]);
	   // CMM: I don't actually know what should be here as the vector of vector of vector
	   //      of double confuses my head...
       value += x * (moments[j][0][0]); 
    }

    return value;

} // Remap::remap

} // namespace Portage

/*-------------------------------------------------------------------------~--*
 * Formatting options for Emacs and vim.
 *
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:t
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *-------------------------------------------------------------------------~--*/
