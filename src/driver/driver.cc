/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/
#include "driver.h"

#include <cstdio>
#include <vector>
#include <numeric>
#include <algorithm>

#include "portage/search/search_simple.h"
#include "portage/intersect/intersectClipper.h"
#include "portage/remap/remap_1st_order.h"
#include "portage/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "portage/wrappers/state/jali/jali_state_wrapper.h"

namespace Portage {


// When we get Portage::Entity_kind, Portage::Entity_ID definitions
// into a header file in the appropriate place, we can make everything
// in the method free of Jali references except the Jali_Mesh_Wrapper
// and Jali_State_Wrapper references

void Driver::run()
{
    std::printf("in Driver::run()...\n");

	const SearchSimple<Jali_Mesh_Wrapper,Jali_Mesh_Wrapper> 
            search(source_mesh_, target_mesh_);

	//Get an instance of the desired intersect algorithm type
	const IntersectClipper<Jali_Mesh_Wrapper, Jali_Mesh_Wrapper> 
            intersect{source_mesh_, target_mesh_};

	// Eventually put this in a loop over remap variable names as well
	// Assume for now that we are only doing cell-based remap
	const Remap_1stOrder<Jali_Mesh_Wrapper,Jali_State_Wrapper,Jali::Entity_kind>
            remap(source_mesh_, source_state_, Jali::CELL, remap_var_names_[0]);

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh "
              << numTargetCells << std::endl;

	// Ask for a StateVector with the name remap_var_names_[0] to be
	// added to the targetState_. If it is already present, the
	// existing StateVector reference is returned. If its not present,
	// it is added. This logic needs to be reversed. The find function
	// should add it if it is not found (if so requested).

    std::vector<double> dummyvals(numTargetCells,0);
    double *target_field = NULL;
    target_state_.get_data(Jali::CELL,remap_var_names_[0],&target_field);

    // Create a cellIndices vector and populates with a sequence of
    // ints starting at 0. 

	composerFunctor<SearchSimple<Jali_Mesh_Wrapper,Jali_Mesh_Wrapper>, 
        IntersectClipper<Jali_Mesh_Wrapper, Jali_Mesh_Wrapper>, 
        Remap_1stOrder<Jali_Mesh_Wrapper,Jali_State_Wrapper,Jali::Entity_kind> >
            composer(&search, &intersect, &remap, remap_var_names_[0]);

    // this populates targetField with the doubles returned from the final remap
    std::transform(target_mesh_.begin(Jali::CELL), target_mesh_.end(Jali::CELL),target_field,
                   composer);

} // Driver::run

} // namespace Portage

/*-------------------------------------------------------------------------~--*
 * Formatting options for Emacs and vim.
 *
 * Local Variables:
 * mode:c++
 * indent-tabs-mode:nil
 * c-basic-offset:4
 * tab-width:4
 * End:
 * vim: set tabstop=4 shiftwidth=4 expandtab :
 *-------------------------------------------------------------------------~--*/
