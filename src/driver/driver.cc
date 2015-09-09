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

#include "Mesh.hh"
#include "MeshFactory.hh"


namespace Portage {

void Driver::run()
{
    std::printf("in Driver::run()...\n");

	const SearchSimple search(&sourceMesh_, &targetMesh_);

	//Get an instance of the desired intersect algorithm type
	IntersectClipper<Jali::Entity_ID> intersect{cellToXY(&sourceMesh_), cellToXY(&targetMesh_)};

	// Eventually put this in a loop over remap variable names as well
	// Assume for now that we are only doing cell-based remap
	const Remap_1stOrder remap(sourceMesh_, sourceState_, 
							   remap_var_names_[0], Jali::CELL);

    int numTargetCells = targetMesh_.num_entities(Jali::CELL,Jali::OWNED);
    std::cout << "Number of target cells in target mesh "
              << numTargetCells << std::endl;

	// Ask for a StateVector with the name remap_var_names_[0] to be added to the targetState_. If it is already present, the existing StateVector reference is returned. If its not present, it is added. This logic needs to be reversed. The find function should add it if it is not found (if so requested).

    std::vector<double> dummyvals(numTargetCells,0);
    Portage::StateVector & targetField = 
            targetState_.add(remap_var_names_[0],Jali::CELL,&(dummyvals[0]));

    // Create a cellIndices vector and populates with a sequence of
    // ints starting at 0. Will go away when Jali has iterators for
    // mesh entities

    std::vector<int> cellIndices(numTargetCells);
    std::iota(cellIndices.begin(), cellIndices.end(), 0);

	composerFunctor<SearchSimple, IntersectClipper<Jali::Entity_ID >, Remap_1stOrder> 
		composer(&search, &intersect, &remap,
				 &sourceMesh_, &targetMesh_,
				 remap_var_names_[0]);

    // this populates targetField with the doubles returned from the final remap
    std::transform(cellIndices.begin(), cellIndices.end(),
                   targetField.begin(),
                   composer);

#ifdef DEBUG_OUTPUT
    Portage::State::const_iterator 
            itc = targetState_.find(remap_var_names_[0], Jali::CELL);
    Portage::StateVector & stateVector =  *itc;
    std::cout << stateVector << std::endl; 
#endif

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
