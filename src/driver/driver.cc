/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "driver.h"

#include <cstdio>
#include <vector>
#include <numeric>
#include <algorithm>

#include "portage/search/search.h"
#include "portage/intersect/intersectClipper.h"
#include "portage/remap/remap.h"

#include "Mesh.hh"
#include "MeshFactory.hh"


namespace Portage {

void Driver::run()
{
    std::printf("in Driver::run()...\n");

	const Search s(&sourceMesh_, &targetMesh_);
	// Use Angela's new Clipper-based intersector
	const IntersectClipper i;
	const Remap r(sourceMesh_, sourceState_, targetMesh_, targetState_);

    int numTargetCells = targetMesh_.num_entities(Jali::CELL, Jali::ALL);

	std::cout << "Number of target cells in target mesh " << numTargetCells << std::endl;

	std::vector<int> cellIndices(numTargetCells);
	// populates cellIndices with a sequence of ints starting at 0...
	std::iota(cellIndices.begin(), cellIndices.end(), 0);

	std::vector<double> newField(numTargetCells);
	// this will populate newField with the doubles returned from the final remap
	std::transform(cellIndices.begin(), cellIndices.end(),
				   newField.begin(),
				   composerFunctor<Search, IntersectClipper, Remap>(&s, &i, &r, 
																	&sourceMesh_, &targetMesh_,
																	remap_var_names_[0]));
	// Add it to the new state
	// CMM: is there a less hacky way of getting vector to array?
	targetState_.add("remapped_data", Jali::CELL, &newField[0]);

#ifdef DEBUG_OUTPUT
    std::vector<StateVector>::const_iterator field = targetState_.find("remapped_data", Jali::CELL);
    Portage::StateVector stateVector =  *field;
    for (unsigned int i=0; i<numTargetCells; i++)
    {
        double x = *(stateVector.begin() + i);
        std::cout << "Remapped field for cell " << i << ": " << x << std::endl; 
    }
#endif

} // Driver::run

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
