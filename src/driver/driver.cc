/*--------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *--------------------------------------------------------------------------~~*/

#include "driver.h"

#include <cstdio>
#include <vector>
//#include <utility>
#include <numeric>
#include <algorithm>

#include "portage/search/search.h"
#include "portage/intersect/intersect.h"
#include "portage/remap/remap.h"

#include "Mesh.hh"
#include "MeshFactory.hh"


namespace Portage {

void Driver::run()
{
    std::printf("in Driver::run()...\n");

	const Search s(&sourceMesh_, &targetMesh_);
	const Intersect i(sourceMesh_, targetMesh_);
	const Remap r(sourceMesh_, sourceState_, targetMesh_, targetState_);

    int numTargetCells = targetMesh_.num_entities(Jali::CELL, Jali::ALL);

	std::vector<int> cellIndices(numTargetCells);
	// populates cellIndices with a sequence of ints starting at 0...
	std::iota(cellIndices.begin(), cellIndices.end(), 0);

	//    double* newField = new double[numTargetCells];
	std::vector<double> newField(numTargetCells);

	// this will populate newField with the doubles returned from the final remap
	std::transform(cellIndices.begin(), cellIndices.end(),
				   newField.begin(),
				   composerFunctor<Search, Intersect, Remap>(&s, &i, &r));
//     for (unsigned int i=0; i<numTargetCells; i++)
//     {
//         // Search for possible intersections.
//         Search s;
//         Jali::Entity_ID_List candidates;
//         s.search(&sourceMesh_, &targetMesh_, i, &candidates);
 
// #ifdef DEBUG_OUTPUT
//         std::cout << std::endl << "Candidates for cell " << i << std::endl;
//         for (unsigned int j=0; j<candidates.size(); j++)
//             std::cout << candidates[j] << " ";
//         std::cout << std::endl;
// #endif

//         // Calculate the overlap of actual intersections.
//         std::vector<float> moments;
//         Intersect isect(sourceMesh_, targetMesh_);
//         isect.intersect(i, &candidates, &moments);

// #ifdef DEBUG_OUTPUT
//         std::cout << std::endl << "Moments for cell " << i << std::endl;
//         for (unsigned int j=0; j<moments.size(); j++)
//             std::cout << moments[j] << " ";
//         std::cout << std::endl;
// #endif

//         // Remap from sourceMesh_ to targetMesh_
//         Remap r(sourceMesh_, sourceState_, targetMesh_, targetState_);
//         newField[i] = r.remap(remap_var_names_[0], i, candidates, moments);
//     }

//     targetState_.add("remapped_data", Jali::CELL, newField);

// #ifdef DEBUG_OUTPUT
//     std::vector<StateVector>::const_iterator field = targetState_.find("remapped_data", Jali::CELL);
//     Portage::StateVector stateVector =  *field;
//     for (unsigned int i=0; i<numTargetCells; i++)
//     {
//         double x = *(stateVector.begin() + i);
//         std::cout << "Remapped field for cell " << i << ": " << x << std::endl; 
//     }
// #endif

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
