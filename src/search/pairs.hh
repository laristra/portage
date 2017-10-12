/*---------------------------------------------------------------------------~
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
 *---------------------------------------------------------------------------~*/

#ifndef pairs_INCLUDED
#define pairs_INCLUDED

#include <vector>
#include <list>
#include <memory>

#include "pile.hh"
#include "lretypes.hh"

namespace Portage {
namespace Meshfree {
namespace Pairs {

  //\///////////////////////////////////////////////////////////////////////////
  // pair finding functions
  //\///////////////////////////////////////////////////////////////////////////

  /// types of pair finder
  enum contain_type{CELLS, SORT, HASHX, HASHY};

  /// Neighbor finding based on containment: driver function
  std::shared_ptr<std::vector<std::list<ulong>>> PairsFind(
      const vpile &x, const vpile &y, const vpile &h,
		  const contain_type type=CELLS, const bool half_pairs=false);

}
}
}

#endif
