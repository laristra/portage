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

  /// data structure
  struct pairs_data_t {
    vpile x, y, h;
    pile hmax;
    vpile yminmax;
    pile delta;
    vulong nsidesm;
    vulong strides;
    vector<vector<ulong>> cells;
  };

  std::list<ulong> PairsContainCellsG(
      const std::shared_ptr<pairs_data_t> pairdata_p,
      const ulong j);

  std::list<ulong> PairsContainCellsS(
      const std::shared_ptr<pairs_data_t> pairdata_p,
      const ulong j);

  /// Neighbor finding based on containment: driver function
  std::shared_ptr<pairs_data_t> PairsFind(
      const vpile &x, const vpile &y, const vpile &h,
		  const bool do_scatter,
      const contain_type type=CELLS, const bool half_pairs=false);

}
}
}

#endif
