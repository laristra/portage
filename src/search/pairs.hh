/*---------------------------------------------------------------------------~
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
 *---------------------------------------------------------------------------~*/

#ifndef pairs_INCLUDED
#define pairs_INCLUDED

#include <vector>
#include <list>

#include "pile.hh"
#include "lretypes.hh"

namespace Portage {
namespace Meshfree {
namespace Pairs {

  //\///////////////////////////////////////////////////////////////////////////
  // pair finding functions
  //\///////////////////////////////////////////////////////////////////////////

  /// search structure
  class CellPairFinder {
   public:

    /// Neighbor finding based on containment: build structure
    CellPairFinder(
      const vpile &x, const vpile &y, const vpile &h,
		  const bool do_scatter);

    /// Neighbor finding based on containment: find for given point
    std::list<ulong> find(const ulong j) const;

   private:
    std::list<ulong> find_gather(const ulong j) const;
    std::list<ulong> find_scatter(const ulong j) const;

    vpile x, y, h;
    size_t dim;
    bool do_scatter;
    pile hmax;
    vpile yminmax;
    pile delta;
    vulong nsidesm;
    vulong strides;
    std::vector<std::vector<ulong>> cells;
  };

}  // namespace Pairs
}  // namespace Meshfree
}  // namespace Portage

#endif
