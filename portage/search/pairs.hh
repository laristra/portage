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

namespace Portage { namespace Meshfree { namespace Pairs {

//\///////////////////////////////////////////////////////////////////////////
// pair finding functions
//\///////////////////////////////////////////////////////////////////////////

/// search structure
class CellPairFinder {
public:

  /**
   * @brief Default constructor
   */
  CellPairFinder() = default;

  /**
   * @brief Build the neighbor finder object.
   *
   * @param in_x: contained points.
   * @param in_y: containing points with attached box.
   * @param in_h: box sizes, centered on y.
   * @param in_do_scatter: specify whether scatter or gather form.
   */
  CellPairFinder(vpile const& in_x, vpile const& in_y,
                 vpile const& in_h, bool in_do_scatter) {
    init(in_x, in_y, in_h, in_do_scatter);
  }

  /**
   * @brief Disabled copy constructor
   */
  CellPairFinder(const CellPairFinder &) = delete;

  /**
   * @brief Disabled assignment constructor
   */
  CellPairFinder &operator=(const CellPairFinder &) = delete;

  /**
   * @brief Destructor
   */
  ~CellPairFinder() = default;

  /**
   * @brief Build internal structures for neighbor finding.
   *
   * @param in_x: contained points.
   * @param in_y: containing points with attached box.
   * @param in_h: box sizes, centered on y.
   * @param in_do_scatter: specify whether scatter or gather form.
   */
  void init(vpile const& in_x, vpile const& in_y,
            vpile const& in_h, bool in_do_scatter);

  /**
   * @brief Find neighbors of a given point based on containment.
   *
   * @param j: current point index.
   * @return a lst of neighbors of the j-th point.
   */
  std::list<ulong> find(ulong j) const {
    return do_scatter ? find_scatter(j) : find_gather(j);
  }

protected:
  std::list<ulong> find_gather(ulong j) const;
  std::list<ulong> find_scatter(ulong j) const;
  vpile PairsMinMax(const vpile &in_y, const vpile &in_h) const;
  vpile PairsMinMax(const vpile &in_c, const pile &in_h) const;
  void PairsIntegize(int in_dim, const pile &in_value,
                     const vpile &in_minmax, const pile &in_delta,
                     const vulong &in_sizes, vulong &in_ivalue) const;
  ulong cellindex(int in_dim, const vulong &in_strides,
                  const vulong &in_indices) const;
  void cellindices(int in_dim, const vulong &in_strides,
                   const ulong &in_index, vulong &in_indices) const;


private:
  int dim = 1;
  bool do_scatter = false;
  vpile x, y, h;
  vpile cminmax;
  pile delta;
  vulong nsidesm;
  vulong strides;
  std::vector<std::vector<ulong>> cells;
};

}}} // namespace Portage::Meshfree::Pairs

#endif
