#ifndef INTERPOLATORS_EAP_SEARCH_CARTESIAN_HH_INCLUDE
#define INTERPOLATORS_EAP_SEARCH_CARTESIAN_HH_INCLUDE

#include <cassert>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "portage/support/portage.h"
#include "portage/support/Point.h"


namespace EAP {

///
/// \class SearchCartesian  Search a Cartesian product mesh

template <int D, typename SourceMeshType, typename TargetMeshType>
class SearchCartesian {

 public:

  //! Default constructor (disabled)
  SearchCartesian() = delete;

  //! Constructor with meshes
  SearchCartesian(const SourceMeshType& source_mesh,
                  const TargetMeshType& target_mesh)
      : sourceMesh_(source_mesh),
        targetMesh_(target_mesh) {}

  /// Assignment operator (disabled)
  SearchCartesian & operator = (const SearchCartesian&) = delete;

  /// Search for source cells that intersect a given target cell
  std::vector<int64_t> operator() (const int tgt_cell) const;

 private:

  template <int N>
  using IntPoint = std::array<int, N>;

  // list overlap cells given bounds in each dimension
  std::vector<int64_t> list_cells(
        const IntPoint<D>& ilo, const IntPoint<D>& ihi) const;

 private:

  const SourceMeshType& sourceMesh_;
  const TargetMeshType& targetMesh_;

};  // class SearchCartesian


template <int D, typename SourceMeshType, typename TargetMeshType>
std::vector<int64_t> 
SearchCartesian<D, SourceMeshType, TargetMeshType>::operator() (const int tgt_cell) const {

  const auto EPSILON = 10. * std::numeric_limits<double>::epsilon();

  Portage::Point<D> tlo, thi;
  targetMesh_.cell_get_bounds(tgt_cell, &tlo, &thi);

  // quick check:  does this cell intersect the source mesh at all?
  // if not, exit early
  Portage::Point<D> sglo, sghi;
  sourceMesh_.get_global_bounds(&sglo, &sghi);
  for (int d = 0; d < D; ++d) {
    // allow for roundoff error
    tlo[d] += EPSILON;
    thi[d] -= EPSILON;
    assert(tlo[d] < thi[d]);
    assert(sglo[d] < sghi[d]);
    if (tlo[d] >= sghi[d] || thi[d] <= sglo[d])
      return std::move(std::vector<int64_t>());
  }

  IntPoint<D> ilo, ihi;

  // find which source cells overlap target cell, in each dimension
  for (int d = 0; d < D; ++d) {
    auto saxis_begin = sourceMesh_.axis_point_begin(d);
    auto saxis_end   = sourceMesh_.axis_point_end(d);

    // find last axis point less than or equal to tlo
    auto is_point_above = [&] (const double x, const int pid) {
      return (x < sourceMesh_.axis_point_coordinate(d, pid));
    };
    auto itrlo = std::upper_bound(saxis_begin, saxis_end, tlo[d],
                                  is_point_above);
    ilo[d] = itrlo - saxis_begin - 1;
    ilo[d] = std::max(ilo[d], 0);

    // find first axis point greater than or equal to thi
    auto is_point_below = [&] (const int pid, const double x) {
      return (sourceMesh_.axis_point_coordinate(d, pid) < x);
    };
    auto itrhi = std::lower_bound(itrlo, saxis_end, thi[d],
                                  is_point_below);
    ihi[d] = itrhi - saxis_begin;
    ihi[d] = std::min(ihi[d], sourceMesh_.axis_num_cells(d));
    assert(ihi[d] > ilo[d]);
  }  // for d

  return std::move(list_cells(ilo, ihi));

}  // operator()


template <int D, typename SourceMeshType, typename TargetMeshType>
std::vector<int64_t>
SearchCartesian<D, SourceMeshType, TargetMeshType>::list_cells(
      const IntPoint<D>& ilo, const IntPoint<D>& ihi) const {

  std::vector<int64_t> list;

  switch (D) {
    case 1:
    // 1D case:  iterate over i only
    {
      for (int i = ilo[0]; i < ihi[0]; ++i) {
        int64_t nglobal = (int64_t) i;
        list.emplace_back(nglobal);
      }
    }  // 1D case
    break;

    case 2:
    // 2D case:  iterate over i, j
    {
      // TODO: don't assume ordering, call a mesh wrapper function
      //       to convert ijk to global
      int imax = sourceMesh_.axis_num_cells(0);
      for (int j = ilo[1]; j < ihi[1]; ++j) {
        for (int i = ilo[0]; i < ihi[0]; ++i) {
          int64_t nglobal = ((int64_t) i) +
                            ((int64_t) j) * imax;
          list.emplace_back(nglobal);
        }
      }
    }  // 2D case
    break;

    case 3:
    // 3D case:  iterate over i, j, k
    {
      int imax = sourceMesh_.axis_num_cells(0);
      int jmax = sourceMesh_.axis_num_cells(1);
      for (int k = ilo[2]; k < ihi[2]; ++k) {
        for (int j = ilo[1]; j < ihi[1]; ++j) {
          for (int i = ilo[0]; i < ihi[0]; ++i) {
            int64_t nglobal = ((int64_t) i) +
                              ((int64_t) j) * imax +
                              ((int64_t) k) * imax * jmax;
            list.emplace_back(nglobal);
          }
        }
      }
    }  // 3D case
    break;

  }  // switch D

  return std::move(list);

}  // list_cells

} // namespace EAP

#endif // INTERPOLATORS_EAP_SEARCH_CARTESIAN_HH_INCLUDE
