#ifndef INTERPOLATORS_EAP_SEARCH_DIRECT_PRODUCT_HH_INCLUDE
#define INTERPOLATORS_EAP_SEARCH_DIRECT_PRODUCT_HH_INCLUDE

#include <cassert>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

/*!
  @file search_direct_product.h
  @brief Definition of the SearchDirectProduct class.

  The SearchDirectProduct mesh does the search process for meshes subject to
  the following assumptions:
  - The source mesh is a static, axis-aligned, logically-rectangular mesh, with
    cell widths that are allowed to vary across the mesh.
  - The target mesh may be unstructured, but its wrapper must be able to
    provide an axis-aligned bounding box for any cell.
*/

namespace EAP {

/*!
  @class SearchDirectProduct "search_direct_product.h"
  @brief Definition of the SearchDirectProduct class.

  The SearchDirectProduct mesh does the search process for meshes subject to
  the following assumptions:
  - The source mesh is a static, axis-aligned, logically-rectangular mesh, with
    cell widths that are allowed to vary across the mesh.
  - The target mesh may be unstructured, but its wrapper must be able to
    provide an axis-aligned bounding box for any cell.
*/

template <int D, typename SourceMeshType, typename TargetMeshType>
class SearchDirectProduct {

 public:

  // ==========================================================================
  // Constructors and destructors

  //! Default constructor (disabled)
  SearchDirectProduct() = delete;

  /*!
    @brief Constructor with meshes.
    @param[in] source_mesh The source mesh wrapper.
    @param[in] target_mesh The target mesh wrapper.
  */
  SearchDirectProduct(const SourceMeshType& source_mesh,
                      const TargetMeshType& target_mesh);

  //! Assignment operator (disabled)
  SearchDirectProduct & operator = (const SearchDirectProduct&) = delete;

  // ==========================================================================
  // Callable operation

  /*!
    @brief Search for source cells that intersect a given target cell.
    @param[in] tgt_cell The cell on the target mesh

    Callable routine that takes a target cell and searches for overlapping
    cells on the source mesh.  Returns the list of overlapping source mesh
    cells.
  */
  std::vector<int64_t> operator() (const int tgt_cell) const;

 private:

  // ==========================================================================
  // Define convenient shorthands

  // An N-dimensional point (of integers rather than doubles)
  template <int N>
  using IntPoint = std::array<int, N>;

  // ==========================================================================
  // Private support methods

  //! List cells given bounds in each dimension
  std::vector<int64_t> list_cells(
        const IntPoint<D>& ilo, const IntPoint<D>& ihi) const;

  // ==========================================================================
  // Class data

  const SourceMeshType& sourceMesh_;
  const TargetMeshType& targetMesh_;

};  // class SearchDirectProduct


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// Constructor with meshes.
SearchDirectProduct(const SourceMeshType& source_mesh,
                    const TargetMeshType& target_mesh)
    : sourceMesh_(source_mesh),
      targetMesh_(target_mesh) {
}

// ============================================================================
// Callable

// Search for source cells that intersect a given target cell.
template <int D, typename SourceMeshType, typename TargetMeshType>
std::vector<int64_t> 
    SearchDirectProduct<D, SourceMeshType, TargetMeshType>::operator() (
    const int tgt_cell) const {

  // Tolerance for floating-point round-off
  const auto EPSILON = 10. * std::numeric_limits<double>::epsilon();

  // Get the bounding box for the target mesh cell
  Portage::Point<D> tlo, thi;
  targetMesh_.cell_get_bounds(tgt_cell, &tlo, &thi);
  for (int d = 0; d < D; ++d) {
    // allow for roundoff error: if the intersection is within epsilon, assume
    // this is a precision issue rather than true physical overlap
    tlo[d] += EPSILON;
    thi[d] -= EPSILON;
  }

  // quick check:  does this cell intersect the source mesh at all?
  // if not, exit early
  Portage::Point<D> sglo, sghi;
  sourceMesh_.get_global_bounds(&sglo, &sghi);
  for (int d = 0; d < D; ++d) {
    assert(tlo[d] < thi[d]);
    assert(sglo[d] < sghi[d]);
    if (tlo[d] >= sghi[d] || thi[d] <= sglo[d])
      return std::move(std::vector<int64_t>());
  }

  // find which source cells overlap target cell, in each dimension
  IntPoint<D> ilo, ihi;
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
    // std::max(ilo[d],0) instead of assert(ilo[d] >= 0) to account for cells
    // that only partially overlap the target cell
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

  // Generate list of cells from lower and upper bounds, return list
  return std::move(list_cells(ilo, ihi));

}  // operator()


// ============================================================================
// Private support methods

// List cells given bounds in each dimension
template <int D, typename SourceMeshType, typename TargetMeshType>
std::vector<int64_t>
    SearchDirectProduct<D, SourceMeshType, TargetMeshType>::list_cells(
    const IntPoint<D>& ilo, const IntPoint<D>& ihi) const {

  // Declare the list of cells to be returned
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

} // list_cells

} // namespace EAP

#endif // INTERPOLATORS_EAP_SEARCH_DIRECT_PRODUCT_HH_INCLUDE
