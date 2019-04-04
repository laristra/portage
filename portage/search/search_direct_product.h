#ifndef PORTAGE_SEARCH_SEARCH_DIRECT_PRODUCT_H_
#define PORTAGE_SEARCH_SEARCH_DIRECT_PRODUCT_H_

#include <cassert>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "wonton/wonton/support/CellID.h"
#include "wonton/wonton/support/Point.h"
#include "portage/support/portage.h"

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

namespace Portage {

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
  std::vector<Wonton::CellID> operator() (const int tgt_cell) const;

 private:

  // ==========================================================================
  // Private support methods

  //! List cells given bounds in each dimension
  std::vector<Wonton::CellID> list_cells(
        const std::array<int,D>& ilo, const std::array<int,D>& ihi) const;

  // ==========================================================================
  // Class data

  const SourceMeshType& sourceMesh_;
  const TargetMeshType& targetMesh_;

};  // class SearchDirectProduct


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// Constructor with meshes.
template <int D, typename SourceMeshType, typename TargetMeshType>
SearchDirectProduct<D,SourceMeshType,TargetMeshType>::SearchDirectProduct(
    const SourceMeshType& source_mesh, const TargetMeshType& target_mesh)
    : sourceMesh_(source_mesh), targetMesh_(target_mesh) {
}

// ============================================================================
// Callable

// Search for source cells that intersect a given target cell.
template <int D, typename SourceMeshType, typename TargetMeshType>
std::vector<Wonton::CellID> 
    SearchDirectProduct<D, SourceMeshType, TargetMeshType>::operator() (
    const int tgt_cell) const {

  // Tolerance for floating-point round-off
  const auto EPSILON = 10. * std::numeric_limits<double>::epsilon();

  // Get the bounding box for the target mesh cell
  Wonton::Point<D> tlo, thi;
  targetMesh_.cell_get_bounds(tgt_cell, &tlo, &thi);
  for (int d = 0; d < D; ++d) {
    // allow for roundoff error: if the intersection is within epsilon, assume
    // this is a precision issue rather than true physical overlap
    auto eps = EPSILON * (thi[d] - tlo[d]);
    tlo[d] += eps;
    thi[d] -= eps;
  }

  // quick check:  does this cell intersect the source mesh at all?
  // if not, exit early
  Wonton::Point<D> sglo, sghi;
  sourceMesh_.get_global_bounds(&sglo, &sghi);
  for (int d = 0; d < D; ++d) {
    assert(tlo[d] < thi[d]);
    assert(sglo[d] < sghi[d]);
    if (tlo[d] >= sghi[d] || thi[d] <= sglo[d])
      return std::move(std::vector<Wonton::CellID>());
  }

  // find which source cells overlap target cell, in each dimension
  std::array<int,D> ilo, ihi;
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
std::vector<Wonton::CellID>
    SearchDirectProduct<D, SourceMeshType, TargetMeshType>::list_cells(
    const std::array<int,D>& ilo, const std::array<int,D>& ihi) const {

  // Declare the list of cells to be returned
  std::vector<Wonton::CellID> list;

  // I think this is less clear than the version below
  /*std::array<int,3> idx_lo = {0,0,0};
    std::array<int,3> idx_hi = {1,1,1};
  for (int d = 0; d < D; ++d) {
    idx_lo[d] = ilo[d];
    idx_hi[d] = ihi[d];
  }
  std::array<int,3> i3;
  for (int i3[2] = idx_lo[2]; i3[2] < idx_hi[2]; ++i3[2]) {
    for (int i3[1] = idx_lo[1]; i3[1] < idx_hi[1]; ++i3[1]) {
      for (int i3[0] = idx_lo[0]; i3[0] < idx_hi[0]; ++i3[0]) {
      std::array<int,D> indices;
        for (int d = 0; d < D; ++d) {
          indices[d] = i3[d];
        }
        list.emplace_back(sourceMesh_.indices_to_cellid<D>(indices));
      }
    }
  }*/

  // TODO: This could be done for any dimensionality using a recursive
  //       function.
  std::array<int,D> idx;
  switch (D) {
    case 1:
      // 1D case:  iterate over i only
      for (idx[0] = ilo[0]; idx[0] < ihi[0]; ++idx[0]) {
        list.emplace_back(sourceMesh_.indices_to_cellid(idx));
      }
      break;
    case 2:
    // 2D case:  iterate over i, j
      for (idx[1] = ilo[1]; idx[1] < ihi[1]; ++idx[1]) {
        for (idx[0] = ilo[0]; idx[0] < ihi[0]; ++idx[0]) {
          list.emplace_back(sourceMesh_.indices_to_cellid(idx));
        }
      }
      break;
    case 3:
      // 3D case:  iterate over i, j, k
      for (idx[2] = ilo[2]; idx[2] < ihi[2]; ++idx[2]) {
        for (idx[1] = ilo[1]; idx[1] < ihi[1]; ++idx[1]) {
          for (idx[0] = ilo[0]; idx[0] < ihi[0]; ++idx[0]) {
            list.emplace_back(sourceMesh_.indices_to_cellid(idx));
          }
        }
      }
      break;
  }  // switch D

  return std::move(list);

} // list_cells

} // namespace Portage

#endif // PORTAGE_SEARCH_SEARCH_DIRECT_PRODUCT_H_
