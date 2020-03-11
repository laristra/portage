#ifndef PORTAGE_SEARCH_SEARCH_DIRECT_PRODUCT_H_
#define PORTAGE_SEARCH_SEARCH_DIRECT_PRODUCT_H_

#include <cassert>
#include <cstdint>
#include <limits>
#include <algorithm>
#include <array>
#include <utility>
#include <vector>

#include "wonton/support/Point.h"
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
  @class SearchDirectProductBase "search_direct_product.h"
  @brief Definition of the SearchDirectProductBase class.

  The SearchDirectProductBase mesh does the search process for meshes subject
  to the following assumptions:
  - The source mesh is a static, axis-aligned, logically-rectangular mesh, with
    cell widths that are allowed to vary across the mesh.
  - The target mesh may be unstructured, but its wrapper must be able to
    provide an axis-aligned bounding box for any cell.

  Unlike most searches, SearchDirectProductBase is additionally tempated on the
  ID type (e.g., global vs local IDs).  Many applications will use a typedef
  that specialized to a specific ID type.
*/

template <
    typename ID_t, int D, typename SourceMeshType, typename TargetMeshType>
class SearchDirectProductBase {

 public:

  // ==========================================================================
  // Constructors and destructors

  //! Default constructor (disabled)
  SearchDirectProductBase() = delete;

  /*!
    @brief Constructor with meshes.
    @param[in] source_mesh The source mesh wrapper.
    @param[in] target_mesh The target mesh wrapper.
  */
  SearchDirectProductBase(const SourceMeshType& source_mesh,
                          const TargetMeshType& target_mesh);

  //! Assignment operator (disabled)
  SearchDirectProductBase & operator= (const SearchDirectProductBase&) =
    delete;

  // ==========================================================================
  // Callable operation

  /*!
    @brief Search for source cells that intersect a given target cell.
    @param[in] tgt_cell The cell on the target mesh

    Callable routine that takes a target cell and searches for overlapping
    cells on the source mesh.  Returns the list of overlapping source mesh
    cells.
  */
  std::vector<int> operator() (const int tgt_cell) const;

 private:

  // ==========================================================================
  // Private support methods

  //! Loop over each dimension recursively and fill the list
  void fill_list_by_dim(
    std::vector<int> &list,
    const int d, const int idx_start, const std::array<int,D> &step_size,
    const std::array<int,D> &ilo, const std::array<int,D> &ihi,
    std::array<int,D> &indices) const;

  //! List cells given index bounds in each dimension
  std::vector<int> list_cells(
        const std::array<int,D>& ilo, const std::array<int,D>& ihi) const;

  // ==========================================================================
  // Class data

  const SourceMeshType& sourceMesh_;
  const TargetMeshType& targetMesh_;

};  // class SearchDirectProductBase


// ============================================================================
// Typedef

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
template<int D, typename SourceMeshType, typename TargetMeshType>
using SearchDirectProduct =
  SearchDirectProductBase<int, D, SourceMeshType, TargetMeshType>;


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// Constructor with meshes.
template <
    typename ID_t, int D, typename SourceMeshType, typename TargetMeshType>
SearchDirectProductBase<ID_t,D,SourceMeshType,TargetMeshType>::
    SearchDirectProductBase(
    const SourceMeshType& source_mesh, const TargetMeshType& target_mesh)
    : sourceMesh_(source_mesh), targetMesh_(target_mesh) {
}

// ============================================================================
// Callable

// Search for source cells that intersect a given target cell.
template <
    typename ID_t, int D, typename SourceMeshType, typename TargetMeshType>
std::vector<int> 
    SearchDirectProductBase<ID_t,D, SourceMeshType, TargetMeshType>::
    operator() ( const int tgt_cell) const {

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
      return std::move(std::vector<int>());
  }

  // find which source cells overlap target cell, in each dimension
  std::array<int,D> ilo, ihi;
  for (int d = 0; d < D; ++d) {
    auto saxis_begin = sourceMesh_.axis_point_begin(d);
    auto saxis_end   = sourceMesh_.axis_point_end(d);

    // find last axis point less than or equal to tlo
    auto is_point_above = [this,d] (const double x, const int pid) {
      return (x < this->sourceMesh_.get_axis_point(d, pid));
    };
    auto itrlo = std::upper_bound(saxis_begin, saxis_end, tlo[d],
                                  is_point_above);
    ilo[d] = itrlo - saxis_begin - 1;
    // std::max(ilo[d],0) instead of assert(ilo[d] >= 0) to account for cells
    // that only partially overlap the target cell
    ilo[d] = std::max(ilo[d], 0);

    // find first axis point greater than or equal to thi
    auto is_point_below = [this,d] (const int pid, const double x) {
      return (this->sourceMesh_.get_axis_point(d, pid) < x);
    };
    auto itrhi = std::lower_bound(itrlo, saxis_end, thi[d],
                                  is_point_below);
    ihi[d] = itrhi - saxis_begin;
    ihi[d] = std::min(ihi[d], sourceMesh_.axis_num_cells(d));
    assert(ihi[d] > ilo[d]);
  }  // for d

  // Generate list of cells from lower and upper bounds, return list
  return(std::move(list_cells(ilo, ihi)));

}  // operator()


// ============================================================================
// Private support methods

// Loop over each dimension recursively and fill the list
template<
    typename ID_t, int D, typename SourceMeshType, typename TargetMeshType>
void SearchDirectProductBase<ID_t,D,SourceMeshType,TargetMeshType>::
    fill_list_by_dim(
    std::vector<int> &list,
    const int d, const int idx_start, const std::array<int,D> &step_size,
    const std::array<int,D> &ilo, const std::array<int,D> &ihi,
    std::array<int,D> &indices) const {
  if (d >= 0) {
    int index = idx_start;
    for (indices[d] = ilo[d]; indices[d] < ihi[d]; ++indices[d]) {
      fill_list_by_dim(list, d-1, index, step_size, ilo, ihi, indices);
      index += step_size[d];
    }
  } else {
    list[idx_start] = sourceMesh_.indices_to_cellid(indices);
  }
}


// List cells given index bounds in each dimension
template <
    typename ID_t, int D, typename SourceMeshType, typename TargetMeshType>
std::vector<int>
    SearchDirectProductBase<ID_t,D, SourceMeshType, TargetMeshType>::
    list_cells(
    const std::array<int,D>& ilo, const std::array<int,D>& ihi) const {

  // Compute step sizes for recursion
  std::array<int,D> stepsize;
  stepsize[0] = 1;
  for(int d = 0; d < D-1; ++d) {
    stepsize[d+1] = stepsize[d] * (ihi[d] - ilo[d]);
  }

  // Declare the list of cells to be returned
  int list_size = stepsize[D-1] * (ihi[D-1] - ilo[D-1]);
  std::vector<int> list(list_size);

  // Recurse across dimensions
  std::array<int,D> indices;
  fill_list_by_dim(list, D-1, 0, stepsize, ilo, ihi, indices);

  return(std::move(list));

} // list_cells

} // namespace Portage

#endif // PORTAGE_SEARCH_SEARCH_DIRECT_PRODUCT_H_
