#ifndef INTERSECT_BOXES_H
#define INTERSECT_BOXES_H

// ============================================================================

#include <cmath>
#include <vector>

#include "wonton/support/CoordinateSystems.h"

#include "portage/support/portage.h"
#include "portage/support/Point.h"

// ============================================================================

namespace Portage {

//! \class IntersectBoxes  Intersection algorithm for axis-aligned boxes
//!
//! It is assumed that both boxes are axis-aligned and in the same geometry.
//! For example, this will work for two axis-aligned cells in 2D axis-aligned
//! cylindrical geometry, but it cannot intersect a 3D axis-aligned box in
//! Cartesian coordinates with a 3D axis-aligned box in spherical coordinates.

template <int D, typename SourceMeshType, typename TargetMeshType,
         class CoordSys = Wonton::DefaultCoordSys>
class IntersectBoxes {

 public:

  //! Default constructor (disabled)
  IntersectBoxes() = delete;

  //! Constructor with meshes
  IntersectBoxes(const SourceMeshType& source_mesh,
                 const TargetMeshType& target_mesh)
      : sourceMeshWrapper(source_mesh),
        targetMeshWrapper(target_mesh) {}

  //! Assignment operator (disabled)
  IntersectBoxes & operator = (const IntersectBoxes&) = delete;

  //! \brief Intersect control volume of a target box with control volumes
  //!        of a set of source boxes.
  //! \param[in] tgt_cell  Entity of target mesh to intersect
  //! \param[in] src_cells Entities of source mesh to intersect against
  //! \return vector of Weights_t structure containing moments of intersection
  std::vector<Portage::Weights_t> operator() (
      const int tgt_cell, const std::vector<int>& src_cells) const;

 private:

  //! Wrapper for the source mesh
  const SourceMeshType& sourceMeshWrapper;

  //! Wrapper for the target mesh
  const TargetMeshType& targetMeshWrapper;

};  // class IntersectBoxes


// ============================================================================
// Callable operator
template <int D, typename SourceMeshType, typename TargetMeshType,
         class CoordSys>
std::vector<Portage::Weights_t>
  IntersectBoxes<D, SourceMeshType, TargetMeshType, CoordSys>::operator() (
  const int tgt_cell, const std::vector<int>& src_cells) const {

  // Get bounding box for target box
  Portage::Point<D> tlo, thi;
  targetMeshWrapper.cell_get_bounds(tgt_cell, &tlo, &thi);

  // Allocate storage for moments
  int nsrc = src_cells.size();
  std::vector<Portage::Weights_t> sources_and_weights(nsrc);

  // Loop over source boxes and compute moments
  int ninserted = 0;
  for (int i = 0; i < nsrc; ++i) {
    int s = src_cells[i];

    // Get source cell bounding box
    Portage::Point<D> slo, shi;
    sourceMeshWrapper.cell_get_bounds(s, &slo, &shi);

    // Compute intersection and volume
    Portage::Point<D> ilo, ihi;   // bounding box of intersection
    double vol0 = 1.;
    for (int d = 0; d < D; ++d) {
      ilo[d] = std::max(slo[d], tlo[d]);
      ihi[d] = std::min(shi[d], thi[d]);
      double delta = std::max(ihi[d] - ilo[d], 0.);
      vol0 *= delta;
    }

    // If the intersection volume is zero, don't bother computing the moments
    // because they will also be zero, and don't bother inserting into the list
    // of actually intersecting boxes.
    if (vol0 <= 0.) continue;

    // At this point the two boxes actually intersect and we need to compute
    // moments and insert into the list.

    // Save the source cell ID
    auto & this_wt = sources_and_weights[ninserted];
    this_wt.entityID = s;

    // Compute and save weights (volume + 1 moment for each dimension)
    auto & weights = this_wt.weights;
    weights.resize(1 + D);
    switch (geometry) {
      weights[0] = CoordSys::modify_volume(vol0, ilo, ihi);
      Wonton::Point<D> mom0;
      for (int d = 0; d < D; ++d) {
        const ibar = 0.5 * (ilo[d] + ihi[d]);
        mom0[d] = ibar * vol0;
      }
      auto first_moments = CoordSys::compute_moments(mom0, ilo, ihi);
      for (int d = 0; d < D; ++d) {
        weights[1+d] = first_moments[d];
      }
    }

    // Increment the count, because we've now inserted a new entry
    ++ninserted;

  }  // for i

  // Not all source cells provided necessary actually intersect.
  sources_and_weights.resize(ninserted);

  return sources_and_weights;

}  // operator()


} // namespace Portage

#endif // INTERSECT_BOXES_H
