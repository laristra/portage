#ifndef INTERSECT_BOXES_H
#define INTERSECT_BOXES_H

// ============================================================================

#include <cmath>
#include <type_traits>
#include <vector>

#include "wonton/support/CoordinateSystem.h"
#include "wonton/support/Point.h"

#include "portage/support/portage.h"

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
        targetMeshWrapper(target_mesh) {
    // TODO: The mesh wrappers need to expose the coordinate systems in some
    //       way.  Then we can assert that the coordinate systems are the same.
    //       The problem is that some mesh wrappers won't have a coordinate
    //       system, so there would then have to be some way to check that and
    //       then declare that Cartesian is the default.
    //static_assert(std::is_same<>::value);
  }

  //! Assignment operator (disabled)
  IntersectBoxes& operator= (const IntersectBoxes&) = delete;

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
  Wonton::Point<D> tlo, thi;
  targetMeshWrapper.cell_get_bounds(tgt_cell, &tlo, &thi);

  // Allocate storage for moments
  int nsrc = src_cells.size();
  std::vector<Portage::Weights_t> sources_and_weights(nsrc);

  // Loop over source boxes and compute moments
  int ninserted = 0;
  for (int i = 0; i < nsrc; ++i) {
    int s = src_cells[i];

    // Get source cell bounding box
    Wonton::Point<D> slo, shi;
    sourceMeshWrapper.cell_get_bounds(s, &slo, &shi);

    // Compute intersection and volume
    Wonton::Point<D> ilo, ihi;   // bounding box of intersection
    double volume = 1.;
    for (int d = 0; d < D; ++d) {
      ilo[d] = std::max(slo[d], tlo[d]);
      ihi[d] = std::min(shi[d], thi[d]);
      double delta = std::max(ihi[d] - ilo[d], 0.);
      volume *= delta;
    }

    // If the intersection volume is zero, don't bother computing the moments
    // because they will also be zero, and don't bother inserting into the list
    // of actually intersecting boxes.
    if (volume <= 0.) continue;

    // At this point the two boxes actually intersect and we need to compute
    // moments and insert into the list.

    // Save the source cell ID
    auto & this_wt = sources_and_weights[ninserted];
    this_wt.entityID = s;

    // Compute moments
    // TODO: How many orders of moments do we need to provide?  For first-order
    //       interpolation, we only need zeroth moments; for second-order
    //       interpolation, we need up through first moments; and so on.
    //       Unfortunately, the intersector doesn't currently know what order
    //       of interpolation is needed.  So we'll just have to hard-code this
    //       for now, assuming second-order interpolation for lack of a better
    //       choice.
    Wonton::Point<D> first_moments;
    for (int d = 0; d < D; ++d) {
      const auto ibar = 0.5 * (ilo[d] + ihi[d]);
      first_moments[d] = ibar * volume;
    }

    // Correct for coordinate system
    // -- Instead of calculating extra moments, use the bounding box to
    //    explicitly update the moments.  But this only works if your cells are
    //    axis-aligned boxes.  In other words: this is an optimization for
    //    intersect_boxes, rather than a general-purpose tool for all
    //    intersectors.  Currently only zeroth and first moments have this
    //    optimization, but they could be computed and implemented for
    //    higher-order moments.
    CoordSys::modify_volume(volume, ilo, ihi);
    CoordSys::modify_first_moments(first_moments, ilo, ihi);

    // Save weights (moments)
    auto & weights = this_wt.weights;
    weights.resize(1+D);
    weights[0] = volume;
    for (int d = 0; d < D; ++d) {
        weights[1+d] = first_moments[d];
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
