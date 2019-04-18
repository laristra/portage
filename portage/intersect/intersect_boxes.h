#ifndef INTERSECT_BOXES_H
#define INTERSECT_BOXES_H

// ============================================================================

#include <cmath>
#include <vector>

#include "portage/support/portage.h"
#include "portage/support/Point.h"

// ============================================================================

namespace Portage {

///
/// \class IntersectBoxes  Intersection algorithm for axis-aligned boxes

template <int D, typename SourceMeshType, typename TargetMeshType>
class IntersectBoxes {

 public:

  //! Default constructor (disabled)
  IntersectBoxes() = delete;

  //! Constructor with meshes
  IntersectBoxes(const SourceMeshType& source_mesh,
                 const TargetMeshType& target_mesh)
      : sourceMeshWrapper(source_mesh),
        targetMeshWrapper(target_mesh),
        geometry(source_mesh.geometry()) {}

  /// Assignment operator (disabled)
  IntersectBoxes & operator = (const IntersectBoxes&) = delete;

  /// \brief Intersect control volume of a target entity with control volumes of a set of source entities
  /// \param[in] tgt_entity  Entity of target mesh to intersect
  /// \param[in] src_entities Entities of source mesh to intersect against
  /// \return vector of Weights_t structure containing moments of intersection
  std::vector<Portage::Weights_t>
  operator() (const int tgt_cell,
              const std::vector<int>& src_cells) const;

 private:
  const SourceMeshType& sourceMeshWrapper;
  const TargetMeshType& targetMeshWrapper;
  const Geometry geometry;

};  // class IntersectBoxes


template <int D, typename SourceMeshType, typename TargetMeshType>
std::vector<Portage::Weights_t>
IntersectBoxes<D, SourceMeshType, TargetMeshType>::operator() (
    const int tgt_cell,
    const std::vector<int>& src_cells) const {

  Portage::Point<D> tlo, thi;
  targetMeshWrapper.cell_get_bounds(tgt_cell, &tlo, &thi);

  int nsrc = src_cells.size();
  std::vector<Portage::Weights_t> sources_and_weights(nsrc);
  int ninserted = 0;
  for (int i = 0; i < nsrc; ++i) {
    int s = src_cells[i];

    Portage::Point<D> slo, shi;
    sourceMeshWrapper.cell_get_bounds(s, &slo, &shi);

    // compute intersection and volume
    Portage::Point<D> ilo, ihi;
    double vol = 1.;
    for (int d = 0; d < D; ++d) {
      ilo[d] = std::max(slo[d], tlo[d]);
      ihi[d] = std::min(shi[d], thi[d]);
      double delta = std::max(ihi[d] - ilo[d], 0.);
      vol *= delta;
    }

    if (vol <= 0.) continue;

    auto & this_wt = sources_and_weights[ninserted];
    this_wt.entityID = s;

    // compute weights (volume + 1 moment for each dimension)
    auto & weights = this_wt.weights;
    weights.resize(1 + D);
    switch (geometry) {
      case Geometry::CARTESIAN:
        weights[0] = vol;
        for (int d = 0; d < D; ++d) {
          weights[1+d] = vol * 0.5 * (ilo[d] + ihi[d]);
        }
        break;
      case Geometry::CYLINDRICAL:
        double xsum = ilo[0] + ihi[0];
        double volrz = vol * 0.5 * xsum;
        weights[0] = volrz;
        weights[1] = vol * (1./3.) * (xsum * xsum - ilo[0] * ihi[0]);
        for (int d = 1; d < D; ++d) {
          weights[1+d] = volrz * 0.5 * (ilo[d] + ihi[d]);
        }
      break;
    }

    ++ninserted;

  }  // for i

  sources_and_weights.resize(ninserted);
  return sources_and_weights;

}  // operator()


} // namespace Portage

#endif // INTERSECT_BOXES_H
