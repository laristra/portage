/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#ifndef PORTAGE_INTERSECT_INTERSECT_R2D_H_
#define PORTAGE_INTERSECT_INTERSECT_R2D_H_

#include <array>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <sstream>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Polytope.h"

// portage includes
extern "C" {
#include "wonton/intersect/r3d/r2d.h"
}
#include "portage/support/portage.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "portage/intersect/intersect_polys_r2d.h"

#ifdef PORTAGE_HAS_TANGRAM
#include "tangram/driver/CellMatPoly.h"
#include "tangram/driver/driver.h"
#include "tangram/support/MatPoly.h"
#endif 


namespace Portage {

namespace {
bool poly2_is_convex(std::vector<Wonton::Point<2>> const& pverts,
                     NumericTolerances_t const& num_tols) {

  // Check if target polygon is convex or at least star-convex
  int npverts = pverts.size();
  for (int i = 0; i < npverts; i++) {
    
    // Compute distance from pverts[i+1] to segment (pverts[i], pverts[i+2])
    int ifv = i, imv = (i+1)%npverts, isv = (i+2)%npverts;
    Wonton::Vector<2> normal( pverts[isv][1] - pverts[ifv][1], 
                             -pverts[isv][0] + pverts[ifv][0]);
    normal.normalize();
    Wonton::Vector<2> fv2mv = pverts[imv] - pverts[ifv];
    double dst = Wonton::dot(fv2mv, normal);
    
    if (dst <= -num_tols.min_absolute_distance)
      return false;
  }

  return true;
}
}


///
/// \class IntersectR2D  2-D intersection algorithm


template <Entity_kind on_what, class SourceMeshType,
          class SourceStateType, class TargetMeshType,
          template<class, int, class, class> class InterfaceReconstructorType =
          DummyInterfaceReconstructor,
          class Matpoly_Splitter = void,
          class Matpoly_Clipper = void>
class IntersectR2D {

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor2D =
      Tangram::Driver<InterfaceReconstructorType, 2, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef PORTAGE_HAS_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               std::shared_ptr<InterfaceReconstructor2D> ir)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir),
        num_tols_(num_tols) {}
#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), num_tols_(num_tols) {}

  /// \brief Set the source mesh material that we have to intersect against

  void set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect control volume of a target entity with control volumes of a set of source entities
  /// \param[in] tgt_entity  Entity of target mesh to intersect
  /// \param[in] src_entities Entities of source mesh to intersect against
  /// \return vector of Weights_t structure containing moments of intersection

  std::vector<Weights_t>
  operator() (const int tgt_entity, std::vector<int> const& src_entities) const {
    std::cerr << "IntersectR3D not implemented for this entity type" << std::endl;
    return std::vector<Weights_t>(0);
  }

  IntersectR2D() = delete;

  /// Assignment operator (disabled)
  IntersectR2D & operator = (const IntersectR2D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
  int matid_ = -1;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  NumericTolerances_t num_tols_;
};  // class IntersectR2D



//////////////////////////////////////////////////////////////////////////////
// Specialization of Intersect2D class for cells

template <class SourceMeshType, class SourceStateType,
          class TargetMeshType,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter, class Matpoly_Clipper>
class IntersectR2D<Entity_kind::CELL, SourceMeshType, SourceStateType, TargetMeshType,
                   InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper> {

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor2D =
      Tangram::Driver<InterfaceReconstructorType, 2, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef PORTAGE_HAS_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               std::shared_ptr<InterfaceReconstructor2D> ir)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir),
        num_tols_(num_tols) {}
#endif

  /// Constructor WITHOUT interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), num_tols_(num_tols) {}

  /// \brief Set the source mesh material that we have to intersect against

  void set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect target cell with a set of source cell
  /// \param[in] tgt_entity  Cell of target mesh to intersect
  /// \param[in] src_entities List of source cells to intersect against
  /// \return vector of Weights_t structure containing moments of intersection

  std::vector<Weights_t> operator() (int tgt_cell, std::vector<int> const& src_cells) const {
    std::vector<Wonton::Point<2>> target_poly;
    targetMeshWrapper.cell_get_coordinates(tgt_cell, &target_poly);

    bool trg_convex = poly2_is_convex(target_poly, num_tols_);

#ifdef DEBUG
    Wonton::Polytope<2> tpoly(target_poly);
    std::vector<double> tmom = tpoly.moments();
    if (tmom[0] < 0.0) {
      std::stringstream sstr;
      sstr << "intersect_r2d.h: Target cell " << tgt_cell << " has a negative volume";
      throw std::runtime_error(sstr.str());
    }
#endif
    
    int nsrc = src_cells.size();
    std::vector<Weights_t> sources_and_weights(nsrc);
    int ninserted = 0;
    for (int i = 0; i < nsrc; i++) {
      int s = src_cells[i];

      Weights_t & this_wt = sources_and_weights[ninserted];
      this_wt.entityID = s;

#ifdef PORTAGE_HAS_TANGRAM
      int nmats = sourceStateWrapper.cell_get_num_mats(s);
      std::vector<int> cellmats;
      sourceStateWrapper.cell_get_mats(s, &cellmats);

      if (!nmats || (matid_ == -1) || (nmats == 1 && cellmats[0] == matid_)) {
        // ---------- Intersection with pure cell ---------------
        // nmats == 0 -- no materials ==> single material
        // matid_ == -1 -- intersect with mesh not a particular material
        // nmats == 1 && cellmats[0] == matid -- intersection with pure cell
        //                                       containing matid

        std::vector<Wonton::Point<2>> source_poly;
        sourceMeshWrapper.cell_get_coordinates(s, &source_poly);

#ifdef DEBUG
        Wonton::Polytope<2> spoly(source_poly);
        std::vector<double> smom = spoly.moments();
        if (smom[0] < 0.0) {
          std::stringstream sstr;
          sstr << "intersect_r2d.h: Source cell " << s << " has a negative volume";
          throw std::runtime_error(sstr.str());
        }
#endif

        // If target polygon is convex we don't care if the source
        // polygon is non-convex because R3D can deal with it.
        if (trg_convex)
          this_wt.weights = intersect_polys_r2d(source_poly, target_poly,
                                                num_tols_);
        else {
          bool src_convex = poly2_is_convex(source_poly, num_tols_);
          if (src_convex)  // flip the order of the polygons
            this_wt.weights = intersect_polys_r2d(target_poly, source_poly,
                                                  num_tols_);
          else  // indicate that target is non-convex so it can be triangulated
            this_wt.weights = intersect_polys_r2d(target_poly, source_poly,
                                                  num_tols_, false);
        }
        
      } else {  // multi-material case
        // How can I check that I didn't get DummyInterfaceReconstructor
        //
        // static_assert(InterfaceReconstructorType<SourceMeshType, 2> !=
        //               DummyInterfaceReconstructor<SourceMeshType, 2>);

        assert(interface_reconstructor != nullptr);  // cannot be nullptr

        if (std::find(cellmats.begin(), cellmats.end(), matid_) !=
            cellmats.end()) {
          // mixed cell containing this material - intersect with
          // polygon approximation of this material in the cell
          // (obtained from interface reconstruction)

          Tangram::CellMatPoly<2> const& cellmatpoly =
              interface_reconstructor->cell_matpoly_data(s);
          std::vector<Tangram::MatPoly<2>> matpolys =
              cellmatpoly.get_matpolys(matid_);

          this_wt.weights.resize(3, 0.0);
          for (auto& matpoly : matpolys) {
            std::vector<Wonton::Point<2>> tpnts = matpoly.points();
            std::vector<Wonton::Point<2>> source_poly;
            source_poly.reserve(tpnts.size());
            for (auto const & p : tpnts) source_poly.push_back(p);

            std::vector<double> momvec;
            
            momvec = intersect_polys_r2d(source_poly, target_poly, num_tols_,
                                         trg_convex);

            for (int k = 0; k < 3; k++)
              this_wt.weights[k] += momvec[k];
          }
        }
      }
#else
      std::vector<Wonton::Point<2>> source_poly;
      sourceMeshWrapper.cell_get_coordinates(s, &source_poly);

#ifdef DEBUG
      Wonton::Polytope<2> spoly(source_poly);
      std::vector<double> smom = spoly.moments();
      if (smom[0] < 0.0) {
        std::stringstream sstr;
        sstr << "intersect_r2d.h: Source cell " << src_cell << " has a negative volume";
        throw std::runtime_error(sstr.str());
      }
#endif

      if (trg_convex)
        this_wt.weights = intersect_polys_r2d(source_poly, target_poly,
                                              num_tols_);
      else {
        src_convex = poly2_is_convex(source_poly, num_tols_);
        if (src_convex) // flip the two order of the polygons
          this_wt.weights = intersect_polys_r2d(target_poly, source_poly,
                                                num_tols_);
        else // indicate that target is non-convex so it can be triangulated
          this_wt.weights = intersect_polys_r2d(target_poly, source_poly,
                                                num_tols_, false);
      }        
#endif

      // Increment if vol of intersection > 0; otherwise, allow overwrite
      if (!this_wt.weights.empty() && this_wt.weights[0] > 0.0)
        ninserted++;
    }

    sources_and_weights.resize(ninserted);
    return sources_and_weights;
  }

  IntersectR2D() = delete;

  /// Assignment operator (disabled)
  IntersectR2D & operator = (const IntersectR2D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
  int matid_ = -1;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  NumericTolerances_t num_tols_;
};  // class IntersectR2D




//////////////////////////////////////////////////////////////////////////////
// Specialization of Intersect2D class for nodes

template <class SourceMeshType, class SourceStateType,
          class TargetMeshType,
          template <class, int, class, class> class InterfaceReconstructorType,
          class Matpoly_Splitter, class Matpoly_Clipper>
class IntersectR2D<Entity_kind::NODE, SourceMeshType, SourceStateType, TargetMeshType,
                   InterfaceReconstructorType, Matpoly_Splitter, Matpoly_Clipper> {

#ifdef PORTAGE_HAS_TANGRAM
  using InterfaceReconstructor2D =
      Tangram::Driver<InterfaceReconstructorType, 2, SourceMeshType,
                      Matpoly_Splitter, Matpoly_Clipper>;
#endif

 public:

#ifdef PORTAGE_HAS_TANGRAM
  /// Constructor with interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols,
               std::shared_ptr<InterfaceReconstructor2D> ir)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), interface_reconstructor(ir),
        num_tols_(num_tols) {}
#endif


  /// Constructor WITHOUT interface reconstructor

  IntersectR2D(SourceMeshType const & source_mesh,
               SourceStateType const & source_state,
               TargetMeshType const & target_mesh,
               NumericTolerances_t num_tols)
      : sourceMeshWrapper(source_mesh), sourceStateWrapper(source_state),
        targetMeshWrapper(target_mesh), num_tols_(num_tols) {}

  /// \brief Set the source mesh material that we have to intersect against

  void set_material(int m) {
    matid_ = m;
  }

  /// \brief Intersect control volume of a target node with control volumes of
  /// a set of source nodes
  /// \param[in] tgt_node  Target mesh node whose control volume we consider
  /// \param[in] src_nodes List of source nodes whose control volumes we will intersect against
  /// \return vector of Weights_t structure containing moments of intersection

  std::vector<Weights_t> operator() (int tgt_node, std::vector<int> const& src_nodes) const {
    std::vector<Wonton::Point<2>> target_poly;
    targetMeshWrapper.dual_cell_get_coordinates(tgt_node, &target_poly);

    bool trg_convex = poly2_is_convex(target_poly, num_tols_);

#ifdef DEBUG
    Wonton::Polytope<2> tpoly(target_poly);
    std::vector<double> tmom = tpoly.moments();
    if (tmom[0] < 0.0) {
      std::stringstream sstr;
      sstr << "intersect_r2d.h: Dual Cell for target node " << tgt_node << " has a negative volume";
      throw std::runtime_error(sstr.str());
    }
#endif
    
    
    int nsrc = src_nodes.size();
    std::vector<Weights_t> sources_and_weights(nsrc);
    int ninserted = 0;
    for (int i = 0; i < nsrc; i++) {
      int s = src_nodes[i];
      std::vector<Wonton::Point<2>> source_poly;
      sourceMeshWrapper.dual_cell_get_coordinates(s, &source_poly);

#ifdef DEBUG
      Wonton::Polytope<2> spoly(source_poly);
      std::vector<double> smom = spoly.moments();
      if (smom[0] < 0.0) {
        std::stringstream sstr;
        sstr << "intersect_r2d.h: Dual cell for source node " << s << " has a negative volume";
        throw std::runtime_error(sstr.str());
      }
#endif

      Weights_t & this_wt = sources_and_weights[ninserted];
      this_wt.entityID = s;
      if (trg_convex)
        this_wt.weights = intersect_polys_r2d(source_poly, target_poly,
                                              num_tols_);
      else {
        bool src_convex = poly2_is_convex(source_poly, num_tols_);
        if (src_convex) // flip the two order of the polygons
          this_wt.weights = intersect_polys_r2d(target_poly, source_poly,
                                                num_tols_);
        else // indicate that target is non-convex so it can be triangulated
          this_wt.weights = intersect_polys_r2d(target_poly, source_poly,
                                                num_tols_, false);
      }        

      // Increment if vol of intersection > 0; otherwise, allow overwrite
      if (!this_wt.weights.empty() && this_wt.weights[0] > 0.0)
        ninserted++;
    }

    sources_and_weights.resize(ninserted);
    return sources_and_weights;
  }

  IntersectR2D() = delete;

  /// Assignment operator (disabled)

  IntersectR2D & operator = (const IntersectR2D &) = delete;

 private:
  SourceMeshType const & sourceMeshWrapper;
  SourceStateType const & sourceStateWrapper;
  TargetMeshType const & targetMeshWrapper;
  int matid_ = -1;
#ifdef PORTAGE_HAS_TANGRAM
  std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  NumericTolerances_t num_tols_;
};  // class IntersectR2D

} // namespace Portage

#endif // PORTAGE_INTERSECT_INTERSECT_R2D_H_
