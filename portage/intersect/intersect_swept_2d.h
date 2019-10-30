/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#pragma once

#include "portage/support/portage.h"
#include "wonton/support/Point.h"

#ifdef HAVE_TANGRAM
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/driver/driver.h"
  #include "tangram/support/MatPoly.h"
#endif
/* -------------------------------------------------------------------------- */
namespace Portage {

  /**
   * @class IntersectSweptFace2D intersect_swept_face_2d.h
   * @brief Kernel to compute swept faces volumes for planar advection-based remap.
   *
   * @tparam on_what: the entity kind we want to remap.
   * @tparam SourceMesh: the source mesh wrapper type.
   * @tparam SourceState: the source state wrapper to query field infos.
   * @tparam TargetMesh: the target mesh wrapper type.
   * @tparam InterfaceReconstructor: materials interface reconstructor type.
   * @tparam Matpoly_Splitter: material polygons splitter type.
   * @tparam Matpoly_Clipper: material polygons clipper type.
   */
  template<
    Entity_kind on_what,
    class SourceMesh,
    class SourceState,
    class TargetMesh,
    template<class, int, class, class>
      class InterfaceReconstructor = DummyInterfaceReconstructor,
    class Matpoly_Splitter = void,
    class Matpoly_Clipper = void
  >
  class IntersectSwept2D {

    // useful aliases
#ifdef HAVE_TANGRAM
    using InterfaceReconstructor2D = Tangram::Driver<
      InterfaceReconstructor, 2, SourceMesh,
      Matpoly_Splitter, Matpoly_Clipper
    >;
#endif

  public:

    /**
     * @brief Default constructor (disabled).
     *
     */
    IntersectSwept2D() = delete;

#ifdef HAVE_TANGRAM

    /**
     * @brief Constructor for multi-material case.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] num_tols: numerical tolerances.
     * @param[in] ir: interface reconstructor for querying matpolys on source mesh.
     */
    IntersectSwept2D(SourceMesh const &source_mesh,
                     SourceState const &source_state,
                     TargetMesh const &target_mesh,
                     NumericTolerances_t num_tols,
                     std::shared_ptr<InterfaceReconstructor2D> ir)
      : source_mesh_(source_mesh),
        source_state_(source_state),
        target_mesh_(target_mesh),
        num_tols_(num_tols),
        interface_reconstructor(ir) {}

#endif

    /**
     * @brief Constructor for single material case.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] num_tols: numerical tolerances.
     */
    IntersectSwept2D(SourceMesh const &source_mesh,
                     SourceState const &source_state,
                     TargetMesh const &target_mesh,
                     NumericTolerances_t num_tols)
      : source_mesh_(source_mesh),
        source_state_(source_state),
        target_mesh_(target_mesh),
        num_tols_(num_tols) {}

    /**
     * @brief Assignment operator (disabled).
     *
     * @param[in] other: the intersector to copy.
     * @return current intersector reference.
     */
    IntersectSwept2D& operator=(IntersectSwept2D const& other) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~IntersectSwept2D() = default;

    /**
     * @brief Set the material we are operating on.
     *
     * @param m: the material ID
     */
    void set_material(int m) { material_id_ = m; }

    /**
     * @brief Perform the actual swept faces volumes computation.
     *
     * @param target_id: the current target cell index.
     * @param source_id: the related source cell index.
     * @param stencil: current source cell and its immediate neighbors.
     * @return: a list of swept faces volume and related source cell pair.
     */
    std::vector<Weights_t> operator()(int target_id,
                                      int source_id,
                                      std::vector<int> const& stencil) const {
      // see specialization for cells
      std::cerr << "Sorry: current entity type not supported." << std::endl;
    }

  private:
    SourceMesh const &source_mesh_;
    TargetMesh const &target_mesh_;
    SourceState const &source_state_;
    int material_id_ = -1;
    NumericTolerances_t num_tols_;

#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  }; // class IntersectSwept2D


  /**
   * @brief Specialization for cell-based remap.
   *
   * @tparam SourceMesh: the source mesh wrapper type.
   * @tparam SourceState: the source state wrapper to query field infos.
   * @tparam TargetMesh: the target mesh wrapper type.
   * @tparam InterfaceReconstructor: materials interface reconstructor type.
   * @tparam Matpoly_Splitter: material polygons splitter type.
   * @tparam Matpoly_Clipper: material polygons clipper type.
   */
  template<
    class SourceMesh, class SourceState, class TargetMesh,
    template<class, int, class, class>
      class InterfaceReconstructor,
    class Matpoly_Splitter, class Matpoly_Clipper
  >
  class IntersectSwept2D<
    Entity_kind::CELL, SourceMesh, SourceState, TargetMesh,
    InterfaceReconstructor, Matpoly_Splitter, Matpoly_Clipper> {

    // useful aliases
#ifdef HAVE_TANGRAM
    using InterfaceReconstructor2D = Tangram::Driver<
      InterfaceReconstructor, 2, SourceMesh,
      Matpoly_Splitter, Matpoly_Clipper
    >;
#endif

  public:

    /**
     * @brief Default constructor (disabled).
     *
     */
    IntersectSwept2D() = delete;

#ifdef HAVE_TANGRAM

    /**
     * @brief Constructor for multi-material case.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] num_tols: numerical tolerances.
     * @param[in] ir: interface reconstructor for querying matpolys on source mesh.
     */
    IntersectSwept2D(SourceMesh const &source_mesh,
                     SourceState const &source_state,
                     TargetMesh const &target_mesh,
                     NumericTolerances_t num_tols,
                     std::shared_ptr<InterfaceReconstructor2D> ir)
      : source_mesh_(source_mesh),
        source_state_(source_state),
        target_mesh_(target_mesh),
        num_tols_(num_tols),
        interface_reconstructor(ir) {}

#endif

    /**
     * @brief Constructor for single material case.
     *
     * @param[in] source_mesh: mesh wrapper used to query source mesh info.
     * @param[in] source_state: state-manager wrapper used to query field info.
     * @param[in] target_mesh: mesh wrapper used to query target mesh info.
     * @param[in] num_tols: numerical tolerances.
     */
    IntersectSwept2D(SourceMesh const &source_mesh,
                     SourceState const &source_state,
                     TargetMesh const &target_mesh,
                     NumericTolerances_t num_tols)
      : source_mesh_(source_mesh),
        source_state_(source_state),
        target_mesh_(target_mesh),
        num_tols_(num_tols) {}

    /**
     * @brief Assignment operator (disabled).
     *
     * @param[in] other: the intersector to copy.
     * @return current intersector reference.
     */
    IntersectSwept2D& operator=(IntersectSwept2D const& other) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~IntersectSwept2D() = default;

    /**
     * @brief Set the material we are operating on.
     *
     * @param m: the material ID
     */
    void set_material(int m) { material_id_ = m; }

    /**
     * @brief Perform the actual swept faces volumes computation.
     *
     * @param target_id: the current target cell index.
     * @param source_id: the related source cell index.
     * @param stencil: current source cell and its immediate neighbors.
     * @return: a list of swept faces volume and related source cell pair.
     */
    std::vector<Weights_t> operator()(int target_id,
                                      int source_id,
                                      std::vector<int> const& stencil) const {
      // TODO
    }

  private:
    SourceMesh const &source_mesh_;
    TargetMesh const &target_mesh_;
    SourceState const &source_state_;
    int material_id_ = -1;
    NumericTolerances_t num_tols_;

#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  }; // class IntersectSwept2D::CELL
/* -------------------------------------------------------------------------- */
} // namespace Portage
