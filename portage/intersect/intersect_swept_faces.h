/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#pragma once

#include "portage/support/portage.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
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
   * @tparam dim: dimension of the problem.
   * @tparam on_what: the entity kind we want to remap.
   * @tparam SourceMesh: the source mesh wrapper type.
   * @tparam SourceState: the source state wrapper to query field infos.
   * @tparam TargetMesh: the target mesh wrapper type.
   * @tparam InterfaceReconstructor: materials interface reconstructor type.
   * @tparam Matpoly_Splitter: material polygons splitter type.
   * @tparam Matpoly_Clipper: material polygons clipper type.
   */
  template<
    int dim,
    Entity_kind on_what,
    class SourceMesh,
    class SourceState,
    class TargetMesh,
    template<class, int, class, class>
      class InterfaceReconstructor = DummyInterfaceReconstructor,
    class Matpoly_Splitter = void,
    class Matpoly_Clipper = void
  >
  class IntersectSwept {

    // useful aliases
#ifdef HAVE_TANGRAM
    using InterfaceReconstructorDriver = Tangram::Driver<
      InterfaceReconstructor, dim, SourceMesh,
      Matpoly_Splitter, Matpoly_Clipper
    >;
#endif

  public:

    /**
     * @brief Default constructor (disabled).
     *
     */
    IntersectSwept() = delete;

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
    IntersectSwept(SourceMesh const &source_mesh,
                   SourceState const &source_state,
                   TargetMesh const &target_mesh,
                   NumericTolerances_t num_tols,
                   std::shared_ptr<InterfaceReconstructorDriver> ir)
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
    IntersectSwept(SourceMesh const &source_mesh,
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
    IntersectSwept& operator=(IntersectSwept const& other) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~IntersectSwept() = default;

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
    std::shared_ptr<InterfaceReconstructorDriver> interface_reconstructor;
#endif
  }; // class IntersectSwept


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
  class IntersectSwept<
    2, Entity_kind::CELL, SourceMesh, SourceState, TargetMesh,
    InterfaceReconstructor, Matpoly_Splitter, Matpoly_Clipper> {

    // useful aliases
#ifdef HAVE_TANGRAM
    using InterfaceReconstructorDriver = Tangram::Driver<
      InterfaceReconstructor, 2, SourceMesh,
      Matpoly_Splitter, Matpoly_Clipper
    >;
#endif

  public:

    /**
     * @brief Default constructor (disabled).
     *
     */
    IntersectSwept() = delete;

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
    IntersectSwept(SourceMesh const &source_mesh,
                   SourceState const &source_state,
                   TargetMesh const &target_mesh,
                   NumericTolerances_t num_tols,
                   std::shared_ptr<InterfaceReconstructorDriver> ir)
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
    IntersectSwept(SourceMesh const &source_mesh,
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
    IntersectSwept& operator=(IntersectSwept const& other) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~IntersectSwept() = default;

    /**
     * @brief Set the material we are operating on.
     *
     * @param m: the material ID
     */
    void set_material(int m) { material_id_ = m; }


    /**
     * @brief Retrieve the cell incident to a given face of a given cell.
     *
     * @param cell: the current cell.
     * @param face: the current face of the given cell.
     * @return the index of the incident cell to the given face or -1.
     */
    int get_face_incident_neigh(int cell, int face) const {
      std::vector<int> facecells;
      source_mesh_.face_get_cells(face, Entity_type::ALL, &facecells);
      int const nb_adj_cells = facecells.size();

      // Interfaces we need are always connected to two cells
      if (nb_adj_cells == 2) {
        int index = (facecells[0] == cell ? 1 : 0);
        return facecells[index];
      }
      return -1;
    }

    /**
     * @brief Perform the actual swept faces volumes computation.
     *
     * @param target_id: the current target cell index.
     * @param source_id: the related source cell index.
     * @param stencil: current source cell and its immediate neighbors.
     * @return: a list of swept faces volume and related source cell pair.
     */
    std::vector<Weights_t> operator()(int target_id,
                                      std::vector<int> const& stencil) const {
      // here the source and target cell have the exact ID.
      int const source_id = target_id;
      int const size = stencil.size();
      assert(std::find(stencil.begin(), stencil.end(), source_id) != stencil.end());

      std::vector<Weights_t> source_weights;

#ifdef HAVE_TANGRAM
      int const nb_mats = source_state_.cell_get_num_mats(source_id);
      bool const single_mat = not nb_mats or nb_mats == 1 or material_id_ == -1;

      if (single_mat) {
#endif
        std::vector<int> edges, dirs, nodes;

        // retrieve current source cell faces and related directions
        source_mesh_.cell_get_faces_and_dirs(source_id, &edges, &dirs);
        int const nb_dirs = dirs.size();
        int const nb_edges = edges.size();

        if (nb_edges != nb_dirs or nb_edges != size - 1) {
          std::cerr << "Error: invalid retrieve edges for cell "<< target_id;
          std::cerr << std::endl;
          source_weights.clear();
          return source_weights;
        }

        #ifdef DEBUG
          // ensure that we have the same face/edge index for source and target.
          std::vector<int> target_edges, target_dirs, target_nodes;
          target_mesh_.cell_get_faces_and_dirs(target_id, &target_edges, &target_dirs);
          int const nb_target_edges = target_edges.size();

          assert(nb_edges == nb_target_edges);
          for (int i = 0; i < nb_edges; ++i) {
            assert(edges[i] == target_edges[i]);
          }
        #endif

        std::map<int, double> swept_volumes;
        for (auto const& i : stencil) {
          swept_volumes[i] = 0.;
        }
        // update source cell area
        swept_volumes[source_id] = source_mesh_.cell_volume(source_id);

        for (int i = 0; i < nb_edges; ++i) {

          // step 0: retrieve nodes and reorder them according to edge direction
          nodes.clear();
          source_mesh_.face_get_nodes(edges[i], &nodes);
          int const nb_nodes = nodes.size();

          #ifdef DEBUG
            // ensure that we have the same nodal indices for source and target.
            target_mesh_.face_get_nodes(target_edges[i], &target_nodes);
            int const nb_target_nodes = target_nodes.size();

            assert(nb_nodes == nb_target_nodes);
            for (int j = 0; j < nb_nodes; ++j) {
              assert(nodes[i] == target_nodes[i]);
            }
          #endif

          // step 1: construct the swept face polygon
          std::vector<Wonton::Point<2>> swept_polygon(nb_nodes * 2);

          if (dirs[i] > 0) {
            // if the edge has the same orientation as the cell, then reverse
            // its nodes order such that we have a positive swept volume on
            // outside and negative swept volume on inside.
            source_mesh_.node_get_coordinates(nodes[1], swept_polygon.data());
            source_mesh_.node_get_coordinates(nodes[0], swept_polygon.data()+1);
            target_mesh_.node_get_coordinates(nodes[0], swept_polygon.data()+2);
            target_mesh_.node_get_coordinates(nodes[1], swept_polygon.data()+3);
          } else {
            // otherwise keep the same nodal order.
            source_mesh_.node_get_coordinates(nodes[0], swept_polygon.data());
            source_mesh_.node_get_coordinates(nodes[1], swept_polygon.data()+1);
            target_mesh_.node_get_coordinates(nodes[1], swept_polygon.data()+2);
            target_mesh_.node_get_coordinates(nodes[0], swept_polygon.data()+3);
          }

          /* step 2: compute its area then:
           * - split face into two triangles.
           * - compute and sum their areas using Heron's formula.
           *   split into two triangles and sum up their signed areas.
           *
           *     3___ c ___2      k:(a,e,d)=(0,1,3)
           *     /         |      k':(b,c,e)=(1,2,3)
           *    d    e     b
           *   /           |      area(k)=sqrt(s(s-a)(s-b)(s-c))
           *  /_____ a ____|      with s=(a+e+d)/2
           * 0             1
           */
          double const a = (swept_polygon[1] - swept_polygon[0]).norm();
          double const b = (swept_polygon[2] - swept_polygon[1]).norm();
          double const c = (swept_polygon[3] - swept_polygon[2]).norm();
          double const d = (swept_polygon[0] - swept_polygon[3]).norm();
          double const e = (swept_polygon[3] - swept_polygon[1]).norm();

          double const s[2] = { 0.5 * (a + e + d), 0.5 * (e + b + c) };
          double const area = std::sqrt(s[0] * (s[0] - a) * (s[0] - e) * (s[0] - d))
                            + std::sqrt(s[1] * (s[1] - e) * (s[1] - b) * (s[1] - c));

          // step 3: check sign and add to corresponding list.
          if (area < 0.) {
            // if negative volume then accumulate to that of the source cell
            // it means that it would be substracted from that source cell.
            swept_volumes[source_id] += area;
          } else {
            // retrieve the cell incident to the current edge.
            int const neigh = get_face_incident_neigh(source_id, edges[i]);

            // just skip in case of a boundary edge
            if (neigh < 0)
              continue;
            // check if incident cell belongs to the current stencil
            else if (not swept_volumes.count(neigh)) {
              std::cerr << "Error: invalid stencil for source cell "<< source_id;
              std::cerr << std::endl;
              source_weights.clear();
              return source_weights;
            } else {
              // update related area if ok
              swept_volumes[neigh] += area;
            }
          }
        }

        // append computed volumes to 'source_weights'
        for (auto&& entry : swept_volumes) {
          int const& id = entry.first;
          std::vector<double> weight { entry.second };
          source_weights.emplace_back(id, weight);
        }
        return source_weights;
#ifdef HAVE_TANGRAM
      } else /* multi-material case */ {
        std::cerr << "Error: multi-material swept face remap not yet supported";
        std::cerr << std::endl;
        source_weights.clear();
        return source_weights;
      }
#endif
    }

  private:
    SourceMesh const &source_mesh_;
    TargetMesh const &target_mesh_;
    SourceState const &source_state_;
    int material_id_ = -1;
    NumericTolerances_t num_tols_;

#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructorDriver> interface_reconstructor;
#endif
  }; // class IntersectSwept::CELL
/* -------------------------------------------------------------------------- */
} // namespace Portage
