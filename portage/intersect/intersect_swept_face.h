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
  #include "tangram/reconstruct/cutting_distance_solver.h"
#endif

/* -------------------------------------------------------------------------- */
namespace Portage {

  /**
   * @brief Kernel to compute interpolation weights for swept-face remap.
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
  class IntersectSweptFace {

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
    IntersectSweptFace() = delete;

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
    IntersectSweptFace(SourceMesh const &source_mesh,
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
    IntersectSweptFace(SourceMesh const &source_mesh,
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
    IntersectSweptFace& operator=(IntersectSweptFace const& other) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~IntersectSweptFace() = default;

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
                                      std::vector<int> const& stencil) const {
      // see specialization for cells
      std::cerr << "Sorry: current entity type not supported." << std::endl;
      return std::vector<Weights_t>();
    }

    #if DEBUG
      void enable_debug_prints() { verbose = true; }
    #endif

  private:
    SourceMesh const &source_mesh_;
    TargetMesh const &target_mesh_;
    SourceState const &source_state_;
    int material_id_ = -1;
    NumericTolerances_t num_tols_;
#if DEBUG
    bool verbose = false;
#endif

#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructorDriver> interface_reconstructor;
#endif
  }; // class IntersectSweptFace


  /**
   * @brief Specialization for 2D cell-based remap.
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
  class IntersectSweptFace<2, Entity_kind::CELL,
                           SourceMesh, SourceState,
                           TargetMesh, InterfaceReconstructor,
                           Matpoly_Splitter, Matpoly_Clipper> {

    // useful aliases
#ifdef HAVE_TANGRAM
    using InterfaceReconstructor2D = Tangram::Driver<
      InterfaceReconstructor, 2, SourceMesh,
      Matpoly_Splitter, Matpoly_Clipper>;
#endif

  public:

    /**
     * @brief Default constructor (disabled).
     *
     */
    IntersectSweptFace() = delete;

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
    IntersectSweptFace(SourceMesh const &source_mesh,
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
    IntersectSweptFace(SourceMesh const &source_mesh,
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
    IntersectSweptFace& operator=(IntersectSweptFace const& other) = delete;

    /**
     * @brief Destructor.
     *
     */
    ~IntersectSweptFace() = default;

    /**
     * @brief Set the material we are operating on.
     *
     * @param m: the material ID
     */
    void set_material(int m) { material_id_ = m; }

  private:

    /**
     * @brief Compute moments using divergence theorem.
     *
     *                       1
     *     p3_____ p2    a = - * sum_i=0^{n-1} det(pi pj) with j=(i+1) % n
     *      /    /           2
     *     /    /
     *    /    /             1 [sum_i=0^{n-1} (xi + xj) * det(pi pj)]
     *   /___ /          c = - [sum_i=0^{n-1} (yi + yj) * det(pi pj)]
     * p0    p1             6*a
     *
     * @param swept_polygon: swept polygon points coordinates.
     * @return swept polygon moments.
     */
    std::vector<double> compute_moments_divergence_theorem
      (std::vector<Wonton::Point<2>> const& swept_polygon) const {

      int const n = swept_polygon.size();
      constexpr double const factor[] = { 0.5, 1./6., 1./6. }; // evaluated at compile-time

      std::vector<double> moments(3);

      for (int i=0; i < n; ++i) {
        auto const& pi = swept_polygon[i];
        auto const& pj = swept_polygon[(i + 1) % n];
        auto const det = (pi[0] * pj[1] - pi[1] * pj[0]);
        moments[0] += det;
        moments[1] += (pi[0] + pj[0]) * det; // manually unrolled for perfs
        moments[2] += (pi[1] + pj[1]) * det;
      }

      for (int i=0; i < 3; ++i) { moments[i] *= factor[i]; }

      return moments;
    }

    /**
     * @brief Retrieve the source cell moments.
     *
     * @param source_id: index of the source cells.
     * @return a list of cell moments.
     */
    std::vector<double> compute_source_moments(int source_id) const {
      double const area = source_mesh_.cell_volume(source_id);
      Wonton::Point<2> centroid;
      source_mesh_.cell_centroid(source_id, &centroid);
      return std::vector<double>{ area, area * centroid[0], area * centroid[1] };
    }

    /**
     * @brief Check that given swept face centroid lies within the given cell.
     *
     * @param cell: given cell index
     * @param moment: swept face moments from which to extract centroid.
     * @return true if centroid is within the cell, false otherwise.
     */
    bool centroid_inside_cell(int cell, std::vector<double> const& moment) const {

      double const cell_area = std::abs(moment[0]);
      if (cell_area < num_tols_.min_absolute_volume)
        return false;

      std::vector<int> edges, dirs, nodes;
      source_mesh_.cell_get_faces_and_dirs(cell, &edges, &dirs);
      int const nb_edges = edges.size();

      Wonton::Point<2> triangle[3];
      triangle[0] = { moment[1] / cell_area, moment[2] / cell_area };

      for (int i = 0; i < nb_edges; ++i) {
        source_mesh_.face_get_nodes(edges[i], &nodes);
        // set triangle according to edge orientation.
        if (dirs[i] > 0) {
          source_mesh_.node_get_coordinates(nodes[0], triangle + 1);
          source_mesh_.node_get_coordinates(nodes[1], triangle + 2);
        } else {
          source_mesh_.node_get_coordinates(nodes[1], triangle + 1);
          source_mesh_.node_get_coordinates(nodes[0], triangle + 2);
        }

        // check determinant sign
        double const& ax = triangle[0][0];
        double const& ay = triangle[0][1];
        double const& bx = triangle[1][0];
        double const& by = triangle[1][1];
        double const& cx = triangle[2][0];
        double const& cy = triangle[2][1];
        double const det = ax * by - ax * cy
                         - bx * ay + bx * cy
                         + cx * ay - cx * by;
        if (det < 0.)
          return false;
      }
      return true;
    }

#ifdef HAVE_TANGRAM
    /**
     * @brief For a given cell and its face finds the moments associated with the
     * intersection of the swept region and MatPoly's with material_id_ that belong
     * to the corresponding face group
     *
     * @param cell_id: index of the cell.
     * @param face_group_id: index of the cell's face group.
     * @param swept_volume: volume to clip off the associated triangle 
     * in the cell's decomposition.
     * @return Intersection moments for the material_id_ 
     */
    std::vector<double> compute_face_group_moments(
      int const cell_id,
      int const face_group_id,
      double const swept_volume) const {
      
      std::vector<int> cfaces, cfdirs;
      source_mesh_.cell_get_faces_and_dirs(cell_id, &cfaces, &cfdirs);
      int nfaces = cfaces.size();
      int cface_id = std::distance(
        cfaces.begin(), std::find(cfaces.begin(), cfaces.end(), face_group_id));
      //Face group should be associated with one of the cell's faces
      assert(cface_id != nfaces);

      //Retrieve tolerance used by the interface reconstructor
      const std::vector<Tangram::IterativeMethodTolerances_t>& ims_tols = 
        interface_reconstructor->iterative_methods_tolerances();
      double dst_tol = ims_tols[0].arg_eps;
      double vol_tol = ims_tols[0].fun_eps;

      //Create a MatPoly for the cell
      Tangram::MatPoly<2> cell_mp;
      Tangram::cell_get_matpoly(source_mesh_, cell_id, &cell_mp, dst_tol);
      //Get the face normal and MatPoly's in the face's group
      std::vector<Tangram::MatPoly<2>> face_group_polys;
      Tangram::Plane_t<2> cutting_plane;
      cutting_plane.normal = 
        cell_mp.face_normal_and_group(cface_id, face_group_id, &face_group_polys);

      //Find the cutting distance for the given swept volume
      Tangram::CuttingDistanceSolver<2, Matpoly_Clipper> cds(face_group_polys,
        cutting_plane.normal, ims_tols[0], true);

      cds.set_target_volume(swept_volume);
      std::vector<double> cds_res = cds();

      //Check if we had enough volume in the face group
      if (cds_res[1] < swept_volume - vol_tol) {
        throw std::runtime_error("Mesh displacement is too big for the implemented swept-face method");
      }

      cutting_plane.dist2origin = cds_res[0];

      //Get the face group MatPoly's with material_id_ from the reconstructor
      Tangram::CellMatPoly<2> const& cellmatpoly =
        interface_reconstructor->cell_matpoly_data(cell_id);
      std::vector<Tangram::MatPoly<2>> group_mat_polys = 
        cellmatpoly.get_face_group_matpolys(face_group_id, material_id_);

      //Clip the MatPoly's with the plane to get moments
      Matpoly_Clipper clip_matpolys(vol_tol);
      clip_matpolys.set_matpolys(group_mat_polys, true);
      clip_matpolys.set_plane(cutting_plane);
      std::vector<double> moments = clip_matpolys();

      return moments;
    }
#endif

  public:
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
      // convenience function to check if a given cell is within the stencil.
      auto in_stencil = [&](int cell) -> bool {
        return std::find(stencil.begin(), stencil.end(), cell) != stencil.end();
      };

      // here the source and target cell have the same ID.
      int const source_id = target_id;

      std::vector<Weights_t> swept_moments;

      // Step 1: Add the source cell moments in the first place. 
      // For the single material case, add the moments of the source cell id.
      // For the multimaterial case, add only the moments for the material
      // the intersector is working on. 
#ifdef HAVE_TANGRAM
      int const nb_mats = source_state_.cell_get_num_mats(source_id);
      std::vector<int> cellmats;
      source_state_.cell_get_mats(source_id, &cellmats);
      // nb_mats == 0 -- no materials ==> single material
      // material_id_ == -1 -- intersect with mesh not a particular material
      // nb_mats == 1 && cellmats[0] == material_id_ -- intersection with pure cell
      //                                                containing material_id_      
      bool const source_cell_mat = 
        (std::find(cellmats.begin(), cellmats.end(), material_id_) != cellmats.end());
      bool const single_mat = !nb_mats || (material_id_ == -1) || 
                              (nb_mats == 1 && source_cell_mat);

      if (single_mat) {
#endif
        // add source cell moments in the first place
        swept_moments.emplace_back(source_id, compute_source_moments(source_id));

#ifdef HAVE_TANGRAM
      } else if (source_cell_mat) {
        // mixed cell should contain this material
        assert(interface_reconstructor != nullptr);

        // add source cell moments in the first place: 
        // obtain the aggregated moments of all MatPoly's with material_id_
        Tangram::CellMatPoly<2> const& cellmatpoly =
          interface_reconstructor->cell_matpoly_data(source_id);

        swept_moments.emplace_back(source_id, cellmatpoly.material_moments(material_id_));
      }
#endif

      // Step 2: Obtain all facets and normals of the source cell
      std::vector<int> edges, dirs, nodes;

      // retrieve current source cell faces/edges and related directions
      source_mesh_.cell_get_faces_and_dirs(source_id, &edges, &dirs);
      int const nb_edges = edges.size();

#if DEBUG
      // ensure that we have the same face/edge index for source and target.
      std::vector<int> target_edges, target_dirs, target_nodes;
      target_mesh_.cell_get_faces_and_dirs(target_id, &target_edges, &target_dirs);
      int const nb_target_edges = target_edges.size();

      assert(nb_edges == nb_target_edges);
      for (int j = 0; j < nb_edges; ++j) {
        assert(edges[j] == target_edges[j]);
      }
#endif

      // Step 3: Main loop over all facets of the source cell
      for (int i = 0; i < nb_edges; ++i) {

        // step 3a: retrieve nodes and reorder them according to edge direction
        nodes.clear();
        source_mesh_.face_get_nodes(edges[i], &nodes);

#if DEBUG
        // ensure that we have the same nodal indices for source and target.
        target_mesh_.face_get_nodes(target_edges[i], &target_nodes);
        int const nb_source_nodes = nodes.size();
        int const nb_target_nodes = target_nodes.size();

        assert(nb_source_nodes == nb_target_nodes);
        for (int j = 0; j < nb_target_nodes; ++j) {
          assert(nodes[j] == target_nodes[j]);
        }
#endif

        // step 3b: construct the swept face polygon
        std::vector<Wonton::Point<2>> swept_polygon(4);

        // if the edge has the same orientation as the cell, then reverse
        // its nodes order such that we have a positive swept volume on
        // outside and negative swept volume on inside.
        // otherwise keep the same nodal order.
        unsigned const j = (dirs[i] > 0 ? 1 : 0);
        unsigned const k = j ^ 1;

        source_mesh_.node_get_coordinates(nodes[j], swept_polygon.data());
        source_mesh_.node_get_coordinates(nodes[k], swept_polygon.data()+1);
        target_mesh_.node_get_coordinates(nodes[k], swept_polygon.data()+2);
        target_mesh_.node_get_coordinates(nodes[j], swept_polygon.data()+3);

        // step 3c: compute swept polygon moment using divergence theorem
        auto moments = compute_moments_divergence_theorem(swept_polygon);

        // step 3d: add the source cells and their correct swept-moments to 
        // the weights vector.  Approaches to compute the amount
        // of swept-moment to be added are different for single and multimaterial cases.

        if (std::fabs(moments[0]) < num_tols_.min_absolute_volume) {
          // skip if the swept polygon is almost empty.
          // it may occur when the cell is shifted only in one direction.
          continue;
        } else if (moments[0] < 0.0) {
          // if the computed swept face area is negative then assign its
          // moments to the source cell: it will be substracted
          // from the source cell area when performing the interpolation.
#ifdef HAVE_TANGRAM
          if (single_mat) {
#endif          
            swept_moments.emplace_back(source_id, moments);
#ifdef HAVE_TANGRAM
          } else if (source_cell_mat) {
            swept_moments.emplace_back(source_id, 
              compute_face_group_moments(source_id, edges[i], std::fabs(moments[0])));
          }
#endif            
        } else {
          // retrieve the cell incident to the current edge.
          int const neigh = source_mesh_.cell_get_face_adj_cell(source_id, edges[i]);

          // just skip in case of a boundary edge
          if (neigh < 0) {
            continue;
          }
          // sanity check: ensure that incident cell belongs to the stencil.
          else if (not in_stencil(neigh)) {
            auto id = std::to_string(source_id);
            throw std::runtime_error("invalid stencil for source cell "+ id);
          }
#if DEBUG            
          // sanity check: ensure that swept face centroid remains
          // inside the neighbor cell.
          else if (not centroid_inside_cell(neigh, moments)) {
            throw std::runtime_error("invalid target mesh for swept face");          
          }
#endif            
          // append to list as current neighbor moment.
          else {
#ifdef HAVE_TANGRAM
            if (single_mat) {
#endif
              swept_moments.emplace_back(neigh, moments);
#ifdef HAVE_TANGRAM
            } else {
              //Skip if the neighboring cell doesn't contain material_id_
              std::vector<int> neighmats;
              source_state_.cell_get_mats(neigh, &neighmats);
              if (std::find(neighmats.begin(), neighmats.end(), material_id_) ==
                  neighmats.end())
                continue;
              
              //Compute and append moments for the neighbor
              swept_moments.emplace_back(neigh, 
                compute_face_group_moments(neigh, edges[i], std::fabs(moments[0])));
            }
#endif  
          }
        }
      } // end for each edge of current cell
     
      return swept_moments;
    }

  private:
    SourceMesh const &source_mesh_;
    TargetMesh const &target_mesh_;
    SourceState const &source_state_;
    int material_id_ = -1;
    NumericTolerances_t num_tols_;
#if DEBUG
    bool verbose = false;
#endif

#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  }; // class IntersectSweptFace::2D::CELL


  /* Define aliases to have the same interface as the other intersectors such
   * as R2D and R3D. This simple workaround allow us to discard the dimension
   * template parameter when instantiating the kernel.
   */
  template<
    Entity_kind entity_kind,
    class SourceMesh, class SourceState, class TargetMesh,
    template<class, int, class, class>
      class InterfaceReconstructor = DummyInterfaceReconstructor,
    class Matpoly_Splitter = void, class Matpoly_Clipper = void
  >
  using IntersectSweptFace2D = IntersectSweptFace<2, entity_kind,
                                                  SourceMesh, SourceState,
                                                  TargetMesh,
                                                  InterfaceReconstructor,
                                                  Matpoly_Splitter,
                                                  Matpoly_Clipper>;

  template<
    Entity_kind entity_kind,
    class SourceMesh, class SourceState, class TargetMesh,
    template<class, int, class, class>
      class InterfaceReconstructor = DummyInterfaceReconstructor,
    class Matpoly_Splitter = void, class Matpoly_Clipper = void
  >
  using IntersectSweptFace3D = IntersectSweptFace<3, entity_kind,
                                                  SourceMesh, SourceState,
                                                  TargetMesh,
                                                  InterfaceReconstructor,
                                                  Matpoly_Splitter,
                                                  Matpoly_Clipper>;

  /* ------------------------------------------------------------------------ */
} // namespace Portage
