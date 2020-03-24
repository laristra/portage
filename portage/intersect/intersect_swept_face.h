/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

#pragma once

#include "portage/support/portage.h"
#include "portage/intersect/dummy_interface_reconstructor.h"
#include "wonton/support/Point.h"
#include "wonton/support/Polytope.h"

#ifdef HAVE_TANGRAM
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/driver/driver.h"
  #include "tangram/support/MatPoly.h"
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

  private:
    SourceMesh const &source_mesh_;
    TargetMesh const &target_mesh_;
    SourceState const &source_state_;
    int material_id_ = -1;
    NumericTolerances_t num_tols_;
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
     * @brief Compute moments using facetized method:
     * - triangulate and compute each simplex orientation using determinant.
     * - compute the intersection point of its couple of diagonals.
     *
     *     d_____c    facets: (a,b,d) and (b,c,d)
     *     /\   /                   1    |ax  ay  1|
     *    / \  /      area(a,b,d) = - det|bx  by  1| > 0 if counterclockwise
     *   /  \ /                     2    |dx  dy  1|
     *  /___\/
     * a     b        centroid(a,b,d):  find (s,t) such that:
     *                                  a + s(c - a) = b + t(d - b)
     *
     *                                resolve the equation: A.X = B
     *                                |dx-bx  ax-cx| |t| |ax-bx|
     *                                |dy-by  ay-cy| |s|=|ay-by|
     *
     * @param swept_polygon: swept polygon points coordinates.
     * @return swept polygon moments.
     */
    std::vector<double> compute_moments_facetized_method
      (std::vector<Wonton::Point<2>> const& swept_polygon) const {

      std::vector<double> moment;

      // retrieve quadrilateral vertices coordinates.
      double const& ax = swept_polygon[0][0];
      double const& ay = swept_polygon[0][1];
      double const& bx = swept_polygon[1][0];
      double const& by = swept_polygon[1][1];
      double const& cx = swept_polygon[2][0];
      double const& cy = swept_polygon[2][1];
      double const& dx = swept_polygon[3][0];
      double const& dy = swept_polygon[3][1];

      /* step 1: compute signed area */
      double const det[] = {
        ax * by - ax * dy - bx * ay + bx * dy + dx * ay - dx * by,
        bx * cy - bx * dy - cx * by + cx * dy + dx * by - dx * cy
      };

      // check that both triangles have the same orientation
      bool const both_positive = (det[0] >= 0 and det[1] >= 0);
      bool const both_negative = (det[0] < 0 and det[1] < 0);

      if (not both_positive and not both_negative) {
        std::cerr << "Error: twisted swept face polygon." << std::endl;
        return moment;
      }

      double const area = 0.5 * (det[0] + det[1]);

      /* step 2: compute centroid */
      std::vector<double> centroid = { 0, 0 };

      // retrieve the determinant of A to compute its inverse A^-1.
      double const denom = (dx - bx) * (ay - cy) - (dy - by) * (ax - cx);

      // check if diagonals are not colinear
      if (std::abs(denom) > 0) {
        // check if given value is in [0,1]
        auto in_range = [](double x) -> bool { return 0. <= x and x <= 1.; };
        // compute the intersection point parameter:
        // compute the inverse of the matrix A and multiply it by the vector B
        double const param[] = {
          std::abs(((by - dy) * (ax - bx) + (dx - bx) * (ay - by)) / denom),
          std::abs(((ay - cy) * (ax - bx) + (ax - cx) * (ay - by)) / denom)
        };
        // check if diagonals intersection lies on their respective segments
        if (in_range(param[0]) and in_range(param[1])) {
          #if DEBUG
            double x[] = { ax + param[0] * (cx - ax), bx + param[1] * (dx - bx) };
            double y[] = { ay + param[0] * (cy - ay), by + param[1] * (dy - by) };
            assert(std::abs(x[0] - x[1]) < num_tols_.polygon_convexity_eps);
            assert(std::abs(y[0] - y[1]) < num_tols_.polygon_convexity_eps);
            centroid = { x[0], y[0] };
          #else
            centroid = { ax + param[0] * (cx - ax), ay + param[0] * (cy - ay) };
          #endif
          moment = { area, area * centroid[0], area * centroid[1] };
        }
      }
      return moment;
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

  public:
    /**
     * @brief Perform the actual swept faces computation.
     *
     * @param target_id: the current target cell index.
     * @param stencil: current source cell and its immediate neighbors.
     * @return: a list of swept faces moments and related source cell pair.
     */
    std::vector<Weights_t> operator()(int target_id,
                                      std::vector<int> const& stencil) const {
      // convenience function to check if a given cell is within the stencil.
      auto in_stencil = [&](int cell) -> bool {
        return std::find(stencil.begin(), stencil.end(), cell) != stencil.end();
      };

      // here the source and target cell have the exact ID.
      int const source_id = target_id;

      std::vector<Weights_t> swept_moments;

#ifdef HAVE_TANGRAM
      int const nb_mats = source_state_.cell_get_num_mats(source_id);
      bool const single_mat = not nb_mats or nb_mats == 1 or material_id_ == -1;

      if (single_mat) {
#endif
        std::vector<int> edges, dirs, nodes;

        // add source cell moments in the first place
        swept_moments.emplace_back(source_id, compute_source_moments(source_id));

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

        for (int i = 0; i < nb_edges; ++i) {

          // step 0: retrieve nodes and reorder them according to edge direction
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

          // step 1: construct the swept face polygon
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

          // compute swept polygon moment using divergence theorem
          auto moments = compute_moments_divergence_theorem(swept_polygon);

          /* step 3: assign the computed moments to the source cell or one
           * of its neighbors according to the sign of the swept face area.
           */
          if (std::abs(moments[0]) < num_tols_.min_absolute_volume) {
            // just skip if the swept polygon is almost flat.
            // it may occur when the cell is shifted only in one direction.
            continue;
          } else if (moments[0] < 0.) {
            // if the computed swept face area is negative then assign its
            // moments to the source cell: it will be substracted
            // from the source cell area when performing the interpolation.
            swept_moments.emplace_back(source_id, moments);
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
            // sanity check: ensure that swept face centroid remains
            // inside the neighbor cell.
            else if (not centroid_inside_cell(neigh, moments)) {
              throw std::runtime_error("invalid target mesh for swept face");
            }
            // append to list as current neighbor moment.
            else {
              swept_moments.emplace_back(neigh, moments);
            }
          }
        } // end for each edge of current cell

        return swept_moments;
#ifdef HAVE_TANGRAM
      } else /* multi-material case */ {
        throw std::runtime_error("multi-material case not yet supported");
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
    std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  }; // class IntersectSweptFace::2D::CELL

  /**
   * @brief Specialization for 3D cell-based remap.
   *
   * @tparam SourceMesh: the source mesh wrapper type.
   * @tparam SourceState: the source state wrapper to query field infos.
   * @tparam TargetMesh: the target mesh wrapper type.
   * @tparam InterfaceReconstructor: materials interface reconstructor type.
   * @tparam Matpoly_Splitter: material polyhedra splitter type.
   * @tparam Matpoly_Clipper: material polyhedra clipper type.
   */
  template<
    class SourceMesh, class SourceState, class TargetMesh,
    template<class, int, class, class>
    class InterfaceReconstructor,
    class Matpoly_Splitter, class Matpoly_Clipper
  >
  class IntersectSweptFace<3, Entity_kind::CELL,
    SourceMesh, SourceState,
    TargetMesh, InterfaceReconstructor,
    Matpoly_Splitter, Matpoly_Clipper> {

    // useful aliases
#ifdef HAVE_TANGRAM
    using InterfaceReconstructor3D = Tangram::Driver<
      InterfaceReconstructor, 3, SourceMesh,
      Matpoly_Splitter, Matpoly_Clipper
    >;
#endif

    using Polyhedron = Wonton::Polytope<3>;

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
    IntersectSweptFace(SourceMesh const& source_mesh,
                       SourceState const& source_state,
                       TargetMesh const& target_mesh,
                       NumericTolerances_t num_tols,
                       std::shared_ptr<InterfaceReconstructor3D> ir)
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
    IntersectSweptFace(SourceMesh const& source_mesh,
                       SourceState const& source_state,
                       TargetMesh const& target_mesh,
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
     * @brief Compute (swept) polyhedron moments using divergence theorem.
     *
     *       ______      - polyhedron faces are split into n triangles.
     *      /|    /|     - their vertices 'pi' are ordered counterclockwise.
     *     / |___/_|     - 'ni' denotes the normal vector of the triangle 'ti'.
     *    /  /  /  /     - 'vi' denotes the area of the triangle 'ti'.
     * ci___/_ /  /      - {e1,e2,e3} denotes the canonical basis of R^3.
     *  |  /  |  /
     *  | /   | /              1
     *  |/____|/           v = - sum_i=1^n vi.ni
     *  ai    bi  e3           6
     *            |             1            [ni.e1([(ai+bi).e1]^2 +[(bi+ci).e1]^2 +[(ci+ai).e1]^2)]
     *            |__ e2   c = --- sum_i=1^n [ni.e2([(ai+bi).e2]^2 +[(bi+ci).e2]^2 +[(ci+ai).e2]^2)]
     *           /             48v           [ni.e3([(ai+bi).e3]^2 +[(bi+ci).e3]^2 +[(ci+ai).e3]^2)]
     *          e1
     * @param coords: the polyhedron vertex coordinates.
     * @param faces: vertices index list of each face of the polyhedron ordered
     *               such that the normal vector of that face points outwards.
     * @return polyhedron moments: [volume, volume * centroid.x,
     *                                      volume * centroid.y,
     *                                      volume * centroid.z].
     */
    std::vector<double> compute_moments
      (std::vector<Wonton::Point<3>> const& coords,
       std::vector<std::vector<int>> const& faces) const {

      // implementation may change in the future
      Polyhedron polyhedron(coords, faces);
      return polyhedron.moments();
    }

    /**
     * @brief Compute the source cell moments.
     *
     * @param source_id: index of the source cell.
     * @return a list of the source cell moments.
     */
    std::vector<double> compute_source_moments(int source_id) const {

      Wonton::Point<3> centroid;
      double const volume = source_mesh_.cell_volume(source_id);
      source_mesh_.cell_centroid(source_id, &centroid);
      return std::vector<double>{volume,
                                 centroid[0] * volume,
                                 centroid[1] * volume,
                                 centroid[2] * volume};
    }

    /**
     * @brief Verify if the displacement is valid or not.
     *
     * Indeed the computed interpolation weights for the swept volume based
     * remap would be completely inaccurate, in case of large displacement due
     * to a large advection timestep for instance.
     * In fact, we should ideally check if the swept volume centroid               
     * related to each face of the source cell would lie in the first layer        
     * of the cell neighbors. However it would be really expensive since           
     * we have to split that face into triangles (even if the cell is convex),     
     * form a tetrahedron from each triangle and the given centroid,               
     * and check the signed volume of the resulting tetrahedron. And we have       
     * to do it for each cell face. Hence the suggested workaround is just to      
     * check that the volume of each swept region associated to one of             
     * the contributing cell does not exceed that of the source cell, and that
     * the absolute swept volume associated to the source cell does not exceed
     * one and a half of the volume of the source cell itself.
     * IMPORTANT: this slightly relaxed check is only valid for single material.
     *
     * @param source_id: the index of the source cell.
     * @param moments: the list of swept region moments.
     * @return true if the displacement is valid, false otherwise.
     */
    bool valid_displacement(int source_id,
                            std::vector<Weights_t> const& moments) const {

      int const nb_moments = moments.size();
      assert(nb_moments > 0);

      // the first entry of the swept region moment list should be
      // the source cell moment, that is its own volume and weighted centroid.
      assert(moments[0].entityID == source_id);
      double const source_cell_volume = moments[0].weights[0];
      double source_swept_volume = 0.;
      assert(source_cell_volume > num_tols_.min_absolute_volume);

      for (int i = 1; i < nb_moments; ++i) {
        double const swept_volume = std::abs(moments[i].weights[0]);
        double const entity_volume = source_mesh_.cell_volume(moments[i].entityID);
        assert(swept_volume > num_tols_.min_absolute_volume);
        // each swept region volume should not exceed that of the attached cell.
        if (swept_volume > entity_volume) {
          return false;
        } else if (moments[i].entityID == source_id) {
          // accumulate the self-contribution for further comparisons.
          source_swept_volume += swept_volume;
        }
      }

      // Normally, we should have p <= 3 swept regions associated with the source
      // cell, with p the number of dimensions involved by the displacement.
      // The volume of each swept region should not normally exceed half the
      // volume of the source cell itself. Hence to ensure that the displacement
      // is not too large, we just check that the total self-contribution
      // is less than twice the source cell volume.
      // Notice that this criterion could be refined later.
      return source_swept_volume <= (1.5 * source_cell_volume) + num_tols_.min_absolute_volume;
    }


  public:
    /**
     * @brief Perform the actual swept moments computation for the given cell.
     *
     * After decomposing the cell into faces, it constructs a swept polyhedron
     * for each face such that we have a positive volume if its centroid
     * lies outside the cell (and hence attached to an incident neighbor) and
     * a negative volume otherwise (and hence attached to the source cell).
     * The sign of the computed moments are actually related to that of the
     * fluxes associated to each face of the cell, and are computed using the
     * divergence theorem. To avoid degenerated cases and to limit the loss
     * of accuracy in case of large displacements, we perform a quick check
     * on the volume of each swept polyhedron such that they don't exceed a
     * certain threshold.
     * Notice that the stencil is not really required but still kept to be
     * consistent with the other intersection functors.
     *
     * @param target_id: the current target cell index.
     * @param stencil: current source cell and its immediate neighbors.
     * @return a list of swept region moments and related source cell pair.
     */
    std::vector<Weights_t> operator()(int target_id,
                                      std::vector<int> const& stencil) const {

      // convenience lambda to check if some cell is within the stencil.
      auto in_stencil = [&](int cell) -> bool {
        return std::find(stencil.begin(), stencil.end(), cell) != stencil.end();
      };

      // here the source and target cell have the exact ID.
      int const source_id = target_id;

      std::vector<Weights_t> swept_moments;

#ifdef HAVE_TANGRAM
      int const nb_mats = source_state_.cell_get_num_mats(source_id);
      bool const single_mat = not nb_mats or nb_mats == 1 or material_id_ == -1;

      if (single_mat) {
#endif
        std::vector<int> faces, dirs, nodes;

        // retrieve current source cell faces/edges and related directions
        source_mesh_.cell_get_faces_and_dirs(source_id, &faces, &dirs);
        int const nb_faces = faces.size();

        // add source cell moments in the first place
        swept_moments.emplace_back(source_id, compute_source_moments(source_id));

        #if DEBUG
          // ensure that we have the same face index for source and target.
          std::vector<int> target_faces, target_dirs, target_nodes;
          target_mesh_.cell_get_faces_and_dirs(target_id, &target_faces, &target_dirs);
          int const nb_target_faces = target_faces.size();

          assert(nb_faces == nb_target_faces);
          for (int j = 0; j < nb_faces; ++j) {
            assert(faces[j] == target_faces[j]);
          }
        #endif

        for (int i = 0; i < nb_faces; ++i) {
          // step 0: retrieve nodes and reorder them according to face orientation
          nodes.clear();
          source_mesh_.face_get_nodes(faces[i], &nodes);

          int const nb_face_nodes = nodes.size();
          int const nb_poly_nodes = 2 * nb_face_nodes;
          int const nb_poly_faces = nb_poly_nodes + 2;

          #if DEBUG
            // ensure that we have the same nodal indices for source and target.
            target_mesh_.face_get_nodes(target_faces[i], &target_nodes);
            int const nb_source_nodes = nodes.size();
            int const nb_target_nodes = target_nodes.size();

            assert(nb_source_nodes == nb_target_nodes);
            for (int j = 0; j < nb_target_nodes; ++j) {
              assert(nodes[j] == target_nodes[j]);
            }
          #endif

          /* step 1: construct the swept volume polyhedron which can be:
           * - a prism for a triangular face,
           * - a hexahedron for a quadrilateral face,
           * - a (n+2)-face polyhedron for an arbitrary n-polygon.
           */
          std::vector<Wonton::Point<3>> swept_poly_coords(nb_poly_nodes);
          std::vector<std::vector<int>> swept_poly_faces(nb_poly_faces);

          /* the swept polyhedron must be formed in a way that the vertices of
           * each of its faces are ordered such that their normals outside that
           * swept volume - this is necessarily true except for the original face
           * inherited from the cell itself.
           * for this face, if the ordering of its vertices is such that the normal
           * points out of the cell, then the normal points into the swept volume.
           * in such a case, the vertex ordering must be reversed.
           *
           *   source hex        target hex       face swept polyhedron:
           *                     7'......6'
           *                      .:    .:             4'......5'
           *    7______6         . :   . :             /:    /:
           *    /|    /|      4'...:..5' :            / :   / :
           *   / |   / |       :   :..:..2'          /  :  /  :  
           * 4 __|__5  |       :  .   :  .         4____:_5...:
           * |   3__|__2       : .    : .          |   /0'|  /1'
           * |  /   |  /       :......:            |  /   | /
           * | /    | /        0'     1'           | /    |/
           * |/_____|/                             |/_____/
           * 0      1                              0      1
           *                                 ∙dirs[f] > 0: [4,5,1,0,4',5',1',0']
           *                                 ∙dirs[f] < 0: [0,1,5,4,0',1',5',4']
           */
          bool const outward_normal = dirs[i] > 0;

          for (int current = 0; current < nb_face_nodes; ++current) {
            int const reverse = (nb_face_nodes - 1) - current;
            int const index   = (outward_normal ? reverse : current);
            int const offset  = current + nb_face_nodes;
            source_mesh_.node_get_coordinates(nodes[index], swept_poly_coords.data() + current);
            target_mesh_.node_get_coordinates(nodes[index], swept_poly_coords.data() + offset);
          }

          /* now build the swept polyhedron faces, which vertices are indexed
           * RELATIVELY to the polyhedron vertices list.
           * - first allocate memory for vertices list of each face.
           * - then add the original face and its twin induced by sweeping.
           * - eventually construct the other faces induced by edge sweeping.
           */
          for (int current = 0; current < nb_poly_faces; ++current) {
            // for each twin face induced by sweeping, its number of vertices is
            // exactly that of the current cell face, whereas the number of
            // vertices of the other faces is exactly 4.
            int const size = (current < 2 ? nb_face_nodes : 4);
            swept_poly_faces[current].resize(size);
          }


          /* swept polyhedron face construction rules:
           *
           *       3'_____2'     n_poly_faces: 2 + n_face_edges = 2 + 4 = 6.
           *       /|    /|      n_poly_nodes: 2 * n_face_edges = 2 * 4 = 8 = n.
           *      / |   / |      ordered vertex list: [3,2,1,0,3',2',1',0']
           *     /  |__/__|
           *    /  /  /  / 1'             absolute         relative
           *  3___/_2/  /        ∙f[0]: (3 |2 |1 |0 )      (0,1,2,3)
           *  |  /  |  /         ∙f[1]: (3'|2'|1'|0')      (4,5,6,7)
           *  | /   | /          ∙f[2]: (3 |3'|2'|2 )  =>  (0,4,5,1)
           *  |/____|/           ∙f[3]: (2 |2'|1'|1 )      (1,5,6,2)
           *  0     1            ∙f[4]: (1 |1'|0'|0 )      (2,6,7,3)
           *                     ∙f[5]: (0 |0'|3'|3 )      (3,7,4,0)
           *
           *  let m = n/2 with n the number of polyhedron vertices.
           *  - twin faces: [0, m-1] and [m-1, n].
           *  - side faces: [i, i+m, ((i+1) % m)+m, (i+1) % m]
           */
          for (int current = 0; current < nb_face_nodes; ++current) {
            // a) set twin faces vertices
            swept_poly_faces[0][current] = current;
            swept_poly_faces[1][current] = nb_poly_nodes - current - 1;

            // b) set side faces vertices while keeping them counterclockwise.
            int const index = current + 2;
            swept_poly_faces[index][0] = current;
            swept_poly_faces[index][3] = (current + 1) % nb_face_nodes;
            swept_poly_faces[index][1] = swept_poly_faces[index][0] + nb_face_nodes;
            swept_poly_faces[index][2] = swept_poly_faces[index][3] + nb_face_nodes;
          }

          /* step 2: compute swept polygon moments using divergence theorem */
          auto moments = compute_moments(swept_poly_coords, swept_poly_faces);

          /* step 3: assign the computed moments to the source cell or one
           * of its neighbors according to the sign of the swept region volume.
           */
          if (std::abs(moments[0]) < num_tols_.min_absolute_volume) {
            // just skip if the swept region is almost flat.
            // it may occur when the cell is shifted only in one direction.
            continue;
          } else if (moments[0] < 0.) {
            // if the computed swept region volume is negative then assign its
            // moments to the source cell: it will be substracted
            // from the source cell area when performing the interpolation.
            swept_moments.emplace_back(source_id, moments);
          } else {
            // retrieve the cell incident to the current edge.
            int const neigh = source_mesh_.cell_get_face_adj_cell(source_id, faces[i]);

            // just skip in case of a boundary edge
            if (neigh < 0) {
              continue;
            }
              // sanity check: ensure that incident cell belongs to the stencil.
            else if (not in_stencil(neigh)) {
              auto id = std::to_string(source_id);
              throw std::runtime_error("invalid stencil for source cell" + id);
            }
              // append to list as current neighbor moment.
            else {
              swept_moments.emplace_back(neigh, moments);
            }
          }
        } // end of for each face of current cell

        if (valid_displacement(source_id, swept_moments)) {
          return swept_moments;
        } else
          throw std::runtime_error("invalid displacement");

#ifdef HAVE_TANGRAM
      } else /* multi-material case */ {
        throw std::runtime_error("multi-material case not yet supported");
      }
#endif
    }

  private:
    SourceMesh const& source_mesh_;
    TargetMesh const& target_mesh_;
    SourceState const& source_state_;
    int material_id_ = -1;
    NumericTolerances_t num_tols_;
#ifdef HAVE_TANGRAM
    std::shared_ptr<InterfaceReconstructor3D> interface_reconstructor;
#endif
  }; // class IntersectSweptFace::3D::CELL

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
