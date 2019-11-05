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
   * @brief Specialization for planar cell-based remap.
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
  class IntersectSweptFace<
    2, Entity_kind::CELL, SourceMesh, SourceState, TargetMesh,
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
     * @brief Retrieve the edges of current cell but skip the boundary ones.
     *
     * @param cell: the current cell.
     * @param filter_edges: list of its edges without the boundary ones.
     * @param filter_dirs: related filtered edge directions.
     */
    void get_filtered_edges(int cell,
                            std::vector<int>* filter_edges,
                            std::vector<int>* filter_dirs) const {
      filter_edges->clear();
      filter_dirs->clear();

      std::vector<int> edges, dirs, adj;
      source_mesh_.cell_get_faces_and_dirs(cell, &edges, &dirs);
      int const nb_edges = edges.size();

      for (int i = 0; i < nb_edges; ++i) {
        adj.clear();
        source_mesh_.face_get_cells(edges[i], Wonton::Entity_type::ALL, &adj);
        if (adj.size() == 2 or (adj.size() == 1 and adj[0] == cell)) {
          filter_edges->emplace_back(edges[i]);
          filter_dirs->emplace_back(dirs[i]);
        }
      }
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
        get_filtered_edges(source_id, &edges, &dirs);
        int const nb_edges = edges.size();

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
        std::map<int, std::pair<double,double>> swept_centroids;

        for (auto const& i : stencil) {
          swept_volumes[i] = 0.;
          swept_centroids[i] = std::make_pair(0.,0.);
        }

        // update source cell area
        swept_volumes[source_id] = source_mesh_.cell_volume(source_id);
        int nb_source_polys = 0;

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
           * - triangulate polygon
           * - compute each simplex orientation using determinant.
           *
           *     d_____c    k1:(a,b,d)
           *     /\   /     k2:(b,c,d)
           *    / \  /
           *   /  \ /                1    |ax  ay  1|
           *  /___\/       area(k1)= - det|bx  by  1| > 0 if counterclockwise
           * a     b                 2    |dx  dy  1|
           */
          double const& ax = swept_polygon[0][0];
          double const& ay = swept_polygon[0][1];
          double const& bx = swept_polygon[1][0];
          double const& by = swept_polygon[1][1];
          double const& cx = swept_polygon[2][0];
          double const& cy = swept_polygon[2][1];
          double const& dx = swept_polygon[3][0];
          double const& dy = swept_polygon[3][1];

          double const det[] = {
            ax * by - ax * dy - bx * ay + bx * dy + dx * ay - dx * by,
            bx * cy - bx * dy - cx * by + cx * dy + dx * by - dx * cy
          };

          // check that both triangles have the same orientation
          bool const both_positive = (det[0] >= 0 and det[1] >= 0);
          bool const both_negative = (det[0] < 0 and det[1] < 0);

          if (not both_positive and not both_negative) {
            std::cerr << "Error: twisted swept face polygon." << std::endl;
            source_weights.clear();
            return source_weights;
          }

          double const signed_area = 0.5 * (det[0] + det[1]);


          /* step 3: compute its centroid.
           * - compute the intersection point of its couple of diagonals.
           *
           *     d_____c    find (s,t) such that:
           *     /\   /     a + s(c - a) = b + t(d - b)
           *    / \  /
           *   /  \ /     resolve the equation: A.X = B
           *  /___\/      |dx-bx  ax-cx| |t| |ax-bx|
           * a     b      |dy-by  ay-cy| |s|=|ay-by|
           */
          Wonton::Point<2> centroid;

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
                assert(std::abs(x[0] - x[1]) < num_tols.polygon_convexity_eps);
                assert(std::abs(y[0] - y[1]) < num_tols.polygon_convexity_eps);
                centroid = Wonton::createP2(x[0], y[0]);

                std::cout << "a: ("<< ax <<" "<< ay <<"), ";
                std::cout << "c: ("<< cx <<" "<< cy <<"), ";
                std::cout << "param: " << param[0] << std::endl;
                std::cout << "bx: "<< bx << ", (dx - bx): "<< dx - bx << std::endl;
                std::cout << "by: "<< by << ", (dy - by): "<< dy - by << std::endl;
              #else
                centroid[0] = ax + param[0] * (cx - ax);
                centroid[1] = ay + param[0] * (cy - ay);
              #endif
            }
          }

          /* step 3: check aree sign, choose the right source cell,
           * and add computed area and centroid to related lists.
           */
          if (signed_area < 0.) {
            // if negative volume then accumulate to that of the source cell
            // it means that it would be substracted from that source cell.
            swept_volumes[source_id] += signed_area;
            swept_centroids[source_id].first  += centroid[0];
            swept_centroids[source_id].second += centroid[1];
            nb_source_polys++;

            #if DEBUG
              std::cout << "source_centroid["<< source_id <<"]: " << centroid;
              std::cout << ", area: "<< swept_volumes[source_id];
              std::cout << ", nb_source_polys: "<< nb_source_polys << std::endl;
            #endif
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
              swept_volumes[neigh] += signed_area;
              swept_centroids[neigh] = std::make_pair(centroid[0], centroid[1]);

              #if DEBUG
                std::cout << "neigh_centroid["<< neigh <<"]: " << centroid;
                std::cout << " of edge: ["<< swept_polygon[0] <<", ("<< bx <<" "<< by<<")]";
                std::cout << std::endl;
              #endif
            }
          }
        } // end for each edge of current cell

        #if DEBUG
          std::cout << " =========== " << std::endl;
        #endif

        // average swept face centroids on source cell.
        // this not accurate if the swept face polygon is not convex.
        if (swept_volumes[source_id] > 0. and nb_source_polys > 1) {
          swept_centroids[source_id].first  /= nb_source_polys;
          swept_centroids[source_id].second /= nb_source_polys;
        }

        // append computed volumes to 'source_weights'
        for (auto&& entry : swept_volumes) {
          auto const& id = entry.first;
          auto const& area = entry.second;
          if (std::abs(area) > 0.) {
            auto const& centroid = swept_centroids[id];
            std::vector<double> weight { area, centroid.first, centroid.second };
            source_weights.emplace_back(id, weight);
          }
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
    std::shared_ptr<InterfaceReconstructor2D> interface_reconstructor;
#endif
  }; // class IntersectSweptFace::CELL
/* -------------------------------------------------------------------------- */
} // namespace Portage
