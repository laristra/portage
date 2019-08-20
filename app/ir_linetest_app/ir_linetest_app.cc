/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>
#include <algorithm>
#include <stdexcept>

// portage includes
#include "portage/search/search_simple.h"
#include "portage/support/portage.h"
extern "C" {
#include "wonton/intersect/r3d/r2d.h"
}
#include "portage/intersect/intersect_r2d.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

// tangram includes
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"

#ifdef PORTAGE_ENABLE_MPI
#include "mpi.h"
#endif

using Wonton::Point;
/*!

 @brief Integration test for Portage and Tangram that uses XMOF2D interface reconstructor.
 The domain is a unit square with two material separated by a linear interface.
 Source and target meshes are rectangular (SimpleMesh2D) with dimensions specified
 by a user.
 For a given material interface the reference volume fractions and centroids are computed
 on the source mesh using R2D. Tangram driver is used to perform interface reconstuction
 and obtain a vector of CellMatPoly objects for multi-material cells. The reconstructor
 used by Tangram is XMOF2D.
 For target cells intersecting with multi-material source cells volume fractions and
 centroids of their intersections with single-material polygons in respective CellMatPoly
 objects are computed, and for both materials their overall volume fractions and centroids
 are found in every target cell.
 The obtained volume fractions and centroids are used to perform interface reconstruction
 on the target mesh using Tangram driver with XMOF2D as a reconstructor.
 For every multi-material cell, the Hausdorff distance between the reconstructed and
 reference material interfaces is computed, the max distance over all multi-material cells
 is then compared to a given tolerance to determine if the interface reconstruction was
 performed correctly.
 Requires Tangram and XMOF2D libraries to be linked.
*/


/* Refence material interface is given in the form
   y = mat_int_a*x + mat_int_b
   Point (Px,Py) is below the line if
   Py < mat_int_a*Px + mat_int_b
*/
const double mat_int_a = 1.5;
const double mat_int_b = -0.25;

/*
   Below the line is material 1, above is material 0
*/
const int mat_id_above = 0;
const int mat_id_below = 1;
const int NUM_MATS = 2;
const int NUM_MOMENTS=3;

const double IR_tol = 1.0e-8;     //Max allowed Hausdorff distance between
                                  //reference and reconstructed material interfaces
const double deps = 1.0e-15;        //Distance tolerance
const double seps = 1.0e-14;        //Size(area) tolerance
const double denom_eps = 1.0e-6;  //Denominator tolerance

/*!
 @brief Convert the line given by y = line_a*x + line_b to the R2D format.
 @param[in] line_a Slope of the line.
 @param[in] line_b y-intercept of the line.
 @param[out] r2d_line Corresponding line in R2D format
*/
void R2DizeLine(const double line_a,
                const double line_b,
                r2d_plane& r2d_line);

/*!
 @brief Convert a mesh cell into an R2D polygon.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] cell_id Index of the cell to convert.
 @param[out] r2d_polygon Corresponding polygon in R2D format
*/
template <class Mesh_Wrapper>
void R2DizeCell(const Mesh_Wrapper& Mesh,
                const int cell_id,
                r2d_poly& r2d_polygon);
/*!
 @brief Determines the position of a mesh cell with respect to a line.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] cell_id Index of the cell to test for position.
 @param[in] line_a Slope of the reference material interface.
 @param[in] line_b y-intercept of the the reference material interface.
 @param[in] yeps Max difference in y-coordinate from the value given by the line equation
                 for points on the line.
 @return -1 if the cell is below the line, 1 if the cell is above the line,
          0 if the line passes through the interior of the cell
*/
template <class Mesh_Wrapper>
int CellPosition(const Mesh_Wrapper& Mesh,
                 const int cell_id,
                 const double line_a,
                 const double line_b,
                 const double yeps = deps);

/*!
 @brief Computes material data for a mesh and a given linear material interface using R2D.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] line_a Slope of the line.
 @param[in] line_b y-intercept of the line.
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
*/
template <class Mesh_Wrapper>
void get_materials_data(const Mesh_Wrapper& Mesh,
                        const double line_a,
                        const double line_b,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Point<2>>& cell_mat_centroids,
                        std::vector<int>& offsets);
/*!
 @brief Determines the position of a point with respect to a line.
 @param[in] pt Point coordinates.
 @param[in] line_nvec Normal vector to the line.
 @param[in] line_pt Point on the line.
 @param[in] eps Max distance between the line and points considered to be on the line.
 @return -1 if the point is below the line, 1 if the point is above the line,
          0 if the point is on the line
*/


int PointPosition(const Wonton::Point2& pt,
                  const Wonton::Vector<2>& line_nvec,
                  const Wonton::Point2& line_pt,
                  const double eps = deps);

/*!
 @brief Computes the intersection of two non-parallel lines.
 @param[in] lines_pts Two pairs of points defining the intersecting lines.
 @param[in] den_eps Minimal value of the denominator for which the analytical formula is
                    to be used, if the denominator is smaller, the bisection algorithm
                    is used instead.
 @param[in] eps Max distance between points that are considered coincident.
 @return Intersection point
*/
Point<2> LinesIntersect(std::vector< std::vector<Point<2>> > lines_pts,
                               double den_eps,
                               double eps = deps);
/*!
 @brief Computes Hausdorff distance between two linear segments.
 @param[in] segments Two pairs of points defining the linear segments.
 @return Hausdorff distance
*/
double SegmentsDistance(const std::vector<std::vector<Point<2>>>& segments);

/*!
 @brief Computes the max Hausdorff distance between the reference and reconstructed
        interfaces over the multi-material cells in a mesh.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] line_a Slope of the reference material interface.
 @param[in] line_b y-intercept of the the reference material interface.
 @param[in] cell_num_mats Reference number of materials for every mesh cell.
 @param[in] cellmatpoly_list Vector of pointers to CellMatPoly objects corresponding to
                             the results of interface reconstruction
 @param[in] eps Max distance between points that are considered coincident.
 @return Max Hausdorff distance
*/
template <class Mesh_Wrapper>
double get_ir_error(const Mesh_Wrapper& Mesh,
                    const double line_a,
                    const double line_b,
                    const std::vector<int>& cell_num_mats,
                    const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
                    cellmatpoly_list,
                    const double eps = deps);

/*!
 @brief Intersects material polygons of source cell with target cell and accumulates
        material moment data.
 @param[in/out] mat_moments 2D vector containing moment data by material (indexed first by
             mat_id and then by moment).
 @param[in] mat_id Material ID corresponding to source_points.
 @param[in] source_points Points corresonding to some region in source mesh (could
            be an entire cell or a matpoly within a cell).
 @param[in] target_points Points corresonding to a cell in the target mesh.
 @param[in/out] tcell_num_ids Indices of materials in each mesh cell, a flat vector, requires
            computations of offsets.
 @param[in/out] tcell_num_mats Reference number of materials for every target mesh cell.
 @param[in] tc Target cell ID.
 @param[in] idx Offset into cell_mat_ids giving start of data for current target cell.
*/
void add_intersect_data(std::vector<std::vector<double> >& mat_moments,
                        int mat_id,
                        const std::vector<Point<2> >& source_points,
                        const std::vector<Point<2> >& target_points,
                        std::vector<int>& tcell_mat_ids,
                        std::vector<int>& tcell_num_mats,
                        int tc,
                        int idx){
  Portage::NumericTolerances_t num_tols;
  num_tols.use_default();
  // Intersect source candidate matpoly with target cell
  std::vector<double> moments =
    Portage::intersect_polys_r2d(source_points, target_points, num_tols);
  // Accumulate moments (if any) from the intersection
  if (moments[0] > seps) {
    for(int moment=0;moment<NUM_MOMENTS;++moment){
      mat_moments[mat_id][moment]+=moments[moment];
    }
    // The same type of material could be contributed by multiple source cells, or
    // multiple matpolys of a single source cell; only record the first time a
    // particular material id is seen for the target cell
    auto begin=&tcell_mat_ids[0]+idx;
    auto end=begin+tcell_num_mats[tc];
    if (std::find(begin, end, mat_id)==end) {
      ++tcell_num_mats[tc];
      tcell_mat_ids.push_back(mat_id);
    }
  }
}

int main(int argc, char** argv) {
#ifdef PORTAGE_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif
  if (argc != 5) {
    std::ostringstream os;
    os << std::endl << "Correct usage: ir_linetest_app  <nx_source> <ny_source> <nx_target> <ny_target>" << std::endl;
    throw std::runtime_error(os.str());
  }

  int s_nx = atoi(argv[1]);
  int s_ny = atoi(argv[2]);
  std::vector<double> xbnds = {0.0, 1.0};
  std::vector<double> ybnds = {0.0, 1.0};

  Wonton::Simple_Mesh source_mesh(xbnds[0], ybnds[0],
                                   xbnds[1], ybnds[1],
                                   s_nx, s_ny);
  Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(source_mesh);

  int t_nx = atoi(argv[3]);
  int t_ny = atoi(argv[4]);
  Wonton::Simple_Mesh target_mesh(xbnds[0], ybnds[0],
                                   xbnds[1], ybnds[1],
                                   t_nx, t_ny);
  Wonton::Simple_Mesh_Wrapper target_mesh_wrapper(target_mesh);

  Portage::SearchSimple<Wonton::Simple_Mesh_Wrapper,
                        Wonton::Simple_Mesh_Wrapper>
    search(source_mesh_wrapper, target_mesh_wrapper);

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Point<2>> cell_mat_centroids;
  std::vector<int> offsets;
  get_materials_data(source_mesh_wrapper, mat_int_a, mat_int_b,
                     cell_num_mats, cell_mat_ids,
                     cell_mat_volfracs, cell_mat_centroids, offsets);

  std::vector<Tangram::IterativeMethodTolerances_t> tols(2,{100, 1e-12, 1e-12});
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Wonton::Simple_Mesh_Wrapper> source_xmof_driver(source_mesh_wrapper, tols);

  source_xmof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids,
                                   cell_mat_volfracs, cell_mat_centroids);
  source_xmof_driver.reconstruct();

  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    cellmatpoly_list = source_xmof_driver.cell_matpoly_ptrs();

  double source_ir_rel_error = get_ir_error(source_mesh_wrapper, mat_int_a, mat_int_b,
                                          cell_num_mats, cellmatpoly_list, deps);
  std::cout << "Reconstruction error for the source mesh -> " << source_ir_rel_error <<
             std::endl;

  Tangram::write_to_gmv(source_mesh_wrapper, 2, cell_num_mats, cell_mat_ids,
                        cellmatpoly_list, "source_mesh.gmv");

  int tcells = target_mesh_wrapper.num_owned_cells();

  std::vector<int> tcell_num_mats=std::vector<int>(tcells, 0);
  std::vector<int> tcell_mat_ids{};
  std::vector<double> tcell_mat_volfracs{};
  std::vector<Point<2>> tcell_mat_centroids{};

  // Target cells loop
  for (int tc = 0; tc < tcells; ++tc) {
    std::vector<Point<2>> target_points;
    target_mesh_wrapper.cell_get_coordinates(tc, &target_points);
    double cell_size = target_mesh_wrapper.cell_volume(tc);

    // Beginning index for current target cell (into flat vector of material ids for all target cells)
    int idx=tcell_mat_ids.size();

    // Accumulator for target cell material moments
    std::vector<std::vector<double> > mat_moments(NUM_MATS, std::vector<double>(NUM_MOMENTS,0));

    // Search for candidate overlapping source cells
    std::vector<int> sc_candidates;
    search(tc, &sc_candidates);
    for (int sc : sc_candidates){

      // Source candidate cellmatpoly
      auto cellmatpoly=cellmatpoly_list[sc];

      if (cellmatpoly){ // Mixed material cell

        // Access all matpolys from source candidate cellmatpoly
        for (int matpoly_id=0;matpoly_id<cellmatpoly->num_matpolys();++matpoly_id){

          // Get the matpoly and determine which material it contains
          auto matpoly = cellmatpoly->get_ith_matpoly(matpoly_id);
          int mat_id=cellmatpoly->matpoly_matid(matpoly_id);

          // Get Tangram points for matpoly, convert to Portage points
          std::vector<Point<2>> source_points;
          for (auto p : matpoly.points()) {
            source_points.push_back(Point<2>(p));
          }
          add_intersect_data(mat_moments,mat_id,
                             source_points,target_points,
                             tcell_mat_ids,tcell_num_mats,tc,idx);
        }

      } else { // Single material cell

        // Get Tangram points for matpoly, convert to Portage points
        std::vector<Point<2>> source_points;
        source_mesh_wrapper.cell_get_coordinates(sc,&source_points);
        int mat_id=cell_mat_ids[offsets[sc]];
        add_intersect_data(mat_moments,mat_id,
                           source_points,target_points,
                           tcell_mat_ids,tcell_num_mats,tc,idx);
      }
    } // Finish gathering material moments for target cell

    // Populate centroid and volfrac data for target cell
    auto begin=&tcell_mat_ids[0]+idx;
    auto end=begin+tcell_num_mats[tc];
    for(auto it=begin; it<end; ++it){
      int mat_id=*it;
      tcell_mat_centroids.push_back({
          mat_moments[mat_id][1]/mat_moments[mat_id][0],
            mat_moments[mat_id][2]/mat_moments[mat_id][0]});
      tcell_mat_volfracs.push_back(mat_moments[mat_id][0]/cell_size);
    }
  }

  // Perform final interface reconstruction on target mesh, then write to file
  // for comparison with source mesh
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Wonton::Simple_Mesh_Wrapper> target_xmof_driver(target_mesh_wrapper, tols);
  target_xmof_driver.set_volume_fractions(tcell_num_mats, tcell_mat_ids,
                                   tcell_mat_volfracs, tcell_mat_centroids);
  target_xmof_driver.reconstruct();  // executor arg defaults to null -> serial
  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    tcellmatpoly_list = target_xmof_driver.cell_matpoly_ptrs();
  double target_ir_rel_error = get_ir_error(target_mesh_wrapper, mat_int_a, mat_int_b,
  tcell_num_mats, tcellmatpoly_list, deps);
  std::cout << "Reconstruction error for the target mesh -> " << target_ir_rel_error <<
           std::endl;
  Tangram::write_to_gmv(target_mesh_wrapper, 2, tcell_num_mats, tcell_mat_ids,
  tcellmatpoly_list, "target_mesh.gmv");

#ifdef PORTAGE_ENABLE_MPI
  MPI_Finalize();
#endif
}

/*!
 @brief Convert the line given by y = line_a*x + line_b to the R2D format.
 @param[in] line_a Slope of the line.
 @param[in] line_b y-intercept of the line.
 @param[out] r2d_line Corresponding line in R2D format
*/
void R2DizeLine(const double line_a,
                const double line_b,
                r2d_plane& r2d_line) {
  double norm = sqrt(1.0 + line_a*line_a);
  r2d_line.n.x = -line_a/norm;
  r2d_line.n.y =  1.0/norm;
  r2d_line.d   = -line_b*r2d_line.n.y;
}

/*!
 @brief Convert a mesh cell into an R2D polygon.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
 implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] cell_id Index of the cell to convert.
 @param[out] r2d_polygon Corresponding polygon in R2D format
*/
template <class Mesh_Wrapper>
void R2DizeCell(const Mesh_Wrapper& Mesh,
                const int cell_id,
                r2d_poly& r2d_polygon) {
  std::vector<Point<2>> cell_points;
  Mesh.cell_get_coordinates(cell_id, &cell_points);
  int npoints = (int) cell_points.size();
  r2d_rvec2* vertices = new r2d_rvec2[npoints];
  for (int ipt = 0; ipt < npoints; ipt++) {
    vertices[ipt].xy[0] = cell_points[ipt][0];
    vertices[ipt].xy[1] = cell_points[ipt][1];
  }

  r2d_init_poly(&r2d_polygon, vertices, npoints);
  delete vertices;
}

/*!
 @brief Determines the position of a mesh cell with respect to a line.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
 implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] cell_id Index of the cell to test for position.
 @param[in] line_a Slope of the reference material interface.
 @param[in] line_b y-intercept of the the reference material interface.
 @param[in] yeps Max difference in y-coordinate from the value given by the line equation
 for points on the line.
 @return -1 if the cell is below the line, 1 if the cell is above the line,
 0 if the line passes through the interior of the cell
*/
template <class Mesh_Wrapper>
int CellPosition(const Mesh_Wrapper& Mesh,
                 const int cell_id,
                 const double line_a,
                 const double line_b,
                 const double yeps) {
  std::vector<Point<2>> cell_points;
  Mesh.cell_get_coordinates(cell_id, &cell_points);
  int npoints = (int) cell_points.size();
  int counter_off = 0, counter_on = 0;
  for (int ipt = 0; ipt < npoints; ipt++) {
    double line_y = line_a*cell_points[ipt][0] + line_b;
    if (cell_points[ipt][1] < line_y - yeps)
      counter_off += -1;
    else if (cell_points[ipt][1] > line_y + yeps)
      counter_off++;
    else counter_on++;
  }
  int result = 0;
  if (counter_on - counter_off == npoints)
    result = -1;
  else if (counter_on + counter_off == npoints)
    result = 1;

  return result;
}

/*!
 @brief Computes material data for a mesh and a given linear material interface using R2D.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
 implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] line_a Slope of the line.
 @param[in] line_b y-intercept of the line.
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
 computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
 vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
 requires computations of offsets
 @param[out] offsets Offset into cell_mat_ids giving the starting section for each cell
*/
template <class Mesh_Wrapper>
void get_materials_data(const Mesh_Wrapper& Mesh,
                        const double line_a,
                        const double line_b,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Point<2>>& cell_mat_centroids,
                        std::vector<int>& offsets) {
  const int POLY_ORDER = 1;  //Max degree of moments to calculate

  cell_num_mats.clear();
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();

  r2d_plane r2d_line;
  R2DizeLine(line_a, line_b, r2d_line);
  int ncells = Mesh.num_owned_cells();
  cell_num_mats.resize(ncells, 1);
  for (int icell = 0; icell < ncells; icell++) {
    double cell_size = Mesh.cell_volume(icell);
    Point<2> cell_centroid;
    Mesh.cell_centroid(icell, &cell_centroid);
    Point<2> cell_cen(cell_centroid);

    int cell_pos = CellPosition(Mesh, icell, line_a, line_b);
    if (cell_pos == 0) {
      r2d_poly cell_r2d_poly;
      R2DizeCell(Mesh, icell, cell_r2d_poly);
      r2d_clip(&cell_r2d_poly, &r2d_line, 1);
      r2d_real moments[R2D_NUM_MOMENTS(POLY_ORDER)];
      r2d_reduce(&cell_r2d_poly, moments, POLY_ORDER);
      if (moments[0] < -seps)
        throw std::runtime_error("Negative area of the clipped polygon");
      if ((moments[0] > seps) && (moments[0] < cell_size - seps)) {
        cell_num_mats[icell] = 2;

        cell_mat_ids.push_back(mat_id_above);
        double poly_above_vfrac = moments[0]/cell_size;
        cell_mat_volfracs.push_back(poly_above_vfrac);
        Point<2> poly_above_cen(
          moments[1]/moments[0], moments[2]/moments[0]);
        cell_mat_centroids.push_back(poly_above_cen);

        cell_mat_ids.push_back(mat_id_below);
        cell_mat_volfracs.push_back(1.0 - poly_above_vfrac);
        Point<2> poly_below_cen(
          (cell_cen - poly_above_vfrac*poly_above_cen)/(1.0 - poly_above_vfrac));
        cell_mat_centroids.push_back(poly_below_cen);
      }
      else {
        cell_mat_volfracs.push_back(1.0);
        cell_mat_centroids.push_back(cell_cen);
        if (moments[0] <= seps)
          cell_mat_ids.push_back(mat_id_below);
        else
          cell_mat_ids.push_back(mat_id_above);
      }
    }
    else {
      cell_mat_volfracs.push_back(1.0);
      cell_mat_centroids.push_back(cell_cen);
      if (cell_pos < 0)
        cell_mat_ids.push_back(mat_id_below);
      else
        cell_mat_ids.push_back(mat_id_above);
    }
    offsets.push_back(cell_mat_ids.size()-1);
  }
}

/*!
 @brief Determines the position of a point with respect to a line.
 @param[in] pt Point coordinates.
 @param[in] line_nvec Normal vector to the line.
 @param[in] line_pt Point on the line.
 @param[in] eps Max distance between the line and points considered to be on the line.
 @return -1 if the point is below the line, 1 if the point is above the line,
 0 if the point is on the line
*/

int PointPosition(const Point<2>& pt,
                  const Wonton::Vector2& line_nvec,
                  const Point<2>& line_pt,
                  const double eps) {
  Wonton::Vector<2> vec2line = line_pt - pt;
  double prj = dot(vec2line, line_nvec);
  int pos;
  if (std::fabs(prj) < eps) pos = 0;
  else pos = std::signbit(prj) ? -1 : 1;

  return pos;
}

/*!
 @brief Computes the intersection of two non-parallel lines.
 @param[in] lines_pts Two pairs of points defining the intersecting lines.
 @param[in] den_eps Minimal value of the denominator for which the analytical formula is
 to be used, if the denominator is smaller, the bisection algorithm
 is used instead.
 @param[in] eps Max distance between points that are considered coincident.
 @return Intersection point
*/
Point<2> LinesIntersect(std::vector< std::vector<Point<2>> > lines_pts,
                               double denom_eps,
                               double eps) {
  Point<2> p_int;

  double denom = (lines_pts[0][0][0] - lines_pts[0][1][0])*(lines_pts[1][0][1] - lines_pts[1][1][1]) -
  (lines_pts[0][0][1] - lines_pts[0][1][1])*(lines_pts[1][0][0] - lines_pts[1][1][0]);

  if (std::fabs(denom) > denom_eps) {
    p_int[0] = (lines_pts[0][0][0]*lines_pts[0][1][1] - lines_pts[0][0][1]*lines_pts[0][1][0])*
    (lines_pts[1][0][0] - lines_pts[1][1][0]) -
    (lines_pts[1][0][0]*lines_pts[1][1][1] - lines_pts[1][0][1]*lines_pts[1][1][0])*
    (lines_pts[0][0][0] - lines_pts[0][1][0]);
    p_int[1] = (lines_pts[0][0][0]*lines_pts[0][1][1] - lines_pts[0][0][1]*lines_pts[0][1][0])*
    (lines_pts[1][0][1] - lines_pts[1][1][1]) -
    (lines_pts[1][0][0]*lines_pts[1][1][1] - lines_pts[1][0][1]*lines_pts[1][1][0])*
    (lines_pts[0][0][1] - lines_pts[0][1][1]);

    p_int /= denom;
  }
  else {
    Wonton::Vector<2> line1_nvec(lines_pts[1][1][1] - lines_pts[1][0][1],
                                lines_pts[1][0][0] - lines_pts[1][1][0]);
    line1_nvec.normalize();
    Point<2> end_pts[2];
    end_pts[0] = lines_pts[0][0]; end_pts[1] = lines_pts[0][1];
    int end_pts_pos[2];
    end_pts_pos[0] = PointPosition(end_pts[0], line1_nvec, lines_pts[1][0], eps);
    end_pts_pos[1] = PointPosition(end_pts[1], line1_nvec, lines_pts[1][0], eps);
    if (end_pts_pos[0] == end_pts_pos[1]) {
      double d2line[2];
      for (int ipt = 0; ipt < 2; ipt++) {
        Wonton::Vector<2> vec2line = lines_pts[1][0] - end_pts[ipt];
        d2line[ipt] = std::fabs(dot(vec2line, line1_nvec));
      }
      int inearest = (d2line[0] < d2line[1]) ? 0 : 1;
      while (end_pts_pos[0] == end_pts_pos[1]) {
        Wonton::Vector<2> dvec = end_pts[inearest] - end_pts[(inearest + 1)%2];
        end_pts[inearest] += 2*dvec;
        end_pts_pos[inearest] = PointPosition(end_pts[inearest], line1_nvec, lines_pts[1][0], eps);
      }
    }
    int mid_pt_pos;
    do {
      Point<2> mid_pt = 0.5*(end_pts[0] + end_pts[1]);
      mid_pt_pos = PointPosition(mid_pt, line1_nvec, lines_pts[1][0], eps);
      if (mid_pt_pos == 0)
        p_int = mid_pt;
      else if (mid_pt_pos == end_pts_pos[0])
        end_pts[0] = mid_pt;
      else
        end_pts[1] = mid_pt;
    } while (mid_pt_pos != 0);
  }

  return p_int;
}

/*!
 @brief Computes Hausdorff distance between two linear segments.
 @param[in] segments Two pairs of points defining the linear segments.
 @return Hausdorff distance
*/
double SegmentsDistance(const std::vector<std::vector<Point<2>>>& segments) {
  double max_dist = 0.0;
  for (int iseg = 0; iseg < 2; iseg++) {
    int iother = (iseg + 1)%2;
    Wonton::Vector<2> other_dvec = segments[iother][1] - segments[iother][0];
    double other_len = other_dvec.norm();
    other_dvec /= other_len;
    for (int ipt = 0; ipt < 2; ipt++) {
      double shift = dot(segments[iseg][ipt] - segments[iother][0], other_dvec);
      double cur_dist;
      if (shift <= 0.0)
        cur_dist = (segments[iseg][ipt] - segments[iother][0]).norm();
      else if (shift >= other_len)
        cur_dist = (segments[iseg][ipt] - segments[iother][1]).norm();
      else {
        Point<2> prj_pt = segments[iother][0] + shift*other_dvec;
        cur_dist = (segments[iseg][ipt] - prj_pt).norm();
      }
      if (cur_dist > max_dist) max_dist = cur_dist;
    }
  }
  return max_dist;
}

/*!
 @brief Computes the max Hausdorff distance between the reference and reconstructed
 interfaces over the multi-material cells in a mesh.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
 implementation that provides required functionality.
 @param[in] Mesh Mesh wrapper.
 @param[in] line_a Slope of the reference material interface.
 @param[in] line_b y-intercept of the the reference material interface.
 @param[in] cell_num_mats Reference number of materials for every mesh cell.
 @param[in] cellmatpoly_list Vector of pointers to CellMatPoly objects corresponding to
 the results of interface reconstruction
 @param[in] eps Max distance between points that are considered coincident.
 @return Max Hausdorff distance
*/
template <class Mesh_Wrapper>
double get_ir_error(const Mesh_Wrapper& Mesh,
                    const double line_a,
                    const double line_b,
                    const std::vector<int>& cell_num_mats,
                    const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
                    cellmatpoly_list,
                    const double eps) {
  int ncells = Mesh.num_owned_cells();
  double max_hdist = 0.0;
  int offset = 0;
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      if (cellmatpoly_list[icell] != nullptr)
        throw std::runtime_error("CellMatPoly was constructed for a single-material cell!");
      continue;
    }
    if (cellmatpoly_list[icell] == nullptr)
      throw std::runtime_error("CellMatPoly was NOT constructed for a multi-material cell!");

    std::vector<int> ifaces, fdirs;
    Mesh.cell_get_faces_and_dirs(icell, &ifaces, &fdirs);
    int nsides = (int) ifaces.size();
    std::vector<Point<2>> ref_int_pts;
    for (int iside = 0; iside < nsides; iside++) {
      std::vector<int> inodes;
      Mesh.face_get_nodes(ifaces[iside], &inodes);
      std::vector<Point<2>> side_pts(2);
      Mesh.node_get_coordinates(inodes[0], &side_pts[0]);
      Mesh.node_get_coordinates(inodes[1], &side_pts[1]);
      int pts_pos[2];
      for (int ipt = 0; ipt < 2; ipt++) {
        double line_y = line_a*side_pts[ipt][0] + line_b;
        if (side_pts[ipt][1] < line_y - eps)
          pts_pos[ipt] = -1;
        else if (side_pts[ipt][1] > line_y + eps)
          pts_pos[ipt] = 1;
        else
          pts_pos[ipt] = 0;
      }

      if (!pts_pos[0] || !pts_pos[1]) {
        Point<2> ref_int = pts_pos[0] ? side_pts[1] : side_pts[0];
        bool new_int_pt = true;
        for (int iip = 0; iip < ref_int_pts.size(); iip++)
          if (approxEq(ref_int, ref_int_pts[iip], eps)) {
            new_int_pt = false;
            break;
          }
        if (new_int_pt)
          ref_int_pts.push_back(ref_int);
      }
      else if (pts_pos[0] != pts_pos[1]) {
        std::vector<Point<2>> line_pts = {
          Point<2>(0.0, line_b), Point<2>(1.0, line_a + line_b) };
        Point<2> ref_int = LinesIntersect({side_pts, line_pts}, denom_eps, eps);
        ref_int_pts.push_back(ref_int);
      }
    }
    if (ref_int_pts.size() != 2)
      throw std::runtime_error("Expected two intersection points!");

    const Tangram::CellMatPoly<2>& cell_mat_poly = *cellmatpoly_list[icell];
    int nfaces = cell_mat_poly.num_matfaces();
    std::vector<int> iintfaces;
    for (int iface = 0; iface < nfaces; iface++)
      if (cell_mat_poly.matface_is_interface(iface))
        iintfaces.push_back(iface);
    if (iintfaces.size() != 1)
      throw std::runtime_error("All CellMatPoly's are expected to have only one interior face!");
    else {
      std::vector<Point<2>> int_face_pts =
        cell_mat_poly.matface_points(iintfaces[0]);
      std::vector<Point<2>> ir_int_pts = { Point<2>(int_face_pts[0]),
                                                  Point<2>(int_face_pts[1]) };
      double cur_mat_int_dist = SegmentsDistance({ ir_int_pts, ref_int_pts });
      if (cur_mat_int_dist > max_hdist)  max_hdist = cur_mat_int_dist;
    }
  }
  return (max_hdist > eps) ? max_hdist : 0.0;
}
