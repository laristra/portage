/*
  This file is part of the Ristra portage project.
  Please see the license file at the root of this repository, or at:
  https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <math.h>
#include "portage/wonton/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "portage/simple_mesh/simple_mesh.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/driver/write_to_gmv.h"

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

extern "C" {
#include "portage/intersect/r2d.h"
}

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

const double deps = 1.0e-15;	 //Distance tolerance
const double seps = 1.0e-14;	 //Size(area) tolerance
const double denom_eps = 1.0e-6; //Denominator tolerance

void R2DizeLine(const double line_a,
                const double line_b,
                r2d_plane& r2d_line);

template <class Mesh_Wrapper>
void R2DizeCell(const Mesh_Wrapper& Mesh,
                const int cell_id,
                r2d_poly& r2d_polygon);

template <class Mesh_Wrapper>
int CellPosition(const Mesh_Wrapper& Mesh,
                 const int cell_id,
                 const double line_a,
                 const double line_b,
                 const double yeps = deps);

template <class Mesh_Wrapper>
void get_materials_data(const Mesh_Wrapper& Mesh,
                        const double line_a,
                        const double line_b,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Tangram::Point2>& cell_mat_centroids);

int PointPosition(const Portage::Point2& pt,
                  const Portage::Vector2& line_nvec,
                  const Portage::Point2& line_pt,
                  const double eps = deps);

Portage::Point2 LinesIntersect(std::vector< std::vector<Portage::Point2> > lines_pts,
                               double den_eps,
                               double eps = deps);

double SegmentsDistance(const std::vector<std::vector<Portage::Point2>>& segments);

template <class Mesh_Wrapper>
double get_ir_error(const Mesh_Wrapper& Mesh,
                    const double line_a,
                    const double line_b,
                    const std::vector<int>& cell_num_mats,
                    const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
                      cellmatpoly_list,
                    const double eps = deps);

int main(int argc, char** argv) {
#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif
  if (argc != 3) {
    std::ostringstream os;
    os << std::endl << "Correct usage: ir_linetest_app <nx> <ny>" << std::endl;
    throw std::runtime_error(os.str());
  }

  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  std::vector<double> xbnds = {0.0, 1.0};
  std::vector<double> ybnds = {0.0, 1.0};
 
  Portage::Simple_Mesh source_mesh(xbnds[0], ybnds[0],
                                   xbnds[1], ybnds[1],
                                   nx, ny);
  Wonton::Simple_Mesh_Wrapper source_mesh_wrapper(source_mesh);

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point2> cell_mat_centroids;
  get_materials_data(source_mesh_wrapper, mat_int_a, mat_int_b,
                     cell_num_mats, cell_mat_ids, 
                     cell_mat_volfracs, cell_mat_centroids);
  
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Wonton::Simple_Mesh_Wrapper> xmof_driver(source_mesh_wrapper);
  
  xmof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                   cell_mat_volfracs, cell_mat_centroids);
  xmof_driver.reconstruct();
  
  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    cellmatpoly_list = xmof_driver.cell_matpoly_ptrs();

  double source_ir_rel_error = get_ir_error(source_mesh_wrapper, mat_int_a, mat_int_b,
                                            cell_num_mats, cellmatpoly_list, deps);
  std::cout << "Reconstruction error for the source mesh -> " << source_ir_rel_error <<
               std::endl;
  
  Tangram::write_to_gmv(source_mesh_wrapper, 2, cell_num_mats, cell_mat_ids,
                        cellmatpoly_list, "source_mesh.gmv");
#ifdef ENABLE_MPI  
  MPI_Finalize();
#endif
}

void R2DizeLine(const double line_a,
                const double line_b,
                r2d_plane& r2d_line) {
  double norm = sqrt(1.0 + line_a*line_a);
  r2d_line.n.x = -line_a/norm;
  r2d_line.n.y =  1.0/norm;
  r2d_line.d   = -line_b*r2d_line.n.y;
}

template <class Mesh_Wrapper>
void R2DizeCell(const Mesh_Wrapper& Mesh,
                const int cell_id,
                r2d_poly& r2d_polygon) {
  std::vector<Portage::Point2> cell_points;
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

template <class Mesh_Wrapper>
int CellPosition(const Mesh_Wrapper& Mesh,
                 const int cell_id,
                 const double line_a,
                 const double line_b,
                 const double yeps) {
  std::vector<Portage::Point2> cell_points;
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

template <class Mesh_Wrapper>
void get_materials_data(const Mesh_Wrapper& Mesh,
                        const double line_a,
                        const double line_b,
                        std::vector<int>& cell_num_mats,
                        std::vector<int>& cell_mat_ids,
                        std::vector<double>& cell_mat_volfracs,
                        std::vector<Tangram::Point2>& cell_mat_centroids) {
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
    Portage::Point2 cell_centroid;
    Mesh.cell_centroid(icell, &cell_centroid);
    Tangram::Point2 cell_cen(cell_centroid);
    
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
        Tangram::Point2 poly_above_cen(
          moments[1]/moments[0], moments[2]/moments[0]);
        cell_mat_centroids.push_back(poly_above_cen);
        
        cell_mat_ids.push_back(mat_id_below);
        cell_mat_volfracs.push_back(1.0 - poly_above_vfrac);
        Tangram::Point2 poly_below_cen(
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
  }
}

int PointPosition(const Portage::Point2& pt,
                  const Portage::Vector2& line_nvec,
                  const Portage::Point2& line_pt,
                  const double eps) {
  Portage::Vector2 vec2line = line_pt - pt;
  double prj = dot(vec2line, line_nvec);
  int pos;
  if (std::fabs(prj) < eps) pos = 0;
  else pos = signbit(prj) ? -1 : 1;
  
  return pos;
}

Portage::Point2 LinesIntersect(std::vector< std::vector<Portage::Point2> > lines_pts,
                               double denom_eps,
                               double eps) {
  Portage::Point2 p_int;

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
    Portage::Vector2 line1_nvec(lines_pts[1][1][1] - lines_pts[1][0][1],
                                lines_pts[1][0][0] - lines_pts[1][1][0]);
    line1_nvec.normalize();
    Portage::Point2 end_pts[2];
    end_pts[0] = lines_pts[0][0]; end_pts[1] = lines_pts[0][1];
    int end_pts_pos[2];
    end_pts_pos[0] = PointPosition(end_pts[0], line1_nvec, lines_pts[1][0], eps);
    end_pts_pos[1] = PointPosition(end_pts[1], line1_nvec, lines_pts[1][0], eps);
    if (end_pts_pos[0] == end_pts_pos[1]) {
      double d2line[2];
      for (int ipt = 0; ipt < 2; ipt++) {
        Portage::Vector2 vec2line = lines_pts[1][0] - end_pts[ipt];
        d2line[ipt] = std::fabs(dot(vec2line, line1_nvec));
      }
      int inearest = (d2line[0] < d2line[1]) ? 0 : 1;
      while (end_pts_pos[0] == end_pts_pos[1]) {
        Portage::Vector2 dvec = end_pts[inearest] - end_pts[(inearest + 1)%2];
        end_pts[inearest] += 2*dvec;
        end_pts_pos[inearest] = PointPosition(end_pts[inearest], line1_nvec, lines_pts[1][0], eps);
      }
    }
    int mid_pt_pos;
    do {
      Portage::Point2 mid_pt = 0.5*(end_pts[0] + end_pts[1]);
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

double SegmentsDistance(const std::vector<std::vector<Portage::Point2>>& segments) {
  double max_dist = 0.0;
  for (int iseg = 0; iseg < 2; iseg++) {
    int iother = (iseg + 1)%2;
    Portage::Vector2 other_dvec = segments[iother][1] - segments[iother][0];
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
        Portage::Point2 prj_pt = segments[iother][0] + shift*other_dvec;
        cur_dist = (segments[iseg][ipt] - prj_pt).norm();
      }
      if (cur_dist > max_dist) max_dist = cur_dist;
    }
  }
  return max_dist;
}

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
    std::vector<Portage::Point2> ref_int_pts;
    for (int iside = 0; iside < nsides; iside++) {
      std::vector<int> inodes;
      Mesh.face_get_nodes(ifaces[iside], &inodes);
      std::vector<Portage::Point2> side_pts(2);
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
        Portage::Point2 ref_int = pts_pos[0] ? side_pts[1] : side_pts[0];
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
        std::vector<Portage::Point2> line_pts = {
          Portage::Point2(0.0, line_b), Portage::Point2(1.0, line_a + line_b) };
        Portage::Point2 ref_int = LinesIntersect({side_pts, line_pts}, denom_eps, eps);
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
      std::vector<Tangram::Point2> int_face_pts =
        cell_mat_poly.matface_points(iintfaces[0]);
      std::vector<Portage::Point2> ir_int_pts = { Portage::Point2(int_face_pts[0]), 
                                                  Portage::Point2(int_face_pts[1]) };
      double cur_mat_int_dist = SegmentsDistance({ ir_int_pts, ref_int_pts });
      if (cur_mat_int_dist > max_hdist)  max_hdist = cur_mat_int_dist;
    }
  }
  return (max_hdist > eps) ? max_hdist : 0.0;
}

