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
#include <iomanip>
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

const double deps = 1.0e-15;	//Distance tolerance
const double seps = 1.0e-14;	//Size(area) tolerance

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

int main(int argc, char** argv) {
#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
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
