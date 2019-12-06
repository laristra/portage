/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <vector>

// portage includes
#include "portage/driver/coredriver.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/interpolate/gradient.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/search/search_kdtree.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/support/Point.h"

// Jali includes
#include "JaliState.h"
#include "JaliStateVector.h"
#include "MeshFactory.hh"

/*!
  @file momentumapp.cc
  @brief An application that drives remap of nodal velocities

  The program is used to showcase the capability of remapping nodal
  velocities such that the total momentum is conserved. The app
  consists of two parts:
    A. Populate input data: corner masses and nodal velocities
    B. Conservative remap of momentum
    C. Verifivation of output data: coner mass and nodal
       velocities on the target mesh
 */


void print_usage() {
  std::printf("usage: momentumapp ncellsx ncellsy [order=1]\n");
}

template<int D>
double dot_product(const Wonton::Point<D>& a, const Wonton::Point<D>& b) {
  double r = 0.0;
  for (int i = 0; i < D; i++) r += a[i] * b[i];
  return r;
}


int main(int argc, char** argv) {
  if (argc <= 2) {
    print_usage();
    return 0;
  }

  // Initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
      MPI_Init(&argc, &argv);
  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // number of CELLS in x and y for input mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  auto order = (argc > 3) ? atoi(argv[3]) : 1;

  // size of computational domain
  double lenx = 1.0;  // [m]
  double leny = 1.0;  // [m]
  
  // physical quantities
  double density = 1.0;  // [kg / m^2]


  // Preliminaries

  // -- setup Jali meshes
  auto srcmesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, lenx, leny, nx, ny);
  auto trgmesh = Jali::MeshFactory(MPI_COMM_WORLD)(0.0, 0.0, lenx, leny, nx + 1, ny + 1);

  // -- setup mesh wrappers
  Wonton::Jali_Mesh_Wrapper srcmesh_wrapper(*srcmesh);
  Wonton::Jali_Mesh_Wrapper trgmesh_wrapper(*trgmesh);

  int ncells_src = srcmesh_wrapper.num_owned_cells();
  int nnodes_src = srcmesh_wrapper.num_owned_nodes();
  int ncorners_src = srcmesh_wrapper.num_owned_corners();

  int ncells_trg = trgmesh_wrapper.num_owned_cells();
  int ncorners_trg = trgmesh_wrapper.num_owned_corners();

  // -- states
  std::shared_ptr<Jali::State> srcstate = Jali::State::create(srcmesh);
  std::shared_ptr<Jali::State> trgstate = Jali::State::create(trgmesh);

  // -- state wrappers
  Wonton::Jali_State_Wrapper srcstate_wrapper(*srcstate);
  Wonton::Jali_State_Wrapper trgstate_wrapper(*trgstate);

  // -- define two constant velocity components in states
  std::vector<double> ux_src(nnodes_src, 1.0);
  std::vector<double> uy_src(nnodes_src, 1.0);

  srcstate->add("velocity_x", srcmesh, Jali::Entity_kind::NODE,
                Jali::Entity_type::ALL, &(ux_src[0]));

  srcstate->add("velocity_y", srcmesh, Jali::Entity_kind::NODE,
                Jali::Entity_type::ALL, &(uy_src[0]));

  trgstate->add<double, Jali::Mesh, Jali::UniStateVector>(
                "velocity_x", trgmesh, Jali::Entity_kind::NODE,
                Jali::Entity_type::ALL);

  trgstate->add<double, Jali::Mesh, Jali::UniStateVector>(
                "velocity_y", trgmesh, Jali::Entity_kind::NODE,
                Jali::Entity_type::ALL);

  // -- define fixed corner masses from constant density
  std::vector<double> mass_cn_src(ncorners_src);

  double total_mass(0.0);
  for (int cn = 0; cn < ncorners_src; ++cn) {
    mass_cn_src[cn] = density * srcmesh_wrapper.corner_volume(cn);
    total_mass += mass_cn_src[cn];
  }

  srcstate->add("mass_corner", srcmesh, Jali::Entity_kind::CORNER,
                Jali::Entity_type::ALL, &(mass_cn_src[0]));

  trgstate->add<double, Jali::Mesh, Jali::UniStateVector>(
                "mass_corner", trgmesh, Jali::Entity_kind::CORNER,
                Jali::Entity_type::ALL);

  // -- summary
  std::cout << "mesh:           " << nx << " x " << ny << std::endl;
  std::cout << "density:        " << density << " kg/m^2" << std::endl;
  std::cout << "total mass:     " << total_mass << " kg" << std::endl;

  //
  // Main six-step remap algorithm
  //

  // -- Step 1: collect cell-centered mass from corner masses
  //    names of auxiliary qunatities carry no suffix _src
  std::vector<double> mass_c(ncells_src, 0.0);
  std::vector<int> corners;
  
  for (int c = 0; c < ncells_src; ++c) {
    srcmesh_wrapper.cell_get_corners(c, &corners);
    for (auto cn : corners) mass_c[c] += mass_cn_src[cn];
  }

  // -- compute density in each mesh cell
  std::vector<double> density_c(ncells_src);
  for (int c = 0; c < ncells_src; ++c) {
    density_c[c] = mass_c[c] / srcmesh_wrapper.cell_volume(c);
  }
  
  // -- Step 2: collect cell-centered momentum from corner contributions
  //    and compute specific momentum
  std::vector<double> momentum_x(ncells_src, 0.0);
  std::vector<double> momentum_y(ncells_src, 0.0);

  Wonton::Point<2> total_momentum = {0.0, 0.0};
  for (int c = 0; c < ncells_src; ++c) {
    srcmesh_wrapper.cell_get_corners(c, &corners);

    double& mx = momentum_x[c];
    double& my = momentum_y[c];
    for (auto cn : corners) { 
      int v = srcmesh_wrapper.corner_get_node(cn);
      mx += mass_cn_src[cn] * ux_src[v];
      my += mass_cn_src[cn] * uy_src[v];
    }
    total_momentum[0] += mx;
    total_momentum[1] += my;

    double volume = srcmesh_wrapper.cell_volume(c);
    mx /= volume;
    my /= volume;
  }

  std::cout << "total momentum: " << total_momentum << " kg m/s" << std::endl;

  // -- Step 3: remap density and specific momentum followin three basic
  //    step: search, intersect, interpolate.
  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      cd(srcmesh_wrapper, srcstate_wrapper,
         trgmesh_wrapper, trgstate_wrapper);

  Portage::NumericTolerances_t default_num_tols;
  default_num_tols.use_default();
  cd.set_num_tols(default_num_tols);

  auto candidates = cd.search<Portage::SearchKDTree>();
  auto srcwts = cd.intersect_meshes<Portage::IntersectR2D>(candidates);
  bool has_mismatch = cd.check_mesh_mismatch(srcwts);

  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  // -- we need register the fields we want to remap
  std::vector<std::string> field_names;
  std::vector<const double*> field_pointers;

  field_names.push_back("density");
  field_pointers.push_back(&(density_c[0]));

  field_names.push_back("momentum_x");
  field_pointers.push_back(&(momentum_x[0]));

  field_names.push_back("momentum_y");
  field_pointers.push_back(&(momentum_y[0]));

  for (int i = 0; i < field_names.size(); ++i) {
    srcstate->add(field_names[i], srcmesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, field_pointers[i]);

    trgstate->add<double, Jali::Mesh, Jali::UniStateVector>(
                  field_names[i], trgmesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL);

    cd.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
        "density", "density", srcwts, dblmin, dblmax);
  }

  // -- Step 4: create linear reconstruction (limited or unlimited) 
  //    of density and specific momentun on the target mesh.
  std::vector<Portage::vector<Wonton::Vector<2>>> gradients;

  for (int i = 0; i < field_names.size(); ++i) {
    Portage::Limited_Gradient<2, Wonton::Entity_kind::CELL, 
                              Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
        gradient_kernel(trgmesh_wrapper, trgstate_wrapper,
                        field_names[i], Portage::NOLIMITER, Portage::BND_NOLIMITER);

    Portage::vector<Wonton::Vector<2>> gradient(ncells_trg);

    Portage::transform(trgmesh_wrapper.begin(Wonton::Entity_kind::CELL),
                       trgmesh_wrapper.end(Wonton::Entity_kind::CELL),
                       gradient.begin(), gradient_kernel);

    gradients.push_back(gradient);  // could be optimized
  }

  // -- Step 5: integrate density and specific momentum at corners
  //    extract primary data
  Jali::UniStateVector<double, Jali::Mesh> mass_cn_trg;
  Jali::UniStateVector<double, Jali::Mesh> ux_trg, uy_trg;

  trgstate->get("mass_corner", trgmesh, Jali::Entity_kind::CORNER,
                Jali::Entity_type::ALL, &mass_cn_trg);

  trgstate->get("velocity_x", trgmesh, Jali::Entity_kind::NODE,
                Jali::Entity_type::ALL, &ux_trg);

  trgstate->get("velocity_y", trgmesh, Jali::Entity_kind::NODE,
                Jali::Entity_type::ALL, &uy_trg);

  //    extract auxiliary (cell-centerd) data
  Jali::UniStateVector<double, Jali::Mesh> density_trg;
  Jali::UniStateVector<double, Jali::Mesh> momentum_x_trg, momentum_y_trg;

  trgstate->get("density", trgmesh, Jali::Entity_kind::CELL,
                Jali::Entity_type::ALL, &density_trg);

  trgstate->get("momentum_x", trgmesh, Jali::Entity_kind::CELL,
                Jali::Entity_type::ALL, &momentum_x_trg);

  trgstate->get("momentum_y", trgmesh, Jali::Entity_kind::CELL,
                Jali::Entity_type::ALL, &momentum_y_trg);

  // -- integrate
  std::vector<double> momentum_cn_x(ncorners_trg);
  std::vector<double> momentum_cn_y(ncorners_trg);

  Wonton::Point<2> xc;
  std::vector<int> wedges;
  std::array<Wonton::Point<2>, 3> wcoords;

  for (int c = 0; c < ncells_trg; ++c) {
    trgmesh_wrapper.cell_centroid(c, &xc);
    trgmesh_wrapper.cell_get_corners(c, &corners);

    for (auto cn : corners) {
      double volume = trgmesh_wrapper.corner_volume(cn); 
      trgmesh_wrapper.corner_get_wedges(c, &wedges);

      Wonton::Point<2> xcn = {0.0, 0.0};  // corner centroid
      for (auto w : wedges) {
        double frac = trgmesh_wrapper.wedge_volume(w) / (3 * volume); 
        trgmesh_wrapper.wedge_get_coordinates(w, &wcoords);
        for (int i = 0; i < 3; ++i) {
          xcn += frac * wcoords[i];
        }
      }

      // integral is the value at centroid
      mass_cn_trg[cn] = density_trg[c] + dot_product(gradients[0][c], xcn - xc);
      // momentum_cn_x[cn] = momentum_x_trg[c] + dot(gradients[1][c], xcn - xc);
      // momentum_cn_y[cn] = momentum_y_trg[c] + dot(gradients[2][c], xcn - xc);
    }
  }

  // -- Step 6: gather data to mesh nodes and compute velocity

  //
  // Verification
  //

  MPI_Finalize();
}
