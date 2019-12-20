/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/


#include <fstream>
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
  @brief An application that remaps velocity and mass for SGH and CCH

  The program is used to showcase the capability of remapping velocity 
  and mass for staggered and cell-centered hydro codes. Velocity remap
  conserves the total momentum. The app is controlled by a few input 
  commands. Unnessesary longer code is used for the implementation clarity.

  SGH:
    A. Populate input data: corner masses and nodal velocities
    B. Conservative remap of momentum
    C. Verifivation of output data: coner masses and nodal 
       velocities on the target mesh

  CCH: 
    A. Populate input data: cell-centered masses and velocities
    B. Conservative remap of momentum
    C. Verifivation of output data: cell-centered massed and velocities

  Assumptions: meshes are matching.
 */


const int SGH = 1;
const int CCH = 2;


void print_usage() {
  std::cout << "Usage: ./momentumapp nx ny method limiter ini_density ini_velocity\n\n"
            << "   source mesh:  nx x ny rectangular cells and covers the unit square\n" 
            << "   target mesh:  (nx + 2) x (ny + 4) rectangular cells\n\n"
            << "   method:       SGH=1, CCH=2\n"
            << "   limiter:      0 - limiter is off, otherwise Barth-Jespersen is used\n\n"
            << "   ini_density:  0 - constant density (1)\n"
            << "                 1 - linear density (1 + x * 2y)\n"
            << "                 2 - quadratic density (1 + x + xy)\n\n"
            << "   ini_velocity: 0 - constant velocity (1, 2)\n"
            << "                 1 - linear velocity (x, 2 y)\n"
            << "                 2 - quadratic velocity (x^2, 2 y^2)\n"
            << "                 3 - discontinuous velocity along central vertical line\n";
}


/* ******************************************************************
* Utility that computes corner centroid
****************************************************************** */
void corner_get_centroid(
    int cn, const Wonton::Jali_Mesh_Wrapper& mesh,
    Wonton::Point<2>& xcn)
{
  std::vector<int> wedges;
  std::array<Wonton::Point<2>, 3> wcoords;

  double volume = mesh.corner_volume(cn); 
  mesh.corner_get_wedges(cn, &wedges);

  xcn = {0.0, 0.0};
  for (auto w : wedges) {
    double frac = mesh.wedge_volume(w) / (3 * volume); 
    mesh.wedge_get_coordinates(w, &wcoords);
    for (int i = 0; i < 3; ++i) {
      xcn += frac * wcoords[i];
    }
  }
}


/* ******************************************************************
* App class that handles initialization and verification of fields
* since it is different in SCH and CCH methods.
****************************************************************** */
class MomentumRemap {
 public:
  MomentumRemap(int method) : method_(method) {};
  ~MomentumRemap() {};

  // initialization
  void InitMass(const Wonton::Jali_Mesh_Wrapper& mesh,
                int ini_method,
                std::vector<double>& mass);

  void InitVelocity(const Wonton::Jali_Mesh_Wrapper& mesh,
                    int ini_method,
                    std::vector<double>& ux,
                    std::vector<double>& uy);

  // field type
  Jali::Entity_kind MassKind() const {
    return (method_ == SGH) ? Jali::Entity_kind::CORNER : Jali::Entity_kind::CELL;
  }
  Jali::Entity_kind VelocityKind() const {
    return (method_ == SGH) ? Jali::Entity_kind::NODE : Jali::Entity_kind::CELL;
  }

  // V&V
  template<class T>
  double TotalMass(const Wonton::Jali_Mesh_Wrapper& mesh, const T& mass) const;

  template<class T> 
  Wonton::Point<2> TotalMomentum(const Wonton::Jali_Mesh_Wrapper& mesh,
                                 const T& mass, const T& ux, const T& uy) const;

  template<class T>
  Wonton::Point<2> VelocityMin(const Wonton::Jali_Mesh_Wrapper& mesh,
                               const T& ux, const T& uy);
  template<class T>
  Wonton::Point<2> VelocityMax(const Wonton::Jali_Mesh_Wrapper& mesh,
                               const T& ux, const T& uy);

  template<class T>
  void ErrorVelocity(const Wonton::Jali_Mesh_Wrapper& mesh,
                     int ini_method,
                     const T& ux, const T& uy,
                     double* l2err, double* l2norm);

 private:
  int method_;
};


/* ******************************************************************
* Initiaization of mass
****************************************************************** */
void MomentumRemap::InitMass(
    const Wonton::Jali_Mesh_Wrapper& mesh,
    int ini_method,
    std::vector<double>& mass)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_corners() + mesh.num_ghost_corners()
                               : mesh.num_owned_cells() + mesh.num_ghost_cells();
  mass.resize(nrows);

  double vol, den;
  Wonton::Point<2> xyz;

  for (int n = 0; n < nrows; ++n) {
    if (method_ == SGH) {
      corner_get_centroid(n, mesh, xyz);
      vol = mesh.corner_volume(n);
    } else {
      mesh.cell_centroid(n, &xyz);
      vol = mesh.cell_volume(n);
    }

    if (ini_method == 0) {
      den = 1.0;
    } else if (ini_method == 1) {
      den = 1.0 + xyz[0] + 2 * xyz[1];
    } else if (ini_method == 2) {
      den = 1.0 + xyz[0] + xyz[0] * xyz[1];
    }
    mass[n] = den * vol;
  }
}


/* ******************************************************************
* Initiaization of velocity
****************************************************************** */
void MomentumRemap::InitVelocity(
    const Wonton::Jali_Mesh_Wrapper& mesh,
    int ini_method,
    std::vector<double>& ux,
    std::vector<double>& uy)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() + mesh.num_ghost_nodes()
                               : mesh.num_owned_cells() + mesh.num_ghost_cells();
  ux.resize(nrows);
  uy.resize(nrows);

  Wonton::Point<2> xyz;

  for (int n = 0; n < nrows; ++n) {
    if (method_ == SGH)
      mesh.node_get_coordinates(n, &xyz);
    else
      mesh.cell_centroid(n, &xyz);

    if (ini_method == 0) {
      ux[n] = 1.0;
      uy[n] = 2.0;
    } else if (ini_method == 1) {
      ux[n] = xyz[0];
      uy[n] = 2.0 * xyz[1];
    } else if (ini_method == 2) {
      ux[n] = xyz[0] * xyz[0];
      uy[n] = 2.0 * xyz[1] * xyz[1];
    } else if (ini_method == 3) {
      ux[n] = (xyz[0] < 0.5) ? 1.0 : 2.0;
      uy[n] = 2.0 * xyz[1] * xyz[1];
    }
  }
}


/* ******************************************************************
* Calculate total momentum 
****************************************************************** */
template<class T>
double MomentumRemap::TotalMass(
    const Wonton::Jali_Mesh_Wrapper& mesh, const T& mass) const
{
  int nrows = (method_ == SGH) ? mesh.num_owned_corners() : mesh.num_owned_cells();
  double sum(0.0), sum_glb;
  for (int n = 0; n < nrows; ++n) sum += mass[n];

  MPI_Reduce(&sum, &sum_glb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  return sum_glb;
}


/* ******************************************************************
* Calculate total momentum 
****************************************************************** */
template<class T> 
Wonton::Point<2> MomentumRemap::TotalMomentum(
    const Wonton::Jali_Mesh_Wrapper& mesh,
    const T& mass, const T& ux, const T& uy) const
{
  double mx(0.0), my(0.0), mx_glb, my_glb;

  if (method_ == SGH) {
    std::vector<int> corners;

    for (int c = 0; c < mesh.num_owned_cells(); ++c) {
      mesh.cell_get_corners(c, &corners);

      for (auto cn : corners) { 
        int v = mesh.corner_get_node(cn);
        mx += mass[cn] * ux[v];
        my += mass[cn] * uy[v];
      }
    }
  }

  else if (method_ == CCH) {
    for (int c = 0; c < mesh.num_owned_cells(); ++c) {
      mx += mass[c] * ux[c];
      my += mass[c] * uy[c];
    }
  }

  MPI_Reduce(&mx, &mx_glb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&my, &my_glb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  Wonton::Point<2> momentum = {mx_glb, my_glb};
  return momentum;
}


/* ******************************************************************
* Velocity bounds
****************************************************************** */
template<class T>
Wonton::Point<2> MomentumRemap::VelocityMin(
    const Wonton::Jali_Mesh_Wrapper& mesh, const T& ux, const T& uy)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  double uxmin(ux[0]), uymin(uy[0]), uxmin_glb, uymin_glb;
 
  for (int n = 1; n < nrows; ++n) {
    uxmin = std::min(uxmin, ux[n]);
    uymin = std::min(uymin, uy[n]);
  }

  MPI_Reduce(&uxmin, &uxmin_glb, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&uymin, &uymin_glb, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

  Wonton::Point<2> umin = { uxmin_glb, uymin_glb };
  return umin;
}


template<class T>
Wonton::Point<2> MomentumRemap::VelocityMax(
    const Wonton::Jali_Mesh_Wrapper& mesh, const T& ux, const T& uy)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  double uxmax(ux[0]), uymax(uy[0]), uxmax_glb, uymax_glb;
 
  for (int n = 1; n < nrows; ++n) {
    uxmax = std::max(uxmax, ux[n]);
    uymax = std::max(uymax, uy[n]);
  }

  MPI_Reduce(&uxmax, &uxmax_glb, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&uymax, &uymax_glb, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  Wonton::Point<2> umax = { uxmax_glb, uymax_glb };
  return umax;
}


/* ******************************************************************
* Error in velocity
****************************************************************** */
template<class T>
void MomentumRemap::ErrorVelocity(
    const Wonton::Jali_Mesh_Wrapper& mesh,
    int ini_method,
    const T& ux, const T& uy,
    double* l2err, double* l2norm)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  double ux_exact, uy_exact;
  Wonton::Point<2> xyz;

  *l2err = 0.0;
  *l2norm = 0.0;

  for (int n = 0; n < nrows; ++n) {
    if (method_ == SGH)
      mesh.node_get_coordinates(n, &xyz);
    else
      mesh.cell_centroid(n, &xyz);

    if (ini_method == 0) {
      ux_exact = 1.0;
      uy_exact = 2.0;
    } else if (ini_method == 1) {
      ux_exact = xyz[0];
      uy_exact = 2.0 * xyz[1];
    } else if (ini_method == 2) {
      ux_exact = xyz[0] * xyz[0];
      uy_exact = 2.0 * xyz[1] * xyz[1];
    } else if (ini_method == 3) {
      ux_exact = (xyz[0] < 0.5) ? 1.0 : 2.0;
      uy_exact = 2.0 * xyz[1] * xyz[1];
    }

    *l2err += (ux_exact - ux[n]) * (ux_exact - ux[n])
            + (uy_exact - uy[n]) * (uy_exact - uy[n]);

    *l2norm += ux_exact * ux_exact + uy_exact * uy_exact;
  }

  int nrows_glb;
  double l2err_glb, l2norm_glb;

  MPI_Reduce(&nrows, &nrows_glb, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(l2err, &l2err_glb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(l2norm, &l2norm_glb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  *l2err = std::sqrt(l2err_glb / nrows_glb);
  *l2norm = std::sqrt(l2norm_glb / nrows_glb);
}


/* ******************************************************************
* Main driver for the momentum remap
****************************************************************** */
int main(int argc, char** argv) {
  // initialize MPI
  int mpi_init_flag;
  MPI_Initialized(&mpi_init_flag);
  if (!mpi_init_flag)
      MPI_Init(&argc, &argv);

  int numpe, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc < 7) {
    if (rank == 0) print_usage();
    MPI_Finalize();
    return 0;
  }

  // number of cells in x and y directions for the source mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  auto method = atoi(argv[3]);
  auto limiter = (atoi(argv[4]) == 0) ? Portage::NOLIMITER: Portage::BARTH_JESPERSEN;
  auto ini_density = atoi(argv[5]);
  auto ini_velocity = atoi(argv[6]);

  if ((method != SGH && method != CCH) || 
      (ini_density < 0 || ini_density > 2) ||
      (ini_velocity < 0 || ini_velocity > 3)) {
    if (rank == 0) {
      std::cout << "=== Input ERROR ===\n";
      print_usage();    
    }
    MPI_Finalize();
    return 0;
  }

  if (numpe > 1 && method == SGH) {
    if (rank == 0) {
      std::cout << "=== Input ERROR ===\n";
      std::cout << "method=SGH runs only in the serial mode, since ghost data\n";
      std::cout << "           on a target mesh cannot be populated easily\n";
    }
    MPI_Finalize();
    return 0;
  }

  // size of computational domain
  double lenx = 1.0;  // [m]
  double leny = 1.0;
  
  //
  // Preliminaries, common for SGH and CCH
  //

  // -- setup Jali meshes
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);
  mesh_factory.included_entities(Jali::Entity_kind::ALL_KIND);

  auto srcmesh = mesh_factory(0.0, 0.0, lenx, leny, nx, ny);
  auto trgmesh = mesh_factory(0.0, 0.0, lenx, leny, nx + 2, ny + 4);

  // -- setup mesh wrappers
  Wonton::Jali_Mesh_Wrapper srcmesh_wrapper(*srcmesh);
  Wonton::Jali_Mesh_Wrapper trgmesh_wrapper(*trgmesh);

  int ncells_src = srcmesh_wrapper.num_owned_cells() + srcmesh_wrapper.num_ghost_cells();
  int nnodes_src = srcmesh_wrapper.num_owned_nodes() + srcmesh_wrapper.num_ghost_nodes();

  int ncells_trg = trgmesh_wrapper.num_owned_cells();
  int nnodes_trg = trgmesh_wrapper.num_owned_nodes();
  int ncorners_trg = trgmesh_wrapper.num_owned_corners();

  // -- states
  std::shared_ptr<Jali::State> srcstate = Jali::State::create(srcmesh);
  std::shared_ptr<Jali::State> trgstate = Jali::State::create(trgmesh);

  // -- state wrappers
  Wonton::Jali_State_Wrapper srcstate_wrapper(*srcstate);
  Wonton::Jali_State_Wrapper trgstate_wrapper(*trgstate);

  // -- register velocity components with the states
  MomentumRemap mr(method);

  Jali::Entity_kind kind;
  Jali::Entity_type type = Jali::Entity_type::ALL;
  std::vector<double> ux_src(nnodes_src);
  std::vector<double> uy_src(nnodes_src);

  kind = mr.VelocityKind();
  mr.InitVelocity(srcmesh_wrapper, ini_velocity, ux_src, uy_src);

  srcstate->add("velocity_x", srcmesh, kind, type, &(ux_src[0]));
  srcstate->add("velocity_y", srcmesh, kind, type, &(uy_src[0]));

  trgstate->add<double, Jali::Mesh, Jali::UniStateVector>("velocity_x", trgmesh, kind, type);
  trgstate->add<double, Jali::Mesh, Jali::UniStateVector>("velocity_y", trgmesh, kind, type);

  // -- register mass with the states
  std::vector<double> mass_src;

  kind = mr.MassKind();
  mr.InitMass(srcmesh_wrapper, ini_density, mass_src);

  srcstate->add("mass", srcmesh, kind, type, &(mass_src[0]));
  trgstate->add<double, Jali::Mesh, Jali::UniStateVector>(
                "mass", trgmesh, kind, type);

  // -- summary
  auto total_mass_src = mr.TotalMass(srcmesh_wrapper, mass_src);
  auto total_momentum_src = mr.TotalMomentum(srcmesh_wrapper, mass_src, ux_src, uy_src);
  auto umin = mr.VelocityMin(srcmesh_wrapper, ux_src, uy_src);
  auto umax = mr.VelocityMax(srcmesh_wrapper, ux_src, uy_src);
  if (rank == 0) {
    std::cout << "=== SOURCE data ===" << std::endl;
    std::cout << "mesh:           " << nx << " x " << ny << std::endl;
    std::cout << "total mass:     " << total_mass_src << " kg" << std::endl;
    std::cout << "total momentum: " << total_momentum_src << " kg m/s" << std::endl;
    std::cout << "limiter:        " << ((limiter == Portage::NOLIMITER) ? "none\n" : "BJ\n");
    std::cout << "velocity bounds," << " min: " << umin << " max: " << umax << std::endl;
  }

  //
  // SEVEN-step REMAP algorithm 
  //

  // Step 1 (SGH only)
  // -- gather cell-centered mass from corner masses
  std::vector<double> mass_c(ncells_src, 0.0);
  std::vector<int> corners;

  if (method == SGH) {
    for (int c = 0; c < ncells_src; ++c) {
      srcmesh_wrapper.cell_get_corners(c, &corners);
      for (auto cn : corners) mass_c[c] += mass_src[cn];
    }
  }

  // Step 2 (SGH and CCH)
  // -- compute density in each mesh cell
  std::vector<double> density(ncells_src);
  for (int c = 0; c < ncells_src; ++c) {
    double tmp = (method == SGH) ? mass_c[c] : mass_src[c];
    density[c] = tmp / srcmesh_wrapper.cell_volume(c);
  }
  
  // Step 3 (SGH and CCH)
  // -- compute cell-centered specific momentum
  std::vector<double> momentum_x_src(ncells_src, 0.0);
  std::vector<double> momentum_y_src(ncells_src, 0.0);

  for (int c = 0; c < ncells_src; ++c) {
    if (method == SGH) {
      srcmesh_wrapper.cell_get_corners(c, &corners);

      for (auto cn : corners) { 
        int v = srcmesh_wrapper.corner_get_node(cn);
        momentum_x_src[c] += mass_src[cn] * ux_src[v];
        momentum_y_src[c] += mass_src[cn] * uy_src[v];
      }
    } else {
      momentum_x_src[c] += mass_src[c] * ux_src[c];
      momentum_y_src[c] += mass_src[c] * uy_src[c];
    }

    double volume = srcmesh_wrapper.cell_volume(c);
    momentum_x_src[c] /= volume;
    momentum_y_src[c] /= volume;
  }

  // Step 4 (SGH and CCH)
  // -- remap density and specific momentum following three basic
  //    steps: search, intersect, and interpolate
  Portage::CoreDriver<2, Wonton::Entity_kind::CELL,
                      Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
      cd(srcmesh_wrapper, srcstate_wrapper,
         trgmesh_wrapper, trgstate_wrapper);

  Portage::NumericTolerances_t default_num_tols;
  default_num_tols.use_default();
  cd.set_num_tols(default_num_tols);

  auto candidates = cd.search<Portage::SearchKDTree>();
  auto srcwts = cd.intersect_meshes<Portage::IntersectR2D>(candidates);

  double dblmin = -std::numeric_limits<double>::max();
  double dblmax =  std::numeric_limits<double>::max();

  // -- we need to register fields that we want to remap
  std::vector<std::string> field_names;
  std::vector<const double*> field_pointers;

  field_names.push_back("density");
  field_pointers.push_back(&(density[0]));

  field_names.push_back("momentum_x");
  field_pointers.push_back(&(momentum_x_src[0]));

  field_names.push_back("momentum_y");
  field_pointers.push_back(&(momentum_y_src[0]));

  for (int i = 0; i < field_names.size(); ++i) {
    srcstate->add(field_names[i], srcmesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL, field_pointers[i]);

    trgstate->add<double, Jali::Mesh, Jali::UniStateVector>(
                  field_names[i], trgmesh, Jali::Entity_kind::CELL,
                  Jali::Entity_type::ALL);

    auto gradients = cd.compute_source_gradient(field_names[i], limiter);

    cd.interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
        field_names[i], field_names[i], srcwts, dblmin, dblmax,
        nullptr, &gradients);
  }

  // Step 5 (SGH only)
  // -- create linear reconstruction (limited or unlimited) 
  //    of density and specific momentum on the target mesh
  std::vector<Portage::vector<Wonton::Vector<2>>> gradients;

  if (method == SGH) {
    for (int i = 0; i < field_names.size(); ++i) {
      Portage::Limited_Gradient<2, Wonton::Entity_kind::CELL, 
                                Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper>
          gradient_kernel(trgmesh_wrapper, trgstate_wrapper,
                          field_names[i], limiter, Portage::BND_NOLIMITER);

      int ncells_all = ncells_trg + trgmesh_wrapper.num_ghost_cells();
      Portage::vector<Wonton::Vector<2>> gradient(ncells_all);

      Portage::transform(trgmesh_wrapper.begin(Wonton::Entity_kind::CELL),
                         trgmesh_wrapper.end(Wonton::Entity_kind::CELL),
                         gradient.begin(), gradient_kernel);

      gradients.push_back(gradient);  // could be optimized
    }
  }

  // Step 6 (SGH and CCH)
  // -- integrate density and specific momentum on the target mesh
  Jali::UniStateVector<double, Jali::Mesh> mass_trg;
  Jali::UniStateVector<double, Jali::Mesh> ux_trg, uy_trg;

  kind = mr.MassKind();
  type = Jali::Entity_type::ALL;
  trgstate->get("mass", trgmesh, kind, type, &mass_trg);

  kind = mr.VelocityKind();
  trgstate->get("velocity_x", trgmesh, kind, type, &ux_trg);
  trgstate->get("velocity_y", trgmesh, kind, type, &uy_trg);

  //    extract auxiliary (cell-centerd) data
  Jali::UniStateVector<double, Jali::Mesh> density_trg;
  Jali::UniStateVector<double, Jali::Mesh> momentum_x_trg, momentum_y_trg;

  kind = Jali::Entity_kind::CELL;
  trgstate->get("density", trgmesh, kind, type, &density_trg);

  trgstate->get("momentum_x", trgmesh, kind, type, &momentum_x_trg);
  trgstate->get("momentum_y", trgmesh, kind, type, &momentum_y_trg);

  //    integrate
  std::vector<double> momentum_cn_x, momentum_cn_y;  // work memory

  if (method == SGH) {
    momentum_cn_x.resize(ncorners_trg);
    momentum_cn_y.resize(ncorners_trg);

    Wonton::Point<2> xc, xcn;
    Wonton::Vector<2> grad;
    std::vector<int> wedges;
    std::array<Wonton::Point<2>, 3> wcoords;

    for (int c = 0; c < ncells_trg; ++c) {
      trgmesh_wrapper.cell_centroid(c, &xc);
      trgmesh_wrapper.cell_get_corners(c, &corners);

      for (auto cn : corners) {
        double volume = trgmesh_wrapper.corner_volume(cn); 
        trgmesh_wrapper.corner_get_wedges(cn, &wedges);

        corner_get_centroid(cn, trgmesh_wrapper, xcn);

        // integral is the value at centroid
        double vol = trgmesh_wrapper.corner_volume(cn);

        grad = gradients[0][c];
        mass_trg[cn] = vol * (density_trg[c] + dot(grad, xcn - xc));

        grad = gradients[1][c];
        momentum_cn_x[cn] = vol * (momentum_x_trg[c] + dot(grad, xcn - xc));

        grad = gradients[2][c];
        momentum_cn_y[cn] = vol * (momentum_y_trg[c] + dot(grad, xcn - xc));
      }
    }
  } 
  else {
    for (int c = 0; c < ncells_trg; ++c) {
      double vol = trgmesh_wrapper.cell_volume(c);
      mass_trg[c] = density_trg[c] * vol;

      ux_trg[c] = momentum_x_trg[c] / density_trg[c];
      uy_trg[c] = momentum_y_trg[c] / density_trg[c];
    }
  }

  // Step 7 (SGH only)
  // -- gather data to mesh nodes and compute velocity
  if (method == SGH) {
    int nnodes_all = nnodes_trg + trgmesh_wrapper.num_ghost_nodes();
    std::vector<double> mass_v(nnodes_all, 0.0);  // work memory
    std::vector<double> momentum_v_x(nnodes_all, 0.0);
    std::vector<double> momentum_v_y(nnodes_all, 0.0);

    for (int cn = 0; cn < ncorners_trg; ++cn) {
      int v = trgmesh_wrapper.corner_get_node(cn);
      mass_v[v] += mass_trg[cn];
      momentum_v_x[v] += momentum_cn_x[cn];
      momentum_v_y[v] += momentum_cn_y[cn];
    }

    for (int v = 0; v < nnodes_trg; ++v) {
      ux_trg[v] += momentum_v_x[v] / mass_v[v];
      uy_trg[v] += momentum_v_y[v] / mass_v[v];
    }
  }

  //
  // Verification
  //
  auto total_mass_trg = mr.TotalMass(trgmesh_wrapper, mass_trg);
  auto total_momentum_trg = mr.TotalMomentum(trgmesh_wrapper, mass_trg, ux_trg, uy_trg);
  umin = mr.VelocityMin(trgmesh_wrapper, ux_trg, uy_trg);
  umax = mr.VelocityMax(trgmesh_wrapper, ux_trg, uy_trg);

  if (rank == 0) {
    std::cout << "\n=== TARGET data ===" << std::endl;
    std::cout << "mesh:           " << nx + 2 << " x " << ny + 4 << std::endl;
    std::cout << "total mass:     " << total_mass_trg << " kg" << std::endl;
    std::cout << "total momentum: " << total_momentum_trg << " kg m/s" << std::endl;
    std::cout << "velocity bounds," << " min: " << umin << " max: " << umax << std::endl;
  }

  auto err = total_momentum_trg - total_momentum_src;
  double cons_law0 = std::fabs(total_mass_trg - total_mass_src);
  double cons_law1 = std::sqrt(err[0] * err[0] + err[1] * err[1]);

  if (rank == 0) {
    std::cout << "\n=== Conservation error ===" << std::endl;
    std::cout << "in total mass:     " << cons_law0 << std::endl;
    std::cout << "in total momentum: " << cons_law1 << std::endl;
  }

  double l2err, l2norm;
  mr.ErrorVelocity(trgmesh_wrapper, ini_velocity, ux_trg, uy_trg, &l2err, &l2norm);

  if (rank == 0) {
    std::cout << "\n=== Remap error ===" << std::endl;
    std::cout << "in velocity: l2-err=" << l2err << " l2-norm=" << l2norm << std::endl;
  }

  // save data
  std::ofstream datafile;
  datafile.open("errors.txt");
  datafile << "0 " << cons_law0 << std::endl;
  datafile << "1 " << cons_law1 << std::endl;
  datafile << "2 " << l2err << std::endl;
  datafile << "3 " << l2norm << std::endl;
  datafile.close();

  MPI_Finalize();
  return 0;
}
