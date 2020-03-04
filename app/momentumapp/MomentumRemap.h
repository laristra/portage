/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef MOMENTUM_REMAP_ND_HH_
#define MOMENTUM_REMAP_ND_HH_

#include <fstream>
#include <sstream>
#include <vector>

// portage includes
#include "portage/driver/coredriver.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/interpolate/gradient.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/search/search_kdtree.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/support/Point.h"

// App includes
#include "corner_get_centroid.h"

const int SGH = 1;
const int CCH = 2;

/* ******************************************************************
* App class that handles initialization and verification of fields
* since it is different in SCH and CCH methods.
****************************************************************** */
template<int D, class Mesh_Wrapper>
class MomentumRemap {
 public:
  MomentumRemap(int method) : method_(method) {};
  ~MomentumRemap() {};

  // initialization
  void InitMass(
      const Mesh_Wrapper& mesh,
      user_field_t& formula, std::vector<double>& mass);

  void InitVelocity(
      const Mesh_Wrapper& mesh,
      user_field_t& formula, std::vector<double>& u);

  // field type
  Wonton::Entity_kind MassKind() const {
    return (method_ == SGH) ? Wonton::Entity_kind::CORNER : Wonton::Entity_kind::CELL;
  }
  Wonton::Entity_kind VelocityKind() const {
    return (method_ == SGH) ? Wonton::Entity_kind::NODE : Wonton::Entity_kind::CELL;
  }

  // main remap method
  template<class State_Wrapper>
  void RemapND(
      const Mesh_Wrapper& srcmesh_wrapper, State_Wrapper& srcstate_wrapper,
      const Mesh_Wrapper& trgmesh_wrapper, State_Wrapper& trgstate_wrapper,
      Portage::Limiter_type limiter);

  // V&V
  template<class T>
  double TotalMass(const Mesh_Wrapper& mesh, const T& mass) const;

  template<class T> 
  Wonton::Point<D> TotalMomentum(
      const Mesh_Wrapper& mesh,
      const T& mass, const T& ux, const T& uy, const T& uz) const;

  template<class T>
  Wonton::Point<D> VelocityMin(
      const Mesh_Wrapper& mesh, const T& ux, const T& uy, const T& uz);
  template<class T>
  Wonton::Point<D> VelocityMax(
      const Mesh_Wrapper& mesh, const T& ux, const T& uy, const T& uz);

  template<class T>
  void ErrorVelocity(
      const Mesh_Wrapper& mesh,
      user_field_t& formula_x, user_field_t& formula_y, user_field_t& formula_z,
      const T& ux, const T& uy, const T& uz,
      double* l2err, double* l2norm);

 private:
  int method_;
};


/* ******************************************************************
* Initiaization of mass
****************************************************************** */
template<int D, class Mesh_Wrapper>
void MomentumRemap<D, Mesh_Wrapper>::InitMass(
    const Mesh_Wrapper& mesh,
    user_field_t& formula, std::vector<double>& mass)
  {
  int nrows = (method_ == SGH) ? mesh.num_owned_corners() + mesh.num_ghost_corners()
                               : mesh.num_owned_cells() + mesh.num_ghost_cells();
  mass.resize(nrows);

  double vol, den;
  Wonton::Point<D> xyz;

  for (int n = 0; n < nrows; ++n) {
    if (method_ == SGH) {
      corner_get_centroid(n, mesh, xyz);
      vol = mesh.corner_volume(n);
    } else {
      mesh.cell_centroid(n, &xyz);
      vol = mesh.cell_volume(n);
    }

    den = formula(xyz);
    mass[n] = den * vol;
  }
}


/* ******************************************************************
* Initiaization of velocity
****************************************************************** */
template<int D, class Mesh_Wrapper>
void MomentumRemap<D, Mesh_Wrapper>::InitVelocity(
    const Mesh_Wrapper& mesh, user_field_t& formula, std::vector<double>& u)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() + mesh.num_ghost_nodes()
                               : mesh.num_owned_cells() + mesh.num_ghost_cells();
  u.resize(nrows);

  Wonton::Point<D> xyz;

  for (int n = 0; n < nrows; ++n) {
    if (method_ == SGH)
      mesh.node_get_coordinates(n, &xyz);
    else
      mesh.cell_centroid(n, &xyz);

    u[n] = formula(xyz);
  }
}


/* ******************************************************************
* Calculate total momentum 
****************************************************************** */
template<int D, class Mesh_Wrapper>
template<class T>
double MomentumRemap<D, Mesh_Wrapper>::TotalMass(
    const Mesh_Wrapper& mesh, const T& mass) const
{
  int nrows = (method_ == SGH) ? mesh.num_owned_corners() : mesh.num_owned_cells();
  double sum(0.0), sum_glb;
  for (int n = 0; n < nrows; ++n) sum += mass[n];

  MPI_Allreduce(&sum, &sum_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sum_glb;
}


/* ******************************************************************
* Calculate total momentum 
****************************************************************** */
template<int D, class Mesh_Wrapper>
template<class T> 
Wonton::Point<D> MomentumRemap<D, Mesh_Wrapper>::TotalMomentum(
    const Mesh_Wrapper& mesh,
    const T& mass, const T& ux, const T& uy, const T& uz) const
{
  double mx(0.0), my(0.0), mz(0.0), mx_glb, my_glb, mz_glb;

  if (method_ == SGH) {
    std::vector<int> corners;

    for (int c = 0; c < mesh.num_owned_cells(); ++c) {
      mesh.cell_get_corners(c, &corners);

      for (auto cn : corners) { 
        int v = mesh.corner_get_node(cn);
        mx += mass[cn] * ux[v];
        my += mass[cn] * uy[v];
        if (D == 3) mz += mass[cn] * uz[v];
      }
    }
  }

  else if (method_ == CCH) {
    for (int c = 0; c < mesh.num_owned_cells(); ++c) {
      mx += mass[c] * ux[c];
      my += mass[c] * uy[c];
      if (D == 3) mz += mass[c] * uz[c];
    }
  }

  Wonton::Point<D> momentum;

  MPI_Allreduce(&mx, &mx_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&my, &my_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (D == 2) {
    momentum = {mx_glb, my_glb};
  } else {
    MPI_Allreduce(&mz, &mz_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    momentum = {mx_glb, my_glb, mz_glb};
  }

  return momentum;
}


/* ******************************************************************
* Velocity bounds
****************************************************************** */
template<int D, class Mesh_Wrapper>
template<class T>
Wonton::Point<D> MomentumRemap<D, Mesh_Wrapper>::VelocityMin(
    const Mesh_Wrapper& mesh,
    const T& ux, const T& uy, const T& uz)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  double uxmin(ux[0]), uymin(uy[0]), uzmin(uz[0]), uxmin_glb, uymin_glb, uzmin_glb;
 
  for (int n = 1; n < nrows; ++n) {
    uxmin = std::min(uxmin, ux[n]);
    uymin = std::min(uymin, uy[n]);
    if (D == 3) uzmin = std::min(uzmin, uz[n]);
  }

  MPI_Allreduce(&uxmin, &uxmin_glb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&uymin, &uymin_glb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  Wonton::Point<D> umin;
  if (D == 2) {
    umin = { uxmin_glb, uymin_glb };
  } else {
    MPI_Allreduce(&uzmin, &uzmin_glb, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    umin = { uxmin_glb, uymin_glb, uzmin_glb };
  }

  return umin;
}


template<int D, class Mesh_Wrapper>
template<class T>
Wonton::Point<D> MomentumRemap<D, Mesh_Wrapper>::VelocityMax(
    const Mesh_Wrapper& mesh, const T& ux, const T& uy, const T& uz)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  double uxmax(ux[0]), uymax(uy[0]), uzmax(uz[0]), uxmax_glb, uymax_glb, uzmax_glb;
 
  for (int n = 1; n < nrows; ++n) {
    uxmax = std::max(uxmax, ux[n]);
    uymax = std::max(uymax, uy[n]);
    if (D == 3) uzmax = std::max(uzmax, uz[n]);
  }

  MPI_Allreduce(&uxmax, &uxmax_glb, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&uymax, &uymax_glb, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  Wonton::Point<D> umax;
  if (D == 2) {
    umax = { uxmax_glb, uymax_glb };
  } else {
    MPI_Allreduce(&uzmax, &uzmax_glb, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    umax = { uxmax_glb, uymax_glb, uzmax_glb };
  }
  return umax;
}


/* ******************************************************************
* Error in velocity
****************************************************************** */
template<int D, class Mesh_Wrapper>
template<class T>
void MomentumRemap<D, Mesh_Wrapper>::ErrorVelocity(
    const Mesh_Wrapper& mesh,
    user_field_t& formula_x, user_field_t& formula_y, user_field_t& formula_z,
    const T& ux, const T& uy, const T& uz,
    double* l2err, double* l2norm)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  double ux_exact, uy_exact, uz_exact;
  Wonton::Point<D> xyz;

  *l2err = 0.0;
  *l2norm = 0.0;

  for (int n = 0; n < nrows; ++n) {
    if (method_ == SGH)
      mesh.node_get_coordinates(n, &xyz);
    else
      mesh.cell_centroid(n, &xyz);

    ux_exact = formula_x(xyz);
    uy_exact = formula_y(xyz);

    *l2err += (ux_exact - ux[n]) * (ux_exact - ux[n])
            + (uy_exact - uy[n]) * (uy_exact - uy[n]);

    *l2norm += ux_exact * ux_exact + uy_exact * uy_exact;

    if (D == 3) {
      uz_exact = formula_z(xyz);
      *l2err += (uz_exact - uz[n]) * (uz_exact - uz[n]);
      *l2norm += uz_exact * uz_exact;
    }
  }

  int nrows_glb;
  double l2err_glb, l2norm_glb;

  MPI_Allreduce(&nrows, &nrows_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l2err, &l2err_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l2norm, &l2norm_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *l2err = std::sqrt(l2err_glb / nrows_glb);
  *l2norm = std::sqrt(l2norm_glb / nrows_glb);
}


/* ******************************************************************
* 2D and 3D remap algorithm
****************************************************************** */
template<int D, class Mesh_Wrapper>
template<class State_Wrapper>
inline
void MomentumRemap<D, Mesh_Wrapper>::RemapND(
    const Mesh_Wrapper& srcmesh_wrapper, State_Wrapper& srcstate_wrapper,
    const Mesh_Wrapper& trgmesh_wrapper, State_Wrapper& trgstate_wrapper,
    Portage::Limiter_type limiter)
{
  // mesh data
  int ncells_src = srcmesh_wrapper.num_owned_cells() + srcmesh_wrapper.num_ghost_cells();
  int ncells_trg = trgmesh_wrapper.num_owned_cells();
  int nnodes_trg = trgmesh_wrapper.num_owned_nodes();
  int ncorners_trg = trgmesh_wrapper.num_owned_corners();

  // state fields
  auto kind = MassKind();

  const double *mass_src, *ux_src, *uy_src, *uz_src;
  srcstate_wrapper.mesh_get_data(kind, "mass", &mass_src);

  kind = VelocityKind();
  srcstate_wrapper.mesh_get_data(kind, "velocity_x", &ux_src);
  srcstate_wrapper.mesh_get_data(kind, "velocity_y", &uy_src);
  if (D == 3) srcstate_wrapper.mesh_get_data(kind, "velocity_z", &uz_src);

  // Step 1 (SGH only)
  // -- gather cell-centered mass from corner masses
  std::vector<double> mass_c(ncells_src, 0.0);
  std::vector<int> corners;

  if (method_ == SGH) {
    for (int c = 0; c < ncells_src; ++c) {
      srcmesh_wrapper.cell_get_corners(c, &corners);
      for (auto cn : corners) mass_c[c] += mass_src[cn];
    }
  }

  // Step 2 (SGH and CCH)
  // -- compute density in each mesh cell
  std::vector<double> density(ncells_src);
  for (int c = 0; c < ncells_src; ++c) {
    double tmp = (method_ == SGH) ? mass_c[c] : mass_src[c];
    density[c] = tmp / srcmesh_wrapper.cell_volume(c);
  }
  
  // Step 3 (SGH and CCH)
  // -- compute cell-centered specific momentum
  std::vector<double> momentum_x_src(ncells_src, 0.0);
  std::vector<double> momentum_y_src(ncells_src, 0.0);
  std::vector<double> momentum_z_src(ncells_src, 0.0);

  for (int c = 0; c < ncells_src; ++c) {
    if (method_ == SGH) {
      srcmesh_wrapper.cell_get_corners(c, &corners);

      for (auto cn : corners) { 
        int v = srcmesh_wrapper.corner_get_node(cn);
        momentum_x_src[c] += mass_src[cn] * ux_src[v];
        momentum_y_src[c] += mass_src[cn] * uy_src[v];
        if (D == 3) momentum_z_src[c] += mass_src[cn] * uz_src[v];
      }
    } else {
      momentum_x_src[c] += mass_src[c] * ux_src[c];
      momentum_y_src[c] += mass_src[c] * uy_src[c];
      if (D == 3) momentum_z_src[c] += mass_src[c] * uz_src[c];
    }

    double volume = srcmesh_wrapper.cell_volume(c);
    momentum_x_src[c] /= volume;
    momentum_y_src[c] /= volume;
    if (D == 3) momentum_z_src[c] /= volume;
  }

  // Step 4 (SGH and CCH)
  // -- remap density and specific momentum following three basic
  //    steps: search, intersect, and interpolate
  Portage::CoreDriver<D, Wonton::Entity_kind::CELL, Mesh_Wrapper, State_Wrapper>
      cd(srcmesh_wrapper, srcstate_wrapper,
         trgmesh_wrapper, trgstate_wrapper);

  Portage::NumericTolerances_t default_num_tols;
  default_num_tols.use_default();
  cd.set_num_tols(default_num_tols);

  auto candidates = cd.template search<Portage::SearchKDTree>();
  auto srcwts = cd.template intersect_meshes<Portage::IntersectRND<D>::template Intersect>(candidates);

  // -- we need to register fields that we want to remap
  std::vector<std::string> field_names;
  std::vector<const double*> field_pointers;

  field_names.push_back("density");
  field_pointers.push_back(&(density[0]));

  field_names.push_back("momentum_x");
  field_pointers.push_back(&(momentum_x_src[0]));

  field_names.push_back("momentum_y");
  field_pointers.push_back(&(momentum_y_src[0]));

  if (D == 3) {
    field_names.push_back("momentum_z");
    field_pointers.push_back(&(momentum_z_src[0]));
  }

  std::vector<double> tmp(ncells_trg);
  for (int i = 0; i < field_names.size(); ++i) {
    srcstate_wrapper.mesh_add_data(Wonton::Entity_kind::CELL, field_names[i], field_pointers[i]);
    trgstate_wrapper.mesh_add_data(Wonton::Entity_kind::CELL, field_names[i], &(tmp[0]));

    auto gradients = cd.compute_source_gradient(field_names[i], limiter);

    cd.template interpolate_mesh_var<double, Portage::Interpolate_2ndOrder>(
        field_names[i], field_names[i], srcwts, &gradients);
  }

  // Step 5 (SGH only)
  // -- create linear reconstruction (limited or unlimited) 
  //    of density and specific momentum on the target mesh
  int ncells_all = ncells_trg + trgmesh_wrapper.num_ghost_cells();
  std::vector<Portage::vector<Wonton::Vector<D>>> gradients(
      field_names.size(), Portage::vector<Wonton::Vector<D>>(ncells_all));

  if (method_ == SGH) {
    for (int i = 0; i < field_names.size(); ++i) {
      Portage::Limited_Gradient<D, Wonton::Entity_kind::CELL, 
                                Mesh_Wrapper, State_Wrapper>
          gradient_kernel(trgmesh_wrapper, trgstate_wrapper,
                          field_names[i], limiter, Portage::BND_NOLIMITER);

      Portage::vector<Wonton::Vector<D>> gradient(ncells_all);

      Portage::transform(trgmesh_wrapper.begin(Wonton::Entity_kind::CELL),
                         trgmesh_wrapper.end(Wonton::Entity_kind::CELL),
                         gradients[i].begin(), gradient_kernel);
    }
  }

  // Step 6 (SGH and CCH)
  // -- integrate density and specific momentum on the target mesh
  double *mass_trg, *ux_trg, *uy_trg, *uz_trg;

  kind = MassKind();
  trgstate_wrapper.mesh_get_data(kind, "mass", &mass_trg);

  kind = VelocityKind();
  trgstate_wrapper.mesh_get_data(kind, "velocity_x", &ux_trg);
  trgstate_wrapper.mesh_get_data(kind, "velocity_y", &uy_trg);
  if (D == 3) trgstate_wrapper.mesh_get_data(kind, "velocity_z", &uz_trg);

  //    extract auxiliary (cell-centerd) data
  double *density_trg, *momentum_x_trg, *momentum_y_trg, *momentum_z_trg;

  kind = Wonton::Entity_kind::CELL;
  trgstate_wrapper.mesh_get_data(kind, "density", &density_trg);
  trgstate_wrapper.mesh_get_data(kind, "momentum_x", &momentum_x_trg);
  trgstate_wrapper.mesh_get_data(kind, "momentum_y", &momentum_y_trg);
  if (D == 3) trgstate_wrapper.mesh_get_data(kind, "momentum_z", &momentum_z_trg);

  //    integrate
  std::vector<double> momentum_cn_x, momentum_cn_y, momentum_cn_z;  // work memory

  if (method_ == SGH) {
    momentum_cn_x.resize(ncorners_trg);
    momentum_cn_y.resize(ncorners_trg);
    if (D == 3) momentum_cn_z.resize(ncorners_trg);

    Wonton::Point<D> xc, xcn;
    Wonton::Vector<D> grad;
    std::vector<int> wedges;
    std::array<Wonton::Point<D>, D + 1> wcoords;

    for (int c = 0; c < ncells_trg; ++c) {
      trgmesh_wrapper.cell_centroid(c, &xc);
      trgmesh_wrapper.cell_get_corners(c, &corners);

      for (auto cn : corners) {
        double volume = trgmesh_wrapper.corner_volume(cn); 
        trgmesh_wrapper.corner_get_wedges(cn, &wedges);

        corner_get_centroid(cn, trgmesh_wrapper, xcn);

        // integral is the value at centroid
        double cnvol = trgmesh_wrapper.corner_volume(cn);

        grad = gradients[0][c];
        mass_trg[cn] = cnvol * (density_trg[c] + dot(grad, xcn - xc));

        grad = gradients[1][c];
        momentum_cn_x[cn] = cnvol * (momentum_x_trg[c] + dot(grad, xcn - xc));

        grad = gradients[2][c];
        momentum_cn_y[cn] = cnvol * (momentum_y_trg[c] + dot(grad, xcn - xc));

        if (D == 3) {
          grad = gradients[3][c];
          momentum_cn_z[cn] = cnvol * (momentum_z_trg[c] + dot(grad, xcn - xc));
        }
      }
    }
  } 
  else {
    for (int c = 0; c < ncells_trg; ++c) {
      double vol = trgmesh_wrapper.cell_volume(c);
      mass_trg[c] = density_trg[c] * vol;

      ux_trg[c] = momentum_x_trg[c] / density_trg[c];
      uy_trg[c] = momentum_y_trg[c] / density_trg[c];
      if (D == 3) uz_trg[c] = momentum_z_trg[c] / density_trg[c];
    }
  }

  // Step 7 (SGH only)
  // -- gather data to mesh nodes and compute velocity
  if (method_ == SGH) {
    int nnodes_all = nnodes_trg + trgmesh_wrapper.num_ghost_nodes();
    std::vector<double> mass_v(nnodes_all, 0.0);  // work memory
    std::vector<double> momentum_v_x(nnodes_all, 0.0);
    std::vector<double> momentum_v_y(nnodes_all, 0.0);
    std::vector<double> momentum_v_z(nnodes_all, 0.0);

    for (int cn = 0; cn < ncorners_trg; ++cn) {
      int v = trgmesh_wrapper.corner_get_node(cn);
      mass_v[v] += mass_trg[cn];
      momentum_v_x[v] += momentum_cn_x[cn];
      momentum_v_y[v] += momentum_cn_y[cn];
      if (D == 3) momentum_v_z[v] += momentum_cn_z[cn];
    }

    for (int v = 0; v < nnodes_trg; ++v) {
      ux_trg[v] = momentum_v_x[v] / mass_v[v];
      uy_trg[v] = momentum_v_y[v] / mass_v[v];
      if (D == 3) uz_trg[v] = momentum_v_z[v] / mass_v[v];
    }
  }
}

#endif
