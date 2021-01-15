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

// wonton includes
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

// portage includes
#include "portage/driver/coredriver.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/interpolate/gradient.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/search/search_kdtree.h"
#include "portage/support/portage.h"

// App includes
#include "corner_get_centroid.h"
#include "MomentumRemapDefs.h"

const int SGH = 1;
const int CCH = 2;

/* ******************************************************************
* App class that handles initialization and verification of fields
* since it is different in SCH and CCH methods.
****************************************************************** */
template<int D, class Mesh_Wrapper>
class MomentumRemap {
 public:
  explicit MomentumRemap(int method) : method_(method) {}
  ~MomentumRemap() = default;

  // initialization using formula for density
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
      const Mesh_Wrapper& mesh, const T& mass, const T u[D]) const;

  template<class T>
  Wonton::Point<D> VelocityMin(const Mesh_Wrapper& mesh, const T u[D]);
  template<class T>
  Wonton::Point<D> VelocityMax(const Mesh_Wrapper& mesh, const T u[D]);

  template<class T>
  void ErrorVelocity(
      const Mesh_Wrapper& mesh,
      user_field_t& formula_x, user_field_t& formula_y, user_field_t& formula_z,
      const T u[D], double* l2err, double* l2norm);

 private:
  int method_;
};


/* ******************************************************************
* Initialization of mass
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
    const Mesh_Wrapper& mesh, const T& mass, const T u[D]) const
{
  std::vector<double> m_loc(D, 0.0), m_glb(D);

  if (method_ == SGH) {
    std::vector<int> corners;

    for (int c = 0; c < mesh.num_owned_cells(); ++c) {
      mesh.cell_get_corners(c, &corners);

      for (auto cn : corners) { 
        int v = mesh.corner_get_node(cn);
        for (int i = 0; i < D; ++i) {
          m_loc[i] += mass[cn] * u[i][v];
        }
      }
    }
  }

  else if (method_ == CCH) {
    for (int c = 0; c < mesh.num_owned_cells(); ++c) {
      for (int i = 0; i < D; ++i) {
        m_loc[i] += mass[c] * u[i][c];
      }
    }
  }


  MPI_Allreduce(&(m_loc[0]), &(m_glb[0]), D, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Wonton::Point<D> momentum(m_glb);

  return momentum;
}


/* ******************************************************************
* Velocity bounds
****************************************************************** */
template<int D, class Mesh_Wrapper>
template<class T>
Wonton::Point<D> MomentumRemap<D, Mesh_Wrapper>::VelocityMin(
    const Mesh_Wrapper& mesh, const T u[D])
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  std::vector<double> umin_loc(D, 1e+99), umin_glb(D);
 
  for (int n = 1; n < nrows; ++n) {
    for (int i = 0; i < D; ++i) {
      umin_loc[i] = std::min(umin_loc[i], u[i][n]);
    }
  }

  MPI_Allreduce(&(umin_loc[0]), &(umin_glb[0]), D, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  Wonton::Point<D> umin(umin_glb);

  return umin;
}


template<int D, class Mesh_Wrapper>
template<class T>
Wonton::Point<D> MomentumRemap<D, Mesh_Wrapper>::VelocityMax(
    const Mesh_Wrapper& mesh, const T u[D])
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();
  std::vector<double> umax_loc(D, -1e+99), umax_glb(D);
 
  for (int n = 1; n < nrows; ++n) {
    for (int i = 0; i < D; ++i) {
      umax_loc[i] = std::max(umax_loc[i], u[i][n]);
    }
  }

  MPI_Allreduce(&(umax_loc[0]), &(umax_glb[0]), D, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  Wonton::Point<D> umax(umax_glb);

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
    const T u[D], double* l2err, double* l2norm)
{
  int nrows = (method_ == SGH) ? mesh.num_owned_nodes() : mesh.num_owned_cells();

  *l2err = 0.0;
  *l2norm = 0.0;

  for (int n = 0; n < nrows; ++n) {
    Wonton::Point<D> xyz {};

    if (method_ == SGH)
      mesh.node_get_coordinates(n, &xyz);
    else
      mesh.cell_centroid(n, &xyz);

    double const ux_exact = formula_x(xyz);
    double const uy_exact = formula_y(xyz);

    *l2err += (ux_exact - u[0][n]) * (ux_exact - u[0][n])
            + (uy_exact - u[1][n]) * (uy_exact - u[1][n]);

    *l2norm += ux_exact * ux_exact + uy_exact * uy_exact;

    if (D == 3) {
      double const uz_exact = formula_z(xyz);
      *l2err += (uz_exact - u[2][n]) * (uz_exact - u[2][n]);
      *l2norm += uz_exact * uz_exact;
    }
  }

  int nrows_glb = 0;
  double l2err_glb = 0.;
  double l2norm_glb = 0.;

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
  std::string velocity[3] = { "velocity_x", "velocity_y", "velocity_z" };

  auto kind = MassKind();
  const double *mass_src;
  srcstate_wrapper.mesh_get_data(kind, "mass", &mass_src);

  kind = VelocityKind();
  const double *u_src[D];
  for (int i = 0; i < D; ++i) {
    srcstate_wrapper.mesh_get_data(kind, velocity[i], &u_src[i]);
  }

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
  std::string momentum[3] = { "momentum_x", "momentum_y", "momentum_z" };
  std::vector<double> momentum_src[D];
  for (int i = 0; i < D; ++i) {
    momentum_src[i].resize(ncells_src, 0.0);
  }

  for (int c = 0; c < ncells_src; ++c) {
    if (method_ == SGH) {
      srcmesh_wrapper.cell_get_corners(c, &corners);

      for (auto cn : corners) { 
        int v = srcmesh_wrapper.corner_get_node(cn);
        for (int i = 0; i < D; ++i) {
          momentum_src[i][c] += mass_src[cn] * u_src[i][v];
        }
      }
    } else {
      for (int i = 0; i < D; ++i) {
        momentum_src[i][c] += mass_src[c] * u_src[i][c];
      }
    }

    double volume = srcmesh_wrapper.cell_volume(c);
    for (int i = 0; i < D; ++i) {
      momentum_src[i][c] /= volume;
    }
  }

  // Step 4 (SGH and CCH)
  // -- remap density and specific momentum following three basic
  //    steps: search, intersect, and interpolate
  Portage::CoreDriver<D, Wonton::Entity_kind::CELL, Mesh_Wrapper, State_Wrapper>
      cd(srcmesh_wrapper, srcstate_wrapper,
         trgmesh_wrapper, trgstate_wrapper);

  auto candidates = cd.template search<Portage::SearchKDTree>();
  auto srcwts = cd.template intersect_meshes<Portage::IntersectRnD>(candidates);

  // -- we need to register fields that we want to remap
  std::vector<std::string> field_names;
  std::vector<const double*> field_pointers;

  field_names.emplace_back("density");
  field_pointers.emplace_back(density.data());

  for (int i = 0; i < D; ++i) {
    field_names.emplace_back(momentum[i]);
    field_pointers.emplace_back(momentum_src[i].data());
  }

  int const num_fields = field_names.size();
  std::vector<double> tmp(ncells_trg);

  for (int i = 0; i < num_fields; ++i) {
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
  std::vector<Wonton::vector<Wonton::Vector<D>>> gradients(
      field_names.size(), Wonton::vector<Wonton::Vector<D>>(ncells_all));

  if (method_ == SGH) {
    using Gradient = Portage::Limited_Gradient<D, Wonton::CELL,
                                               Mesh_Wrapper, State_Wrapper>;
    Gradient gradient_kernel(trgmesh_wrapper, trgstate_wrapper);

    for (int i = 0; i < num_fields; ++i) {
      Wonton::vector<Wonton::Vector<D>> gradient(ncells_all);
      gradient_kernel.set_interpolation_variable(field_names[i], limiter, Portage::BND_NOLIMITER);
      Wonton::transform(trgmesh_wrapper.begin(Wonton::Entity_kind::CELL),
                         trgmesh_wrapper.end(Wonton::Entity_kind::CELL),
                         gradients[i].begin(), gradient_kernel);
    }
  }

  // Step 6 (SGH and CCH)
  // -- integrate density and specific momentum on the target mesh
  double *mass_trg;
  kind = MassKind();
  trgstate_wrapper.mesh_get_data(kind, "mass", &mass_trg);

  kind = VelocityKind();
  double *u_trg[D];
  for (int i = 0; i < D; ++i) {
    trgstate_wrapper.mesh_get_data(kind, velocity[i], &u_trg[i]);
  }

  //    extract auxiliary (cell-centerd) data
  double *density_trg;
  const double *momentum_trg[D];

  kind = Wonton::Entity_kind::CELL;
  trgstate_wrapper.mesh_get_data(kind, "density", &density_trg);
  for (int i = 0; i < D; ++i) {
    trgstate_wrapper.mesh_get_data(kind, momentum[i], &momentum_trg[i]);
  }

  //    integrate
  std::vector<double> momentum_cn[D]; // work memory

  if (method_ == SGH) {
    for (int i = 0; i < D; ++i) {
      momentum_cn[i].resize(ncorners_trg);
    }

    Wonton::Point<D> xc, xcn;
    Wonton::Vector<D> grad;
    std::vector<int> wedges;
    std::array<Wonton::Point<D>, D + 1> wcoords;

    for (int c = 0; c < ncells_trg; ++c) {
      trgmesh_wrapper.cell_centroid(c, &xc);
      trgmesh_wrapper.cell_get_corners(c, &corners);

      for (auto&& cn : corners) {
        trgmesh_wrapper.corner_get_wedges(cn, &wedges);

        corner_get_centroid(cn, trgmesh_wrapper, xcn);

        // integral is the value at centroid
        double cnvol = trgmesh_wrapper.corner_volume(cn);

        grad = gradients[0][c];
        mass_trg[cn] = cnvol * (density_trg[c] + dot(grad, xcn - xc));

        for (int i = 0; i < D; ++i) {
          grad = gradients[i + 1][c];
          momentum_cn[i][cn] = cnvol * (momentum_trg[i][c] + dot(grad, xcn - xc));
        }
      }
    }
  } 
  else {
    for (int c = 0; c < ncells_trg; ++c) {
      double vol = trgmesh_wrapper.cell_volume(c);
      mass_trg[c] = density_trg[c] * vol;
      for (int i = 0; i < D; ++i) {
        u_trg[i][c] = momentum_trg[i][c] / density_trg[c];
      }
    }
  }

  // Step 7 (SGH only)
  // -- gather data to mesh nodes and compute velocity
  if (method_ == SGH) {
    int nnodes_all = nnodes_trg + trgmesh_wrapper.num_ghost_nodes();
    std::vector<double> mass_v(nnodes_all, 0.0);  // work memory
    std::vector<double> momentum_v[D];
    for (int i = 0; i < D; ++i) {
      momentum_v[i].resize(nnodes_all, 0.0);
    }

    for (int cn = 0; cn < ncorners_trg; ++cn) {
      int v = trgmesh_wrapper.corner_get_node(cn);
      mass_v[v] += mass_trg[cn];
      for (int i = 0; i < D; ++i) {
        momentum_v[i][v] += momentum_cn[i][cn];
      }
    }

    for (int v = 0; v < nnodes_trg; ++v) {
      for (int i = 0; i < D; ++i) {
        u_trg[i][v] = momentum_v[i][v] / mass_v[v];
      }
    }
  }
}

#endif
