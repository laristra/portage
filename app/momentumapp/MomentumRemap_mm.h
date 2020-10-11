/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef MOMENTUM_REMAP_ND_HH_
#define MOMENTUM_REMAP_ND_HH_

#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

// portage includes
#include "portage/driver/mmdriver.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/interpolate/gradient.h"
#include "portage/intersect/intersect_rNd.h"
#include "portage/intersect/simple_intersect_for_tests.h"
#include "portage/search/search_kdtree.h"
#include "portage/support/portage.h"

// wonton includes
#include "wonton/support/Point.h"

// tangram includes
#include "tangram/utility/get_material_moments.h"
#include "tangram/reconstruct/MOF.h"

#ifdef HAVE_XMOF2D
#include "tangram/reconstruct/xmof2D_wrapper.h"
#define IR_2D XMOF2D_Wrapper
#else
#define IR_2D MOF
#endif

#include "MomentumRemapDefs.h"

/* ******************************************************************
* App class that handles initialization and verification of fields.
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
class MomentumRemap_mm {
 public:
  MomentumRemap_mm(const Method& method)
    : method_(method),
      mass_name_("mass"),
      velocity_name_({ "velocity_x", "velocity_y", "velocity_z" }) {};

  void InitMaterials(const Mesh_Wrapper& mesh, State_Wrapper& state);

  // initialization using formula for density
  void SetMass(
      const Mesh_Wrapper& mesh, State_Wrapper& state,
      const std::string& name, std::vector<user_field_t>& formula);

  void SetVelocity(
      const Mesh_Wrapper& mesh, State_Wrapper& state,
      const std::string& name, std::vector<user_field_t>& formula);

  // field type
  Wonton::Entity_kind MassKind() const {
    return (method_ == Method::SGH) ? Wonton::Entity_kind::CORNER : Wonton::Entity_kind::CELL;
  }
  Wonton::Entity_kind VelocityKind() const {
    return (method_ == Method::SGH) ? Wonton::Entity_kind::NODE : Wonton::Entity_kind::CELL;
  }

  // main remap method
  void RemapND(
      const Mesh_Wrapper& srcmesh_wrapper, State_Wrapper& srcstate_wrapper,
      const Mesh_Wrapper& trgmesh_wrapper, State_Wrapper& trgstate_wrapper,
      Portage::Limiter_type limiter);

  // access
  const std::vector<std::string>& mat_names() const { return mat_names_; }

  // V&V
  double TotalMass(const Mesh_Wrapper& mesh, const State_Wrapper& state) const;

  Wonton::Point<D> TotalMomentum(
      const Mesh_Wrapper& mesh, const State_Wrapper& state) const;

  Wonton::Point<D> VelocityMin(
      const Mesh_Wrapper& mesh, const State_Wrapper& state);

  Wonton::Point<D> VelocityMax(
      const Mesh_Wrapper& mesh, const State_Wrapper& state);

  void ErrorVelocity(
      const Mesh_Wrapper& mesh, const State_Wrapper& state,
      std::vector<user_field_t>& formula_x, 
      std::vector<user_field_t>& formula_y, std::vector<user_field_t>& formula_z,
      double* l2err, double* l2norm);

  void ErrorDensity(
      const Mesh_Wrapper& mesh, const State_Wrapper& state,
      std::vector<user_field_t>& formula_rho, 
      double* l2err, double* l2norm);

 private:
  Method method_;
  std::vector<std::string> mat_names_;

  std::string mass_name_;
  std::vector<std::string> velocity_name_;
};


/* ******************************************************************
* Initialization of materials
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
void MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::InitMaterials(
   const Mesh_Wrapper& mesh, State_Wrapper& state)
{
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();

  // create cell-centered multi-material data
  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;  // flattened 2D array
  std::vector<double> cell_mat_volfracs;  // flattened 2D array
  std::vector<Wonton::Point<D> > cell_mat_centroids;  // flattened 2D array

  std::vector<int> material_IDs({1, 0});
  Wonton::Point<2> circle_cen(1.0, 1.0);
  double circle_rad = 0.5;
  int nquadrant_samples = 90;

  double dst_tol = 100 * std::numeric_limits<double>::epsilon();
  double vol_tol = 100 * std::numeric_limits<double>::epsilon();

  get_material_moments<Mesh_Wrapper>(mesh, material_IDs, 
    circle_cen, circle_rad, nquadrant_samples,
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids, 
    vol_tol, dst_tol, false);  

  // Count the number of local materials and gather their IDs
  std::set<int> mat_ids;
  for (auto id : cell_mat_ids) mat_ids.insert(id);
  std::vector<int> local_material_IDs(mat_ids.begin(), mat_ids.end());
 
  // Compute offsets into flattened arrays based on cell_num_mats
  std::vector<int> offsets(ncells, 0);
  for (int i = 0; i < ncells - 1; i++)
    offsets[i + 1] = offsets[i] + cell_num_mats[i];

  // Convert data from cell-centric to material-centric form as we
  // will need it for adding it to the state manager
  std::unordered_map<int, std::vector<int>> matcells;
  std::unordered_map<int, std::vector<double>> mat_volfracs;
  std::unordered_map<int, std::vector<Wonton::Point<D> > > mat_centroids;

  for (int c = 0; c < ncells; c++) {
    int ibeg = offsets[c];
    for (int j = 0; j < cell_num_mats[c]; j++) {
      int mat_id = cell_mat_ids[ibeg + j];
      matcells[mat_id].push_back(c);
      mat_volfracs[mat_id].push_back(cell_mat_volfracs[ibeg + j]);
      mat_centroids[mat_id].push_back(cell_mat_centroids[ibeg + j]);
    }
  }

  // Add materials, volume fractions, and centroids to the source state
  // If the material isn't found in the partition, add empty data to the source state
  // Note names are dimensioned on all materials whether in the partition or not
  int nmats = material_IDs.size();
  mat_names_.resize(nmats);

  for (int m = 0; m < nmats; m++) {
    mat_names_[m] = "mat" + std::to_string(m);
    if (matcells.find(m) != matcells.end()) {
      state.add_material(mat_names_[m], matcells[m]);
      state.mat_add_celldata("mat_volfracs", m, mat_volfracs[m].data());
      state.mat_add_celldata("mat_centroids", m, mat_centroids[m].data());
    }
    else {
      state.add_material(mat_names_[m], {});
      state.mat_add_celldata("mat_volfracs", m, static_cast<double*>(nullptr));
      state.mat_add_celldata("mat_centroids", m, static_cast<Wonton::Point<D>* >(nullptr));
    }
  }
}


/* ******************************************************************
* Initialization of mass
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
void MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::SetMass(
    const Mesh_Wrapper& mesh, State_Wrapper& state,
    const std::string& name, std::vector<user_field_t>& formula)
{
  Wonton::Point<D> xyz;

  int nmat_names = mat_names_.size();
  for (int m = 0; m < nmat_names; ++m) {
    std::vector<int> matcells;
    state.mat_get_cells(m, &matcells);
    int nmatcells = matcells.size();

    double const *matvf;
    state.mat_get_celldata("mat_volfracs", m, &matvf);

    Wonton::Point<D> const *matcen;
    state.mat_get_celldata("mat_centroids", m, &matcen);

    std::vector<double> mass(nmatcells);

    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcells[ic];
      double vol = mesh.cell_volume(c);
      mass[ic] = formula[m](matcen[ic]) * vol * matvf[ic];
    }

    state.mat_add_celldata(name, m, &(mass[0]));
  }
}


/* ******************************************************************
* Initiaization of cell-centered velocity
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
void MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::SetVelocity(
    const Mesh_Wrapper& mesh, State_Wrapper& state,
    const std::string& name, std::vector<user_field_t>& formula)
{
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  std::vector<double> u(ncells);

  Wonton::Point<D> xyz;

  for (int c = 0; c < ncells; ++c) {
    mesh.cell_centroid(c, &xyz);
    u[c] = formula[0](xyz);
  }

  state.mesh_add_data(Wonton::Entity_kind::CELL, name, &(u[0]));
}


/* ******************************************************************
* Calculate total momentum 
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
double MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::TotalMass(
    const Mesh_Wrapper& mesh, const State_Wrapper& state) const
{
  double sum(0.0), sum_glb;

  int nmat_names = mat_names_.size();
  for (int m = 0; m < nmat_names; ++m) {
    std::vector<int> matcells;
    state.mat_get_cells(m, &matcells);
    int nmatcells = matcells.size();

    double const *mass;
    state.mat_get_celldata(mass_name_, m, &mass);

    for (int ic = 0; ic < nmatcells; ic++) sum += mass[ic];
  }

  // global communications
  MPI_Allreduce(&sum, &sum_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return sum_glb;
}


/* ******************************************************************
* Calculate total momentum 
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
Wonton::Point<D> MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::TotalMomentum(
    const Mesh_Wrapper& mesh, const State_Wrapper& state) const
{
  std::vector<double> m_loc(D, 0.0), m_glb(D);

  double const *vel[D];
  for (int i = 0; i < D; ++i) {
    state.mesh_get_data(Wonton::Entity_kind::CELL, velocity_name_[i], &vel[i]);
  }

  int nmat_names = mat_names_.size();
  for (int m = 0; m < nmat_names; ++m) {
    std::vector<int> matcells;
    state.mat_get_cells(m, &matcells);
    int nmatcells = matcells.size();

    double const *mass;
    state.mat_get_celldata(mass_name_, m, &mass);

    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcells[ic];
      for (int i = 0; i < D; ++i) {
        m_loc[i] += mass[ic] * vel[i][c];
      }
    }
  }

  // global communications
  MPI_Allreduce(&(m_loc[0]), &(m_glb[0]), D, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  Wonton::Point<D> momentum(m_glb);

  return momentum;
}


/* ******************************************************************
* Velocity bounds
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
Wonton::Point<D> MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::VelocityMin(
    const Mesh_Wrapper& mesh, const State_Wrapper& state)
{
  int ncells = mesh.num_owned_cells();
  std::vector<double> umin_loc(D, 1e+99), umin_glb(D);
 
  for (int i = 0; i < D; ++i) {
    double const *vel;
    state.mesh_get_data(Wonton::Entity_kind::CELL, velocity_name_[i], &vel);

    for (int c = 0; c < ncells; ++c) {
      umin_loc[i] = std::min(umin_loc[i], vel[c]);
    }
  }

  MPI_Allreduce(&(umin_loc[0]), &(umin_glb[0]), D, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  Wonton::Point<D> umin(umin_glb);

  return umin;
}


template<int D, class Mesh_Wrapper, class State_Wrapper>
Wonton::Point<D> MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::VelocityMax(
    const Mesh_Wrapper& mesh, const State_Wrapper& state)
{
  int ncells = mesh.num_owned_cells();
  std::vector<double> umax_loc(D, -1e+99), umax_glb(D);
 
  for (int i = 0; i < D; ++i) {
    double const *vel;
    state.mesh_get_data(Wonton::Entity_kind::CELL, velocity_name_[i], &vel);

    for (int c = 0; c < ncells; ++c) {
      umax_loc[i] = std::max(umax_loc[i], vel[c]);
    }
  }

  MPI_Allreduce(&(umax_loc[0]), &(umax_glb[0]), D, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  Wonton::Point<D> umax(umax_glb);

  return umax;
}


/* ******************************************************************
* Error in cell-centered velocity
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
void MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::ErrorVelocity(
    const Mesh_Wrapper& mesh, const State_Wrapper& state,
    std::vector<user_field_t>& formula_x,
    std::vector<user_field_t>& formula_y, std::vector<user_field_t>& formula_z,
    double* l2err, double* l2norm)
{
  int ncells = mesh.num_owned_cells();
  double ux_exact, uy_exact, uz_exact;
  Wonton::Point<D> xyz;

  *l2err = 0.0;
  *l2norm = 0.0;

  double const *vel[D];
  for (int i = 0; i < D; ++i) {
    state.mesh_get_data(Wonton::Entity_kind::CELL, velocity_name_[i], &vel[i]);
  }

  for (int c = 0; c < ncells; ++c) {
    mesh.cell_centroid(c, &xyz);

    ux_exact = formula_x[0](xyz);
    uy_exact = formula_y[0](xyz);

    *l2err += (ux_exact - vel[0][c]) * (ux_exact - vel[0][c])
            + (uy_exact - vel[1][c]) * (uy_exact - vel[1][c]);

    *l2norm += ux_exact * ux_exact + uy_exact * uy_exact;

    if (D == 3) {
      uz_exact = formula_z[0](xyz);
      *l2err += (uz_exact - vel[2][c]) * (uz_exact - vel[2][c]);
      *l2norm += uz_exact * uz_exact;
    }
  }

  int ncells_glb;
  double l2err_glb, l2norm_glb;

  MPI_Allreduce(&ncells, &ncells_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l2err, &l2err_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l2norm, &l2norm_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *l2err = std::sqrt(l2err_glb / ncells_glb);
  *l2norm = std::sqrt(l2norm_glb / ncells_glb);
}


/* ******************************************************************
* Error in material-centered density 
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
void MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::ErrorDensity(
    const Mesh_Wrapper& mesh, const State_Wrapper& state,
    std::vector<user_field_t>& formula_rho,
    double* l2err, double* l2norm)
{
  int ndofs = 0;
  *l2err = 0.0;
  *l2norm = 0.0;

  int nmat_names = mat_names_.size();
  for (int m = 0; m < nmat_names; ++m) {
    std::vector<int> matcells;
    state.mat_get_cells(m, &matcells);
    int nmatcells = matcells.size();
    ndofs += nmatcells;

    double const *density;
    state.mat_get_celldata("density", m, &density);

    Wonton::Point<D> const *matcen;
    state.mat_get_celldata("mat_centroids", m, &matcen);

    for (int ic = 0; ic < nmatcells; ic++) {
      double den_exact = formula_rho[m](matcen[ic]);
      *l2err += (den_exact - density[ic]) * (den_exact - density[ic]);
      *l2norm += den_exact * den_exact;
    }
  }

  int ndofs_glb;
  double l2err_glb, l2norm_glb;

  MPI_Allreduce(&ndofs, &ndofs_glb, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l2err, &l2err_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l2norm, &l2norm_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  *l2err = std::sqrt(l2err_glb / ndofs_glb);
  *l2norm = std::sqrt(l2norm_glb / ndofs_glb);
}


/* ******************************************************************
* 2D and 3D remap algorithm
****************************************************************** */
template<int D, class Mesh_Wrapper, class State_Wrapper>
inline
void MomentumRemap_mm<D, Mesh_Wrapper, State_Wrapper>::RemapND(
    const Mesh_Wrapper& srcmesh, State_Wrapper& srcstate,
    const Mesh_Wrapper& trgmesh, State_Wrapper& trgstate,
    Portage::Limiter_type limiter)
{
  // mesh data
  int ncells_trg = trgmesh.num_owned_cells();

  // state fields
  // -- materials and mass
  int nmats = mat_names_.size();
  std::vector<int> matcells_src[nmats];
  const double *matvf_src[nmats], *mass_src[nmats];

  for (int m = 0; m < nmats; ++m) {
    srcstate.mat_get_cells(m, &matcells_src[m]);
    srcstate.mat_get_celldata(mass_name_, m, &mass_src[m]);
    srcstate.mat_get_celldata("mat_volfracs", m, &matvf_src[m]);
  }

  // -- velocity
  auto kind = VelocityKind();
  const double *u_src[D];
  for (int i = 0; i < D; ++i) {
    srcstate.mesh_get_data(kind, velocity_name_[i], &u_src[i]);
  }

  // === Step 1 (SGH only)
  // -- gather cell-centered mass from corner masses
  //    will be implemented later

  // === Step 2 (SGH and CCH)
  // -- compute material-centered density and save it to the state
  for (int m = 0; m < nmats; ++m) {
    int nmatcells = matcells_src[m].size();
    std::vector<double> density(nmatcells);

    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcells_src[m][ic];
      double vol = srcmesh.cell_volume(c);
      density[ic] = mass_src[m][ic] / (vol * matvf_src[m][ic]);
    }
    srcstate.mat_add_celldata("density", m, &(density[0]));
  }

  // === Step 3 (SGH and CCH)
  // -- compute material-centered specific momentum and save it to the state
  //    SGH will be implemented later
  std::vector<std::string> field_names({ "density", "momentum_x", "momentum_y"});
  if (D == 3) field_names.push_back("momentum_z");

  int nfield_names = field_names.size();
  for (int i = 0; i < nfield_names; ++i) {
    trgstate.template mat_add_celldata<double>(field_names[i]);
  }

  int nmat_names = mat_names_.size();
  for (int m = 0; m < nmat_names; ++m) {
    int nmatcells = matcells_src[m].size();
    std::vector<double> momentum(nmatcells);

    const double* density;
    srcstate.mat_get_celldata("density", m, &density);

    for (int i = 0; i < D; ++i) {
      for (int ic = 0; ic < nmatcells; ic++) {
        int c = matcells_src[m][ic];
        momentum[ic] = density[ic] * u_src[i][c];
      }
      srcstate.mat_add_celldata(field_names[i + 1], m, &(momentum[0]));
    }
  }

  // === Step 4 (SGH and CCH)
  // -- remap density and specific momentum
  Portage::MMDriver<Portage::SearchKDTree, Portage::IntersectR2D,
                    Portage::Interpolate_2ndOrder, 2,
                    Mesh_Wrapper, State_Wrapper,
                    Mesh_Wrapper, State_Wrapper,
                    Tangram::IR_2D,
                    Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>
      mmd(srcmesh, srcstate,
          trgmesh, trgstate);

  // -- register fields that we want to remap
  mmd.set_remap_var_names(field_names);

  // -- select limiters
  mmd.set_limiter(limiter);
  mmd.set_bnd_limiter(Portage::Boundary_Limiter_type::BND_NOLIMITER);

  // -- specify tolerances
  double dst_tol = 100 * std::numeric_limits<double>::epsilon();
  double vol_tol = std::numeric_limits<double>::epsilon();

  std::vector<Tangram::IterativeMethodTolerances_t> ims_tols(2);
  ims_tols[0] = {1000, dst_tol, vol_tol};
  ims_tols[1] = {100, 1.0e-15, 1.0e-15};

  mmd.set_reconstructor_options(true, ims_tols);

  // -- execute
  int numpe;
  MPI_Comm_size(MPI_COMM_WORLD, &numpe);
  Wonton::MPIExecutor_type mpiexecutor(MPI_COMM_WORLD);
  Wonton::Executor_type *executor = (numpe > 1) ? &mpiexecutor : nullptr;

  mmd.run(executor);

  // -- extract multi-material data
  std::vector<int> matcells_trg[nmats];
  const double *matvf_trg[nmats];

  for (int m = 0; m < nmats; ++m) {
    trgstate.mat_get_cells(m, &matcells_trg[m]);
    trgstate.mat_get_celldata("mat_volfracs", m, &matvf_trg[m]);
  }

  // === Step 5 (SGH only)
  // -- create linear reconstruction (limited or unlimited) 
  //    of density and specific momentum on the target mesh
  //    SGH will be implemented later

  // === Step 6 (SGH and CCH)
  // -- calculate material-centered mass and consistent velocity 
  //    on the target mesh
  std::vector<double> mass_trg(ncells_trg, 0.0);

  kind = VelocityKind();
  double *u_trg[D];
  for (int i = 0; i < D; ++i) {
    trgstate.mesh_get_data(kind, velocity_name_[i], &u_trg[i]);
  }

  for (int c = 0; c < ncells_trg; ++c) {
    for (int i = 0; i < D; ++i) u_trg[i][c] = 0.0;
  }

  for (int m = 0; m < nmat_names; ++m) {
    const double *density, *momentum[D];

    trgstate.mat_get_celldata("density", m, &density);
    for (int i = 0; i < D; ++i) {
      trgstate.mat_get_celldata(field_names[i + 1], m, &momentum[i]);
    }

    int nmatcells = matcells_trg[m].size();
    std::vector<double> mass(nmatcells);

    for (int ic = 0; ic < nmatcells; ic++) {
      int c = matcells_trg[m][ic];
      double vol = trgmesh.cell_volume(c);

      mass[ic] = density[ic] * vol * matvf_trg[m][ic];
      mass_trg[c] += mass[ic];
      for (int i = 0; i < D; ++i) {
        u_trg[i][c] += momentum[i][ic] * vol * matvf_trg[m][ic];
      }
    }

    trgstate.mat_add_celldata(mass_name_, m, &(mass[0]));
  }

  for (int c = 0; c < ncells_trg; ++c) {
    for (int i = 0; i < D; ++i) {
      u_trg[i][c] /= mass_trg[c];
    }
  }

  // Step 7 (SGH only)
  // -- gather data to mesh nodes and compute velocity
  //    SGH will be implemented later
}

#endif
