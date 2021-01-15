/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <fstream>
#include <sstream>
#include <vector>

// portage includes
#include "portage/support/portage.h"

// app includes
#include "user_field.h"
#include "MomentumRemap_mm.h"

// wonton includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/support/Point.h"

// Jali includes
#include "MeshFactory.hh"

/*!
  @file momentumapp2D_mm.cc
  @brief A 2D application that remaps velocity for CCH and two materials

  The program is used to showcase the capability of remapping velocity 
  and mass for staggered (not yet) and cell-centered hydro codes. Velocity
  remap conserves the total momentum. The app is controlled by a few input 
  parameters.

  Assumptions:
    - meshes occupy the same domain;
    - no limiters on the boundary.
*/


void print_usage() {
  std::cout << "Usage: ./momentumapp2D_mm nx ny limiter \"density formula\" \"velx formula\" \"vely formula\"\n\n"
            << "   source mesh:     nx x ny rectangular cells inside unit square\n" 
            << "   target mesh:     max(nx + 2, 1.1 * nx) x max(ny + 3, 1.2 * ny) rectangular cells\n\n"
            << "   limiter:         0 - limiter is off, otherwise Barth-Jespersen is used\n\n"
            << "   density formula: mathematical expression for density, a comma-separated list\n"
            << "   velx formula:    mathematical expression for x-component of velocity\n"
            << "   vely formula:    mathematical expression for y-component of velocity\n\n"
            << "Example: ./momentumapp2D_mm 10 11  1  \"1+x, 1+y\" \"x\" \"y\"\n";
}


void StringToStrings(int D, const std::string& in, std::vector<std::string>& out)
{
  size_t pos;
  std::string tmp(in);

  while((pos = tmp.find_first_of(',')) != std::string::npos) {
    out.push_back(tmp.substr(0, pos));
    tmp.erase(0, pos + 1);
  }

  // populating missing formulas for the one-material case
  for (int i = out.size(); i < D; ++i) out.push_back(tmp);
}


/* ******************************************************************
* Supporting routines
****************************************************************** */
template<typename Mesh_Wrapper>
double EvaluateCurl(Mesh_Wrapper const& mesh,
                    const std::vector<Wonton::vector<Wonton::Vector<2>>>& gradient)
{
  int ncells = mesh.num_owned_cells();

  double L2norm_loc(0.0), L2norm_glb;
  for (int c = 0; c < ncells; ++c) {
    const Wonton::Vector<2>& a = gradient[0][c];
    const Wonton::Vector<2>& b = gradient[1][c]; 
    L2norm_loc += std::pow(a[1] - b[0], 2) * mesh.cell_volume(c); 
  }
 
  MPI_Allreduce(&L2norm_loc, &L2norm_glb, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return std::sqrt(L2norm_glb);
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

  if (argc != 7) {
    if (rank == 0) print_usage();
    MPI_Finalize();
    return 0;
  }

  // number of cells in x and y directions for the source mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  auto limiter = (atoi(argv[3]) == 0) ? Portage::NOLIMITER : Portage::BARTH_JESPERSEN;
  auto formula_rho = std::string(argv[4]);
  auto formula_velx = std::string(argv[5]);
  auto formula_vely = std::string(argv[6]);

  std::vector<std::string> vformula_rho, vformula_velx, vformula_vely;
  StringToStrings(2, formula_rho, vformula_rho);
  StringToStrings(2, formula_velx, vformula_velx);
  StringToStrings(2, formula_vely, vformula_vely);

  int n = vformula_rho.size();
  std::vector<user_field_t> ini_rho(n), ini_velx(n), ini_vely(n);

  for (int i = 0; i < n; ++i) {
    if ((!ini_rho[i].initialize(2, vformula_rho[i])) ||
        (!ini_velx[i].initialize(2, vformula_velx[i])) ||
        (!ini_vely[i].initialize(2, vformula_vely[i]))) {
      if (rank == 0) {
        std::cout << "=== Input ERROR ===\n";
        print_usage();    
      }
      MPI_Finalize();
      return 0;
    }
  }

  // size of computational domain
  double lenx = 1.0;  // [m]
  double leny = 1.0;
  
  //
  // ======= Preliminaries, common for SGH and CCH =======
  //

  // -- setup Jali meshes
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);
  mesh_factory.included_entities(Jali::Entity_kind::ALL_KIND);

  int nxx = std::max(nx + 2, (int)(1.1 * nx));
  int nxy = std::max(ny + 3, (int)(1.2 * ny));
  auto srcmesh = mesh_factory(0.0, 0.0, lenx, leny, nx, ny);
  auto trgmesh = mesh_factory(0.0, 0.0, lenx, leny, nxx, nxy);

  // -- setup mesh wrappers
  Wonton::Jali_Mesh_Wrapper srcmesh_wrapper(*srcmesh);
  Wonton::Jali_Mesh_Wrapper trgmesh_wrapper(*trgmesh);

  // -- states
  std::shared_ptr<Jali::State> srcstate = Jali::State::create(srcmesh);
  std::shared_ptr<Jali::State> trgstate = Jali::State::create(trgmesh);

  // -- state wrappers
  Wonton::Jali_State_Wrapper srcstate_wrapper(*srcstate);
  Wonton::Jali_State_Wrapper trgstate_wrapper(*trgstate);

  // -- register materials with the states
  MomentumRemap_mm<2, Wonton::Jali_Mesh_Wrapper, Wonton::Jali_State_Wrapper> mr(Method::CCH);
  mr.InitMaterials(srcmesh_wrapper, srcstate_wrapper);

  int nmat_names = mr.mat_names().size();
  for (int m = 0; m < nmat_names; ++m)
    trgstate_wrapper.add_material(mr.mat_names()[m], {});
    
  trgstate_wrapper.mat_add_celldata<double>("mat_volfracs");
  trgstate_wrapper.mat_add_celldata<Wonton::Point<2> >("mat_centroids");

  // -- register mass and velocity with the states
  mr.SetMass(srcmesh_wrapper, srcstate_wrapper, "mass", ini_rho);
  mr.SetVelocity(srcmesh_wrapper, srcstate_wrapper, "velocity_x", ini_velx);
  mr.SetVelocity(srcmesh_wrapper, srcstate_wrapper, "velocity_y", ini_vely);

  // -- we re-use the code to allocate memory
  mr.SetMass(trgmesh_wrapper, trgstate_wrapper, "mass", ini_rho);
  mr.SetVelocity(trgmesh_wrapper, trgstate_wrapper, "velocity_x", ini_velx);
  mr.SetVelocity(trgmesh_wrapper, trgstate_wrapper, "velocity_y", ini_vely);

  Wonton::Point<2> xyz;
  
  // -- summary
  auto total_mass_src = mr.TotalMass(srcmesh_wrapper, srcstate_wrapper);
  auto total_momentum_src = mr.TotalMomentum(srcmesh_wrapper, srcstate_wrapper);
  auto umin = mr.VelocityMin(srcmesh_wrapper, srcstate_wrapper);
  auto umax = mr.VelocityMax(srcmesh_wrapper, srcstate_wrapper);
  if (rank == 0) {
    std::cout << "=== SOURCE data ===" << std::endl;
    std::cout << "mesh:           " << nx << " x " << ny << std::endl;
    std::cout << "total mass:     " << total_mass_src << " kg" << std::endl;
    std::cout << "total momentum: " << total_momentum_src << " kg m/s" << std::endl;
    std::cout << "limiter:        " << ((limiter == Portage::NOLIMITER) ? "none\n" : "BJ\n");
    std::cout << "velocity bounds," << " min: " << umin << " max: " << umax << std::endl;
  }

  //
  // ======= SEVEN-step REMAP algorithm =======
  //
  mr.RemapND(srcmesh_wrapper, srcstate_wrapper,
             trgmesh_wrapper, trgstate_wrapper,
             limiter);

  //
  // Verification 
  //
  auto total_mass_trg = mr.TotalMass(trgmesh_wrapper, trgstate_wrapper);
  auto total_momentum_trg = mr.TotalMomentum(trgmesh_wrapper, trgstate_wrapper);
  umin = mr.VelocityMin(trgmesh_wrapper, trgstate_wrapper);
  umax = mr.VelocityMax(trgmesh_wrapper, trgstate_wrapper);

  if (rank == 0) {
    std::cout << "\n=== TARGET data ===" << std::endl;
    std::cout << "mesh:           " << nxx << " x " << nxy << std::endl;
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

  double ul2err, ul2norm, dl2err, dl2norm;
  mr.ErrorVelocity(trgmesh_wrapper, trgstate_wrapper,
                   ini_velx, ini_vely, ini_vely,
                   &ul2err, &ul2norm);

  mr.ErrorDensity(trgmesh_wrapper, trgstate_wrapper,
                  ini_rho, &dl2err, &dl2norm);

  if (rank == 0) {
    std::cout << "\n=== Remap error ===" << std::endl;
    std::cout << "in velocity: l2-err=" << ul2err << " l2-norm=" << ul2norm << std::endl;
    std::cout << "in dencity:  l2-err=" << dl2err << " l2-norm=" << dl2norm << std::endl;
  }

  //
  // Additional analysis: angular momentum
  //
  double curl_src(0.0), curl_trg(0.0);
  std::vector<std::string> vel_names = { "velocity_x", "velocity_y" };

  using Gradient = Portage::Limited_Gradient<2, Wonton::CELL,
                                             Wonton::Jali_Mesh_Wrapper,
                                             Wonton::Jali_State_Wrapper>;

  {
    // gradient field computation is only supported on owned cells,
    // besides the curl operator uses the gradient on owned cells (line 75).
    int ncells_src = srcmesh_wrapper.num_owned_cells();
    std::vector<Wonton::vector<Wonton::Vector<2>>> grads_src(2, Portage::vector<Wonton::Vector<2>>(ncells_src));
    Gradient gradient_kernel(srcmesh_wrapper, srcstate_wrapper);

    for (int i = 0; i < 2; ++i) {
      gradient_kernel.set_interpolation_variable(vel_names[i], limiter, Portage::BND_NOLIMITER);
      Wonton::transform(srcmesh_wrapper.begin(Wonton::CELL, Wonton::PARALLEL_OWNED),
                        srcmesh_wrapper.end(Wonton::CELL, Wonton::PARALLEL_OWNED),
                        grads_src[i].begin(), gradient_kernel);
    }

    curl_src = EvaluateCurl(srcmesh_wrapper, grads_src);
  }

  {
    // gradient field computation is only supported on owned cells,
    // besides the curl operator uses the gradient on owned cells (line 75).
    int ncells_trg = trgmesh_wrapper.num_owned_cells();
    std::vector<Portage::vector<Wonton::Vector<2>>> grads_trg(2, Portage::vector<Wonton::Vector<2>>(ncells_trg));
    Gradient gradient_kernel(trgmesh_wrapper, trgstate_wrapper);

    for (int i = 0; i < 2; ++i) {
      gradient_kernel.set_interpolation_variable(vel_names[i], limiter, Portage::BND_NOLIMITER);
      Wonton::transform(trgmesh_wrapper.begin(Wonton::CELL, Wonton::PARALLEL_OWNED),
                        trgmesh_wrapper.end(Wonton::CELL, Wonton::PARALLEL_OWNED),
                        grads_trg[i].begin(), gradient_kernel);
    }

    curl_trg = EvaluateCurl(trgmesh_wrapper, grads_trg);
  }

  if (rank == 0) {
    std::cout << "\n=== Rotational component ===" << std::endl;
    std::cout << "inp velocity: L2-norm=" << curl_src << std::endl;
    std::cout << "out velocity: L2-norm=" << curl_trg << std::endl;
  }

  // save data
  if (rank == 0) {
    std::ofstream datafile;
    std::stringstream ss;
    ss << "errorsMM_0" << ".txt"; 
    datafile.open(ss.str());
    datafile << "0 " << cons_law0 << std::endl;
    datafile << "1 " << cons_law1 << std::endl;
    datafile << "2 " << ul2err << std::endl;
    datafile << "3 " << ul2norm << std::endl;
    datafile << "4 " << curl_src << std::endl;
    datafile << "5 " << curl_trg << std::endl;
    datafile.close();
  }

  MPI_Finalize();
  return 0;
}



