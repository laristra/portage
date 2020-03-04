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
#include "MomentumRemap.h"

// wonton includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"
#include "wonton/support/Point.h"

// Jali includes
#include "MeshFactory.hh"

/*!
  @file momentumapp3D.cc
  @brief A 3D application that remaps velocity and mass for SGH and CCH

  The program is used to showcase the capability of remapping velocity 
  and mass for staggered and cell-centered hydro codes. Velocity remap
  conserves the total momentum. The app is controlled by a few input 
  commands. Unnecessary longer code is used for the implementation clarity.

  SGH:
    A. Populate input data: corner masses and nodal velocities
    B. Conservative remap of momentum
    C. Verification of output data: corner masses and nodal 
       velocities on the target mesh

  CCH: 
    A. Populate input data: cell-centered masses and velocities
    B. Conservative remap of momentum
    C. Verification of output data: cell-centered massed and velocities

  Assumptions: meshes are matching.
*/


void print_usage() {
  std::cout << "Usage: ./momentumapp3D nx ny nz method limiter \"rho(x,y,z)\" \"vx(x,y,z)\" \"vy(x,y,z)\"  \"vz(x,y,z)\"\n\n"
            << "   source mesh:  nx x ny x nz bricks inside the unit cube\n" 
            << "   target mesh:  (nx + 2) x (ny + 3) x (nz + 4) bricks\n\n"
            << "   method:       SGH=1, CCH=2\n"
            << "   limiter:      0 - limiter is off, otherwise Barth-Jespersen is used\n\n"
            << "   rho(x,y,z):   mathematical expression for density\n"
            << "   vx(x,y,z):    mathematical expression for x-component of velocity\n"
            << "   vy(x,y,z):    mathematical expression for y-component of velocity\n\n"
            << "   vz(x,y,z):    mathematical expression for z-component of velocity\n\n"
            << "Example: ./momentumapp3D 10 10 10 2 1  \"1+x\" \"x\" \"x*z\" \"if((x < 0.5), 1+x, 2+z)\"\n";
}


/* ******************************************************************
* Utility that computes corner centroid
****************************************************************** */
template<>
void corner_get_centroid(
    int cn, const Wonton::Jali_Mesh_Wrapper& mesh,
    Wonton::Point<3>& xcn)
{
  std::vector<int> wedges;
  std::array<Wonton::Point<3>, 4> wcoords;

  double volume = mesh.corner_volume(cn); 
  mesh.corner_get_wedges(cn, &wedges);

  xcn = {0.0, 0.0, 0.0};
  for (auto w : wedges) {
    double frac = mesh.wedge_volume(w) / volume; 
    mesh.wedge_get_coordinates(w, &wcoords);
    for (int i = 0; i < 4; ++i) {
      xcn += frac * wcoords[i] / 4;
    }
  }
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

  if (argc < 10) {
    if (rank == 0) print_usage();
    MPI_Finalize();
    return 0;
  }

  // number of cells in x and y directions for the source mesh
  auto nx = atoi(argv[1]);
  auto ny = atoi(argv[2]);
  auto nz = atoi(argv[3]);
  auto method = atoi(argv[4]);
  auto limiter = (atoi(argv[5]) == 0) ? Portage::NOLIMITER : Portage::BARTH_JESPERSEN;
  auto formula_rho = std::string(argv[6]);
  auto formula_velx = std::string(argv[7]);
  auto formula_vely = std::string(argv[8]);
  auto formula_velz = std::string(argv[9]);

  user_field_t ini_rho, ini_velx, ini_vely, ini_velz;

  if ((method != SGH && method != CCH) ||
      (!ini_rho.initialize(3, formula_rho)) ||
      (!ini_velx.initialize(3, formula_velx)) ||
      (!ini_vely.initialize(3, formula_vely)) ||
      (!ini_velz.initialize(3, formula_velz))) {
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
  double lenz = 1.0;

  //
  // Preliminaries, common for SGH and CCH
  //

  // -- setup Jali meshes
  Jali::MeshFactory mesh_factory(MPI_COMM_WORLD);
  mesh_factory.included_entities(Jali::Entity_kind::ALL_KIND);

  auto srcmesh = mesh_factory(0.0, 0.0, 0.0, lenx, leny, lenz, nx, ny, nz);
  auto trgmesh = mesh_factory(0.0, 0.0, 0.0, lenx, leny, lenz, nx + 2, ny + 3, nz + 4);
  // auto trgmesh = mesh_factory("test/prism40.exo");

  // -- setup mesh wrappers
  Wonton::Jali_Mesh_Wrapper srcmesh_wrapper(*srcmesh);
  Wonton::Jali_Mesh_Wrapper trgmesh_wrapper(*trgmesh);

  // -- states
  std::shared_ptr<Jali::State> srcstate = Jali::State::create(srcmesh);
  std::shared_ptr<Jali::State> trgstate = Jali::State::create(trgmesh);

  // -- state wrappers
  Wonton::Jali_State_Wrapper srcstate_wrapper(*srcstate);
  Wonton::Jali_State_Wrapper trgstate_wrapper(*trgstate);

  // -- register velocity components with the states
  //    target state does not need data, but code re-use does the task easier
  MomentumRemap<3, Wonton::Jali_Mesh_Wrapper> mr(method);

  std::vector<double> u_src[3], tmp;

  auto kind = mr.VelocityKind();
  mr.InitVelocity(srcmesh_wrapper, ini_velx, u_src[0]);
  mr.InitVelocity(srcmesh_wrapper, ini_vely, u_src[1]);
  mr.InitVelocity(srcmesh_wrapper, ini_velz, u_src[2]);

  srcstate_wrapper.mesh_add_data(kind, "velocity_x", &(u_src[0][0]));
  srcstate_wrapper.mesh_add_data(kind, "velocity_y", &(u_src[1][0]));
  srcstate_wrapper.mesh_add_data(kind, "velocity_z", &(u_src[2][0]));

  mr.InitVelocity(trgmesh_wrapper, ini_velx, tmp);
  trgstate_wrapper.mesh_add_data(kind, "velocity_x", &(tmp[0]));
  trgstate_wrapper.mesh_add_data(kind, "velocity_y", &(tmp[0]));
  trgstate_wrapper.mesh_add_data(kind, "velocity_z", &(tmp[0]));

  // -- register mass with the states
  std::vector<double> mass_src;

  kind = mr.MassKind();
  mr.InitMass(srcmesh_wrapper, ini_rho, mass_src);

  srcstate_wrapper.mesh_add_data(kind, "mass", &(mass_src[0]));
  mr.InitMass(trgmesh_wrapper, ini_rho, tmp);
  trgstate_wrapper.mesh_add_data(kind, "mass", &(tmp[0]));

  // -- summary
  auto total_mass_src = mr.TotalMass(srcmesh_wrapper, mass_src);
  auto total_momentum_src = mr.TotalMomentum(srcmesh_wrapper, mass_src, u_src);
  auto umin = mr.VelocityMin(srcmesh_wrapper, u_src);
  auto umax = mr.VelocityMax(srcmesh_wrapper, u_src);
  if (rank == 0) {
    std::cout << "=== SOURCE data ===" << std::endl;
    std::cout << "mesh:           " << nx << " x " << ny << " x " << nz << std::endl;
    std::cout << "total mass:     " << total_mass_src << " kg" << std::endl;
    std::cout << "total momentum: " << total_momentum_src << " kg m/s" << std::endl;
    std::cout << "limiter:        " << ((limiter == Portage::NOLIMITER) ? "none\n" : "BJ\n");
    std::cout << "velocity bounds," << " min: " << umin << " max: " << umax << std::endl;
  }

  //
  // SEVEN-step REMAP algorithm 
  //
  mr.RemapND<Wonton::Jali_State_Wrapper>(srcmesh_wrapper, srcstate_wrapper,
                                         trgmesh_wrapper, trgstate_wrapper,
                                         limiter);

  //
  // Verification 
  //
  const double *mass_trg;
  const double *u_trg[3];

  kind = mr.MassKind();
  trgstate_wrapper.mesh_get_data(kind, "mass", &mass_trg);

  kind = mr.VelocityKind();
  trgstate_wrapper.mesh_get_data(kind, "velocity_x", &u_trg[0]);
  trgstate_wrapper.mesh_get_data(kind, "velocity_y", &u_trg[1]);
  trgstate_wrapper.mesh_get_data(kind, "velocity_z", &u_trg[2]);

  // use 2D/3D routines with dummy parameters 
  auto total_mass_trg = mr.TotalMass(trgmesh_wrapper, mass_trg);
  auto total_momentum_trg = mr.TotalMomentum(trgmesh_wrapper, mass_trg, u_trg);
  umin = mr.VelocityMin(trgmesh_wrapper, u_trg);
  umax = mr.VelocityMax(trgmesh_wrapper, u_trg);

  if (rank == 0) {
    std::cout << "\n=== TARGET data ===" << std::endl;
    std::cout << "mesh:           " << nx + 2 << " x " << ny + 3 << " x " << nz + 4 << std::endl;
    std::cout << "total mass:     " << total_mass_trg << " kg" << std::endl;
    std::cout << "total momentum: " << total_momentum_trg << " kg m/s" << std::endl;
    std::cout << "velocity bounds," << " min: " << umin << " max: " << umax << std::endl;
  }

  auto err = total_momentum_trg - total_momentum_src;
  double cons_law0 = std::fabs(total_mass_trg - total_mass_src);
  double cons_law1 = std::sqrt(err[0] * err[0] + err[1] * err[1] + err[2] * err[2]);

  if (rank == 0) {
    std::cout << "\n=== Conservation error ===" << std::endl;
    std::cout << "in total mass:     " << cons_law0 << std::endl;
    std::cout << "in total momentum: " << cons_law1 << std::endl;
  }

  double l2err, l2norm;
  mr.ErrorVelocity(trgmesh_wrapper, ini_velx, ini_vely, ini_velz,
                   u_trg, &l2err, &l2norm);

  if (rank == 0) {
    std::cout << "\n=== Remap error ===" << std::endl;
    std::cout << "in velocity: l2-err=" << l2err << " l2-norm=" << l2norm << std::endl;
  }

  // save data
  if (rank == 0) {
    std::ofstream datafile;
    std::stringstream ss;
    ss << "errors3D_" << method - 1 << ".txt"; 
    datafile.open(ss.str());
    datafile << "0 " << cons_law0 << std::endl;
    datafile << "1 " << cons_law1 << std::endl;
    datafile << "2 " << l2err << std::endl;
    datafile << "3 " << l2norm << std::endl;
    datafile.close();
  }

  /*
  srcstate->export_to_mesh();
  srcmesh->write_to_exodus_file("input.exo");
  trgstate->export_to_mesh();
  trgmesh->write_to_exodus_file("output.exo");
  */

  MPI_Finalize();
  return 0;
}
