/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_MMDRIVER_H_
#define PORTAGE_MMDRIVER_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>

#ifdef HAVE_TANGRAM
#include "tangram/driver/driver.h"
#endif

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/search/search_kdtree.h"
#include "portage/intersect/intersect_r2d.h"
#include "portage/intersect/intersect_r3d.h"
#include "portage/interpolate/interpolate_1st_order.h"
#include "portage/interpolate/interpolate_2nd_order.h"
#include "portage/wonton/mesh/flat/flat_mesh_wrapper.h"
#include "portage/wonton/state/flat/flat_state_wrapper.h"

#ifdef ENABLE_MPI
#include "portage/distributed/mpi_bounding_boxes.h"
#endif

/*!
  @file driver.h
  @brief Example driver for mapping between two meshes.

  This should serve as a good example for how to write your own driver routine
  and datastructures.
*/

namespace Portage {

using namespace Wonton;


// Dummy interface reconstruction class, if external interface
// reconstruction methods are not found
template<class Mesh_Wrapper, int Dim> class DummyInterfaceReconstructor {
 public:
  DummyInterfaceReconstructor(Mesh_Wrapper const & mesh) {};
  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs,
                            std::vector<Tangram::Point<Dim>> const& cell_mat_centroids) {}; 
  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs) {};

  void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {}
  std::shared_ptr<Tangram::CellMatPoly<Dim>> operator()(const int cell_op_ID) const {}
};

/*!
  @class MMDriver "mmdriver.h"
  @brief MMDriver provides the API to mapping multi-material data from one mesh to another.

  @tparam Search  A search method that takes the dimension, source mesh class
  and target mesh class as template parameters

  @tparam Intersect  A polyhedron-polyhedron intersection class that takes
  the source and taget mesh classes as template parameters

  @tparam Interpolate An interpolation class that takes the source and
  target mesh classes, the source state class (that stores source
  field values), the kind of entity the interpolation is on and the
  dimension of the problem as template parameters

  @tparam SourceMesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.

  @tparam SourceState_Wrapper A lightweight wrapper to a specific input state
  manager implementation that provides certain functionality.

  @tparam TargetMesh_Wrapper A lightweight wrapper to a specific target mesh
  implementation that provides certain functionality.

  @tparam TargetState_Wrapper A lightweight wrapper to a specific target state
  manager implementation that provides certain functionality.

  @tparam InterfaceReconstructor An interface reconstruction class
  that takes the raw interface reconstruction method, the dimension of
  the problem and the source mesh class as template parameters
*/
template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class> class Intersect,
          template<int, Entity_kind, class, class, class> class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper = SourceMesh_Wrapper,
          class TargetState_Wrapper = SourceState_Wrapper,
          template <class, int> class InterfaceReconstructor = DummyInterfaceReconstructor>
class MMDriver {

  // Something like this would be very helpful to users
  // static_assert(
  //   D == Interpolate::D,
  //   "The dimension of Driver and Interpolate do not match!"
  // );


 public:
  /*!
    @brief Constructor for running the interpolation driver.
    @param[in] sourceMesh A @c SourceMesh_Wrapper to the source mesh.
    @param[in] sourceState A @c SourceState_Wrapperfor the data that lives on the
    source mesh.
    @param[in] targetMesh A @c TargetMesh_Wrapper to the target mesh.
    @param[in,out] targetState A @c TargetState_Wrapper for the data that will
    be mapped to the target mesh.
  */
  MMDriver(SourceMesh_Wrapper const& sourceMesh,
           SourceState_Wrapper const& sourceState,
           TargetMesh_Wrapper const& targetMesh,
           TargetState_Wrapper& targetState)
      : source_mesh_(sourceMesh), source_state_(sourceState),
        target_mesh_(targetMesh), target_state_(targetState),
        dim_(sourceMesh.space_dimension()) {
    assert(sourceMesh.space_dimension() == targetMesh.space_dimension());
  }

  /// Copy constructor (disabled)
  MMDriver(const MMDriver &) = delete;

  /// Assignment operator (disabled)
  MMDriver & operator = (const MMDriver &) = delete;

  /// Destructor
  ~MMDriver() {}

  /*!
    @brief Specify the names of the variables to be interpolated along with the
    limiter to use for all of them
    @param[in] remap_var_names A list of variable names of the variables to
    interpolate from the source mesh to the target mesh.  This variable must
    exist in both meshes' state manager
    @param[in] limiter_type The type of limiter to use (NOLIMITER, BARTH_JESPERSEN)
  */
  void set_remap_var_names(std::vector<std::string> const &remap_var_names,
                           LimiterType limiter_type = NOLIMITER) {
    // All variables use the same type of limiters
    std::vector<LimiterType> limiters(remap_var_names.size(), limiter_type);

    // remap variable names same in source and target mesh
    set_remap_var_names(remap_var_names, remap_var_names, limiters);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_remap_var_names A list of the variables names of the
    variables to interpolate from the source mesh.
    @param[in] target_remap_var_names  A list of the variables names of the
    variables to interpolate to the target mesh.
    @param[in] limiter_type The limiter to use for higher order remaps (NOLIMITER, BARTH_JESPERSEN)
  */
  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names,
      LimiterType limiter_type = NOLIMITER) {
    std::vector<LimiterType> limiters(source_remap_var_names.size(),
                                      limiter_type);

    set_remap_var_names(source_remap_var_names, target_remap_var_names,
                        limiters);
  }

  /*!
    @brief Specify the names of the variables to be interpolated
    @param[in] source_remap_var_names A list of the variables names of the
    variables to interpolate from the source mesh.
    @param[in] target_remap_var_names  A list of the variables names of the
    variables to interpolate to the target mesh.
    @param[in] limiter_types Limiters to use for each remapped variable (NOLIMITER, BARTH_JESPERSEN)
  */

  void set_remap_var_names(
      std::vector<std::string> const & source_remap_var_names,
      std::vector<std::string> const & target_remap_var_names,
      std::vector<LimiterType> const & limiter_types) {
    assert(source_remap_var_names.size() == target_remap_var_names.size());

    int nvars = source_remap_var_names.size();
    for (int i = 0; i < nvars; ++i)
      assert(source_state_.get_entity(source_remap_var_names[i]) ==
             target_state_.get_entity(target_remap_var_names[i]));

    source_remap_var_names_ = source_remap_var_names;
    target_remap_var_names_ = target_remap_var_names;
    limiters_ = limiter_types;
  }

  /*!
    @brief Get the names of the variables to be remapped from the
    source mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> source_remap_var_names() const {
    return source_remap_var_names_;
  }

  /*!
    @brief Get the names of the variables to be remapped to the
    target mesh.
    @return A vector of variable names to be remapped.
  */
  std::vector<std::string> target_remap_var_names() const {
    return target_remap_var_names_;
  }

  /*!
    @brief Get the dimensionality of the meshes.
    @return The dimensionality of the meshes.
  */
  unsigned int dim() const {
    return dim_;
  }


  /*!
    @brief remap for a given set of variables on a given entity kind
    @tparam entity_kind  Kind of entity that variables live on
    @param source_remap_var_names  names of remap variables on source mesh (MUST ALL BE DEFINED ON THE SAME Entity_kind - NODE, CELL)
    @param target_remap_var_names  names of remap variables on target mesh (MUST ALL BE DEFINED ON THE SAME Entity_kind - NODE, CELL)
  */

  template<Entity_kind onwhat>
  void remap(std::vector<std::string> const &source_var_names,
             std::vector<std::string> const &target_var_names);

#ifdef ENABLE_MPI
  /*!
    @brief remap for a given set of variables on a given entity kind in a distributed setting with redistribution of data if needed
    @tparam entity_kind  Kind of entity that variables live on
    @param source_remap_var_names  names of remap variables on source mesh (MUST ALL BE DEFINED ON THE SAME Entity_kind - NODE, CELL)
    @param target_remap_var_names  names of remap variables on target mesh (MUST ALL BE DEFINED ON THE SAME Entity_kind - NODE, CELL)
  */

  template<Entity_kind onwhat>
  void remap_distributed(std::vector<std::string> const &source_var_names,
                         std::vector<std::string> const &target_var_names);
#endif



  /*!
    @brief Execute the remapping process
  */
  void run(bool distributed) {

#ifndef ENABLE_MPI
    if (distributed) {
      std::cout << "Request is for a parallel run but Portage is compiled " <<
          "for serial runs only\n";
      return;
    }
#endif

    int comm_rank = 0;
#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

    if (comm_rank == 0)
      std::printf("in MMDriver::run()...\n");

    int numTargetCells = target_mesh_.num_owned_cells();
    std::cout << "Number of target cells in target mesh on rank "
              << comm_rank << ": "
              << numTargetCells << std::endl;


    int nvars = source_remap_var_names_.size();

    // Collect all cell based variables and remap them

    std::vector<std::string> source_cellvar_names;
    std::vector<std::string> target_cellvar_names;
    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat = source_state_.get_entity(source_remap_var_names_[i]);
      if (onwhat == CELL) {
        source_cellvar_names.emplace_back(source_remap_var_names_[i]);
        target_cellvar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    if (source_cellvar_names.size()) {
#ifdef ENABLE_MPI
      if (distributed)
        remap_distributed<CELL>(source_cellvar_names, target_cellvar_names);
      else
#endif
        remap<CELL>(source_cellvar_names, target_cellvar_names);
    }


    // Collect all node based variables and remap them

    std::vector<std::string> source_nodevar_names;
    std::vector<std::string> target_nodevar_names;

    for (int i = 0; i < nvars; ++i) {
      Entity_kind onwhat = source_state_.get_entity(source_remap_var_names_[i]);
      if (onwhat == NODE) {
        source_nodevar_names.emplace_back(source_remap_var_names_[i]);
        target_nodevar_names.emplace_back(target_remap_var_names_[i]);
      }
    }

    if (source_nodevar_names.size()) {
#ifdef ENABLE_MPI
      if (distributed)
        remap_distributed<NODE>(source_nodevar_names, target_nodevar_names);
      else
#endif
        remap<NODE>(source_nodevar_names, target_nodevar_names);
    }
  }  // run

 private:
  SourceMesh_Wrapper const& source_mesh_;
  TargetMesh_Wrapper const& target_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetState_Wrapper& target_state_;
  std::vector<std::string> source_remap_var_names_;
  std::vector<std::string> target_remap_var_names_;
  std::vector<LimiterType> limiters_;
  unsigned int dim_;
};  // class MMDriver



// Serial remap or Distributed remap with no redistributon of data

template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class> class Intersect,
          template<int, Entity_kind, class, class, class> class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper,
          class TargetState_Wrapper,
          template<class, int> class InterfaceReconstructor>
template<Entity_kind onwhat>
void MMDriver<Search, Intersect, Interpolate, D,
              SourceMesh_Wrapper, SourceState_Wrapper,
              TargetMesh_Wrapper, TargetState_Wrapper,
              InterfaceReconstructor
            >::remap(std::vector<std::string> const &src_varnames,
                     std::vector<std::string> const &trg_varnames) {

  static_assert(onwhat == NODE || onwhat == CELL,
                "Remap implemented only for CELL and NODE variables");

  int comm_rank = 0;
#ifdef ENABLE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

  int ntarget_ents_owned = target_mesh_.num_entities(onwhat, PARALLEL_OWNED);
  std::cout << "Number of target entities of kind " << onwhat <<
      " in target mesh on rank " << comm_rank << ": " <<
      ntarget_ents_owned << std::endl;

  int ntarget_ents = target_mesh_.num_entities(onwhat, ALL);

  float tot_seconds = 0.0, tot_seconds_srch = 0.0,
      tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;


  // SEARCH

  Portage::vector<std::vector<int>> candidates(ntarget_ents);
  Portage::vector<std::vector<Weights_t>> source_ents_and_weights(ntarget_ents);

  // Get an instance of the desired search algorithm type
  gettimeofday(&begin_timeval, 0);
  const Search<D, onwhat, SourceMesh_Wrapper, TargetMesh_Wrapper>
      search(source_mesh_, target_mesh_);

  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(), search);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

#ifdef HAVE_TANGRAM
  // Call interface reconstruction only if we got a method from the
  // calling app
  
  if (typeid(InterfaceReconstructor<SourceMesh_Wrapper, D>) !=
      typeid(DummyInterfaceReconstructor<SourceMesh_Wrapper, D>)) {
    // INTERFACE RECONSTRUCTION (BUT RIGHT NOW WE ARE NOT USING IT)
    
    Tangram::Driver<InterfaceReconstructor, D, SourceMesh_Wrapper>
        interface_reconstructor(source_mesh_);
    
    // Since we don't really have multiple materials and volume
    // fractions of those materials in cells, we will pretend that we
    // have two materials in each cell and that the volume fraction of
    // each material in each cell is 0.5
    
    int nmats = 2;
    int nsourcecells = source_mesh_.num_entities(CELL, ALL);
    std::vector<int> cell_num_mats(nsourcecells, nmats);
    std::vector<int> cell_mat_ids(nsourcecells*nmats);
    std::vector<double> cell_mat_volfracs(nsourcecells*nmats, 0.5);
    std::vector<Tangram::Point<D>> cell_mat_centroids(nsourcecells*nmats);
    interface_reconstructor.set_volume_fractions(cell_num_mats,
                                                 cell_mat_ids,
                                                 cell_mat_volfracs,
                                                 cell_mat_centroids);
    interface_reconstructor.reconstruct();
  }
#endif

  // INTERSECT

  gettimeofday(&begin_timeval, 0);

  // Get an instance of the desired intersect algorithm type
  const Intersect<onwhat, SourceMesh_Wrapper, TargetMesh_Wrapper>
      intersect(source_mesh_, target_mesh_);

  // For each cell in the target mesh get a list of candidate-weight
  // pairings (in a traditional mesh, not particle mesh, the weights
  // are moments). Note that this candidate list is different from the
  // search candidate list in that it may not include some of the
  // search candidates. Also, note that for 2nd order and higher
  // remaps, we get multiple moments (0th, 1st, etc) for each
  // target-source cell intersection

  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(),
                     source_ents_and_weights.begin(),
                     intersect);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_xsect = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  // INTERPOLATE (one variable at a time)

  gettimeofday(&begin_timeval, 0);

  int nvars = src_varnames.size();
  if (comm_rank == 0)
    std::cout << "Number of variables on entity kind " << onwhat <<
        " to remap is " << nvars << std::endl;

  // Get an instance of the desired interpolate algorithm type
  Interpolate<D, onwhat, SourceMesh_Wrapper, TargetMesh_Wrapper,
              SourceState_Wrapper>
      interpolate(source_mesh_, target_mesh_, source_state_);

  for (int i = 0; i < nvars; ++i) {
    interpolate.set_interpolation_variable(src_varnames[i],
                                           limiters_[i]);

    double *target_field_raw = nullptr;
    target_state_.get_data(onwhat, trg_varnames[i], &target_field_raw);
    Portage::pointer<double> target_field(target_field_raw);

    Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                       target_mesh_.end(onwhat, PARALLEL_OWNED),
                       source_ents_and_weights.begin(),
                       target_field, interpolate);
  }

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

  std::cout << "Transform Time for Entity Kind " << onwhat << " on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;
  std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;
}



#ifdef ENABLE_MPI

// Distributed Remap with redistribution of mesh and data

template <template <int, Entity_kind, class, class> class Search,
          template <Entity_kind, class, class> class Intersect,
          template<int, Entity_kind, class, class, class> class Interpolate,
          int D,
          class SourceMesh_Wrapper,
          class SourceState_Wrapper,
          class TargetMesh_Wrapper,
          class TargetState_Wrapper,
          template <class, int> class InterfaceReconstructor>
template<Entity_kind onwhat>
void MMDriver<Search, Intersect, Interpolate, D,
              SourceMesh_Wrapper, SourceState_Wrapper,
              TargetMesh_Wrapper, TargetState_Wrapper,
              InterfaceReconstructor
              >::remap_distributed(std::vector<std::string> const &src_varnames,
                                   std::vector<std::string> const &trg_varnames) {

  static_assert(onwhat == NODE || onwhat == CELL,
                "Remap implemented only for CELL and NODE variables");

  int comm_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  int ntarget_ents_owned = target_mesh_.num_entities(onwhat, PARALLEL_OWNED);
  std::cout << "Number of target entities of kind " << onwhat <<
      " in target mesh on rank " << comm_rank << ": " <<
      ntarget_ents_owned << std::endl;

  int ntarget_ents = target_mesh_.num_entities(onwhat, ALL);

  Flat_Mesh_Wrapper<> source_mesh_flat;
  Flat_State_Wrapper<> source_state_flat;

  float tot_seconds = 0.0, tot_seconds_srch = 0.0,
      tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;


  // SEARCH

  Portage::vector<std::vector<int>> candidates(ntarget_ents);
  Portage::vector<std::vector<Weights_t>> source_ents_and_weights(ntarget_ents);

  // Create flat wrappers to distribute source cells
  gettimeofday(&begin_timeval, 0);

  source_mesh_flat.initialize(source_mesh_);
  source_state_flat.initialize(source_state_, source_remap_var_names_);
  MPI_Bounding_Boxes distributor;
  distributor.distribute(source_mesh_flat, source_state_flat,
                         target_mesh_, target_state_);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  float tot_seconds_flat = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;
  std::cout << "Redistribution Time Rank " << comm_rank << " (s): " <<
      tot_seconds_flat << std::endl;

  // Get an instance of the desired search algorithm type
  gettimeofday(&begin_timeval, 0);
  const Search<D, onwhat, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>
      search(source_mesh_flat, target_mesh_);

  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(), search);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_srch = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

#ifdef HAVE_TANGRAM
  /* WE ARE UNABLE TO CALL TANGRAM DRIVER FROM DISTRIBUTED REMAP
   * BECAUSE TANGRAM DRIVER CALLS get_global_id WHICH DOES NOT EXIST
   * FOR FLAT_MESH_WRAPPER. WE CAN UNCOMMENT THIS WHEN TICKET LNK-781
   * IS RESOLVED */

  // // Call interface reconstruction only if we got a method from the
  // // calling app

  // if (typeid(InterfaceReconstructor<Flat_Mesh_Wrapper<>, D>) !=
  //     typeid(DummyInterfaceReconstructor<Flat_Mesh_Wrapper<>, D>)) {
  //   // INTERFACE RECONSTRUCTION (BUT RIGHT NOW WE ARE NOT USING IT)
    
  //   Tangram::Driver<InterfaceReconstructor, D, Flat_Mesh_Wrapper<>>
  //       interface_reconstructor(source_mesh_flat);

  //   // Since we don't really have multiple materials and volume
  //   // fractions of those materials in cells, we will pretend that we
  //   // have two materials in each cell and that the volume fraction of
  //   // each material in each cell is 0.5
    
  //   int nmats = 2;
  //   int nsourcecells = source_mesh_.num_entities(CELL, ALL);
  //   std::vector<int> cell_num_mats(nsourcecells, nmats);
  //   std::vector<int> cell_mat_ids(nsourcecells*nmats);
  //   std::vector<double> cell_mat_volfracs(nsourcecells*nmats, 0.5);
  //   std::vector<Tangram::Point<D>> cell_mat_centroids(nsourcecells*nmats);
  //   interface_reconstructor.set_volume_fractions(cell_num_mats,
  //                                                cell_mat_ids,
  //                                                cell_mat_volfracs,
  //                                                cell_mat_centroids);
  //   interface_reconstructor.reconstruct();
  // }
#endif

  // INTERSECT

  gettimeofday(&begin_timeval, 0);

  // Get an instance of the desired intersect algorithm type
  const Intersect<onwhat, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper>
      intersect(source_mesh_flat, target_mesh_);

  // For each cell in the target mesh get a list of candidate-weight
  // pairings (in a traditional mesh, not particle mesh, the weights
  // are moments). Note that this candidate list is different from the
  // search candidate list in that it may not include some of the
  // search candidates. Also, note that for 2nd order and higher
  // remaps, we get multiple moments (0th, 1st, etc) for each
  // target-source cell intersection

  Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                     target_mesh_.end(onwhat, PARALLEL_OWNED),
                     candidates.begin(),
                     source_ents_and_weights.begin(),
                     intersect);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_xsect = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  // INTERPOLATE (one variable at a time)

  gettimeofday(&begin_timeval, 0);

  int nvars = src_varnames.size();
  if (comm_rank == 0)
    std::cout << "Number of variables on entity kind " << onwhat <<
        " to remap is " << nvars << std::endl;

  // Get an instance of the desired interpolate algorithm type
  Interpolate<D, onwhat, Flat_Mesh_Wrapper<>, TargetMesh_Wrapper,
              Flat_State_Wrapper<>>
      interpolate(source_mesh_flat, target_mesh_, source_state_flat);

  for (int i = 0; i < nvars; ++i) {
    interpolate.set_interpolation_variable(src_varnames[i],
                                           limiters_[i]);

    double *target_field_raw = nullptr;
    target_state_.get_data(onwhat, trg_varnames[i], &target_field_raw);
    Portage::pointer<double> target_field(target_field_raw);

    Portage::transform(target_mesh_.begin(onwhat, PARALLEL_OWNED),
                       target_mesh_.end(onwhat, PARALLEL_OWNED),
                       source_ents_and_weights.begin(),
                       target_field, interpolate);
  }

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  tot_seconds_interp = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  tot_seconds = tot_seconds_srch + tot_seconds_xsect + tot_seconds_interp;

  std::cout << "Transform Time for Entity Kind " << onwhat << " on Rank " <<
      comm_rank << " (s): " << tot_seconds << std::endl;
  std::cout << "   Search Time Rank " << comm_rank << " (s): " <<
      tot_seconds_srch << std::endl;
  std::cout << "   Intersect Time Rank " << comm_rank << " (s): " <<
      tot_seconds_xsect << std::endl;
  std::cout << "   Interpolate Time Rank " << comm_rank << " (s): " <<
      tot_seconds_interp << std::endl;
}
#endif  // ENABLE_MPI

}  // namespace Portage

#endif  // PORTAGE_DRIVER_H_
