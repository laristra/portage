/*
This is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_DRIVER_FIX_MISMATCH_H_
#define PORTAGE_DRIVER_FIX_MISMATCH_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>
#include <limits>

#include "portage/support/portage.h"

/*!
  @file mismatch_fixup.h
  @brief Repair step for remap when the boundaries of the source and target mesh are not exactly coincident
*/

namespace Portage {

using Wonton::Entity_kind;
using Wonton::Entity_type;
using Wonton::Weights_t;

// Check if we have mismatch of mesh boundaries
// T is the target mesh, S is the source mesh
// T_i is the i'th cell of the target mesh and
// |*| is signifies the extent/volume of an entity)
//
// If sum_i(|T_i intersect S|) NE |S|, some source cells are not
// completely covered by the target mesh
//
// If sum_i(|T_i intersect S|) NE |T|, some target cells are not
// completely covered by the source mesh
//
// If there is a mismatch, adjust values for fields to account so that
// integral quantities are conserved and a constant field is readjusted
// to a different constant

template<int D, Entity_kind onwhat,
         class SourceMesh_Wrapper, class SourceState_Wrapper,
         class TargetMesh_Wrapper, class TargetState_Wrapper>
bool fix_mismatch(SourceMesh_Wrapper const & source_mesh,
                  SourceState_Wrapper const & source_state,
                  TargetMesh_Wrapper const & target_mesh,
                  TargetState_Wrapper & target_state,
                  Portage::vector<std::vector<Weights_t>> const & source_ents_and_weights,
                  std::vector<std::string> const & src_var_names,
                  std::vector<std::string> const & trg_var_names) {

  int rank=0;
#ifdef ENABLE_MPI
  int nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif
  

  bool mismatch = false;


  // collect source and target cell volumes because we will need it
  // a few times for each variable we process

  int nsourceents = (onwhat == Entity_kind::CELL) ?
      source_mesh.num_owned_cells() : source_mesh.num_owned_nodes();

  // We will use a mask to avoid computing with source entities that occur on
  // other ranks, thereby avoiding if statements
  
  std::vector<double> source_mask(nsourceents, 1.0);

#ifdef ENABLE_MPI
  if (nprocs > 1) {
    std::vector<int> nsourceents_all(nprocs, 0);
    MPI_Allgather(&nsourceents, 1, MPI_INT, &(nsourceents_all[0]), 1, MPI_INT,
                  MPI_COMM_WORLD);
    int maxsourceents = *std::max_element(nsourceents_all.begin(),
                                          nsourceents_all.end());
  
    std::vector<int> source_gids(maxsourceents, -1);
    for (int s = 0; s < nsourceents; s++)
      source_gids[s] = source_mesh.get_global_id(s, onwhat);
    
    std::vector<int> source_gids_all(nprocs*maxsourceents);
    MPI_Allgather(&(source_gids[0]), maxsourceents, MPI_INT,
                  &(source_gids_all[0]), maxsourceents, MPI_INT,
                  MPI_COMM_WORLD);

    // Add source gids from all lower ranks into a set (no duplicates)
    
    std::unordered_set<int> unique_source_gids;
    for (int p = 0; p < rank; p++)
      for (int s = 0; s < nsourceents_all[p]; s++) {
        int gid = source_gids_all[p*maxsourceents+s];
        unique_source_gids.insert(gid);
      }
    
    // Now go through gids of this rank and mask the ones that are already
    // occurring in other ranks.
    for (int s = 0; s < nsourceents; s++) {
      auto itpair = unique_source_gids.insert(source_gids[s]);
      if (!itpair.second)
        source_mask[s] = 0.0;  // ent already in set; mask this instance
    }
  }
#endif

  // collect volumes of entities that are not masked out and sum them up
  
  std::vector<double> source_cell_volume(nsourceents, 0.0);
  for (int s = 0; s < nsourceents; s++)
    if (source_mask[s])
      source_cell_volume[s] = (onwhat == Entity_kind::CELL) ?
          source_mesh.cell_volume(s) : source_mesh.dual_cell_volume(s);

  double source_volume = std::accumulate(source_cell_volume.begin(),
                                         source_cell_volume.end(), 0.0);
#ifdef ENABLE_MPI
  double global_source_volume = 0.0;
  MPI_Allreduce(&source_volume, &global_source_volume, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#else
  double global_source_volume = source_volume;
#endif


  int ntargetents = (onwhat == Entity_kind::CELL) ?
      target_mesh.num_owned_cells() : target_mesh.num_owned_nodes();

  std::vector<double> target_cell_volume(ntargetents);
  for (int t = 0; t < ntargetents; t++)
    target_cell_volume[t] = (onwhat == Entity_kind::CELL) ?
        target_mesh.cell_volume(t) : target_mesh.dual_cell_volume(t);

  double target_volume = std::accumulate(target_cell_volume.begin(),
                                         target_cell_volume.end(), 0.0);
#ifdef ENABLE_MPI
  double global_target_volume = 0.0;
  MPI_Allreduce(&target_volume, &global_target_volume, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
#else
  double global_target_volume = target_volume;
#endif



  // Sum up intersection volumes and along the way let's detect if any
  // target mesh entities are not at all covered by source mesh
  // entities. In our initial redistribution phase, we will move as
  // many source cells as needed to from different partitions to cover
  // the target partition, UNLESS NO SOURCE CELLS COVERING THE TARGET
  // CELLS EXIST GLOBALLY. So, at this stage, we can just check
  // on-rank intersections only

  std::vector<double> xsect_volumes(ntargetents);
  std::vector<bool> is_empty(ntargetents, false);
  std::vector<int> emptyents;
  for (int t = 0; t < ntargetents; t++) {
    std::vector<Weights_t> const& sw_vec = source_ents_and_weights[t];

    xsect_volumes[t] = 0.0;
    for (auto const& sw : sw_vec)
      xsect_volumes[t] += sw.weights[0];

    if (fabs(xsect_volumes[t]) < 1.0e-16) {
      emptyents.push_back(t);
      is_empty[t] = true;
    }
  }
  int nempty = emptyents.size();

  double xsect_volume = std::accumulate(xsect_volumes.begin(),
                                        xsect_volumes.end(), 0.0);
#ifdef ENABLE_MPI
  double global_xsect_volume = 0.0;
  MPI_Allreduce(&xsect_volume, &global_xsect_volume, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  double global_xsect_volume = xsect_volume;
#endif




  // Are some source cells not fully covered by target cells?
  
  if (fabs((global_xsect_volume-global_source_volume)/global_source_volume)
      > 1.0e-12) {
    mismatch = true;
    std::cerr << "\n** MESH MISMATCH -" <<
        " some source cells are not fully covered by the target mesh\n";

#ifdef DEBUG
    // Find one source cell (or dual cell) that is not fully covered
    // by the target mesh and output its ID. Unfortunately, that means
    // processing all source cells We initialize each source cell to
    // its volume and subtract any intersection volume we find between
    // it and a target cell

    // Also, it will likely give a false positive in distributed
    // meshes because a source cell may be covered by target cells
    // from multiple processors

    std::vector<double> source_covered_vol(source_cell_volume);
    for (auto it = target_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
         it != target_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++) {
      auto const& sw_vec = source_ents_and_weights[*it];
      for (auto const& sw : sw_vec)
        source_covered_vol[sw.entityID] -= sw.weights[0];
    }

    for (auto it = source_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
         it != source_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++)
      if (source_covered_vol[*it] > 1.0e-12) {
        if (onwhat == Entity_kind::CELL)
          std::cerr << "Source cell " << *it << " not fully covered by target cells \n";
        else
          std::cerr << "Source dual cell " << *it << " not fully covered by target dual cells \n";
        break;
      }
#endif

    std::cerr << "\n";
  }


  // Are some target cells not fully coverd by source cells?
  
  if (fabs((global_xsect_volume-global_target_volume)/global_target_volume)
      > 1.0e-12) {
    mismatch = true;
    std::cerr << "\n** MESH MISMATCH -" <<
        " some target cells are not fully covered by the source mesh\n";

#ifdef DEBUG
    // Find one target cell that is not fully covered by the source mesh and
    // output its ID
    for (auto it = target_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
         it != target_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++) {
      int t = *it;
      double covered_vol = 0.0;
      std::vector<Weights_t> const& sw_vec = source_ents_and_weights[t];
      for (auto const& sw : sw_vec)
        covered_vol += sw.weights[0];
      if (fabs(covered_vol-target_cell_volume[t])/target_cell_volume[t] > 1.0e-12) {
        if (onwhat == Entity_kind::CELL)
          std::cerr << "Target cell " << *it << " not fully covered by source cells \n";
        else
          std::cerr << "Target dual cell " << *it << " not fully covered by source dual cells \n";
        break;
      }
    }
#endif

    std::cerr << "\n";
  }



  if (!mismatch) return true;
  
  

  // Collect the empty target cells in layers starting from the ones
  // next to partially or fully covered cells. At the end of this
  // section, partially or fully covered cells will all have a layer
  // number of 0 and empty cell layers will have positive layer
  // numbers starting from 1

  std::vector<int> layernum(ntargetents, 0);
  std::vector<std::vector<int>> emptylayers;
  if (nempty) {
    if (onwhat == Entity_kind::CELL)
      std::cerr << "One or more target cells are not covered by " <<
          "any source cells.\n" <<
          " Will assign values based on their neighborhood\n";
    else
      std::cerr << "One or more target dual cells are not covered by " <<
          "any source dual cells.\n" <<
          " Will assign values based on their neighborhood\n";
     
    int nlayers = 0;
    int ntagged = 0;
    while (ntagged < nempty) {
      std::vector<int> curlayerents;

      for (int ent : emptyents) {
        if (layernum[ent] != 0) continue;

        std::vector<int> nbrs;
        if (onwhat == Entity_kind::CELL)
          target_mesh.cell_get_node_adj_cells(ent, Entity_type::ALL, &nbrs);
        else
          target_mesh.node_get_cell_adj_nodes(ent, Entity_type::ALL, &nbrs);
        for (int nbr : nbrs)
          if (!is_empty[nbr] || layernum[nbr] != 0) {
            // At least one neighbor has some material or will
            // receive some material (indicated by having a +ve
            // layer number)

            curlayerents.push_back(ent);
            break;
          }
      }  // for (ent in emptyents)

      // Tag the current layer cells with the next layer number
      for (int ent : curlayerents)
        layernum[ent] = nlayers+1;
      ntagged += curlayerents.size();

      emptylayers.push_back(curlayerents);
      nlayers++;
    }
  }  // if nempty



  // Now process remap variables

  int nvars = src_var_names.size();

  for (int i = 0; i < nvars; i++) {
    std::string const& src_var_name = src_var_names[i];
    double const *source_data;
    source_state.mesh_get_data(onwhat, src_var_names[i], &source_data);

    // for now get lower and upper bounds from source mesh itself
    // later we may need to get them from the calling application
    double lower_bound = *std::min_element(source_data, source_data+nsourceents);
    double upper_bound = *std::max_element(source_data, source_data+nsourceents);
    
    std::string const& trg_var_name = trg_var_names[i];
    double *target_data;
    target_state.mesh_get_data(onwhat, trg_var_names[i], &target_data);

    // Do something here to populate fully uncovered target cells. We
    // have layers of empty cells starting out from fully or partially
    // populated cells. We will assign every empty cell in a layer the
    // average value of all its populated neighbors.  IN A DISTRIBUTED
    // MESH IT _IS_ POSSIBLE THAT AN EMPTY ENTITY WILL NOT NOT HAVE ANY
    // OWNED NEIGHBOR IN THIS PARTITION THAT HAS MEANINGFUL DATA TO
    // EXTRAPOLATE FROM (we remap data only to owned entities)

    int curlayernum = 1;
    for (std::vector<int> const& curlayer : emptylayers) {
      for (int ent : curlayer) {
        std::vector<int> nbrs;
        if (onwhat == Entity_kind::CELL)
          target_mesh.cell_get_node_adj_cells(ent, Entity_type::PARALLEL_OWNED,
                                              &nbrs);
        else
          target_mesh.node_get_cell_adj_nodes(ent, Entity_type::PARALLEL_OWNED,
                                              &nbrs);

        double aveval = 0.0;
        int nave = 0;
        for (int nbr : nbrs)
          if (layernum[nbr] < curlayernum) {
            aveval += target_data[nbr];
            nave++;
          }
        if (nave)
          aveval /= nave;
#ifdef DEBUG
        else
          std::cerr <<
              "No owned neighbors of empty entity to extrapolate data from\n";
#endif

        target_data[ent] = aveval;
      }
      curlayernum++;
    }


    // At this point assume that all cells have some value in them
    // for each of the variables

    // Now compute the net discrepancy between integrals over source
    // and target. Excess comes from target cells are not fully covered
    // by source cells and deficit from source cells not fully covered
    // by target cells

    double source_sum =
        std::inner_product(source_data, source_data + nsourceents,
                           source_cell_volume.begin(), 0.0);

    double global_source_sum = source_sum;
#ifdef ENABLE_MPI
    MPI_Allreduce(&source_sum, &global_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    double target_sum =
        std::inner_product(target_data, target_data + ntargetents,
                           target_cell_volume.begin(), 0.0);

    double global_target_sum = target_sum;
#ifdef ENABLE_MPI
    MPI_Allreduce(&target_sum, &global_target_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
#endif

    double global_diff = global_target_sum - global_source_sum;
    double reldiff = global_diff/global_source_sum;

    if (fabs(reldiff) < 1.0e-14)
      return true;  // discrepancy is too small - nothing to do

    // Now redistribute the discrepancy among cells in proportion to
    // their volume. This will restore conservation and if the
    // original distribution was a constant it will make the field a
    // slightly different constant. Do multiple iterations because
    // we may have some leftover quantity that we were not able to
    // add/subtract from some cells in any given iteration. We don't
    // expect that process to take more than two iterations at the
    // most.

    double adj_target_volume = target_volume;
    double global_adj_target_volume = global_target_volume;

    // sort of a "unit" discrepancy or "specific" discrepancy
    double udiff = global_diff/global_adj_target_volume;

    int iter = 0;
    while (fabs(reldiff) > 1.0e-14 && iter < 5) {
      for (auto it = target_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
           it != target_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++) {
        int t = *it;

        if ((target_data[t]-udiff) < lower_bound) {
          // Subtracting the full excess will make this cell violate the
          // lower bound. So subtract only as much as will put this cell
          // exactly at the lower bound

          target_data[t] = lower_bound;

          std::cerr << "Hit lower bound for cell " << t << " on rank " << rank
                    << "\n";

          // this cell is no longer in play for adjustment - so remove its
          // volume from the adjusted target_volume

          adj_target_volume -= target_cell_volume[t];

        } else if ((target_data[t]-udiff) > upper_bound) {  // udiff < 0
          // Adding the full deficit will make this cell violate the
          // upper bound. So add only as much as will put this cell
          // exactly at the upper bound

          target_data[t] = upper_bound;

          std::cerr << "Hit upper bound for cell " << t << " on rank " << rank
                    << "\n";

          // this cell is no longer in play for adjustment - so remove its
          // volume from the adjusted target_volume

          adj_target_volume -= target_cell_volume[t];

        } else {
          // This is the equivalent of
          //           [curval*cellvol - diff*cellvol/meshvol]
          // curval = ---------------------------------------
          //                       cellvol

          target_data[t] -= udiff;

        }
      }  // iterate through mesh cells

      // Compute the new integral over all processors
        
      target_sum = std::inner_product(target_data, target_data + ntargetents,
                                      target_cell_volume.begin(), 0.0);

      global_target_sum = target_sum;
#ifdef ENABLE_MPI
      MPI_Allreduce(&target_sum, &global_target_sum, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      // If we did not hit lower or upper bounds, this should be
      // zero after the first iteration.  If we did hit some bounds,
      // then recalculate the discrepancy and discrepancy per unit
      // volume, but only taking into account volume of cells that
      // are not already at the bounds - if we use the entire mesh
      // volume for the recalculation, the convergence slows down

      global_diff = global_target_sum - global_source_sum;

      global_adj_target_volume = adj_target_volume;
#ifdef ENABLE_MPI
      MPI_Allreduce(&adj_target_volume, &global_adj_target_volume, 1,
                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

      udiff = global_diff/global_adj_target_volume;
      reldiff = global_diff/global_source_sum;


      // Now reset adjusted target mesh volume to be the full volume
      // in preparation for the next iteration
      global_adj_target_volume = global_target_volume;

      iter++;
    }  // while leftover is not zero
    
    if (fabs(reldiff) > 1.0e-14) {
      if (rank == 0) {
        std::cerr << "Redistribution not entirely successfully for variable " <<
            src_var_name << "\n";
        std::cerr << "Relative conservation error is " << reldiff << "\n";
        std::cerr << "Absolute conservation error is " << global_diff << "\n";
        return false;
      }
    }

  }  // for each variable

}  // fix_mismatch

}  // namespace Portage

#endif  // PORTAGE_DRIVER_FIX_MISMATCH_H_
