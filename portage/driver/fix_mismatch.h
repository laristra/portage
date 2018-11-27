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

template<int D, Entity_kind onwhat,
         class SourceMesh_Wrapper, class TargetMesh_Wrapper>
bool detect_mismatch(SourceMesh_Wrapper const & source_mesh,
                     TargetMesh_Wrapper const & target_mesh,
                     Portage::vector<std::vector<Weights_t>> const & source_ents_and_weights) {

  bool mismatch = false;

  double xsect_volume = 0.0;
  for (auto const& sw_vec : source_ents_and_weights) {
    for (auto const& sw : sw_vec)
      xsect_volume += sw.weights[0];
  }
#ifdef MPI
  double global_xsect_volume = 0.0;
  MPI_Allreduce(&xsect_volume, &global_xsect_volume, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  double global_xsect_volume = xsect_volume;
#endif

  double source_volume = 0.0;
  for (auto it = source_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
       it != source_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++)
    source_volume += (onwhat == Entity_kind::CELL) ?
        source_mesh.cell_volume(*it) : source_mesh.dual_cell_volume(*it);

#ifdef MPI
  double global_source_volume = 0.0;
  MPI_Allreduce(&source_volume, &global_source_volume, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
  double global_source_volume = source_volume;
#endif


  double target_volume = 0.0;
  for (auto it = target_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
       it != target_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++)
    target_volume += (onwhat == Entity_kind::CELL) ?
        target_mesh.cell_volume(*it) : target_mesh.dual_cell_volume(*it);

#ifdef MPI
  double global_target_volume = 0.0;
  MPI_Allreduce(&target_volume, &global_target_volume, 1,
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
  double global_target_volume = target_volume;
#endif


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

    std::vector<double> source_covered_vol;
    if (onwhat == Entity_kind::CELL)
      source_covered_vol.resize(source_mesh.num_owned_cells());
    else
      source_covered_vol.resize(source_mesh.num_owned_nodes());
    
    for (auto it = source_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
         it != source_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++)
      source_covered_vol[*it] = (onwhat == Entity_kind::CELL) ?
          source_mesh.cell_volume(*it) : source_mesh.dual_cell_volume(*it);

    for (auto it = target_mesh.begin(onwhat, Entity_type::PARALLEL_OWNED);
         it != target_mesh.end(onwhat, Entity_type::PARALLEL_OWNED); it++) {
      auto const& sw_vec = source_ents_and_weights[*it];
      for (auto const& sw : sw_vec)
        source_covered_vol[sw.entityID] -= sw.weights[0];
      break;
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
      double covered_vol = 0.0;
      auto const& sw_vec = source_ents_and_weights[*it];
      for (auto const& sw : sw_vec)
        covered_vol += sw.weights[0];
      double target_cell_vol = (onwhat == Entity_kind::CELL) ?
          target_mesh.cell_volume(*it) : target_mesh.dual_cell_volume(*it);
      if (fabs(covered_vol-target_cell_vol)/target_cell_vol > 1.0e-12) {
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

  return mismatch;
}


// Adjust values for fields to account for boundary mismatch between
// the meshes. We adjust so that integral quantities are conserved
// and a constant field field remains constant (BUT THE VALUE WILL
// BE DIFFERENT)

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

    // For now peg the lower bound at 0.0 and the upper bound at
    // (practically) infinity. Eventually, we have to get this per
    // variable from the application

    double lower_bound = 0.0, upper_bound = std::numeric_limits<double>::max();

    // collect source and target cell volumes because we will need it
    // a few times for each variable we process

    int nsourceents = source_mesh.num_owned_cells();
    Portage::vector<double> source_cell_volumes(nsourceents);
    for (int s = 0; s < nsourceents; s++) 
      source_cell_volumes[s] = source_mesh.cell_volume(s);

    int ntargetents = target_mesh.num_owned_cells();
    Portage::vector<double> target_cell_volumes(ntargetents);
    for (int s = 0; s < ntargetents; s++) 
      target_cell_volumes[s] = target_mesh.cell_volume(s);

    double target_mesh_volume = std::accumulate(target_cell_volumes.begin(),
                                                target_cell_volumes.end(),
                                                0.0);

    // First let's detect if any target mesh entities are not at all
    // covered by source mesh entities. In our initial redistribution
    // phase, we will move as many source cells as needed to from
    // different partitions to cover the target partition, UNLESS NO
    // SOURCE CELLS COVERING THE TARGET CELLS EXIST GLOBALLY. So, at
    // this stage, we can just check on-rank intersections only

    std::vector<double> xsect_volumes(ntargetents);
    std::vector<bool> cell_is_empty(ntargetents, false);
    std::vector<int> emptycells;
    for (int i = 0; i < ntargetents; i++) {
      auto const& sw_vec = source_ents_and_weights[i];

      xsect_volumes[i] = 0.0;
      for (auto const& sw : sw_vec)
        xsect_volumes[i] += sw.weights[0];

      if (fabs(xsect_volumes[i]) < 1.0e-16) {
        emptycells.push_back(i);
        cell_is_empty[i] = true;
      }
    }
    int nempty = emptycells.size();
    
    // Collect the empty target cells in layers starting from the ones
    // next to partially or fully covered cells. At the end of this
    // section, partially or fully covered cells will all have a layer
    // number of 0 and empty cell layers will have positive layer
    // numbers starting from 1

    std::vector<int> layernum(ntargetents, 0);
    std::vector<std::vector<int>> emptylayers;
    if (nempty) {
      std::cerr << "One or more target cells are not covered by " <<
          "any source cells.\n" <<
          " Will assign values based on their neighborhood\n";

      int nlayers = 0;
      int ntagged = 0;
      while (ntagged < nempty) {
        std::vector<int> curlayercells;

        for (int c : emptycells) {
          if (layernum[c] != 0) continue;

          std::vector<int> cnbrs;
          target_mesh.cell_get_node_adj_cells(c, Entity_type::ALL, &cnbrs);
          for (int cnbr : cnbrs)
            if (!cell_is_empty[cnbr] || layernum[cnbr] != 0) {
              // At least one neighbor has some material or will
              // receive some material (indicated by having a +ve
              // layer number)

              curlayercells.push_back(c);
              break;
            }
        }  // for (c in emptcells)
        
        // Tag the current layer cells with the next layer number
        for (int c : curlayercells)
          layernum[c] = nlayers+1;
        ntagged += curlayercells.size();
        
        emptylayers.push_back(curlayercells);
        nlayers++;
      }
    }  // if nempty
      
    // Now process remap variables

    int nvars = src_var_names.size();

    for (int i = 0; i < nvars; i++) {
      std::string const& src_var_name = src_var_names[i];
      double const *source_data_raw;
      source_state.mesh_get_data(onwhat, src_var_names[i], &source_data_raw);
      Portage::pointer<const double> source_data(source_data_raw);

      std::string const& trg_var_name = trg_var_names[i];
      double *target_data_raw;
      target_state.mesh_get_data(onwhat, trg_var_names[i], &target_data_raw);
      Portage::pointer<double> target_data(target_data_raw);

      // Do something here to populate fully uncovered target
      // cells. We have layers of empty cells starting out from fully
      // or partially populated cells. We will assign every empty cell
      // in a layer 75% of the minimum value of all its populated
      // neighbors.

      int curlayernum = 1;
      for (std::vector<int> const& curlayer : emptylayers) {
        for (int c : curlayer) {
          std::vector<int> cnbrs;
          target_mesh.cell_get_node_adj_cells(c, Entity_type::ALL, &cnbrs);

          double aveval = 0.0;
          int nave = 0;
          for (int cnbr : cnbrs)
            if (layernum[cnbr] < curlayernum) {
              aveval += target_data[cnbr];
              nave++;
            }
          aveval /= nave;

          target_data[c] = aveval;
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
                             source_cell_volumes.begin(), 0.0);
      double target_sum =
          std::inner_product(target_data, target_data + ntargetents,
                             target_cell_volumes.begin(), 0.0);

      double diff = target_sum - source_sum;
      double reldiff = diff/source_sum;

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

      double adj_target_mesh_volume = target_mesh_volume;
      double udiff = diff/adj_target_mesh_volume;

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

            // this cell is no longer in play for adjustment - so remove its
            // volume from the target_meshvolume

            adj_target_mesh_volume -= target_cell_volumes[t];

          } else if ((target_data[t]-udiff) > upper_bound) {  // udiff < 0
            // Adding the full deficit will make this cell violate the
            // upper bound. So add only as much as will put this cell
            // exactly at the upper bound

            target_data[t] = upper_bound;

            // this cell is no longer in play for adjustment - so remove its
            // volume from the target_mesh_volume

            adj_target_mesh_volume -= target_cell_volumes[t];

          } else {
            // This is the equivalent of
            //           [curval*cellvol - diff*cellvol/meshvol]
            // curval = ---------------------------------------
            //                       cellvol

            target_data[t] -= udiff;

          }
        }  // iterate through mesh cells

        target_sum =
            std::inner_product(target_data, target_data + ntargetents,
                               target_cell_volumes.begin(), 0.0);

        // If we did not hit lower or upper bounds, this should be
        // zero after the first iteration.  If we did hit some bounds,
        // then recalculate the discrepancy and discrepancy per unit
        // volume, but only taking into account volume of cells that
        // are not already at the bounds - if we use the entire mesh
        // volume for the recalculation, the convergence slows down

        diff = target_sum - source_sum;
        udiff = diff/adj_target_mesh_volume;
        reldiff = diff/source_sum;

        adj_target_mesh_volume = target_mesh_volume;

        iter++;
      }  // while leftover is not zero

      if (fabs(reldiff) < 1.0e-12 && fabs(reldiff) > 1.0e-14) {
        // just subtract the remaining from the first cell
        target_data[0] -= diff/target_cell_volumes[0];
      } else {
        std::cerr << "Redistribution not entirely successfully for variable " <<
            src_var_name << "\n";
        std::cerr << "Relative conservation error is " << diff/source_sum <<
            "\n";
        std::cerr << "Absolute conservation error is " << diff << "\n";
        return false;
      }

    }  // for each variable

  }  // fix_mismatch

}  // namespace Portage

#endif  // PORTAGE_DRIVER_FIX_MISMATCH_H_
