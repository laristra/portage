/*
This is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#ifndef PORTAGE_DRIVER_DETECT_MISMATCH_H_
#define PORTAGE_DRIVER_DETECT_MISMATCH_H_

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
  @file detect_mismatch.h
  @brief Detect if the boundaries of the source and target mesh are not exactly coincident
*/

namespace Portage {

using Wonton::Entity_kind;
using Wonton::Entity_type;
using Wonton::Weights_t;


// Helper function
//
// Get a mask indicating which entities on this processor have not
// been encountered on previous processors. Mask of 1 means this
// entity is being seen for the first time and should figure in any
// calculation, 0 means this entity has been encountered on a previous
// processor. This is useful for meshes where the partitioning of cells
// on ranks is not mutually exclusive.

template<Entity_kind onwhat, class Mesh_Wrapper>
bool get_unique_entity_masks(Mesh_Wrapper const &mesh,
                             std::vector<int> *unique_mask) {
  int rank = 0;
  int nprocs = 1;
#ifdef ENABLE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif


  // collect source and target cell volumes because we will need it
  // a few times for each variable we process

  int nents = (onwhat == Entity_kind::CELL) ?
      mesh.num_owned_cells() : mesh.num_owned_nodes();

  unique_mask->resize(nents, 1);

#ifdef ENABLE_MPI
  if (nprocs > 1) {
    std::vector<int> nents_all(nprocs, 0);
    MPI_Allgather(&nents, 1, MPI_INT, &(nents_all[0]), 1, MPI_INT,
                  MPI_COMM_WORLD);
    int maxents = *std::max_element(nents_all.begin(), nents_all.end());

    std::vector<int> gids(maxents, -1);
    for (int e = 0; e < nents; e++)
      gids[e] = mesh.get_global_id(e, onwhat);

    std::vector<int> gids_all(nprocs*maxents);
    MPI_Allgather(&(gids[0]), maxents, MPI_INT, &(gids_all[0]), maxents,
                  MPI_INT, MPI_COMM_WORLD);

    // Add gids from all lower ranks into a set (no duplicates)

    std::unordered_set<int> unique_gids;
    for (int p = 0; p < rank; p++)
      for (int e = 0; e < nents_all[p]; e++) {
        int gid = gids_all[p*maxents+e];
        unique_gids.insert(gid);
      }

    // Now go through gids of this rank and mask the ones that are already
    // occurring in other ranks.
    for (int e = 0; e < nents; e++) {
      auto itpair = unique_gids.insert(gids[e]);
      if (!itpair.second)
        (*unique_mask)[e] = 0;  // ent already in set; mask this instance
    }
  }
#endif
}  // get_unique_entity_masks



// Check if we have mismatch of mesh boundaries
// T is the target mesh, S is the source mesh
// T_i is the i'th cell of the target mesh and
// |*| signifies the extent/volume of an entity
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
class MismatchFixer {
 public:

  MismatchFixer(SourceMesh_Wrapper const& source_mesh,
                SourceState_Wrapper const& source_state,
                TargetMesh_Wrapper const& target_mesh,
                TargetState_Wrapper & target_state,
                Portage::vector<std::vector<Weights_t>> const & source_ents_and_weights) :
      source_mesh_(source_mesh), source_state_(source_state),
      target_mesh_(target_mesh), target_state_(target_state) {

#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_);
#endif

    // GLOBAL SOURCE VOLUME
    nsourceents_ = (onwhat == Entity_kind::CELL) ?
        source_mesh_.num_owned_cells() : source_mesh_.num_owned_nodes();

    // Get info about which entities on this processor should be
    // masked out and not accounted for in calculations because they
    // were already encountered on a lower rank processor. We have to
    // do this because in our distributed runs, our source partitions
    // don't form a strict tiling (no overlaps) after redistribution

    std::vector<int> source_ent_masks(nsourceents_, 1);
    get_unique_entity_masks<onwhat, SourceMesh_Wrapper>(source_mesh_,
                                                        &source_ent_masks);

    // collect volumes of entities that are not masked out and sum them up

    source_ent_volumes_.resize(nsourceents_, 0.0);
    for (int s = 0; s < nsourceents_; s++)
      source_ent_volumes_[s] = (onwhat == Entity_kind::CELL) ?
          source_ent_masks[s]*source_mesh_.cell_volume(s) :
          source_ent_masks[s]*source_mesh_.dual_cell_volume(s);

    source_volume_ =
        std::accumulate(source_ent_volumes_.begin(), source_ent_volumes_.end(),
                        0.0);

#ifdef ENABLE_MPI
    global_source_volume_ = 0.0;
    MPI_Allreduce(&source_volume_, &global_source_volume_, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
#else
    global_source_volume_ = source_volume_;
#endif


    // GLOBAL TARGET VOLUME
    ntargetents_ = (onwhat == Entity_kind::CELL) ?
        target_mesh_.num_owned_cells() : target_mesh_.num_owned_nodes();

    target_ent_volumes_.resize(ntargetents_, 0.0);
    for (int t = 0; t < ntargetents_; t++)
      target_ent_volumes_[t] = (onwhat == Entity_kind::CELL) ?
          target_mesh_.cell_volume(t) : target_mesh_.dual_cell_volume(t);

    target_volume_ = std::accumulate(target_ent_volumes_.begin(),
                                     target_ent_volumes_.end(), 0.0);
#ifdef ENABLE_MPI
    global_target_volume_ = 0.0;
    MPI_Allreduce(&target_volume_, &global_target_volume_, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);
#else
    global_target_volume_ = target_volume_;
#endif


    // GLOBAL INTERSECTION VOLUME
    // In our initial redistribution phase, we will move as many
    // source cells as needed from different partitions to cover the
    // target partition, UNLESS NO SOURCE CELLS COVERING THE TARGET
    // CELLS EXIST GLOBALLY. So, at this stage, we can just check
    // on-rank intersections only

    xsect_volumes_.resize(ntargetents_, 0.0);
    for (int t = 0; t < ntargetents_; t++) {
      std::vector<Weights_t> const& sw_vec = source_ents_and_weights[t];
      for (auto const& sw : sw_vec)
        xsect_volumes_[t] += sw.weights[0];
    }

    double xsect_volume = std::accumulate(xsect_volumes_.begin(),
                                          xsect_volumes_.end(), 0.0);
#ifdef ENABLE_MPI
    global_xsect_volume_ = 0.0;
    MPI_Allreduce(&xsect_volume, &global_xsect_volume_, 1,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
    global_xsect_volume_ = xsect_volume;
#endif




    // Are some source cells not fully covered by target cells?
    relvoldiff_source_ =
        fabs(global_xsect_volume_-global_source_volume_)/global_source_volume_;
    if (relvoldiff_source_ > voldifftol_) {

      mismatch_ = true;

      if (rank_ == 0) 
        std::cerr << "\n** MESH MISMATCH -" <<
            " some source cells are not fully covered by the target mesh\n";

#ifdef DEBUG
      // Find one source cell (or dual cell) that is not fully covered
      // by the target mesh and output its ID. Unfortunately, that means
      // processing all source cells. We initialize each source cell to
      // its volume and subtract any intersection volume we find between
      // it and a target cell

      // Also, it will likely give a false positive in distributed
      // meshes because a source cell may be covered by target cells
      // from multiple processors

      std::vector<double> source_covered_vol(source_ent_volumes_);
      for (auto it = target_mesh_.begin(onwhat, Entity_type::PARALLEL_OWNED);
           it != target_mesh_.end(onwhat, Entity_type::PARALLEL_OWNED); it++) {
        auto const& sw_vec = source_ents_and_weights[*it];
        for (auto const& sw : sw_vec)
          source_covered_vol[sw.entityID] -= sw.weights[0];
      }

      for (auto it = source_mesh_.begin(onwhat, Entity_type::PARALLEL_OWNED);
           it != source_mesh_.end(onwhat, Entity_type::PARALLEL_OWNED); it++)
        if (source_covered_vol[*it] > voldifftol_) {
          if (onwhat == Entity_kind::CELL)
            std::cerr << "Source cell " << *it <<
                " not fully covered by target cells \n";
          else
            std::cerr << "Source dual cell " << *it <<
                " not fully covered by target dual cells \n";
          break;
        }
      std::cerr << "\n";
#endif
    }

    // Are some target cells not fully covered by source cells?

    relvoldiff_target_ =
        fabs(global_xsect_volume_-global_target_volume_)/global_target_volume_;
    if (relvoldiff_target_ > voldifftol_) {

      mismatch_ = true;

      if (rank_ == 0)
        std::cerr << "\n** MESH MISMATCH -" <<
            " some target cells are not fully covered by the source mesh\n";

#ifdef DEBUG
      // Find one target cell that is not fully covered by the source mesh and
      // output its ID
      for (auto it = target_mesh_.begin(onwhat, Entity_type::PARALLEL_OWNED);
           it != target_mesh_.end(onwhat, Entity_type::PARALLEL_OWNED); it++) {
        int t = *it;
        double covered_vol = 0.0;
        std::vector<Weights_t> const& sw_vec = source_ents_and_weights[t];
        for (auto const& sw : sw_vec)
          covered_vol += sw.weights[0];
        if (fabs(covered_vol-target_ent_volumes_[t])/target_ent_volumes_[t] > voldifftol_) {
          if (onwhat == Entity_kind::CELL)
            std::cerr << "Target cell " << *it << " on rank " << rank_ <<
                " not fully covered by source cells \n";
          else
            std::cerr << "Target dual cell " << *it << " on rank " << rank_ <<
                " not fully covered by source dual cells \n";
          break;
        }
      }
      std::cerr << "\n";
#endif

    }

    if (!mismatch_) return;


    // Discrepancy between intersection volume and source mesh volume PLUS
    // Discrepancy between intersection volume and target mesh volume
    relvoldiff_ = relvoldiff_source_ + relvoldiff_target_;


    // Collect the empty target cells in layers starting from the ones
    // next to partially or fully covered cells. At the end of this
    // section, partially or fully covered cells will all have a layer
    // number of 0 and empty cell layers will have positive layer
    // numbers starting from 1

    std::vector<bool> is_empty(ntargetents_, false);
    std::vector<int> emptyents;
    for (int t = 0; t < ntargetents_; t++) {
      if (fabs(xsect_volumes_[t]) < std::numeric_limits<double>::epsilon()) {
        emptyents.push_back(t);
        is_empty[t] = true;
      }
    }
    int nempty = emptyents.size();

    int global_nempty = nempty;
#ifdef ENABLE_MPI    
    int *nempty_all = new int[nprocs_];
    MPI_Gather(&nempty, 1, MPI_INT, nempty_all, 1, MPI_INT, 0, MPI_COMM_WORLD);
    global_nempty = std::accumulate(nempty_all, nempty_all+nprocs_, 0.0);
#endif
    
    if (global_nempty && rank_ == 0) {
      if (onwhat == Entity_kind::CELL)
        std::cerr << "One or more target cells are not covered by " <<
            "ANY source cells.\n" <<
            " Will assign values based on their neighborhood\n";
      else
        std::cerr << "One or more target dual cells are not covered by " <<
            "ANY source dual cells.\n" <<
            " Will assign values based on their neighborhood\n";
    }
    
    if (nempty) {
      layernum_.resize(ntargetents_, 0);

      int nlayers = 0;
      int ntagged = 0;
      while (ntagged < nempty) {
        std::vector<int> curlayerents;

        for (int ent : emptyents) {
          if (layernum_[ent] != 0) continue;

          std::vector<int> nbrs;
          if (onwhat == Entity_kind::CELL)
            target_mesh_.cell_get_node_adj_cells(ent, Entity_type::ALL, &nbrs);
          else
            target_mesh_.node_get_cell_adj_nodes(ent, Entity_type::ALL, &nbrs);
          for (int nbr : nbrs)
            if (!is_empty[nbr] || layernum_[nbr] != 0) {
              // At least one neighbor has some material or will
              // receive some material (indicated by having a +ve
              // layer number)

              curlayerents.push_back(ent);
              break;
            }
        }  // for (ent in emptyents)

        // Tag the current layer cells with the next layer number
        for (int ent : curlayerents)
          layernum_[ent] = nlayers+1;
        ntagged += curlayerents.size();

        emptylayers_.push_back(curlayerents);
        nlayers++;
      }
    }  // if nempty

    delete [] nempty_all;

  }  // MismatchFixer




  // has this problem been found to have mismatched mesh boundaries?

  bool has_mismatch() const {return mismatch_;}




  /// @brief Repair the remapped field to account for boundary mismatch
  /// @param src_var_name        field variable on source mesh
  /// @param trg_var_name        field variable on target mesh
  /// @param global_lower_bound  lower limit on variable
  /// @param global_upper_bound  upper limit on variable
  /// @param partial_fixup_type  type of fixup in case of partial mismatch
  /// @param empty_fixup_type    type of fixup in empty target entities
  ///
  /// partial_fixup_type can be one of three types:
  ///
  /// CONSTANT - Fields will see no perturbations BUT REMAP WILL BE
  ///            NON-CONSERVATIVE (constant preserving, not linearity
  ///            preserving)
  /// CONSERVATIVE - REMAP WILL BE CONSERVATIVE but perturbations will
  ///                occur in the field (constant fields may not stay
  ///                constant if there is mismatch)
  /// SHIFTED_CONSERVATIVE - REMAP WILL BE CONSERVATIVE and field
  ///                        perturbations will be minimum but field
  ///                        values may be shifted (Constant fields
  ///                        will be shifted to different constant; no
  ///                        guarantees on linearity preservation)
  ///
  /// empty_fixup_type can be one of two types:
  ///
  /// LEAVE_EMPTY - Leave empty cells as is
  /// EXTRAPOLATE - Fill empty cells with extrapolated values
  /// FILL        - Fill empty cells with specified values (not yet implemented)

  bool fix_mismatch(std::string const & src_var_name,
                    std::string const & trg_var_name,
                    double global_lower_bound = -std::numeric_limits<double>::max(),
                    double global_upper_bound = std::numeric_limits<double>::max(),
                    double conservation_tol = 1e2*std::numeric_limits<double>::epsilon(),
                    Partial_fixup_type partial_fixup_type =
                    Partial_fixup_type::SHIFTED_CONSERVATIVE,
                    Empty_fixup_type empty_fixup_type =
                    Empty_fixup_type::EXTRAPOLATE) {


    if (source_state_.field_type(onwhat, src_var_name) ==
        Field_type::MESH_FIELD)
      return fix_mismatch_meshvar(src_var_name, trg_var_name,
                                  global_lower_bound, global_upper_bound,
                                  conservation_tol,
                                  partial_fixup_type, empty_fixup_type);
  }


  /// @brief Repair a remapped mesh field to account for boundary mismatch

  bool fix_mismatch_meshvar(std::string const & src_var_name,
                            std::string const & trg_var_name,
                            double global_lower_bound,
                            double global_upper_bound,
                            double conservation_tol = 1e2*std::numeric_limits<double>::epsilon(),
                            Partial_fixup_type partial_fixup_type =
                            Partial_fixup_type::SHIFTED_CONSERVATIVE,
                            Empty_fixup_type empty_fixup_type =
                            Empty_fixup_type::EXTRAPOLATE) {

    // Now process remap variables
    double const *source_data;
    source_state_.mesh_get_data(onwhat, src_var_name, &source_data);

    double *target_data;
    target_state_.mesh_get_data(onwhat, trg_var_name, &target_data);

    if (empty_fixup_type != Empty_fixup_type::LEAVE_EMPTY) {
      // Do something here to populate fully uncovered target
      // cells. We have layers of empty cells starting out from fully
      // or partially populated cells. We will assign every empty cell
      // in a layer the average value of all its populated neighbors.
      // IN A DISTRIBUTED MESH IT _IS_ POSSIBLE THAT AN EMPTY ENTITY
      // WILL NOT HAVE ANY OWNED NEIGHBOR IN THIS PARTITION THAT HAS
      // MEANINGFUL DATA TO EXTRAPOLATE FROM (we remap data only to
      // owned entities)

      int curlayernum = 1;
      for (std::vector<int> const& curlayer : emptylayers_) {
        for (int ent : curlayer) {
          std::vector<int> nbrs;
          if (onwhat == Entity_kind::CELL)
            target_mesh_.cell_get_node_adj_cells(ent, Entity_type::PARALLEL_OWNED,
                                                 &nbrs);
          else
            target_mesh_.node_get_cell_adj_nodes(ent, Entity_type::PARALLEL_OWNED,
                                                 &nbrs);
          
          double aveval = 0.0;
          int nave = 0;
          for (int nbr : nbrs) {
            if (layernum_[nbr] < curlayernum) {
              aveval += target_data[nbr];
              nave++;
            }
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
    }


    if (partial_fixup_type == Partial_fixup_type::CONSTANT) {
      // In interpolate step, we divided the accumulated integral in a
      // target cell by the intersection volume to preserve a constant
      // and violate conservation. So nothing to do here.
    
      return true;

    } else if (partial_fixup_type == Partial_fixup_type::CONSERVATIVE) {
      // In interpolate step, we divided the accumulated integral (U)
      // in a target cell by the intersection volume (v_i) instead of
      // the target cell volume (v_c) to give a target field of u_t =
      // U/v_i. In partially filled cells, this will preserve a
      // constant source field but fill the cell with too much
      // material. To restore conservation, we undo the division by
      // the intersection volume and then divide by the cell volume
      // (u'_t = U/v_c = u_t*v_i/v_c). This does not affect the values
      // in fully filled cells
      
      for (int t = 0; t < ntargetents_; t++) {
        if (fabs(xsect_volumes_[t]-target_ent_volumes_[t])/target_ent_volumes_[t] > voldifftol_)
          target_data[t] *= xsect_volumes_[t]/target_ent_volumes_[t];
      }

      return true;

    } else if (partial_fixup_type == Partial_fixup_type::SHIFTED_CONSERVATIVE) {
    
      // At this point assume that all cells have some value in them
      // for the variable

      // Now compute the net discrepancy between integrals over source
      // and target. Excess comes from target cells not fully covered by
      // source cells and deficit from source cells not fully covered by
      // target cells

      double source_sum =
          std::inner_product(source_data, source_data + nsourceents_,
                             source_ent_volumes_.begin(), 0.0);

      double global_source_sum = source_sum;
#ifdef ENABLE_MPI
      MPI_Allreduce(&source_sum, &global_source_sum, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      double target_sum =
          std::inner_product(target_data, target_data + ntargetents_,
                             target_ent_volumes_.begin(), 0.0);

      double global_target_sum = target_sum;
#ifdef ENABLE_MPI
      MPI_Allreduce(&target_sum, &global_target_sum, 1, MPI_DOUBLE, MPI_SUM,
                    MPI_COMM_WORLD);
#endif

      double global_diff = global_target_sum - global_source_sum;
      double reldiff = global_diff/global_source_sum;

      if (fabs(reldiff) < conservation_tol)
        return true;  // discrepancy is too small - nothing to do

      // Now redistribute the discrepancy among cells in proportion to
      // their volume. This will restore conservation and if the
      // original distribution was a constant it will make the field a
      // slightly different constant. Do multiple iterations because
      // we may have some leftover quantity that we were not able to
      // add/subtract from some cells in any given iteration. We don't
      // expect that process to take more than two iterations at the
      // most.

      double adj_target_volume = target_volume_;
      double global_adj_target_volume = global_target_volume_;

      // sort of a "unit" discrepancy or difference per unit volume
      double udiff = global_diff/global_adj_target_volume;

      int iter = 0;
      while (fabs(reldiff) > conservation_tol && iter < 5) {
        for (auto it = target_mesh_.begin(onwhat, Entity_type::PARALLEL_OWNED);
             it != target_mesh_.end(onwhat, Entity_type::PARALLEL_OWNED); it++) {
          int t = *it;

          if ((target_data[t]-udiff) < global_lower_bound) {
            // Subtracting the full excess will make this cell violate the
            // lower bound. So subtract only as much as will put this cell
            // exactly at the lower bound

            target_data[t] = global_lower_bound;

#ifdef DEBUG            
            std::cerr << "Hit lower bound for cell " << t << " on rank " <<
                rank_ << "\n";
#endif
            
            // this cell is no longer in play for adjustment - so remove its
            // volume from the adjusted target_volume

            adj_target_volume -= target_ent_volumes_[t];

          } else if ((target_data[t]-udiff) > global_upper_bound) {  // udiff < 0
            // Adding the full deficit will make this cell violate the
            // upper bound. So add only as much as will put this cell
            // exactly at the upper bound

            target_data[t] = global_upper_bound;

#ifdef DEBUG            
            std::cerr << "Hit upper bound for cell " << t << " on rank " <<
                rank_ << "\n";
#endif
            
            // this cell is no longer in play for adjustment - so remove its
            // volume from the adjusted target_volume

            adj_target_volume -= target_ent_volumes_[t];

          } else {
            // This is the equivalent of
            //           [curval*cellvol - diff*cellvol/meshvol]
            // curval = ---------------------------------------
            //                       cellvol

            target_data[t] -= udiff;

          }
        }  // iterate through mesh cells

        // Compute the new integral over all processors

        target_sum = std::inner_product(target_data, target_data + ntargetents_,
                                        target_ent_volumes_.begin(), 0.0);

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
        global_adj_target_volume = global_target_volume_;

        iter++;
      }  // while leftover is not zero

      if (fabs(reldiff) > 1.0e-14) {
        if (rank_ == 0) {
          std::cerr << "Redistribution not entirely successfully for variable " <<
              src_var_name << "\n";
          std::cerr << "Relative conservation error is " << reldiff << "\n";
          std::cerr << "Absolute conservation error is " << global_diff << "\n";
          return false;
        }
      }

      return true;
    } else {
      std::cerr << "Unknown Partial fixup type\n";
      return false;
    }

  }  // fix_mismatch_meshvar

 private:
  SourceMesh_Wrapper const& source_mesh_;
  SourceState_Wrapper const& source_state_;
  TargetMesh_Wrapper const& target_mesh_;
  TargetState_Wrapper & target_state_;
  int nsourceents_, ntargetents_;
  std::vector<double> source_ent_volumes_, target_ent_volumes_, xsect_volumes_;
  double source_volume_, target_volume_;
  double global_source_volume_, global_target_volume_, global_xsect_volume_;
  double relvoldiff_source_, relvoldiff_target_, relvoldiff_;
  std::vector<int> is_empty_, layernum_;
  std::vector<std::vector<int>> emptylayers_;
  bool mismatch_ = false;
  int rank_ = 0, nprocs_ = 1;
  double voldifftol_ = 1e2*std::numeric_limits<double>::epsilon();
};  // MismatchFixer

}  // namespace Portage

#endif  // PORTAGE_DRIVER_FIX_MISMATCH_H_
