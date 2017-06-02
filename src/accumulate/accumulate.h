/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/

#ifndef ACCUMULATE_H_INC_
#define ACCUMULATE_H_INC_

#include <vector>
#include <memory>
#include <cassert>

#include "portage/support/Point.h"
#include "portage/support/weight.h"
#include "portage/support/basis.h"
#include "portage/swarm/swarm.h"
#include "portage/swarm/swarm_state.h"
#include "portage/support/Matrix.h"

namespace Portage {
namespace Meshfree {

enum EstimateType {
  KernelDensity,
  LocalRegression
};

/// Denote which points smoothing length is centered on
enum WeightCenter {
  Gather,  ///< target
  Scatter  ///< source
};

using std::vector;
using std::shared_ptr;

/*!
 * @brief Compute the local regression estimator corrected weights
 */
template<size_t dim>
class Accumulate {
 public:
  Accumulate(
      shared_ptr<Swarm<dim>> source,
      shared_ptr<Swarm<dim>> target,
      EstimateType estimate,
      WeightCenter center,
      shared_ptr<vector<Weight::Kernel>> kernels,
      shared_ptr<vector<Weight::Geometry>> geometries,
      shared_ptr<vector<vector<vector<double>>>> smoothing,
      Basis::Type basis):
   source_(source),
   target_(target),
   estimate_(estimate),
   center_(center),
   kernels_(kernels),
   geometries_(geometries),
   smoothing_(smoothing),
   basis_(basis)
  {
    // check sizes of inputs are consistent
    size_t n_particles;
    if (center == Gather) {
      n_particles = target_->num_owned_cells();
    } else if (center == Scatter) {
      n_particles = source_->num_owned_cells();
    }
    assert(n_particles == kernels_->size());
    assert(n_particles == geometries_->size());
    assert(n_particles == smoothing_->size());
  }

  /** @brief Evaluate meshfree weight function
   * @param particleA source index
   * @param particleB target index
   */
  double weight(const size_t particleA, const size_t particleB)
  {
    double result;
    Point<dim> x = target_->get_particle_coordinates(particleB);
    Point<dim> y = source_->get_particle_coordinates(particleA);
    if (center_ == Gather) {
      result = Weight::eval<dim>((*geometries_)[particleB],
                                 (*kernels_)[particleB],
                                 x,y,
                                 (*smoothing_)[particleB]);
    } else if (center_ == Scatter) {
      result = Weight::eval<dim>((*geometries_)[particleA],
                                 (*kernels_)[particleA],
                                 x,y,
                                 (*smoothing_)[particleA]);
    }
    return result;
  }

  vector<vector<double>>
  operator() (size_t const particleB, vector<size_t> const& source_particles) {
    vector<vector<double>> result;
    result.reserve(source_particles.size());
    
    switch (estimate_) {
      case KernelDensity:  {
        for (auto const& particleA : source_particles) {
          double weight_val = weight(particleA, particleB);
          vector<double> pair_result(1, weight_val);
          result.push_back(pair_result);
        }
        break;
      }
      case LocalRegression: {
        size_t nbasis = Basis::function_size<dim>(basis_);
        Point<dim> x = target_->get_particle_coordinates(particleB);
        
        // Calculate weights and moment matrix (P*W*transpose(P))
        vector<double> weight_val(source_particles.size());
	Matrix moment(nbasis,nbasis,0.);
	size_t iA = 0;
        for (auto const& particleA : source_particles) {
          weight_val[iA] = weight(particleA, particleB); // save weights for later
          Point<dim> y = source_->get_particle_coordinates(particleA);
          auto basis = Basis::shift<dim>(basis_,x,y);
          for (size_t i=0; i<nbasis; i++) {
            for (size_t j=0; j<nbasis; j++) {
              moment[i][j] += basis[i]*basis[j]*weight_val[iA];
            }
          }
	  iA++;
        }

        auto inverse_moment = moment.inverse();
        
        // Calculate inverse(P*W*transpose(P))*P*W
	iA = 0;
        for (auto const& particleA : source_particles) {
          vector<double> pair_result(nbasis);
          Point<dim> y = source_->get_particle_coordinates(particleA);
          auto basis = Basis::shift<dim>(basis_,x,y);
          pair_result = inverse_moment*basis;
          for (size_t i=0; i<nbasis; i++) pair_result[i] *= weight_val[iA];
          result.push_back(pair_result);
	  iA++;
        }
        break;
      }
      default:  assert(false);
    }
    return result;
  }
  
 private:
  shared_ptr<Swarm<dim>> source_;
  shared_ptr<Swarm<dim>> target_;
  EstimateType estimate_;
  WeightCenter center_;
  shared_ptr<vector<Weight::Kernel>> kernels_;
  shared_ptr<vector<Weight::Geometry>> geometries_;
  shared_ptr<vector<vector<vector<double>>>> smoothing_;
  Basis::Type basis_;
};

}}

#endif // ACCUMULATE_H_INC
