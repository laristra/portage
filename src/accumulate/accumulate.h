/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef ACCUMULATE_H_INC_
#define ACCUMULATE_H_INC_

#include <vector>
#include <memory>
#include <cassert>

#include "portage/support/portage.h"
#include "portage/support/Point.h"
#include "portage/support/weight.h"
#include "portage/support/basis.h"
#include "portage/swarm/swarm.h"
#include "portage/support/Matrix.h"

namespace Portage {
namespace Meshfree {

using std::vector;
using std::shared_ptr;

/// Different kinds of estimates to do
enum EstimateType {
  KernelDensity,
  LocalRegression
};

/// Denote which points smoothing length is centered on
enum WeightCenter {
  Gather,  ///< centered on target points
  Scatter  ///< centered on source points
};

/**
 * @class Accumulate portage/accumulate/accumulate.h
 * @brief Compute the local regression estimator corrected weights
 * @class Accumulate portage/accumulate/accumulate.h
 * 
 * This class does the meat of local regression. It computes weight functions, 
 * and the local regression corrections to those weights. 
 */
template<size_t dim,
         class SourceSwarm,
         class TargetSwarm>
class Accumulate {
 public:

  /**
   * @brief Constructor
   * @param source source swarm
   * @param target target swarm
   * @param estimate what type of estimate to do: currently kernel density or local regression
   * @param center where the smoothing paramaters are attached: Gather->target, Scatter->source
   * @param kernels kernel specifiers
   * @param geometries geometry specifiers
   * @param smoothing smoothing lengths (or bandwidths)
   * 
   * The parameters @code kernels@endcode, @code geometries@endcode and @code smoothing@endcode 
   * all must have the same length. If @code center@endcode is @code Gather@endcode, then 
   * the length is the size of the target swarm. If @code center@endcode is @code Scatter@endcode, 
   * the length is the size of the source swarm.
   */
  Accumulate(
      SourceSwarm const& source,
      TargetSwarm const& target,
      EstimateType estimate,
      WeightCenter center,
      vector<Weight::Kernel> const& kernels,
      vector<Weight::Geometry> const& geometries,
      vector<vector<vector<double>>> const& smoothing,
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
      n_particles = target_.num_owned_particles();
    } else if (center == Scatter) {
      n_particles = source_.num_owned_particles();
    }
    assert(n_particles == kernels_.size());
    assert(n_particles == geometries_.size());
    assert(n_particles == smoothing_.size());
  }

  /** 
   * @brief Evaluate meshfree weight function
   * @param particleA source index
   * @param particleB target index
   * @return value of weight function
   *
   * Information in constructor arguments decides the details of the weight function.
   */
  double weight(const size_t particleA, const size_t particleB)
  {
    double result;
    Point<dim> x = target_.get_particle_coordinates(particleB);
    Point<dim> y = source_.get_particle_coordinates(particleA);
    if (center_ == Gather) {
      result = Weight::eval<dim>(geometries_[particleB],
                                 kernels_[particleB],
                                 x,y,
                                 smoothing_[particleB]);
    } else if (center_ == Scatter) {
      result = Weight::eval<dim>(geometries_[particleA],
                                 kernels_[particleA],
                                 x,y,
                                 smoothing_[particleA]);
    }
    return result;
  }

  /** 
   * @brief Compute local regression correction to weight function
   * @param particleB target particle index in target swarm
   * @param source_particles list of source particle neighbors of target particle
   * @return the weight or the corrected weight function according to @code estimate_@endcode
   *
   * The return matrix is of size n x m  where n is the length of source_particles, 
   * and m is the size of the basis. 
   */
  vector<Weights_t>
  operator() (size_t const particleB, vector<unsigned int> const& source_particles) {
    vector<Weights_t> result;
    result.reserve(source_particles.size());
    
    switch (estimate_) {
      case KernelDensity:  {
        for (auto const& particleA : source_particles) {
          double weight_val = weight(particleA, particleB);
          vector<double> pair_result(1, weight_val);
          result.emplace_back(particleA, pair_result);
        }
        break;
      }
      case LocalRegression: {
        size_t nbasis = Basis::function_size<dim>(basis_);
        Point<dim> x = target_.get_particle_coordinates(particleB);
        
        // Calculate weights and moment matrix (P*W*transpose(P))
        vector<double> weight_val(source_particles.size());
        Matrix moment(nbasis,nbasis,0.);
        size_t iA = 0;
        for (auto const& particleA : source_particles) {
          weight_val[iA] = weight(particleA, particleB); // save weights for later
          Point<dim> y = source_.get_particle_coordinates(particleA);
          auto basis = Basis::shift<dim>(basis_,x,y);
          for (size_t i=0; i<nbasis; i++) {
            for (size_t j=0; j<nbasis; j++) {
              moment[i][j] += basis[i]*basis[j]*weight_val[iA];
            }
          }
          iA++;
        }
        
        // Calculate inverse(P*W*transpose(P))*P*W
        iA = 0;
        for (auto const& particleA : source_particles) {
          vector<double> pair_result(nbasis);
          Point<dim> y = source_.get_particle_coordinates(particleA);
          vector<double> basis = Basis::shift<dim>(basis_,x,y);

          // recast as a Portage::Matrix
          Matrix basis_matrix(nbasis,1);
          for (size_t i=0; i<nbasis; i++) basis_matrix[i][0] = basis[i];

          // solve the linear system
#ifdef HAVE_LAPACKE 
          Matrix pair_result_matrix = moment.solve(basis_matrix, "lapack-posv");
#else
          Matrix pair_result_matrix = moment.solve(basis_matrix);
#endif

          for (size_t i=0; i<nbasis; i++) pair_result[i] = pair_result_matrix[i][0]*weight_val[iA];
          result.emplace_back(particleA, pair_result);
          iA++;
        }
        break;
      }
      default:  assert(false);
    }
    return result;
  }
  
 private:
  SourceSwarm const& source_;
  TargetSwarm const& target_;
  EstimateType estimate_;
  WeightCenter center_;
  vector<Weight::Kernel> const& kernels_;
  vector<Weight::Geometry> const& geometries_;
  vector<vector<vector<double>>> const& smoothing_;
  Basis::Type basis_;
};

}}

#endif // ACCUMULATE_H_INC
