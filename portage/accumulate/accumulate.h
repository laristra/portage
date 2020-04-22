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

// portage includes
#include "portage/support/portage.h"
#include "portage/support/weight.h"
#include "portage/support/basis.h"
#include "portage/support/operator.h"
#include "portage/swarm/swarm.h"

// wonton includes
#include "wonton/support/Matrix.h"

namespace Portage { namespace Meshfree {

/**
 * @brief Different kinds of estimates to do
 */
enum EstimateType {
  KernelDensity,
  LocalRegression,
  OperatorRegression
};

/**
 * @brief Denote which points smoothing length is centered on
 */
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
template<int dim,
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
   * @param geometries geometry specifiers for weight functions
   * @param smoothing smoothing lengths (or bandwidths)
   * 
   * The parameters @code kernels@endcode, @code geometries@endcode and @code smoothing@endcode 
   * all must have the same length. If @code center@endcode is @code Gather@endcode, then 
   * the length is the size of the target swarm. If @code center@endcode is @code Scatter@endcode, 
   * the length is the size of the source swarm.
   */
  Accumulate(SourceSwarm const& source, TargetSwarm const& target,
             EstimateType estimate, WeightCenter center,
             Wonton::vector<Weight::Kernel> const& kernels,
             Wonton::vector<Weight::Geometry> const& geometries,
             Wonton::vector<std::vector<std::vector<double>>> const& smoothing,
             basis::Type basis, oper::Type operator_spec = oper::LastOperator,
             Wonton::vector<oper::Domain> const& operator_domain = {},
             Wonton::vector<std::vector<Wonton::Point<dim>>> const& operator_data = {})
    : source_(source),
      target_(target),
      estimate_(estimate),
      center_(center),
      kernels_(kernels),
      geometries_(geometries),
      smoothing_(smoothing),
      basis_(basis),
      operator_spec_(operator_spec),
      operator_domain_(operator_domain),
      operator_data_(operator_data)
  {
#ifdef DEBUG
    // check sizes of inputs are consistent
    size_t n_particles = (center == Gather ? target_.num_owned_particles()
                                           : source_.num_particles());

    assert(n_particles == kernels_.size());
    assert(n_particles == geometries_.size());
    assert(n_particles == smoothing_.size());
    if (operator_spec_ != oper::LastOperator) {
      unsigned const num_target_particles = target_.num_owned_particles();
      assert(operator_data_.size()   == num_target_particles);
      assert(operator_domain_.size() == num_target_particles);
    }
#endif
  }

  /** 
   * @brief Evaluate meshfree weight function
   * @param particleA target index
   * @param particleB source index
   * @return value of weight function
   *
   * Information in constructor arguments decides the details of the weight function.
   */
  double weight(const size_t particleA, const size_t particleB) {
    double result = 0.0;
    Wonton::Point<dim> x = target_.get_particle_coordinates(particleA);
    Wonton::Point<dim> y = source_.get_particle_coordinates(particleB);
    if (center_ == Gather) {
      result = Weight::eval<dim>(geometries_[particleA],
                                 kernels_[particleA],
                                 x,y,
                                 smoothing_[particleA]);
    } else if (center_ == Scatter) {
      result = Weight::eval<dim>(geometries_[particleB],
                                 kernels_[particleB],
                                 y,x, // faceted weights are asymmetric
                                 smoothing_[particleB]);
    }
    return result;
  }

  /** 
   * @brief Compute local regression correction to weight function
   * @param particleA target particle index in target swarm
   * @param source_particles list of source particle neighbors of target particle
   * @return the weight or the corrected weight function according to @code estimate_@endcode
   *
   * The return matrix is of size n x m  where n is the length of source_particles, 
   * and m is the size of the basis. 
   */
  std::vector<Weights_t>
  operator() (size_t const particleA, std::vector<int> const& source_particles) {
    std::vector<Weights_t> result;
    result.reserve(source_particles.size());

    switch (estimate_) {
      case KernelDensity:  {
        for (auto const& particleB : source_particles) {
          double weight_val = weight(particleA, particleB);
          std::vector<double> pair_result(1, weight_val);
          result.emplace_back(particleB, pair_result);
        }
        break;
      }
      case OperatorRegression:
      case LocalRegression: {
        size_t nbasis = basis::function_size<dim>(basis_);
        Wonton::Point<dim> x = target_.get_particle_coordinates(particleA);

	      // If too few particles, set estimate to zero for this target
	      bool zilchit = false;
        if (source_particles.size() < nbasis) {
          zilchit = true;
        }

        // Calculate weights and moment matrix (P*W*transpose(P))
	      std::vector<double> weight_val(source_particles.size());
        Wonton::Matrix moment(nbasis,nbasis,0.);
        size_t iB = 0;
        if (not zilchit) {
          for (auto const& particleB : source_particles) {
            weight_val[iB] = weight(particleA, particleB); // save weights for later
            Wonton::Point<dim> y = source_.get_particle_coordinates(particleB);
            auto basis = basis::shift<dim>(basis_, x, y);
            for (size_t i=0; i<nbasis; i++) {
              for (size_t j=0; j<nbasis; j++) {
                moment[i][j] += basis[i]*basis[j]*weight_val[iB];
              }
            }
            iB++;
          }
        }

        // Calculate inverse(P*W*transpose(P))*P*W
        iB = 0;

        for (auto const& particleB : source_particles) {
	        std::vector<double> pair_result(nbasis);
          Wonton::Point<dim> y = source_.get_particle_coordinates(particleB);
	        std::vector<double> basis = basis::shift<dim>(basis_, x, y);

          // recast as a Portage::Matrix
          Wonton::Matrix basis_matrix(nbasis,1);
          for (size_t i=0; i<nbasis; i++) basis_matrix[i][0] = basis[i];

          // solve the linear system
          if (not zilchit) {
            std::string error="check";
#ifdef HAVE_LAPACKE
            Wonton::Matrix pair_result_matrix = moment.solve(basis_matrix, "lapack-sytr", error);
#else
            Wonton::Matrix pair_result_matrix = moment.solve(basis_matrix, "inverse", error);
#endif
            for (size_t i=0; i<nbasis; i++)
              pair_result[i] = pair_result_matrix[i][0]*weight_val[iB];
              // NOTE: THIS CODE HAS TO WAIT FOR AN UPDATED WONTON
              //if (basis_matrix.is_singular() == 2 or error != "none") {
              //  bad_count++;
              //}
          } else {
            for (size_t i=0; i<nbasis; i++)
              pair_result[i] = 0.;
          }

          // If an operator is being applied, adjust final weights.
          if (estimate_ == OperatorRegression) {
            auto ijet = basis::inverse_jet<dim>(basis_, x);
            std::vector<std::vector<double>> basisop;
            oper::apply<dim>(operator_spec_, basis_,
                             operator_domain_[particleA],
                             operator_data_[particleA], basisop);
            int const num_basis = nbasis;
            int const opsize = oper::size_info(operator_spec_, basis_,
                                            operator_domain_[particleA])[0];
            std::vector<double> operator_result(opsize, 0.);

            for (int j = 0; j < opsize; j++) {
              for (int k = 0; k < num_basis; k++) {
                for (int m = 0; m < num_basis; m++) {
                  operator_result[j] += pair_result[k]*ijet[k][m]*basisop[m][j];
                }
              }
            }
            for (int j = 0; j < num_basis; j++)
              pair_result[j] = operator_result[j];
          }
          result.emplace_back(particleB, pair_result);
          iB++;
        }
	      break;
      }
      default:  // invalid estimate
	      assert(false);
    }
    return result;
  }

 private:
  SourceSwarm const& source_;
  TargetSwarm const& target_;
  EstimateType estimate_;
  WeightCenter center_;
  Wonton::vector<Weight::Kernel> const& kernels_;
  Wonton::vector<Weight::Geometry> const& geometries_;
  Wonton::vector<std::vector<std::vector<double>>> const& smoothing_;
  basis::Type basis_;
  oper::Type operator_spec_;
  Wonton::vector<oper::Domain> operator_domain_;
  Wonton::vector<std::vector<Wonton::Point<dim>>> operator_data_;
};

}}

#endif // ACCUMULATE_H_INC
