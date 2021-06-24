//
// Created by nils-christian on 09.06.21.
//

#include "IpoptRegularization.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include "pipsdef.h"
#include "PIPSIPMppOptions.h"

IpoptRegularization::IpoptRegularization(unsigned int positive_eigenvalues_expected,
   unsigned int negaitve_eigenvalues_expected, MPI_Comm mpi_comm) :
   RegularizationStrategy(positive_eigenvalues_expected, negaitve_eigenvalues_expected, mpi_comm),
   primal_regularization_absolute_minimum{pipsipmpp_options::get_double_parameter("IPOPT_REGULARIZATION_MIN_PRIMAL")},
   primal_regularization_absolute_maximum{pipsipmpp_options::get_double_parameter("IPOPT_REGULARIZATION_MAX_PRIMAL")} {
}

void IpoptRegularization::notify_new_step() {
   RegularizationStrategy::notify_new_step();
   if (primal_regularization_current != -1.0) {
      primal_regularization_last = primal_regularization_current;
   }
}

Regularization IpoptRegularization::get_default_regularization() {
   return {0.0, 0.0, 0.0};
}

Regularization IpoptRegularization::get_regularization_parameters(const Inertia& inertia, double barrier_parameter) {
   if (new_factorization) {
      new_factorization = false;
      return get_regularization_new_matrix(inertia, barrier_parameter);
   }
   else {
      return get_regularization_nth_try(inertia, barrier_parameter);
   }
}

Regularization IpoptRegularization::get_regularization_new_matrix(const Inertia& inertia, double barrier_parameter) {
   const auto[positive_eigenvalues, negative_eigenvalues, zero_eigenvalues] = inertia;
   (void) positive_eigenvalues;
   (void) negative_eigenvalues;

   assert(positive_eigenvalues != positive_eigenvalues_expected || negative_eigenvalues != negative_eigenvalues_expected);

   if (zero_eigenvalues) {
      dual_equality_regularization_current = std::pow(barrier_parameter, barrier_exponent_dual);
   }
   else {
      dual_equality_regularization_current = 1e-4;
   }

   if (primal_regularization_last == 0.0) {
      primal_regularization_current = primal_regularization_initial;
   }
   else {
      primal_regularization_current = std::max(primal_regularization_absolute_minimum,
         primal_regularization_decrease_factor * primal_regularization_last);
   }

   dual_inequality_regularization_current = dual_equality_regularization_current;

   if (pipsipmpp_options::get_bool_parameter("REGULARIZATION_VERBOSE") && PIPS_MPIgetRank(mpi_comm) == 0) {
      std::cout << "IPOPT regularization: returning suggested regularization : " << primal_regularization_current
         << " " << dual_equality_regularization_current << " " << dual_inequality_regularization_current << "\n";
   }
   return {primal_regularization_current, dual_equality_regularization_current, dual_inequality_regularization_current};
}

Regularization IpoptRegularization::get_regularization_nth_try(const Inertia& inertia, double barrier_parameter) {
   // TODO : unused variable in release mode
   const auto[positive_eigenvalues, negative_eigenvalues, zero_eigenvalues] = inertia;
   (void) positive_eigenvalues; (void) negative_eigenvalues; (void) zero_eigenvalues;
   assert(positive_eigenvalues != positive_eigenvalues_expected || negative_eigenvalues != negative_eigenvalues_expected);

   if (primal_regularization_last == 0.0) {
      primal_regularization_current *= primal_regularization_increase_factor_initial;
   }
   else {
      primal_regularization_current *= primal_regularization_increase_factor;
   }

   PIPS_MPIabortIf(primal_regularization_current > primal_regularization_absolute_maximum,
      "ERROR : cannot factorize matrix after excessive error correction");

   if (pipsipmpp_options::get_bool_parameter("REGULARIZATION_VERBOSE") && PIPS_MPIgetRank(mpi_comm) == 0) {
      std::cout << "IPOPT regularization: returning suggested regularization : " << primal_regularization_current
         << " " << dual_equality_regularization_current << " " << dual_inequality_regularization_current << "\n";
   }

   return {primal_regularization_current, dual_equality_regularization_current, dual_inequality_regularization_current};
}
