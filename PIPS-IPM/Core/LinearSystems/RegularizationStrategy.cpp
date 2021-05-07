//
// Created by bzfkempk on 27.04.21.
//

#include "RegularizationStrategy.h"
#include "DistributedOptions.h"
#include "pipsdef.h"
#include<cmath>

RegularizationStrategy::RegularizationStrategy(unsigned int positive_eigenvalues_expected_, unsigned int negative_eigenvalues_expected_)
      : positive_eigenvalues_expected{positive_eigenvalues_expected_}, negative_eigenvalues_expected{negative_eigenvalues_expected_},
      primal_regularization_absolute_minimum{pips_options::get_double_parameter("REGULARIZATION_MIN_PRIMAL")},
      primal_regularization_absolute_maximum{pips_options::get_double_parameter("REGULARIZATION_MAX_PRIMAL")} {
}

bool RegularizationStrategy::is_inertia_correct(const Inertia& inertia) const {
   auto[positive_eigenvalues, negative_eigenvalues, zero_eigenvalues] = inertia;

   std::cout << "comparing inertia of (" << positive_eigenvalues << "," << negative_eigenvalues << "," << zero_eigenvalues << ") against expected: (" <<
      positive_eigenvalues_expected << "," << negative_eigenvalues_expected << "," << "0)\n";
   return positive_eigenvalues_expected == positive_eigenvalues && negative_eigenvalues_expected == negative_eigenvalues && zero_eigenvalues == 0;
}

void RegularizationStrategy::notify_new_step() {
   new_factorization = true;

   if (primal_regularization_current != -1.0)
      primal_regularization_last = primal_regularization_current;

   primal_regularization_current = -1.0;
   dual_equality_regularization_current = -1.0;
   dual_inequality_regularization_current = -1.0;
};

Regularization RegularizationStrategy::get_regularization_parameters(const Inertia& inertia, double barrier_parameter) {
   if (new_factorization) {
      new_factorization = false;
      return get_regularization_new_matrix(inertia, barrier_parameter);
   }
   else {
      return get_regularization_nth_try(inertia, barrier_parameter);
   }
}

Regularization RegularizationStrategy::get_regularization_new_matrix(const Inertia& inertia, double barrier_parameter) {
   auto[positive_eigenvalues, negative_eigenvalues, zero_eigenvalues] = inertia;
   assert(positive_eigenvalues != positive_eigenvalues_expected || negative_eigenvalues != negative_eigenvalues_expected);
   assert(primal_regularization_current == -1.0 && dual_equality_regularization_current == -1.0 && dual_inequality_regularization_current == -1.0);

   if (zero_eigenvalues) {
      dual_equality_regularization_current = std::pow(barrier_parameter, barrier_exponent_dual);
   }
   else {
      dual_equality_regularization_current = 0.0;
   }

   if (primal_regularization_last == 0.0) {
      primal_regularization_current = primal_regularization_initial;
   }
   else {
      primal_regularization_current = std::max(primal_regularization_absolute_minimum,
            primal_regularization_decrease_factor * primal_regularization_last);
   }

   dual_inequality_regularization_current = dual_equality_regularization_current;

   std::cout << "returning suggested regularization : " << primal_regularization_current << " " << dual_equality_regularization_current << " " << dual_inequality_regularization_current << std::endl;
   return {primal_regularization_current, dual_equality_regularization_current, dual_inequality_regularization_current};
}

Regularization RegularizationStrategy::get_regularization_nth_try(const Inertia& inertia, double barrier_parameter) {
   // TODO : unused variable in release mode
   auto[positive_eigenvalues, negative_eigenvalues, zero_eigenvalues] = inertia;
   assert(positive_eigenvalues != positive_eigenvalues_expected || negative_eigenvalues != negative_eigenvalues_expected);

   if (primal_regularization_last == 0.0) {
      primal_regularization_current *= primal_regularization_increase_factor_initial;
   }
   else {
      primal_regularization_current *= primal_regularization_increase_factor;
   }

   PIPS_MPIabortIf(primal_regularization_current < primal_regularization_absolute_maximum,
         "ERROR : cannot factorize matrix after excessive error correction");

   return {primal_regularization_current, dual_equality_regularization_current, dual_inequality_regularization_current};
}

