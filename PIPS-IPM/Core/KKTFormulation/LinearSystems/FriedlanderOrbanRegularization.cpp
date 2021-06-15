//
// Created by nils-christian on 09.06.21.
//

#include "FriedlanderOrbanRegularization.hpp"
#include "PIPSIPMppOptions.h"

FriedlanderOrbanRegularization::FriedlanderOrbanRegularization(unsigned int positive_eigenvalues_expected,
   unsigned int negaitve_eigenvalues_expected, MPI_Comm mpi_comm) :
   RegularizationStrategy(positive_eigenvalues_expected, negaitve_eigenvalues_expected, mpi_comm),
   primal_regularization_initial{
      pipsipmpp_options::get_double_parameter("FRIEDLANDER_ORBAN_REGULARIZATION_INITIAL_PRIMAL")},
   dual_equality_regularization_initial{
      pipsipmpp_options::get_double_parameter("FRIEDLANDER_ORBAN_REGULARIZATION_INITIAL_DUAL_Y")},
   dual_inequality_regularization_initial{
      pipsipmpp_options::get_double_parameter("FRIEDLANDER_ORBAN_REGULARIZATION_INITIAL_DUAL_Z")} {
   primal_regularization_current = primal_regularization_initial / primal_decrease_factor;
   dual_equality_regularization_current = dual_equality_regularization_initial / dual_decrease_factor;
   dual_inequality_regularization_current = dual_inequality_regularization_initial / dual_decrease_factor;
}

Regularization FriedlanderOrbanRegularization::get_regularization_parameters(const Inertia&, double) {
   return get_regularization_parameters();
}

Regularization FriedlanderOrbanRegularization::get_regularization_parameters() {
   if (new_factorization) {
      new_factorization = false;
      primal_regularization_current *= primal_decrease_factor;
      primal_regularization_current = std::max(primal_regularization_minimum, primal_regularization_current);

      dual_equality_regularization_current *= dual_decrease_factor;
      dual_equality_regularization_current = std::max(dual_regularization_minimum,
         dual_equality_regularization_current);

      dual_inequality_regularization_current *= dual_decrease_factor;
      dual_inequality_regularization_current = std::max(dual_regularization_minimum,
         dual_inequality_regularization_current);

      return {primal_regularization_current, dual_equality_regularization_current, dual_inequality_regularization_current};
   } else {
      primal_regularization_current *= primal_increase_factor;
      dual_equality_regularization_current *= dual_increase_factor;
      dual_inequality_regularization_current *= dual_increase_factor;

      return {primal_regularization_current, dual_equality_regularization_current, dual_inequality_regularization_current};
   }
}

Regularization FriedlanderOrbanRegularization::get_default_regularization() {
   return get_regularization_parameters();
}
