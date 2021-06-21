//
// Created by nils-christian on 18.06.21.
//

#include "CurtisReidScaler.h"
#include "ProblemFactory.h"
#include "AbstractMatrix.h"

CurtisReidScaler::CurtisReidScaler(const ProblemFactory& problem_factory, const Problem& problem, bool bitshifting) : Scaler(problem_factory, problem, bitshifting) {
   if (PIPS_MPIgetRank() == 0 && scaling_output)
      std::cout << "Creating CurtisReidScaler... bitshifting=" << bitshifting << "\n";
}

void CurtisReidScaler::scale() {
   create_scaling_vectors();

   const auto [sum_non_zeros_columns, sum_non_zeros_equalities, sum_non_zeros_inequalities] = get_nonzero_vectors();

   const auto [log_sum_columns, log_sum_equalities, log_sum_inequalities] = get_log_sum_vectors();

   /// initialize
   scaling_factors_columns->setToZero();

   scaling_factors_equalities->copyFrom(*log_sum_equalities);
   scaling_factors_equalities->componentDiv(*sum_non_zeros_equalities);

   scaling_factors_inequalitites->copyFrom(*log_sum_inequalities);
   scaling_factors_equalities->componentDiv(*sum_non_zeros_inequalities);

   /// calculate initial residuals
   auto least_squares_residuals = problem_factory.make_primal_vector();

   least_squares_residuals->copyFrom(*log_sum_columns);

};


std::tuple<std::unique_ptr<Vector<double>>, std::unique_ptr<Vector<double>>, std::unique_ptr<Vector<double>>> CurtisReidScaler::get_nonzero_vectors() const{
   std::unique_ptr<Vector<double>> sum_non_zeros_equalities = problem_factory.make_equalities_dual_vector();
   std::unique_ptr<Vector<double>> sum_non_zeros_inequalities = problem_factory.make_inequalities_dual_vector();

   std::unique_ptr<Vector<double>> sum_non_zeros_columns = problem_factory.make_primal_vector();

   sum_non_zeros_columns->setToZero();
   sum_non_zeros_equalities->setToZero();
   sum_non_zeros_inequalities->setToZero();

   auto to_one = [](const double& val) {
      return val != 0.0 ? 1 : 0;
   };

   this->A->sum_transform_rows(*sum_non_zeros_equalities, to_one);
   this->C->sum_transform_rows(*sum_non_zeros_inequalities, to_one);

   this->A->sum_transform_columns(*sum_non_zeros_columns, to_one);
   this->C->sum_transform_columns(*sum_non_zeros_columns, to_one);

   return {std::move(sum_non_zeros_columns), std::move(sum_non_zeros_equalities), std::move(sum_non_zeros_inequalities)};
}

std::tuple<std::unique_ptr<Vector<double>>, std::unique_ptr<Vector<double>>, std::unique_ptr<Vector<double>>> CurtisReidScaler::get_log_sum_vectors() const {
   std::unique_ptr<Vector<double>> log_sum_equalities = problem_factory.make_equalities_dual_vector();
   std::unique_ptr<Vector<double>> log_sum_inequalities = problem_factory.make_inequalities_dual_vector();

   std::unique_ptr<Vector<double>> log_sum_columns = problem_factory.make_primal_vector();

   log_sum_columns->setToZero();
   log_sum_equalities->setToZero();
   log_sum_inequalities->setToZero();

   auto two_log_if_nonzero = [](const double& val) {
      return val == 0.0 ? 0.0 : std::log2(val);
   };

   this->A->sum_transform_rows(*log_sum_equalities, two_log_if_nonzero);
   this->C->sum_transform_rows(*log_sum_inequalities, two_log_if_nonzero);

   this->A->sum_transform_columns(*log_sum_columns, two_log_if_nonzero);
   this->C->sum_transform_columns(*log_sum_columns, two_log_if_nonzero);

   return {std::move(log_sum_columns), std::move(log_sum_equalities), std::move(log_sum_inequalities)};
}

