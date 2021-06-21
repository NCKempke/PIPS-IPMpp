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

/** comments indicate notation in Curtis Reid paper (_1 and _2 for equality and inequality contribution)
 * A.R. Curtis, J.K. Reid, On the automatic scaling of matrices for Gaussian elimination (1971)
 */
void CurtisReidScaler::scale() {
   create_scaling_vectors();

   const auto [sum_non_zeros_columns, sum_non_zeros_equalities, sum_non_zeros_inequalities] = get_nonzero_vectors(); // [N, M_1, M_2]

   const auto [log_sum_columns, log_sum_equalities, log_sum_inequalities] = get_log_sum_vectors(); // [tau, sigma_1, sigma_2]

   /// initialize
   scaling_factors_columns->setToZero(); // c

   scaling_factors_equalities->copyFrom(*log_sum_equalities); // p_1
   scaling_factors_equalities->componentDiv(*sum_non_zeros_equalities);

   scaling_factors_inequalities->copyFrom(*log_sum_inequalities); // p_2
   scaling_factors_inequalities->componentDiv(*sum_non_zeros_inequalities);

   /// calculate initial residuals
   auto least_squares_residuals = problem_factory.make_primal_vector();
   auto tmp_vec_equalities = problem_factory.make_equalities_dual_vector();
   auto tmp_vec_inequalities = problem_factory.make_inequalities_dual_vector();

   tmp_vec_equalities->copyFrom(*log_sum_equalities);
   tmp_vec_equalities->componentDiv(*sum_non_zeros_equalities);

   tmp_vec_inequalities->copyFrom(*log_sum_inequalities);
   tmp_vec_inequalities->componentDiv(*sum_non_zeros_inequalities);

   least_squares_residuals->copyFrom(*log_sum_columns);

   auto transform_to_one = [](const auto& val) {
      return val == 0.0 ? 0.0 : 1.0;
   };

   A->transpose_mult_transform(1.0, *least_squares_residuals, -1.0, *tmp_vec_equalities, transform_to_one);
   C->transpose_mult_transform(1.0, *least_squares_residuals, -1.0, *tmp_vec_inequalities, transform_to_one);
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

#ifndef NDEBUG
   auto positive_predicate = [](const auto& val) {
      return val > 0.0;
   };
   sum_non_zeros_columns->all_of(positive_predicate);
   sum_non_zeros_equalities->all_of(positive_predicate);
   sum_non_zeros_inequalities->all_of(positive_predicate);
#endif

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

