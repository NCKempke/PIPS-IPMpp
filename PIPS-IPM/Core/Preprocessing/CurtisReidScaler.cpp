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
   auto least_squares_primal_residuals = problem_factory.make_primal_vector();
   auto least_squares_equality_residuals = problem_factory.make_equalities_dual_vector();
   auto least_squares_inequality_residuals = problem_factory.make_inequalities_dual_vector();

   auto tmp_vec_primal = problem_factory.make_primal_vector();
   auto tmp_vec_equalities = problem_factory.make_equalities_dual_vector();
   auto tmp_vec_inequalities = problem_factory.make_inequalities_dual_vector();

   tmp_vec_equalities->copyFrom(*log_sum_equalities);
   tmp_vec_equalities->componentDiv(*sum_non_zeros_equalities);

   tmp_vec_inequalities->copyFrom(*log_sum_inequalities);
   tmp_vec_inequalities->componentDiv(*sum_non_zeros_inequalities);

   least_squares_primal_residuals->copyFrom(*log_sum_columns);
   least_squares_equality_residuals->setToZero();
   least_squares_inequality_residuals->setToZero();

   auto transform_to_one = [](const auto& val) {
      return val == 0.0 ? 0.0 : 1.0;
   };

   A->transpose_mult_transform(1.0, *least_squares_primal_residuals, -1.0, *tmp_vec_equalities, transform_to_one);
   C->transpose_mult_transform(1.0, *least_squares_primal_residuals, -1.0, *tmp_vec_inequalities, transform_to_one);

   bool done{false};

   // e_{k-1}
   double e{0.0};
   // q_k
   double q{1.0};
   // s_k
   double s{0.0};

   s = least_squares_primal_residuals->scaled_dot_product_self(*sum_non_zeros_columns);
   int k;
   for (k = 0; !done && k < max_iter; ++k) {
      const bool even_iter = (k % 2 == 0);

      // compute r_{k+1}
      const double scale_last = -1.0 / q;
      const double scale_first = - e / q;

      if (even_iter) {
         // r_{k+1} = -1/q_k E N^-1 r_k - e_{k-1}/q_k r_{k-1}
         tmp_vec_primal->copyFrom(*least_squares_primal_residuals);
         tmp_vec_primal->componentDiv(*sum_non_zeros_columns);

         A->mult_transform(scale_first, *least_squares_equality_residuals, scale_last, *tmp_vec_primal, transform_to_one);
         C->mult_transform(scale_first, *least_squares_inequality_residuals, scale_last, *tmp_vec_primal, transform_to_one);
      } else {
         // r_{k+1} = -1/q_k E^T M^-1 r_k - e_{k-1}/q_k r_{k-1}
         tmp_vec_equalities->copyFrom(*least_squares_equality_residuals);
         tmp_vec_inequalities->copyFrom(*least_squares_inequality_residuals);

         tmp_vec_equalities->componentDiv(*sum_non_zeros_equalities);
         tmp_vec_inequalities->componentDiv(*sum_non_zeros_inequalities);

         A->transpose_mult_transform(scale_first, *least_squares_primal_residuals, scale_last, *tmp_vec_equalities, transform_to_one);
         C->transpose_mult_transform(1.0, *least_squares_primal_residuals, scale_last, *tmp_vec_inequalities, transform_to_one);
      }

      const double s_old = s;
      // compute s_{k+1}
      if (even_iter) {
         // s_{k+1} = r_{k+1}^T M^-1 r_{k+1}
         s = least_squares_equality_residuals->scaled_dot_product_self(*sum_non_zeros_equalities) +
            least_squares_inequality_residuals->scaled_dot_product_self(*sum_non_zeros_inequalities);
      } else {
         // s_k{+1} = r_{k+1}^T N^-1 r_{k+1}
         s = least_squares_primal_residuals->scaled_dot_product_self(*sum_non_zeros_columns);
      }

      // compute e_k = q_k s_{k+1} / s_{k}
      e = q * s / s_old;

      // q_{k+1}
      const double q_old = q;
      q = 1 - e;

      // update factors
      if (even_iter) {
         // c_{2m + 2} = c_{2m} + 1/(q_{2m} q{2m+1}) [ N^-1 r_{2m} + e_{2m-1} e_{2m-2}(c_{2m} - c_{2m-2} ]
      }
   }

   // final iteration
   if (k % 2 == 0) {

   } else {

   }
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

