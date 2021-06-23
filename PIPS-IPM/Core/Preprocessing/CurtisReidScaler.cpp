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
   initialize_temp_vectors();

   auto [last_scaling_factors_columns, last_scaling_factors_equalities, last_scaling_factors_inequalities] =
      create_primal_dual_vector_triplet();

   const auto [sum_non_zeros_columns, sum_non_zeros_equalities, sum_non_zeros_inequalities] = get_nonzero_vectors(); // [N, M_1, M_2]
   const auto [log_sum_columns, log_sum_equalities, log_sum_inequalities] = get_log_sum_vectors(); // [tau, sigma_1, sigma_2]
   const int sum_non_zeros = sum_non_zeros_columns->inf_norm();

   /// initialize
   set_initial_scaling_factors(*log_sum_equalities, *log_sum_inequalities, *sum_non_zeros_equalities, *sum_non_zeros_inequalities);
   last_scaling_factors_equalities->copyFrom(*scaling_factors_equalities);
   last_scaling_factors_inequalities->copyFrom(*scaling_factors_inequalities);

   /// calculate initial residuals
   auto [least_squares_primal_residuals, least_squares_equality_residuals, least_squares_inequality_residuals] =
      get_and_calculate_initial_residuals(*log_sum_columns, *log_sum_equalities, *log_sum_inequalities, *sum_non_zeros_equalities,
         *sum_non_zeros_inequalities);

   auto transform_to_one = [](const auto& val) {
      return val == 0.0 ? 0.0 : 1.0;
   };

   bool done{false};

   // e_{k-1}, e_{k-2}
   double e_curr{0.0}, e_last{0.0}, e_lastlast{0.0};
   // q_k, q_{k+1}
   double q_curr{1.0}, q_next{0.0};
   // s_k+1, s_{k}
   double s_next{0.0}, s_curr{0.0};

   // s_0 = r_0^T N^-1 r_0
   s_curr = least_squares_primal_residuals->scaled_dot_product_self(*sum_non_zeros_columns);

   const double conv_tol = convergence_constant * sum_non_zeros;
   int k;
   for (k = 0; !done && k < max_iter; ++k) {
      const bool even_iter = (k % 2 == 0);

      // compute r_{k+1}
      const double scale_last = -1.0 / q_curr;
      const double scale_first = -e_last / q_curr;

      if (even_iter) {
         // N^-1 r_k
         temp_primal->copyFrom(*least_squares_primal_residuals);
         temp_primal->componentDiv(*sum_non_zeros_columns);

         // r_{k+1} = -e_{k-1}/q_k r_{k-1} -1/q_k E N^-1 r_k
         // least_squares_equality_residuals, least_squares_inequality_residuals contain r_{k-1}
         A->mult_transform(scale_first, *least_squares_equality_residuals, scale_last, *temp_primal, transform_to_one);
         C->mult_transform(scale_first, *least_squares_inequality_residuals, scale_last, *temp_primal, transform_to_one);
      } else {
         // M^-1 r_k
         temp_dual_equalities->copyFrom(*least_squares_equality_residuals);
         temp_dual_inequalities->copyFrom(*least_squares_inequality_residuals);

         temp_dual_equalities->componentDiv(*sum_non_zeros_equalities);
         temp_dual_inequalities->componentDiv(*sum_non_zeros_inequalities);

         // r_{k+1} = -e_{k-1}/q_k r_{k-1} -1/q_k E^T M^-1 r_k
         // least_squares_primal_residuals contains r_{k-1}
         A->transpose_mult_transform(scale_first, *least_squares_primal_residuals, scale_last, *temp_dual_equalities, transform_to_one);
         C->transpose_mult_transform(1.0, *least_squares_primal_residuals, scale_last, *temp_dual_inequalities, transform_to_one);
      }

      // compute s_{k+1}
      if (even_iter) {
         // s_{k+1} = r_{k+1}^T M^-1 r_{k+1}
         s_next = least_squares_equality_residuals->scaled_dot_product_self(*sum_non_zeros_equalities) +
            least_squares_inequality_residuals->scaled_dot_product_self(*sum_non_zeros_inequalities);
      } else {
         // s_k{+1} = r_{k+1}^T N^-1 r_{k+1}
         s_next = least_squares_primal_residuals->scaled_dot_product_self(*sum_non_zeros_columns);
      }

      // compute e_k = q_k s_{k+1} / s_{k}
      e_curr = q_curr * s_next / s_curr;

      // q_{k+1} = 1 - e_curr
      q_next = 1 - e_curr;

      // update one of the factors
      const double qq = 1.0 / q_curr * q_next;
      const double ee_qq = e_last * e_lastlast * qq;

      if (even_iter) {
         // k = 2m
         // c_{2m + 2} = c_{2m} + 1/(q_{2m} q{2m+1}) [ N^-1 r_{2m} + e_{2m-1} e_{2m-2}(c_{2m} - c_{2m-2}) ]
         //       = -ee_qq c_{2m-2} + (1 + ee_qq) c_{2m} + qq N^-1 r_{2m}

         last_scaling_factors_columns->scale(-1.0 * ee_qq);
         last_scaling_factors_columns->axpy(1 + ee_qq, *scaling_factors_columns);
         last_scaling_factors_columns->axpy(qq, *temp_primal);

         std::swap(last_scaling_factors_columns, scaling_factors_columns);
      } else {
         // k = 2m + 1
         // p_{2m + 3} = p_{2m + 1} + 1/(q_{2m + 1} q_{2m + 2}) [ M^-1 r_{2m+1} + e_{2m -1} e_{2m} (p_{2m + 1} - p_{2m - 1}) ]
         //       = -ee_qq p_{2m - 1} + (1 + ee_qq) p_{2m + 1} + qq M^-1 r_{2m+1}

         last_scaling_factors_equalities->scale(-1.0 * ee_qq);
         last_scaling_factors_inequalities->scale(-1.0 * ee_qq);

         last_scaling_factors_equalities->axpy(1.0 + ee_qq, *scaling_factors_equalities);
         last_scaling_factors_inequalities->axpy(1.0 + ee_qq, *scaling_factors_inequalities);

         last_scaling_factors_equalities->axpy(qq, *temp_dual_equalities);
         last_scaling_factors_inequalities->axpy(qq, *temp_dual_inequalities);

         std::swap(last_scaling_factors_equalities, scaling_factors_equalities);
         std::swap(last_scaling_factors_inequalities, scaling_factors_inequalities);
      }

      // update for next iteration
      q_curr = q_next;
      e_lastlast = e_last;
      e_last = e_curr;
      s_curr = s_next;

      done = s_next <= conv_tol;
   }

   // final iteration - getting the factors into phase
   // update one of the factors
   const double qq = 1.0 / q_curr;
   const double ee_qq = e_last * e_lastlast * qq;

   if (k % 2 == 0) {
      // N^-1 r_k
      temp_primal->copyFrom(*least_squares_primal_residuals);
      temp_primal->componentDiv(*sum_non_zeros_columns);

      last_scaling_factors_columns->scale(-1.0 * ee_qq);
      last_scaling_factors_columns->axpy(1 + ee_qq, *scaling_factors_columns);
      last_scaling_factors_columns->axpy(qq, *temp_primal);

      std::swap(last_scaling_factors_columns, scaling_factors_columns);
   } else {
      // M^-1 r_k
      temp_dual_equalities->copyFrom(*least_squares_equality_residuals);
      temp_dual_inequalities->copyFrom(*least_squares_inequality_residuals);

      last_scaling_factors_equalities->scale(-1.0 * ee_qq);
      last_scaling_factors_inequalities->scale(-1.0 * ee_qq);

      last_scaling_factors_equalities->axpy(1.0 + ee_qq, *scaling_factors_equalities);
      last_scaling_factors_inequalities->axpy(1.0 + ee_qq, *scaling_factors_inequalities);

      last_scaling_factors_equalities->axpy(qq, *temp_dual_equalities);
      last_scaling_factors_inequalities->axpy(qq, *temp_dual_inequalities);

      std::swap(last_scaling_factors_equalities, scaling_factors_equalities);
      std::swap(last_scaling_factors_inequalities, scaling_factors_inequalities);
   }

   free_temp_vectors();
};

void CurtisReidScaler::initialize_temp_vectors() {
   std::tie(temp_primal, temp_dual_equalities, temp_dual_inequalities) = create_primal_dual_vector_triplet();
}

void CurtisReidScaler::free_temp_vectors() {
   temp_primal = nullptr;
   temp_dual_equalities = nullptr;
   temp_dual_inequalities = nullptr;
}

void CurtisReidScaler::set_initial_scaling_factors(const Vector<double>& log_sum_equalities, const Vector<double>& log_sum_inequalities,
   const Vector<double>& sum_non_zeros_equalities, const Vector<double>& sum_non_zeros_inequalities) {

   scaling_factors_columns->setToZero(); // c

   // p_0 = M^{-1} \sigma
   scaling_factors_equalities->copyFrom(log_sum_equalities); // p_1
   scaling_factors_equalities->componentDiv(sum_non_zeros_equalities);

   scaling_factors_inequalities->copyFrom(log_sum_inequalities); // p_2
   scaling_factors_inequalities->componentDiv(sum_non_zeros_inequalities);
}

PrimalDualTriplet CurtisReidScaler::get_and_calculate_initial_residuals(const Vector<double>& log_sum_columns,
   const Vector<double>& log_sum_equalities, const Vector<double>& log_sum_inequalities,
   const Vector<double>& sum_non_zeros_equalities, const Vector<double>& sum_non_zeros_inequalities) const {

   auto [least_squares_primal_residuals, least_squares_equality_residuals, least_squares_inequality_residuals] =
      create_primal_dual_vector_triplet();

   // M^{-1} \sigma
   temp_dual_equalities->copyFrom(log_sum_equalities);
   temp_dual_equalities->componentDiv(sum_non_zeros_equalities);

   temp_dual_inequalities->copyFrom(log_sum_inequalities);
   temp_dual_inequalities->componentDiv(sum_non_zeros_inequalities);

   least_squares_primal_residuals->copyFrom(log_sum_columns);

   auto transform_to_one = [](const auto& val) {
      return val == 0.0 ? 0.0 : 1.0;
   };

   // r_0 = \tau - E^T M^-1 \sigma
   A->transpose_mult_transform(1.0, *least_squares_primal_residuals, -1.0, *temp_dual_equalities, transform_to_one);
   C->transpose_mult_transform(1.0, *least_squares_primal_residuals, -1.0, *temp_dual_inequalities, transform_to_one);

   return {std::move(least_squares_primal_residuals), std::move(least_squares_equality_residuals), std::move(least_squares_inequality_residuals)};
}

PrimalDualTriplet CurtisReidScaler::get_nonzero_vectors() const{
   auto [sum_non_zeros_columns, sum_non_zeros_equalities, sum_non_zeros_inequalities] = create_primal_dual_vector_triplet();

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

PrimalDualTriplet CurtisReidScaler::get_log_sum_vectors() const {
   auto [log_sum_columns, log_sum_equalities, log_sum_inequalities] = create_primal_dual_vector_triplet();

   auto two_log_if_nonzero = [](const double& val) {
      return val == 0.0 ? 0.0 : std::log2(val);
   };

   this->A->sum_transform_rows(*log_sum_equalities, two_log_if_nonzero);
   this->C->sum_transform_rows(*log_sum_inequalities, two_log_if_nonzero);

   this->A->sum_transform_columns(*log_sum_columns, two_log_if_nonzero);
   this->C->sum_transform_columns(*log_sum_columns, two_log_if_nonzero);

   return {std::move(log_sum_columns), std::move(log_sum_equalities), std::move(log_sum_inequalities)};
}

void CurtisReidScaler::scale_objective() const {
   assert(scaling_factors_columns);
   // no "real" scaling for the objective implemented right now
   obj->componentMult(*scaling_factors_columns);
}