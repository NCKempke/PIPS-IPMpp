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

   auto [sum_non_zeros_columns, sum_non_zeros_equalities, sum_non_zeros_inequalities] = get_nonzero_vectors();

   auto [log_sum_columns, log_sum_equalities, log_sum_inequalities] = get_log_sum_vectors();


};


std::tuple<std::unique_ptr<Vector<int>>, std::unique_ptr<Vector<int>>, std::unique_ptr<Vector<int>>> CurtisReidScaler::get_nonzero_vectors() const{
   std::unique_ptr<Vector<int>> sum_non_zeros_equalities = problem_factory.make_equalities_dual_integral_vector();
   std::unique_ptr<Vector<int>> sum_non_zeros_inequalities = problem_factory.make_inequalities_dual_integral_vector();

   std::unique_ptr<Vector<int>> sum_non_zeros_columns = problem_factory.make_primal_integral_vector();

   sum_non_zeros_columns->setToZero();
   sum_non_zeros_equalities->setToZero();
   sum_non_zeros_inequalities->setToZero();

   this->A->getNnzPerRow(*sum_non_zeros_equalities);
   this->C->getNnzPerRow(*sum_non_zeros_inequalities);

   this->A->getNnzPerCol(*sum_non_zeros_columns);
   this->C->getNnzPerCol(*sum_non_zeros_columns);

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

