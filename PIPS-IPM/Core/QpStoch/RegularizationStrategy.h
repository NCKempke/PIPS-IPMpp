//
// Created by bzfkempk on 27.04.21.
//

#ifndef PIPSIPMPP_REGULARIZATIONSTRATEGY_H
#define PIPSIPMPP_REGULARIZATIONSTRATEGY_H

#include <tuple>

using Inertia = std::tuple<unsigned int,unsigned int,unsigned int>;
using Regularization = std::tuple<double,double,double>;

class RegularizationStrategy {

public:
   RegularizationStrategy(unsigned int positive_eigenvalues_expected, unsigned int negaitve_eigenvalues_expected);

   bool is_inertia_correct(const Inertia& inertia) const;
   void notify_new_step();
   Regularization get_regularization_parameters(const Inertia& inertia, double barrier_parameter);

private:
   Regularization get_regularization_new_matrix(const Inertia& inertia, double barrier_parameter);
   Regularization get_regularization_nth_try(const Inertia& inertia, double barrier_parameter);

   constexpr static double barrier_exponent_dual{0.25};

   constexpr static double primal_regularization_initial{1e-4};
   constexpr static double primal_regularization_decrease_factor{1.0/3.0};
   constexpr static double primal_regularization_increase_factor_initial{100.0};
   constexpr static double primal_regularization_increase_factor{8.0};

   const unsigned int positive_eigenvalues_expected{};
   const unsigned int negative_eigenvalues_expected{};

   const double primal_regularization_absolute_minimum{};

   const double primal_regularization_absolute_maximum{};

   double primal_regularization_current{-1.0};
   double dual_equality_regularization_current{-1.0};
   double dual_inequality_regularization_current{-1.0};

   double primal_regularization_last{0.0};

   bool new_factorization{true};
};

#endif //PIPSIPMPP_REGULARIZATIONSTRATEGY_H
