//
// Created by bzfkempk on 27.04.21.
//

#ifndef PIPSIPMPP_REGULARIZATIONSTRATEGY_H
#define PIPSIPMPP_REGULARIZATIONSTRATEGY_H

#include <tuple>

using Inertia = std::tuple<unsigned int, unsigned int, unsigned int>;
using Regularization = std::tuple<double, double, double>;

class RegularizationStrategy {

public:
   RegularizationStrategy(unsigned int positive_eigenvalues_expected, unsigned int negaitve_eigenvalues_expected);

   [[nodiscard]] bool is_inertia_correct(const Inertia& inertia) const;
   virtual void notify_new_step();
   [[nodiscard]] virtual Regularization get_default_regularization() = 0;
   virtual Regularization get_regularization_parameters(const Inertia& inertia, double barrier_parameter) = 0;

   virtual ~RegularizationStrategy() = default;
protected:
   double primal_regularization_current{-1.0};
   double dual_equality_regularization_current{-1.0};
   double dual_inequality_regularization_current{-1.0};

   const unsigned int positive_eigenvalues_expected{};
   const unsigned int negative_eigenvalues_expected{};

   bool new_factorization{true};
};

#endif //PIPSIPMPP_REGULARIZATIONSTRATEGY_H
