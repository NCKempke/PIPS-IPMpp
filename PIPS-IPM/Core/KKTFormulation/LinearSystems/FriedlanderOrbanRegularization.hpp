//
// Created by nils-christian on 09.06.21.
//

#ifndef PIPSIPMPP_FRIEDLANDERORBANREGULARIZATION_HPP
#define PIPSIPMPP_FRIEDLANDERORBANREGULARIZATION_HPP

#include "RegularizationStrategy.h"

class FriedlanderOrbanRegularization : public RegularizationStrategy {
public:
   FriedlanderOrbanRegularization(unsigned int positive_eigenvalues_expected, unsigned int negaitve_eigenvalues_expected);

   [[nodiscard]] Regularization get_default_regularization() override;
   Regularization get_regularization_parameters(const Inertia&, double) override;

   ~FriedlanderOrbanRegularization() override = default;
private:
   [[nodiscard]] Regularization get_regularization_parameters();

   const double primal_increase_factor{100};
   const double primal_decrease_factor{1.0/10.0};
   const double dual_increase_factor{100};
   const double dual_decrease_factor{1.0/10.0};

   const double primal_regularization_minimum{1e-8};
   const double dual_regularization_minimum{1e-8};

   const double primal_regularization_initial{1.0};
   const double dual_equality_regularization_initial{1.0};
   const double dual_inequality_regularization_initial{1.0};
};

#endif //PIPSIPMPP_FRIEDLANDERORBANREGULARIZATION_HPP