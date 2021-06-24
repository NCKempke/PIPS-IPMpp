//
// Created by nils-christian on 09.06.21.
//

#ifndef PIPSIPMPP_IPOPTREGULARIZATION_HPP
#define PIPSIPMPP_IPOPTREGULARIZATION_HPP

#include "RegularizationStrategy.h"

class IpoptRegularization : public RegularizationStrategy {
public:
   IpoptRegularization(unsigned int positive_eigenvalues_expected, unsigned int negaitve_eigenvalues_expected, MPI_Comm mpi_comm = MPI_COMM_WORLD);
   ~IpoptRegularization() override = default;

   [[nodiscard]] Regularization get_default_regularization() override;
   Regularization get_regularization_parameters(const Inertia& inertia, double barrier_parameter) override;
   void notify_new_step() override;
private:
   Regularization get_regularization_new_matrix(const Inertia& inertia, double barrier_parameter);
   Regularization get_regularization_nth_try(const Inertia& inertia, double barrier_parameter);

   constexpr static double barrier_exponent_dual{0.25};

   constexpr static double primal_regularization_initial{1e-4};
   constexpr static double primal_regularization_decrease_factor{1.0 / 3.0};
   constexpr static double primal_regularization_increase_factor_initial{100.0};
   constexpr static double primal_regularization_increase_factor{8.0};

   const double primal_regularization_absolute_minimum{};
   const double primal_regularization_absolute_maximum{};

   double primal_regularization_last{0.0};
};

#endif //PIPSIPMPP_IPOPTREGULARIZATION_HPP

