//
// Created by nils-christian on 18.06.21.
//

#ifndef PIPSIPMPP_CURTISREIDSCALER_H
#define PIPSIPMPP_CURTISREIDSCALER_H

#include "Scaler.hpp"
#include "DistributedVector.h"

class Problem;


class CurtisReidScaler : public Scaler {
private:
   std::unique_ptr<Vector<double>> temp_primal;
   std::unique_ptr<Vector<double>> temp_dual_equalities;
   std::unique_ptr<Vector<double>> temp_dual_inequalities;

   const double convergence_constant{1e-2};
   const int max_iter{10};

   PrimalDualTriplet get_nonzero_vectors() const;
   PrimalDualTriplet get_log_sum_vectors() const;

   void initialize_temp_vectors();
   void free_temp_vectors();
   void set_initial_scaling_factors(const Vector<double>& log_sum_equalities, const Vector<double>& log_sum_inequalities,
      const Vector<double>& sum_non_zeros_equalities, const Vector<double>& sum_non_zeros_inequalities);
   PrimalDualTriplet get_and_calculate_initial_residuals(const Vector<double>& log_sum_columns, const Vector<double>& log_sum_equalities,
      const Vector<double>& log_sum_inequalities, const Vector<double>& sum_non_zeros_equalities, const Vector<double>& sum_non_zeros_inequalities) const;
   void scaling_factors_to_power2();

protected:
   void scale_objective() const override;

public:

   CurtisReidScaler(const ProblemFactory& problem_factory, const Problem& problem, bool bitshifting = false);
   ~CurtisReidScaler() override = default;

   /** scale */
   void scale() override;
};

#endif //PIPSIPMPP_CURTISREIDSCALER_H
