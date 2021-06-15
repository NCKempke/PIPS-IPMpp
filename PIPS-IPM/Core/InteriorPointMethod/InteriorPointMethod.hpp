#ifndef INTERIORPOINTMETHOD_H
#define INTERIORPOINTMETHOD_H

#include "Solver.hpp"
#include "MehrotraStrategy.hpp"

class Problem;

class Variables;

class DistributedFactory;

class Scaler;

class InteriorPointMethod : public Solver {
public:
   InteriorPointMethod(DistributedFactory& factory, Problem& problem, MehrotraStrategyType mehrotra_strategy_type, const Scaler* scaler = nullptr);
   TerminationStatus solve(Problem& problem, Variables& iterate, Residuals& residuals) override;
   static double predicted_reduction(Problem& problem, Variables& direction, double step_length);
   ~InteriorPointMethod() override;

protected:
   bool verbose{false};
   const Scaler* scaler{};
   int max_iterations{0};

   /** norm of problem data */
   double dnorm{0.};
   /** norm of original unscaled problem */
   double dnorm_orig{0.};

   std::unique_ptr<MehrotraStrategy> mehrotra_strategy;
   std::unique_ptr<Residuals> residuals_unscaled;

   /** history of values of mu obtained on all iterations to date */
   double* mu_history{};
   /** history of values of residual norm obtained on all iterations to date */
   double* residual_norm_history{};
   /** history of values of phi obtained on all iterations to date */
   double* phi_history{};
   /** the i-th entry of this array contains the minimum value of phi encountered by the algorithm on or before iteration i */
   double* phi_min_history{};

   bool print_timestamp{true};
   double start_time{-1.};

   /** termination parameters */
   double mutol{1.e-6};
   double artol{1.e-4};

   std::pair<double, double> compute_unscaled_gap_and_residual_norm(const Residuals& residuals);
   void update_history(double duality_gap, double residual_norm, int iteration, double mu);
   TerminationStatus compute_status(double duality_gap, double residual_norm, int iteration, double mu);
};

#endif /* INTERIORPOINTMETHOD_H */
