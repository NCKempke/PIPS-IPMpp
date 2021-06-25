#ifndef PIPSIPMPPSOLVER_H
#define PIPSIPMPPSOLVER_H

#include "Solver.hpp"
#include "InteriorPointMethodType.hpp"
#include "FilterLineSearch.hpp"

class Problem;

class Variables;

class DistributedFactory;

class Scaler;

class PIPSIPMppSolver : public Solver {
public:
   PIPSIPMppSolver(DistributedFactory& factory, Problem& problem, InteriorPointMethodType interior_point_method_type, const Scaler* scaler = nullptr);
   ~PIPSIPMppSolver() override = default;

   TerminationStatus solve(Problem& problem, Variables& iterate, Residuals& residuals) override;

   static double predicted_reduction(Problem& problem, Variables& iterate, Variables& direction, double step_length);

   [[nodiscard]] int n_iterations() const { return iteration; };
protected:
   bool verbose{false};
   const Scaler* scaler{};
   const int max_iterations{0};

   /** norm of problem data */
   const double dnorm{0.};
   /** norm of original unscaled problem */
   const double dnorm_orig{0.};

   FilterLineSearch filter_line_search;
   std::unique_ptr<Residuals> residuals_unscaled;

   /** history of values of mu obtained on all iterations to date */
   std::vector<double> mu_history{};
   /** history of values of residual norm obtained on all iterations to date */
   std::vector<double> residual_norm_history{};
   /** history of values of phi obtained on all iterations to date */
   std::vector<double> phi_history{};
   /** the i-th entry of this array contains the minimum value of phi encountered by the algorithm on or before iteration i */
   std::vector<double> phi_min_history{};

   /** iterations in last run */
   int iteration{-1};


   bool print_timestamp{true};
   double start_time{-1.};

   /** termination parameters */
   double mutol{1.e-6};
   double artol{1.e-4};

   std::pair<double, double> compute_unscaled_gap_and_residual_norm(const Problem& problem, const Residuals& residuals);
   void update_history(double duality_gap, double residual_norm, int iteration, double mu);
   TerminationStatus compute_status(double duality_gap, double residual_norm, int iteration, double mu);
   static double barrier_directional_derivative(Problem& problem, Variables& iterate, Variables& direction);
};

#endif /* PIPSIPMPPSOLVER_H */
