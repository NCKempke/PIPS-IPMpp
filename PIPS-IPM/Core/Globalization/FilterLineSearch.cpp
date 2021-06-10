#include <InteriorPointMethod.hpp>
#include "FilterLineSearch.hpp"
#include "FilterStrategy.hpp"
#include "Problem.h"
#include "Residuals.h"
#include "Variables.h"

FilterLineSearch::FilterLineSearch(const Scaler* scaler, int max_iterations, double backtracking_ratio, double min_step_length) :
filter_strategy(FilterStrategy()), backtracking_ratio(backtracking_ratio), scaler(scaler), min_step_length(min_step_length),
max_iterations(max_iterations) {
}

void FilterLineSearch::initialize(Residuals& initial_residuals) {
   /* set the filter upper bound */
   this->filter_strategy.initialize(initial_residuals);
}

void FilterLineSearch::compute_acceptable_iterate(Problem& problem, Variables& current_iterate, Variables& direction, Residuals& current_residuals,
      double primal_step_length, double dual_step_length) {
   bool is_accepted = false;
   this->number_iterations = 0;

   while (!this->termination_(is_accepted)) {
      this->number_iterations++;
      if (verbose) std::cout << "Line search current step lengths: " << primal_step_length << ", " << dual_step_length << "\n";
      // compute the trial iterate
      Variables trial_iterate(current_iterate);
      trial_iterate.saxpy_pd(direction, primal_step_length, dual_step_length);

      // evaluate the residuals at the trial iterate
      Residuals trial_residuals(current_residuals);
      trial_residuals.evaluate(problem, trial_iterate);
      trial_residuals.recompute_residual_norm();

      double predicted_reduction = InteriorPointMethod::predicted_reduction(problem, direction, primal_step_length);

      /* check whether the trial step is accepted */
      is_accepted = this->filter_strategy.check_acceptance(current_iterate, current_residuals, trial_iterate, trial_residuals, predicted_reduction,
            primal_step_length);
      // if the trial iterate was accepted, current_iterate was overwritten
      if (is_accepted) {
         current_iterate.copy(trial_iterate);
      }
      else {
         /* decrease the step length */
         primal_step_length *= this->backtracking_ratio;
         dual_step_length *= this->backtracking_ratio;
         if (verbose) std::cout << "Trial iterate rejected, decreasing the step length to " << primal_step_length << ", " << dual_step_length << "\n";
      }
   }
   if (!is_accepted) {
      std::cout << "Enter restoration phase (not implemented yet)\n";
      assert(false);
   }
}

bool FilterLineSearch::termination_(bool is_accepted) {
   if (is_accepted) {
      return true;
   }
   else if (this->max_iterations < this->number_iterations) {
      return true;
   }
   return false;
}
