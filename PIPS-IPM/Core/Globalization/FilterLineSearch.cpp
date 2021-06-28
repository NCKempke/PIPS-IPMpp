#include <PIPSIPMppSolver.hpp>
#include "FilterLineSearch.hpp"
#include "FilterStrategy.hpp"
#include "Problem.hpp"
#include "Residuals.h"
#include "Variables.h"
#include "PIPSIPMppOptions.h"

FilterLineSearch::FilterLineSearch(DistributedFactory& factory, Problem& problem, double dnorm, InteriorPointMethodType interior_point_method_type,
      const Scaler* scaler) :
      filter_strategy(FilterStrategy()),
      interior_point_method(MehrotraFactory::create(factory, problem, dnorm, interior_point_method_type, scaler)),
      scaler(scaler), verbose{PIPS_MPIgetRank() == 0 && pipsipmpp_options::get_bool_parameter("FILTER_VERBOSE")} {
}

void FilterLineSearch::initialize(Residuals& initial_residuals) {
   /* set the filter upper bound */
   this->filter_strategy.initialize(initial_residuals);
}

void FilterLineSearch::register_observer(AbstractLinearSystem* linear_system) {
   this->interior_point_method->register_observer(linear_system);
}

void FilterLineSearch::compute_acceptable_iterate(Problem& problem, Variables& iterate, Residuals& residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration) {
   // compute a direction
   //this->interior_point_method->corrector_predictor(problem, iterate, residuals, step, linear_system, iteration);
   bool small_corr = this->interior_point_method->compute_predictor_step(problem, iterate, residuals, step, linear_system, iteration);
   this->interior_point_method->compute_corrector_step(problem, iterate, residuals, step, linear_system, iteration, small_corr);
}

void FilterLineSearch::compute_acceptable_iterate(Problem& problem, Variables& current_iterate, Variables& direction, Residuals& current_residuals) {
   bool is_accepted = false;
   this->number_iterations = 0;
   double step_length = 1.;

   while (!this->termination_(is_accepted)) {
      this->number_iterations++;
      if (verbose) std::cout << "Line search current step length: " << step_length << "\n";
      // compute the trial iterate
      std::unique_ptr<Variables> trial_iterate = current_iterate.cloneFull();
      trial_iterate->saxpy(direction, step_length);

      // evaluate the residuals at the trial iterate
      std::unique_ptr<Residuals> trial_residuals = current_residuals.cloneFull();
      trial_residuals->evaluate(problem, *trial_iterate);
      trial_residuals->compute_residual_norm();

      const double predicted_reduction = PIPSIPMppSolver::predicted_reduction(problem, current_iterate, direction, step_length);
      if (verbose) std::cout << "Predicted reduction: " << predicted_reduction << "\n";
      /* check whether the trial step is accepted */
      is_accepted = this->filter_strategy.check_acceptance(current_residuals, *trial_residuals, predicted_reduction);
      if (is_accepted) {
         // if the trial iterate was accepted, overwrite current_iterate
         current_iterate.copy(*trial_iterate);
      }
      else {
         // decrease the step length
         step_length *= this->backtracking_ratio;
         if (verbose) std::cout << "LS trial iterate rejected\n";
      }
   }
   if (!is_accepted) {
      assert(false && "Enter restoration phase (not implemented yet)");
   }
}

bool FilterLineSearch::termination_(bool is_accepted) const {
   if (is_accepted) {
      return true;
   }
   else if (this->max_iterations < this->number_iterations) {
      return true;
   }
   return false;
}
