#include <PIPSIPMppSolver.hpp>
#include "FilterLineSearch.hpp"
#include "FilterStrategy.hpp"
#include "Problem.hpp"
#include "Residuals.h"
#include "Variables.h"
#include "PIPSIPMppOptions.h"

extern int print_level;

FilterLineSearch::FilterLineSearch(DistributedFactory& factory, Problem& problem, double dnorm, InteriorPointMethodType interior_point_method_type,
      const Scaler* scaler) :
      filter_strategy(FilterStrategy()),
      interior_point_method(MehrotraFactory::create(factory, problem, dnorm, interior_point_method_type, scaler)),
      verbose{PIPS_MPIgetRank() == 0 && pipsipmpp_options::get_bool_parameter("FILTER_VERBOSE")} {
}

void FilterLineSearch::initialize(Residuals& initial_residuals) {
   /* set the filter upper bound */
   this->filter_strategy.initialize(initial_residuals);
}

void FilterLineSearch::register_observer(AbstractLinearSystem* linear_system) {
   this->interior_point_method->register_observer(linear_system);
}

void FilterLineSearch::compute_acceptable_iterate(Problem& problem, Variables& current_iterate, Residuals& current_residuals, Variables& step,
      AbstractLinearSystem& linear_system, int iteration) {
   double mu = current_iterate.mu();
   // compute the predictor direction
   if (print_level >= 10) {
      this->print_statistics(problem, current_iterate, current_residuals, iteration, mu, TerminationStatus::NOT_FINISHED, 0);
   }
   bool small_corr = this->interior_point_method->compute_predictor_step(current_iterate, current_residuals, step, linear_system, iteration);
   // compute the correctors
   if (print_level >= 10) {
      this->print_statistics(problem, current_iterate, current_residuals, iteration, mu, TerminationStatus::NOT_FINISHED, 2);
   }
   this->interior_point_method->compute_corrector_step(problem, current_iterate, current_residuals, step, linear_system, iteration, small_corr);
   this->interior_point_method->take_step(current_iterate, step, 1.);

//   bool is_accepted = false;
//   this->number_iterations = 0;
//   double step_length = 1.;
//   while (!this->termination(is_accepted)) {
//      this->number_iterations++;
//      if (verbose) std::cout << "Line search current step length: " << step_length << "\n";
//      // compute the trial iterate
//      std::unique_ptr<Variables> trial_iterate = current_iterate.cloneFull();
//      this->interior_point_method->take_step(*trial_iterate, step, step_length);
//
//      // evaluate the residuals at the trial iterate
//      std::unique_ptr<Residuals> trial_residuals = current_residuals.cloneFull();
//      trial_residuals->evaluate(problem, *trial_iterate);
//      trial_residuals->compute_residual_norm();
//
//      const double predicted_reduction = PIPSIPMppSolver::predicted_reduction(problem, current_iterate, step, step_length);
//      if (verbose) std::cout << "Predicted reduction: " << predicted_reduction << "\n";
//
//      /* check whether the trial step is accepted */
//      is_accepted = this->filter_strategy.check_acceptance(current_residuals, *trial_residuals, predicted_reduction);
//      if (is_accepted) {
//         // if the trial iterate was accepted, overwrite current_iterate
//         current_iterate.copy(*trial_iterate);
//      }
//      else {
//         // decrease the step length
//         step_length *= this->backtracking_ratio;
//         if (verbose) std::cout << "LS trial iterate rejected\n";
//      }
//   }
//   if (!is_accepted) {
//      assert(false && "Enter restoration phase (not implemented yet)");
//   }
}

bool FilterLineSearch::termination(bool is_accepted) const {
   if (is_accepted) {
      return true;
   }
   else if (this->max_iterations < this->number_iterations) {
      return true;
   }
   return false;
}

void FilterLineSearch::print_statistics(const Problem& problem, const Variables& iterate, const Residuals& residuals, int i, double mu,
      TerminationStatus stop_code, int level) {
   this->interior_point_method->print_statistics(problem, iterate, residuals, i, mu, stop_code, level);
}