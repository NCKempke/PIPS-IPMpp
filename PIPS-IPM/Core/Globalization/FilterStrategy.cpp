#include <iostream>
#include <cmath>
#include "FilterStrategy.hpp"
#include "Filter.hpp"
#include "Residuals.h"
#include "Variables.h"
#include "PIPSIPMppOptions.h"

FilterStrategy::FilterStrategy(FilterStrategyParameters& filter_strategy_parameters, FilterParameters& filter_parameters) :
filter(filter_parameters), parameters(filter_strategy_parameters), verbose{PIPS_MPIgetRank() == 0 && pipsipmpp_options::get_bool_parameter
      ("FILTER_VERBOSE")} {
}

FilterStrategy::FilterStrategy() : filter(), parameters({0.1, 0.999, 1e2, 1.25}), verbose{PIPS_MPIgetRank() == 0 && pipsipmpp_options::get_bool_parameter("FILTER_VERBOSE")} {
}

void FilterStrategy::initialize(Residuals& initial_residuals) {
   /* set the filter upper bound */
   double upper_bound = std::max(this->parameters.ubd, this->parameters.fact * initial_residuals.get_residual_norm());
   this->filter.upper_bound = upper_bound;
}

/* check acceptability of step(s) (filter & sufficient reduction)
 * precondition: feasible step
 * */
bool FilterStrategy::check_acceptance(Residuals& current_residuals, Residuals& trial_residuals, double predicted_reduction) {
   const double current_feasibility = current_residuals.feasibility_measure();
   const double current_optimality = current_residuals.optimality_measure();
   const double trial_feasibility = trial_residuals.feasibility_measure();
   const double trial_optimality = trial_residuals.optimality_measure();

   verbose = true;
   if (verbose) std::cout << "Filter strategy: feasibility " << trial_feasibility << " vs " << current_feasibility << "\n";
   if (verbose) std::cout << "Filter strategy: objective " << trial_optimality << " vs " << current_optimality << "\n";
   if (verbose) std::cout << this->filter << "\n";

   bool accept = false;
   /* check acceptance */
   bool filter_acceptable = this->filter.accept(trial_feasibility, trial_optimality);
   if (filter_acceptable) {
      if (verbose) std::cout << "Filter acceptable\n";
      // check acceptance wrt current x (h,f)
      filter_acceptable = this->filter.improves_current_iterate(current_feasibility, current_optimality, trial_feasibility, trial_optimality);
      if (filter_acceptable) {
         if (verbose) std::cout << "Current-iterate acceptable\n";
         double actual_reduction = current_optimality - trial_optimality;
         if (verbose) std::cout << "Filter predicted reduction: " << predicted_reduction << "\n";
         if (verbose) std::cout << "Filter actual reduction: " << actual_reduction << "\n";

         if (verbose) std::cout << "Switching condition: " << predicted_reduction << " >= " << this->parameters.Delta * std::pow(current_feasibility, 2) << " ?\n";
         if (verbose) std::cout << "Armijo condition: " << actual_reduction << " >= "
                   << this->parameters.Sigma * std::max(0., predicted_reduction - 1e-9) << " ?\n";

         /* switching condition: predicted reduction is not promising, accept */
         if (predicted_reduction < this->parameters.Delta * std::pow(current_feasibility, 2)) {
            this->filter.add(current_feasibility, current_optimality);
            accept = true;
         }
            /* Armijo sufficient decrease condition: predicted_reduction should be positive */
         else if (actual_reduction >= this->parameters.Sigma * std::max(0., predicted_reduction - 1e-9)) {
            accept = true;
         }
         else {
            if (verbose) std::cout << "Armijo condition not satisfied\n";
         }
      }
   }
   // TODO remove
   accept = true;
   return accept;
}
