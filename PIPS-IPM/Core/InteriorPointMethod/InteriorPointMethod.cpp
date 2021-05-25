/*
 * InteriorPointMethod.C
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#include "InteriorPointMethod.hpp"
#include "Variables.h"
#include "Residuals.h"
#include "AbstractLinearSystem.h"
#include "Problem.h"
#include "DistributedFactory.h"

InteriorPointMethod::InteriorPointMethod(DistributedFactory& factory, Problem& problem, MehrotraHeuristic mehrotra_heuristic, const Scaler* scaler)
: Solver(factory, problem),
mehrotra_strategy(MehrotraFactory::create(factory, problem, mehrotra_heuristic, scaler)), filter_line_search() {}

TerminationStatus InteriorPointMethod::solve(Problem& problem, Variables& iterate, Residuals& residuals) {
   // initialization of (x,y,z) and factorization routine.
   std::unique_ptr<AbstractLinearSystem> linear_system = factory.make_linear_system(problem);

   // solve the augmented linear system
   factory.iterate_started();
   this->solve_linear_system(iterate, problem, residuals, *step, *linear_system);
   factory.iterate_ended();
   // register the linear system to the step computation strategy
   mehrotra_strategy->register_observer(linear_system.get());

   TerminationStatus status_code;
   int iteration = 0;
   bool termination = false;
   while (!termination) {
      iteration++;
      // run Gondzio's multiple corrector scheme
      status_code = mehrotra_strategy->corrector_predictor(factory, problem, iterate, residuals, *step, *linear_system, iteration);
      if (status_code != NOT_FINISHED) {
         termination = true;
      }
   }
   return status_code;
}