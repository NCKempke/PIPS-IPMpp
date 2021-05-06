/*
 * GondzioStochSolver.C
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#include "GondzioStochSolver.h"
#include "Variables.h"
#include "Residuals.h"
#include "AbstractLinearSystem.h"
#include "Problem.h"
#include "DistributedFactory.h"

GondzioStochSolver::GondzioStochSolver(DistributedFactory& factory, Problem& problem, const Scaler* scaler) : Solver(factory, problem),
      mehrotra_heuristic(factory, problem, scaler) {}

TerminationCode GondzioStochSolver::solve(Problem& problem, Variables& iterate, Residuals& residuals) {
   // initialization of (x,y,z) and factorization routine.
   linear_system = factory.make_linear_system(problem);

   // solve the augmented linear system
   factory.iterate_started();
   solve_linear_system(iterate, problem, residuals, *step);
   factory.iterate_ended();

   // run Mehrotra's corrector predictor scheme
   TerminationCode status_code = mehrotra_heuristic.corrector_predictor(factory, problem, iterate, residuals, *step, *linear_system);
   return status_code;
}

void GondzioStochSolver::add_monitor(Monitor* monitor) {
   this->mehrotra_heuristic.add_monitor(monitor);
}