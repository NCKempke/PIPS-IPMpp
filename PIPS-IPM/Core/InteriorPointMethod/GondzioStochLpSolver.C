/*
 * GondzioStochLpSolver.C
 *
 *  Created on: 20.12.2017
 *      Author: Daniel Rehfeldt, Svenja Uslu
 */

#include "GondzioStochLpSolver.h"
#include "Variables.h"
#include "Residuals.h"
#include "AbstractLinearSystem.h"
#include "Problem.h"
#include "DistributedFactory.h"

GondzioStochLpSolver::GondzioStochLpSolver(DistributedFactory& factory, Problem& problem, const Scaler* scaler) : Solver(factory, problem),
      mehrotra_heuristic(factory, problem, scaler) {}

TerminationCode GondzioStochLpSolver::solve(Problem& problem, Variables& iterate, Residuals& residuals) {
   // initialization of (x,y,z) and factorization routine.
   linear_system = factory.make_linear_system(problem);

   // solve the augmented linear system
   factory.iterate_started();
   solve_linear_system(iterate, problem, residuals, *step);
   factory.iterate_ended();

   // run Mehrotra's corrector predictor scheme
   TerminationCode status_code = mehrotra_heuristic.corrector_predictor_pd(factory, problem, iterate, residuals, *step, *linear_system);
   return status_code;
}

void GondzioStochLpSolver::add_monitor(Monitor* monitor) {
   this->mehrotra_heuristic.add_monitor(monitor);
}