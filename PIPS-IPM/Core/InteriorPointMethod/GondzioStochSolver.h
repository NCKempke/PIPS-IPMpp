/*
 * GondzioStochSolver.h
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef GONDZIOSTOCHSOLVER_H
#define GONDZIOSTOCHSOLVER_H

#include "Solver.h"
#include "MehrotraHeuristic.hpp"

class Problem;

class Variables;

class DistributedFactory;

class Scaler;

class GondzioStochSolver : public Solver {
public:
   GondzioStochSolver(DistributedFactory& factory, Problem& problem, const Scaler* scaler = nullptr);
   TerminationCode solve(Problem& problem, Variables& iterate, Residuals& residuals) override;
   virtual void add_monitor(Monitor* monitor) override;
   ~GondzioStochSolver() = default;

protected:
   MehrotraHeuristic mehrotra_heuristic;
};

#endif /* GONDZIOSTOCHSOLVER_H */
