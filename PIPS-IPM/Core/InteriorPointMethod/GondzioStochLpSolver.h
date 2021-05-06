/*
 * GondzioStochLpSolver.h
 *
 *  Created on: 20.12.2017
 *      Author: Svenja Uslu
 */

#ifndef GONDZIOSTOCHLPSOLVER_H
#define GONDZIOSTOCHLPSOLVER_H

#include "GondzioStochSolver.h"

class Problem;

class Variables;

class DistributedFactory;

class GondzioStochLpSolver : public Solver {
public:
   GondzioStochLpSolver(DistributedFactory& factory, Problem& problem, const Scaler* scaler = nullptr);
   TerminationCode solve(Problem& problem, Variables& iterate, Residuals& residuals) override;
   virtual void add_monitor(Monitor* monitor) override;
   ~GondzioStochLpSolver() override = default;

protected:
   MehrotraHeuristic mehrotra_heuristic;
};


#endif /* GONDZIOSTOCHLPSOLVER_H */
