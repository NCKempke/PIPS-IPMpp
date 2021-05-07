/*
 * InteriorPointMethod.h
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef INTERIORPOINTMETHOD_H
#define INTERIORPOINTMETHOD_H

#include "Solver.h"
#include "MehrotraStrategy.hpp"

class Problem;

class Variables;

class DistributedFactory;

class Scaler;

class InteriorPointMethod : public Solver {
public:
   InteriorPointMethod(DistributedFactory& factory, Problem& problem, MehrotraHeuristic mehrotra_heuristic, const Scaler* scaler = nullptr);
   TerminationCode solve(Problem& problem, Variables& iterate, Residuals& residuals) override;
   virtual void add_monitor(Monitor* monitor) override;
   ~InteriorPointMethod() = default;

protected:
   MehrotraStrategy mehrotra_strategy;
};

#endif /* INTERIORPOINTMETHOD_H */
