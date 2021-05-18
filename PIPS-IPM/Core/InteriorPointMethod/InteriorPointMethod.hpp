/*
 * InteriorPointMethod.h
 *
 *  Created on: Dec 7, 2017
 *      Author: Daniel Rehfeldt
 */

#ifndef INTERIORPOINTMETHOD_H
#define INTERIORPOINTMETHOD_H

#include "Solver.hpp"
#include "MehrotraStrategy.hpp"
#include "FilterLineSearch.hpp"

class Problem;

class Variables;

class DistributedFactory;

class Scaler;

class InteriorPointMethod : public Solver {
public:
   InteriorPointMethod(DistributedFactory& factory, Problem& problem, MehrotraHeuristic mehrotra_heuristic, const Scaler* scaler = nullptr);
   TerminationStatus solve(Problem& problem, Variables& iterate, Residuals& residuals) override;
   virtual ~InteriorPointMethod() = default;

protected:
   std::unique_ptr<MehrotraStrategy> mehrotra_strategy;
   FilterLineSearch filter_line_search;
};

#endif /* INTERIORPOINTMETHOD_H */
