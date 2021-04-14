#ifndef MEHALGORITHMSTOCH_H
#define MEHALGORITHMSTOCH_H

#include "MehrotraSolver.h"

class Problem;
class Variables;
class ProblemFormulation;

/** Derived class of Solver implementing the original Mehrotra
 *  predictor-corrector algorithm 
 * @ingroup QpSolvers
 */
class MehrotraStochSolver : public MehrotraSolver
{

public:
  MehrotraStochSolver( ProblemFormulation * opt, Problem * problem, const Scaler * scaler = nullptr );

  ~MehrotraStochSolver();

  int solve( Problem *problem, Variables *iterate, Residuals * resids ) override;

};

#endif

