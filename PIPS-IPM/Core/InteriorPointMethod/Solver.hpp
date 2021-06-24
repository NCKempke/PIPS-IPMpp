#ifndef SOLVER_H
#define SOLVER_H

#include "TerminationStatus.hpp"
#include <memory>

class Problem;

class Variables;

class Residuals;

class AbstractLinearSystem;

class DistributedFactory;

class Statistics;

/**
 * Interior-point solver
 * @{
 */
class Solver {
public:
   Solver(DistributedFactory& factory, Problem& problem);
   virtual ~Solver() = default;

   /** implements the interior-point method for solving the subproblem */
   virtual TerminationStatus solve(Problem& problem, Variables& iterate, Residuals& residuals) = 0;
   /** solve the IPM system */
   virtual void solve_linear_system(Variables& iterate, Problem& problem, Residuals& residuals, Variables& step, AbstractLinearSystem& linear_system);

protected:
   DistributedFactory& factory;
   /**  storage for step vector */
   std::unique_ptr<Variables> step;
   //AbstractLinearSystem* linear_system{};
};

#endif

