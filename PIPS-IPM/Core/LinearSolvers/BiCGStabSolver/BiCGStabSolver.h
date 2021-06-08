#ifndef BICGSTABLINSYS_H
#define BICGSTABLINSYS_H

#include "DoubleLinearSolver.h"
#include "AbstractMatrix.h"
#include "SimpleVector.h"

class MatTimesVec;

/** implements the linear solver class using an implementation of the stabilized BiCG
 *
 * @ingroup LinearSolvers 
 */
class BiCGStabSolver : public DoubleIterativeLinearSolver {
public:
   /** build the class for the linear system
    *  for Ax=b preconditioned with M1 (or M1 and M2)
    */
   BiCGStabSolver(MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2 = nullptr);

   ~BiCGStabSolver() override = default;;

   /** version of the main solve routine that takes argument as an
    * Vector<double>
    *
    * @param drhs on input contains the right-hand side; on output
    * contains the solution
    */
    using DoubleLinearSolver::solve;
   void solve(Vector<double>& rhs);

protected:
   BiCGStabSolver() = default;

   double tol{};
   double iter{};
   int maxit{};
   int flag{};
};


#endif



