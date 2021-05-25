#ifndef CGLINSYS_H
#define CGLINSYS_H

#include "DoubleLinearSolver.h"
#include "AbstractMatrix.h"
#include "SimpleVector.h"

class MatTimesVec;

/** implements the linear solver class using an implementation of the stabilized BiCG
 *
 * @ingroup LinearSolvers 
 */
class CGSolver : public DoubleIterativeLinearSolver {
public:
   /** build the class for the linear system
    *  for Ax=b preconditioned with M1 (or M1 and M2)
    */
   CGSolver(MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2 = nullptr);

   ~CGSolver() override;

   /** version of the main solve routine that takes argument as an
    * Vector<double>
    *
    * @param drhs on input contains the right-hand side; on output
    * contains the solution
    */
    using DoubleLinearSolver::solve;
   void solve(Vector<double>& rhs);

protected:
   double tol{};
   double iter{};
   int maxit{};
   int flag{};

   double* tmpVec1{};
   double* tmpVec2{};
   double* tmpVec3{};
   double* tmpVec4{};
   double* tmpVec5{};
   double* tmpVec6{};
   //int firstSolve;
};


#endif



