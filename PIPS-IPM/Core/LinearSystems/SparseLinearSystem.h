/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENSPARSELINSYS
#define QPGENSPARSELINSYS

#include <ProblemFactory.h>
#include "LinearSystem.h"
#include "SparseSymMatrixHandle.h"

class DoubleLinearSolver;

/** 
 * implements the aspects of the solvers for sparse general QP
 * formulation that are specific to the sparse case.
 *
 * @ingroup QpGen 
 */
class SparseLinearSystem : public LinearSystem {
protected:
   SparseSymMatrixHandle Mat;
   DoubleLinearSolver* solver;
public:
   SparseLinearSystem(ProblemFactory* factory, Problem* problem, SparseSymMatrix* Mat, DoubleLinearSolver* solver);

   /** perform the actual solve using the factors produced in factor.
    *
    * @param rhs on input contains the aggregated right-hand side of
    * the augmented system; on output contains the solution in
    * aggregated form
    */
   void solveCompressed(OoqpVector& rhs) override;

   void putXDiagonal(const OoqpVector& xdiag) override;
   void putZDiagonal(const OoqpVector& zdiag) override;

   /** calls QpGenLinsys::factor to assemble the augmented system
    * matrix, then calls matrixChanged to factor it
    *
    * @see QpGenLinsys::factor
    */
   void factorize(Problem* problem, Variables* vars) override;

   virtual ~SparseLinearSystem();
};

#endif
