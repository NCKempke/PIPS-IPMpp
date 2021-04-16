/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENSPARSELINSYS
#define QPGENSPARSELINSYS

#include "QpGenLinsys.h"
#include "SparseSymMatrixHandle.h"

class DoubleLinearSolver;

/** 
 * implements the aspects of the solvers for sparse general QP
 * formulation that are specific to the sparse case.
 *
 * @ingroup QpGen 
 */
class QpGenSparseLinsys : public QpGenLinsys {
protected:
  SparseSymMatrixHandle Mat;
  DoubleLinearSolver * solver;
public:
  QpGenSparseLinsys(  QpGen * factory,
		QP * data,
		SparseSymMatrix * Mat,
		DoubleLinearSolver * solver );

  /** perform the actual solve using the factors produced in factor.
   *
   * @param rhs on input contains the aggregated right-hand side of
   * the augmented system; on output contains the solution in
   * aggregated form
   */
  void solveCompressed( OoqpVector& rhs ) override;

  void putXDiagonal( OoqpVector& xdiag ) override;
  void putZDiagonal( OoqpVector& zdiag ) override;

  /** calls QpGenLinsys::factor to assemble the augmented system
   * matrix, then calls matrixChanged to factor it
   *
   * @see QpGenLinsys::factor 
   */
  void factor(Problem *prob, Variables *vars) override;
  
  virtual ~QpGenSparseLinsys();
};
#endif
