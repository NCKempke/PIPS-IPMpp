/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef QPGENDELINSYS
#define QPGENDELINSYS

#include "QpGenLinsys.h"
#include "DenseSymMatrixHandle.h"
#include "OoqpVectorHandle.h"

class QpGen;
class QuadraticProblem;
class DoubleLinearSolver;

/** 
 * implements the aspects of the solvers for dense general QP
 * formulation that are specific to the dense case.
 *
 * @ingroup QpGen 
 */

class QpGenDenseLinsys : public QpGenLinsys {
protected:
  DenseSymMatrixHandle kkt;
  DoubleLinearSolver * solver;
public:
  QpGenDenseLinsys( QpGen * factory,
		 QuadraticProblem * data,
		 DenseSymMatrix * Mat,
		 DoubleLinearSolver * solver );

  /** perform the actual solve using the factors produced in factor.
   *
   * @param rhs on input contains the aggregated right-hand side of
   * the augmented system; on output contains the solution in
   * aggregated form
   */
  virtual void solveCompressed( OoqpVector& rhs );

  virtual void putXDiagonal( OoqpVector& xdiag );
  virtual void putZDiagonal( OoqpVector& zdiag );

  /** calls QpGenLinsys::factor to assemble the augmented system
   * matrix, then calls matrixChanged to factor it
   *
   * @see QpGenLinsys::factor 
   */
  virtual void factor(Problem *prob, Variables *vars);

  virtual ~QpGenDenseLinsys();
};

#endif










