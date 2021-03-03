/* PIPS                                                               *
 * Authors: Miles Lubin                                               *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DESYMINDEFSOLVER2_H
#define DESYMINDEFSOLVER2_H

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "DenseStorageHandle.h"
#include "pipsport.h"

/** Specialized LDL^T solver for saddle point systems
 * @ingroup DenseLinearAlgebra
 * @ingroup LinearSolvers
 */
class DeSymIndefSolver2 : public DoubleLinearSolver {
public:
  DenseStorageHandle mStorage;
protected:
  int nx, ny, n;
public:
  DeSymIndefSolver2( const DenseSymMatrix * storage, int nx );
  void diagonalChanged( int idiag, int extent ) override;
  void matrixChanged() override;
  using DoubleLinearSolver::solve;
  void solve ( OoqpVector& vec ) override;
  virtual ~DeSymIndefSolver2();
};

#endif
