/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MA57SOLVER_MA57SOLVERROOT_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MA57SOLVER_MA57SOLVERROOT_H_

#include "DoubleLinearSolver.h"
#include "Ma57Solver.h"

#include "pipsport.h"
#include "SparseSymMatrixHandle.h"
#include "OoqpVectorHandle.h"

#include "mpi.h"


/** implements linear solver class for root nodes that uses the MA57 solver
 */

class Ma57SolverRoot : public Ma57Solver
{

 public:
  Ma57SolverRoot( SparseSymMatrix * sgm, MPI_Comm mpiComm = MPI_COMM_WORLD );

  ~Ma57SolverRoot();

  void matrixRebuild( DoubleMatrix& matrixNew ) override;
  void matrixChanged() override;

  using Ma57Solver::solve;

  void solve( OoqpVector& rhs ) override;

 private:
  const MPI_Comm comm;

  void factorize();
};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MA57SOLVER_MA57SOLVERROOT_H_ */
