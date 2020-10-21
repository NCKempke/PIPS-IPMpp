/*
 * Ma27SolverRoot.h
 *
 *  Created on: 08.09.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MA27SOLVER_MA27SOLVERROOT_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MA27SOLVER_MA27SOLVERROOT_H_

#include "Ma27Solver.h"
#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "OoqpVectorHandle.h"
#include "pipsport.h"

#include "mpi.h"


/** implements linear solver class for root nodes that uses the MA57 solver
 */

class Ma27SolverRoot : public Ma27Solver
{

 public:
  Ma27SolverRoot( SparseSymMatrix * sgm, bool solve_in_parallel, MPI_Comm mpiComm = MPI_COMM_WORLD );

  ~Ma27SolverRoot();

  void matrixRebuild( DoubleMatrix& matrixNew ) override;
  void matrixChanged() override;

  using Ma27Solver::solve;

  void solve( OoqpVector& rhs ) override;

 private:
  const bool solve_in_parallel;
  const MPI_Comm comm;

  void factorize();
};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MA27SOLVER_MA27SOLVERROOT_H_ */
