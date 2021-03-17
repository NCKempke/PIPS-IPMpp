/*
 * Ma57SolverRoot.h
 *
 *  Created on: 07.09.2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MA57SOLVER_MA57SOLVERROOT_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MA57SOLVER_MA57SOLVERROOT_H_

#include "Ma57Solver.h"
#include "DoubleLinearSolver.h"
#include "SparseSymMatrixHandle.h"
#include "OoqpVectorHandle.h"
#include "pipsport.h"

#include "mpi.h"


/** implements linear solver class for root nodes that uses the MA57 solver
 */

class Ma57SolverRoot : public Ma57Solver
{

 public:
  Ma57SolverRoot( SparseSymMatrix * sgm, bool solve_in_parallel, MPI_Comm mpiComm = MPI_COMM_WORLD, const std::string& name = "root" );

  ~Ma57SolverRoot() = default;

  void matrixRebuild( DoubleMatrix& matrixNew ) override;
  void matrixChanged() override;

  using Ma57Solver::solve;

  void solve( OoqpVector& rhs ) override;

 private:
  const bool solve_in_parallel;
  const MPI_Comm comm;

  void factorize();
};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MA57SOLVER_MA57SOLVERROOT_H_ */
