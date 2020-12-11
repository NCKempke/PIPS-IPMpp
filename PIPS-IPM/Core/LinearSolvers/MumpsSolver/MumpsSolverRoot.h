/*
 * MumpsSolverRoot.h
 *
 *  Created on: 25.06.2019
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_

#include "dmumps_c.h"

#include "mpi.h"
#include "MumpsSolverBase.h"
#include "SparseSymMatrix.h"
#include "OoqpVector.h"
#include "pipsport.h"


/** implements linear solver class for root nodes that uses the MUMPS solver
 */

class MumpsSolverRoot : public MumpsSolverBase {

 public:
  MumpsSolverRoot( SparseSymMatrix * sgm, bool solve_in_parallel );
  MumpsSolverRoot( MPI_Comm mpiComm, SparseSymMatrix * sgm, bool solve_in_parallel );

  ~MumpsSolverRoot();

  void matrixRebuild( DoubleMatrix& matrixNew ) override;
  void matrixChanged() override;

  using DoubleLinearSolver::solve;
  void solve( OoqpVector& rhs ) override;

 private:
  /* indicated wether every process solves or only rank 0 */
  const bool solve_in_parallel;
  void factorize();
};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_ */
