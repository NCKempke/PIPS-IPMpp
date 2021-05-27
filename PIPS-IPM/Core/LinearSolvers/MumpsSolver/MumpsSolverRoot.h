/*
 * MumpsSolverRoot.h
 *
 *  Created on: 25.06.2019
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_

#include "mpi.h"
#include "MumpsSolverBase.h"
#include "SparseSymmetricMatrix.h"
#include "Vector.hpp"
#include "pipsport.h"


/** implements linear solver class for root nodes that uses the MUMPS solver
 */

class MumpsSolverRoot : public MumpsSolverBase {

public:
   MumpsSolverRoot(const SparseSymmetricMatrix* sgm, bool solve_in_parallel);
   MumpsSolverRoot(MPI_Comm mpiComm, const SparseSymmetricMatrix* sgm, bool solve_in_parallel);

   ~MumpsSolverRoot() = default;

   void matrixRebuild( AbstractMatrix& matrixNew) override;
   void matrixChanged() override;

   using DoubleLinearSolver::solve;
   void solve(Vector<double>& rhs) override;

private:
   /* indicated whether every process solves or only rank 0 */
   const bool solve_in_parallel;
   void factorize();
};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERROOT_H_ */
