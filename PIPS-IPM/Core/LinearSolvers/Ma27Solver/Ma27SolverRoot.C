/*
 * Ma27SolverRoot.C
 *
 *  Created on: 08.09.2020
 *      Author: bzfkempk
 */


#include "Ma27SolverRoot.h"

#include "SimpleVector.h"
#include "SparseSymMatrix.h"

Ma27SolverRoot::Ma27SolverRoot(const SparseSymMatrix* sgm, bool solve_in_parallel, MPI_Comm mpiComm, const std::string& name) : Ma27Solver(sgm, name),
      solve_in_parallel(solve_in_parallel), comm(mpiComm) {
   threshold_pivoting_max = 0.5;
   precision = 1e-7;
   assert(mpiComm != MPI_COMM_NULL);
}

void Ma27SolverRoot::matrixRebuild(DoubleMatrix& matrixNew) {
   const int my_rank = PIPS_MPIgetRank(comm);

   if (solve_in_parallel || my_rank == 0) {
      SparseSymMatrix& matrixNewSym = dynamic_cast<SparseSymMatrix&>(matrixNew);

      assert(matrixNewSym.getStorageRef().fortranIndexed());

      freeWorkingArrays();

      mat_storage = matrixNewSym.getStorageHandle();
      nnz = matrixNewSym.numberOfNonZeros();

      init();

      fact.clear();
      matrixChanged();
   }
}

void Ma27SolverRoot::matrixChanged() {
   if (solve_in_parallel || PIPS_MPIgetRank(comm) == 0)
      Ma27Solver::matrixChanged();
}

void Ma27SolverRoot::solve(OoqpVector& rhs) {
   SimpleVector<double>& sv = dynamic_cast<SimpleVector<double>&>(rhs);

   assert(n == rhs.length());

   if (solve_in_parallel || PIPS_MPIgetRank(comm) == 0)
      Ma27Solver::solve(sv);

   if (!solve_in_parallel && PIPS_MPIgetSize(comm) > 0)
      MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, comm);
}
