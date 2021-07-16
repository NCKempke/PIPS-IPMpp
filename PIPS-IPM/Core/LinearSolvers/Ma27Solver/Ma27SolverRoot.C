/*
 * Ma27SolverRoot.C
 *
 *  Created on: 08.09.2020
 *      Author: bzfkempk
 */


#include "Ma27SolverRoot.h"

Ma27SolverRoot::Ma27SolverRoot(const SparseSymmetricMatrix& sgm, bool solve_in_parallel, MPI_Comm mpiComm, const std::string& name) : Ma27Solver(sgm, name),
      solve_in_parallel(solve_in_parallel), comm(mpiComm) {
   threshold_pivoting_max = 0.5;
   precision = 1e-7;
   assert(mpiComm != MPI_COMM_NULL);
}

void Ma27SolverRoot::matrixRebuild(const AbstractMatrix& matrixNew) {
   const int my_rank = PIPS_MPIgetRank(comm);

   if (solve_in_parallel || my_rank == 0) {
      auto& matrixNewSym = dynamic_cast<const SparseSymmetricMatrix&>(matrixNew);

      assert(matrixNewSym.getStorage().fortranIndexed());

      freeWorkingArrays();

      mat_storage = &matrixNewSym.getStorage();
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

void Ma27SolverRoot::solve(Vector<double>& rhs) {
   auto& sv = dynamic_cast<DenseVector<double>&>(rhs);

   assert(n == rhs.length());

   if (solve_in_parallel || PIPS_MPIgetRank(comm) == 0)
      Ma27Solver::solve(sv);

   if (!solve_in_parallel && PIPS_MPIgetSize(comm) > 0)
      MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, comm);
}
