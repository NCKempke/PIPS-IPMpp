/*
 * Ma57SolverRoot.C
 *
 *  Created on: 07.09.2020
 *      Author: bzfkempk
 */

#include "Ma57SolverRoot.h"

Ma57SolverRoot::Ma57SolverRoot(SparseSymmetricMatrix& sgm, bool solve_in_parallel, MPI_Comm mpiComm, std::string name) : Ma57Solver(sgm, std::move(name)),
      solve_in_parallel(solve_in_parallel), comm(mpiComm) {
   assert(mpiComm != MPI_COMM_NULL);
}

void Ma57SolverRoot::matrixRebuild(const AbstractMatrix& matrixNew) {
   assert(omp_get_thread_num() == 0);
   const int my_rank = PIPS_MPIgetRank(comm);

   if (solve_in_parallel || my_rank == 0) {
      const auto& matrixNewSym = dynamic_cast<const SparseSymmetricMatrix&>(matrixNew);

      assert(matrixNewSym.getStorage().fortranIndexed());

      mat_storage = &matrixNewSym.getStorage();
      nnz = matrixNewSym.numberOfNonZeros();
      lkeep = 7 * n + nnz + 2 * std::max(n, nnz) + 42;

      init();

      keep.clear();
      matrixChanged();
   }
}

void Ma57SolverRoot::matrixChanged() {
   assert(omp_get_thread_num() == 0);

   if (solve_in_parallel || PIPS_MPIgetRank(comm) == 0)
      Ma57Solver::matrixChanged();
}

void Ma57SolverRoot::solve(Vector<double>& rhs) {
   auto& sv = dynamic_cast<SimpleVector<double>&>(rhs);

   assert(n == rhs.length());

   if (solve_in_parallel || PIPS_MPIgetRank(comm) == 0)
      Ma57Solver::solve(sv);

   if (!solve_in_parallel && PIPS_MPIgetSize(comm) > 0)
      MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, comm);
}
