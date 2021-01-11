/*
 * Ma27SolverRoot.C
 *
 *  Created on: 08.09.2020
 *      Author: bzfkempk
 */


#include "Ma27SolverRoot.h"

#include "SimpleVector.h"
#include "SparseSymMatrix.h"

Ma27SolverRoot::Ma27SolverRoot( SparseSymMatrix * sgm, OoqpVector* regularization, const std::string& name, bool solve_in_parallel, MPI_Comm mpiComm )
 : Ma27Solver(sgm, regularization, name), solve_in_parallel(solve_in_parallel), comm(mpiComm)
{
   threshold_pivoting_max = 0.1;
   precision = 1e-10;
   assert(mpiComm != MPI_COMM_NULL);
}

Ma27SolverRoot::~Ma27SolverRoot()
{
}

void Ma27SolverRoot::matrixRebuild( DoubleMatrix& matrixNew )
{
   const int rank = PIPS_MPIgetRank( comm );

   if( solve_in_parallel || rank == 0 )
   {
      SparseSymMatrix& matrixNewSym = dynamic_cast<SparseSymMatrix&>(matrixNew);

      assert( matrixNewSym.getStorageRef().fortranIndexed() );

      freeWorkingArrays();
      mat_storage = matrixNewSym.getStorageHandle();

      init();
      matrixChanged();

   }
}

void Ma27SolverRoot::matrixChanged()
{
   if( PIPS_MPIgetRank() == 0 && omp_get_thread_num() == 0 )
      printf("\n Schur complement factorization is starting ...\n ");

   if( solve_in_parallel || PIPS_MPIgetRank(comm) == 0 )
      Ma27Solver::matrixChanged();

   if( PIPS_MPIgetRank() == 0 && omp_get_thread_num() == 0 )
      printf("\n Schur complement factorization completed \n ");

}

void Ma27SolverRoot::solve(OoqpVector& rhs)
{
   SimpleVector& sv = dynamic_cast<SimpleVector &>(rhs);

   assert(n == rhs.length());

   if( solve_in_parallel || PIPS_MPIgetRank(comm) == 0 )
      Ma27Solver::solve(sv);

   if( !solve_in_parallel && PIPS_MPIgetSize() > 0 )
      MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, comm);
}
