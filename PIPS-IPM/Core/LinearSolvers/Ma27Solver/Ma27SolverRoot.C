/*
 * Ma27SolverRoot.C
 *
 *  Created on: 08.09.2020
 *      Author: bzfkempk
 */


#include "Ma27SolverRoot.h"

#include "SimpleVector.h"
#include "SparseSymMatrix.h"

Ma27SolverRoot::Ma27SolverRoot( SparseSymMatrix * sgm, bool solve_in_parallel, MPI_Comm mpiComm )
 : Ma27Solver(sgm), solve_in_parallel(solve_in_parallel), comm(mpiComm)
{
   threshold_pivoting_max = 0.1;
   max_n_iter_refinement = 10;
   precision = 1e-7;
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

      if( rank == 0 )
         printf("\n Schur complement factorization is starting ...\n ");

      freeWorkingArrays();
      mat_storage = matrixNewSym.getStorageHandle();

      init();
      matrixChanged();

      if( rank == 0 )
         printf("\n Schur complement factorization completed \n ");
   }
}

void Ma27SolverRoot::matrixChanged()
{
   if( solve_in_parallel || PIPS_MPIgetRank(comm) == 0 )
      Ma27Solver::matrixChanged();
}

void Ma27SolverRoot::solve(OoqpVector& rhs)
{
   PIPSdebugMessage("MA27 solver: solve (single rhs) \n");

   SimpleVector& sv = dynamic_cast<SimpleVector &>(rhs);

   assert(n == rhs.length());

   if( solve_in_parallel || PIPS_MPIgetRank(comm) == 0 )
      Ma27Solver::solve(sv);

   if( !solve_in_parallel && PIPS_MPIgetSize() > 0 )
      MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, comm);
}
