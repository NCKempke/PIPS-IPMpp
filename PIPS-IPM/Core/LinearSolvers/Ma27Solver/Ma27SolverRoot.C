/*
 * Ma27SolverRoot.C
 *
 *  Created on: 08.09.2020
 *      Author: bzfkempk
 */


#include "Ma27SolverRoot.h"

#include "SimpleVector.h"
#include "SparseSymMatrix.h"

Ma27SolverRoot::Ma27SolverRoot( SparseSymMatrix * sgm, MPI_Comm mpiComm )
 : Ma27Solver(sgm), comm(mpiComm)
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
   if( PIPS_MPIgetRank() == 0 )
   {
      SparseSymMatrix& matrixNewSym = dynamic_cast<SparseSymMatrix&>(matrixNew);

      assert( matrixNewSym.getStorageRef().fortranIndexed() );
      printf("\n Schur complement factorization is starting ...\n ");

      freeWorkingArrays();
      mat_storage = matrixNewSym.getStorageHandle();

      init();
      matrixChanged();

      printf("\n Schur complement factorization completed \n ");
   }
}

void Ma27SolverRoot::matrixChanged()
{
   if( PIPS_MPIgetRank() == 0 )
      Ma27Solver::matrixChanged();
}

void Ma27SolverRoot::solve(OoqpVector& rhs)
{
   PIPSdebugMessage("MA27 solver: solve (single rhs) \n");

   SimpleVector& sv = dynamic_cast<SimpleVector &>(rhs);

   assert(n == rhs.length());

   if( PIPS_MPIgetRank() == 0 )
   {
      std::cout << "SchurSolve" << std::endl;
      Ma27Solver::solve(sv);
      std::cout << "SchurSolveDone" << std::endl;
   }

   if( PIPS_MPIgetSize() > 0 )
      MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, comm);
}
