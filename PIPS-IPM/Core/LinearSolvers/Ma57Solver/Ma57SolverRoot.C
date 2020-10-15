/*
 * Ma57SolverRoot.C
 *
 *  Created on: 07.09.2020
 *      Author: bzfkempk
 */


#include "Ma57SolverRoot.h"

#include "SimpleVector.h"
#include "SparseSymMatrix.h"

Ma57SolverRoot::Ma57SolverRoot( SparseSymMatrix * sgm, MPI_Comm mpiComm )
 : Ma57Solver(sgm), comm(mpiComm)
{
   assert(mpiComm != MPI_COMM_NULL);
}

Ma57SolverRoot::~Ma57SolverRoot()
{
}

void Ma57SolverRoot::matrixRebuild( DoubleMatrix& matrixNew )
{
   if( PIPS_MPIgetRank() == 0 )
   {
      SparseSymMatrix& matrixNewSym = dynamic_cast<SparseSymMatrix&>(matrixNew);

      assert( matrixNewSym.getStorageRef().fortranIndexed() );
      printf("\n Schur complement factorization is starting ...\n ");

      freeWorkingArrays();
      mStorage = matrixNewSym.getStorageHandle();

      init();
      matrixChanged();

      printf("\n Schur complement factorization completed \n ");
   }
}

void Ma57SolverRoot::matrixChanged()
{
   if( PIPS_MPIgetRank() == 0 )
      Ma57Solver::matrixChanged();
}

void Ma57SolverRoot::solve(OoqpVector& rhs)
{
   PIPSdebugMessage("MA57 solver: solve (single rhs) \n");

   SimpleVector& sv = dynamic_cast<SimpleVector &>(rhs);

   assert(n == rhs.length());

   if( PIPS_MPIgetRank() == 0 )
      Ma57Solver::solve(sv);

   if( PIPS_MPIgetSize() > 0 )
      MPI_Bcast(sv.elements(), n, MPI_DOUBLE, 0, comm);
}