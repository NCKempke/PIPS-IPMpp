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
}

Ma57SolverRoot::~Ma57SolverRoot()
{
}

void Ma57SolverRoot::matrixRebuild( DoubleMatrix& matrixNew )
{
   const int myrank = PIPS_MPIgetRank();

   if( myrank == 0 )
   {
      SparseSymMatrix& matrixNewSym = dynamic_cast<SparseSymMatrix&>(matrixNew);

      assert( matrixNewSym.getStorageRef().fortranIndexed() );
      printf("\n Schur complement factorization is starting ...\n ");

      mStorage = matrixNewSym.getStorageHandle();
      n = mStorage->n;
      M = mStorage->M;

      nnz = mStorage->krowM[n];

      if(keep)
         delete[] keep;
      if(irowM)
         delete [] irowM;
      if(jcolM)
         delete [] jcolM;
      if(fact)
         delete [] fact;
      if(ifact)
         delete [] ifact;

      ifact = jcolM = irowM = keep = nullptr;
      fact = nullptr;

      FNAME(ma57id)( cntl, icntl );
      //icntl[1] = -1; // don't print warning messages
      icntl[8] = 10; // up to 10 steps of iterative refinement
      icntl[5] = 5; // 4 use Metis; 5 automatic choice(MA47 or Metis); 3 min
             // degree ordering as in MA27; 2 use MC47;
      icntl[15] = 1;

      kTreatAsZero = 1.e-10; this->setTreatAsZero();

      // set initial value of Threshold parameter
      kThresholdPivoting = 1.e-5;
      this->setThresholdPivoting();

      matrixChanged();

      printf("\n Schur complement factorization completed \n ");
   }
}

void Ma57SolverRoot::matrixChanged()
{
   if( PIPS_MPIgetRank() == 0 )
   {
      if( !keep )
         this->firstCall();

      assert(mStorage->n == mStorage->m);
      assert( n == mStorage->n);
      iworkn = new_iworkn(n);

      FNAME(ma57bd)( &n, &nnz, M, fact, &lfact, ifact,
            &lifact, &lkeep, keep, iworkn, icntl, cntl, info, rinfo );
   }
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
