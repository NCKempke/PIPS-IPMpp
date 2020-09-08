/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP
 * Modified by Cosmin Petra to perform solves with the factors.
 */

#include "Ma57Solver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"
#include "DenseGenMatrix.h"
#include "pipsport.h"

#include <mpi.h>

void dumpdata(int* irow, int* jcol, double*M, int n, int nnz)
{
  printf("======================================================\n");
  for(int i = 0; i < nnz; i++)
     printf("%6d %6d %10.2f\n", irow[i], jcol[i], M[i]);
  printf("\n");
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}

Ma57Solver::Ma57Solver( SparseSymMatrix * sgm )
{
  irowM = 0;
  jcolM = 0;
  fact  = 0;
  ifact = 0;
  keep  = 0;

  iworkn  = nullptr; dworkn = nullptr;
  niworkn = ndworkn = 0;

  mStorage = sgm->getStorageHandle();
  init();
}

void Ma57Solver::init()
{
   freshFactor = false;

   ipessimism = 2;
   rpessimism = 2;

   assert( mStorage->n == mStorage->m );
   n = mStorage->n;
   M = mStorage->M;

   nnz = mStorage->numberOfNonZeros();

   FNAME(ma57id)( cntl, icntl );
   //icntl[1] = -1; // don't print warning messages
   icntl[8] = 10; // up to 10 steps of iterative refinement
   icntl[5] = 5; // 4 use Metis; 5 automatic choice(MA47 or Metis); 3 min
                 // degree ordering as in MA27; 2 use MC47;
   icntl[15] = 1;

   // set initial value of "Treat As Zero" parameter
   kTreatAsZero = 1.e-10; this->setTreatAsZero();

   // set initial value of Threshold parameter
   kThresholdPivoting = 1.e-5; this->setThresholdPivoting();

   // set the largest value of ThresholdPivoting parameter we are
   // willing to tolerate.
   kThresholdPivotingMax = 1.e-1;

   // set the increase factor for ThresholdPivoting parameter
   kThresholdPivotingFactor = 10.0;

   // set the required precision for each linear system solve
   kPrecision = 1.e-9;
}

void Ma57Solver::firstCall()
{
  irowM = new int[nnz];
  jcolM = new int[nnz];

  const int* krowM = mStorage->krowM;
  for( int i = 0; i < mStorage->n; i++ )
  {
     if( mStorage->fortranIndexed() )
     {
        assert(krowM[i] - 1 >= 0);
        for( int k = krowM[i] - 1; k < krowM[i + 1] - 1; k++ )
           irowM[k] = i + 1;
     }
     else
        for( int k = krowM[i]; k < krowM[i + 1]; k++ )
           irowM[k] = i + 1;
  }

  for( int k = 0; k < nnz; k++ )
  {
     if( !mStorage->fortranIndexed() )
        jcolM[k] = mStorage->jcolM[k] + 1;
     else
        jcolM[k] = mStorage->jcolM[k];
  }

  lkeep = 5 * n + nnz + 2 * std::max(n, nnz) + 42;
  keep = new int[lkeep];


  iworkn = new_iworkn( 5 * n );

  FNAME(ma57ad)( &n, &nnz, irowM, jcolM, &lkeep, keep, iworkn, icntl, info, rinfo );
  assert(info[0] >= 0);

  lfact = info[8];
  lfact = 2 * (int) (2 * lfact);
  fact  = new double[lfact];
  lifact = info[9];
  lifact = (int) (2 * lifact);
  ifact  = new int[lifact];
}

void Ma57Solver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void Ma57Solver::matrixChanged()
{
  if( !keep ) this->firstCall();

  assert(mStorage->n == mStorage->m);
  assert( n == mStorage->n);
  iworkn = new_iworkn(n);

#if 0
  do {
#endif

    FNAME(ma57bd)( &n, &nnz, M, fact, &lfact, ifact,
	     &lifact, &lkeep, keep, iworkn, icntl, cntl, info, rinfo );
    freshFactor = true;
#if 0
     int done = 0, tries = 0;;
    if( info[0] != 0 )
       std::cout << "ma57bd: Factorization: info[0]=: " << info[0] << std::endl;
    //assert(false);
    switch( info[0] ) {
    case 0: done = 1;
      break;
    case -3: {
      int ic = 0;
      int lnfact = (int) (info[16] * rpessimism);
      double * newfact = new double[lnfact];
      FNAME(ma57ed)( &n, &ic, keep, fact, &lfact, newfact, &lnfact,
	       ifact, &lifact, ifact, &lifact, info );
      delete [] fact;
      fact = newfact;
      lfact = lnfact;
      rpessimism *= 1.1;
      cout << "Resizing real part. pessimism = " << rpessimism << endl;
    }; break;
    case -4: {
      int ic = 1;
      int lnifact = (int) (info[17] * ipessimism);
      int * nifact = new int[ lnifact ];
      FNAME(ma57ed)( &n, &ic, keep, fact, &lfact, fact, &lfact,
	       ifact, &lifact, nifact, &lnifact, info );
      delete [] ifact;
      ifact = nifact;
      lifact = lnifact;
      ipessimism *= 1.1;
      cout << "Resizing int part. pessimism = " << ipessimism << endl;
    }; break;
    default:
      if( info[0] >= 0 ) done = 1;
      assert( info[0] >= 0 );
    } // end switch
    tries++;
  } while( !done );
  freshFactor = true;

  //delete [] iwork;
#endif
}

void Ma57Solver::solve( OoqpVector& rhs_in )
{
	int job = 0;
	if( freshFactor )
		icntl[8] = 1; // No iterative refinement
	else
		icntl[8] = 10; // Iterative refinement

	SimpleVectorHandle x( new SimpleVector(n) );
	SimpleVectorHandle resid( new SimpleVector(n) );
	SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);

	double * drhs = rhs.elements();
	double * dx   = x->elements();
	double * dresid = resid->elements();

	dworkn = new_dworkn(n * 5);
	iworkn = new_iworkn(n);

	int done = 0;
	int refactorizations = 0;
	int dontRefactor =  (kThresholdPivoting > kThresholdPivotingMax);
	while( !done && refactorizations < 10 )
	{
	   FNAME(ma57dd)( &job,       &n,        &nnz,   M,        irowM,   jcolM,
	         fact,       &lfact,    ifact,  &lifact,  drhs,    dx,
        dresid,      dworkn,    iworkn,  icntl,    cntl,    info,
        rinfo );

    done = 1;
    continue;
    if( resid->infnorm() < kPrecision*( 1 + rhs.infnorm() ) ) {
      // resids are fine, use them
      done = 1;
    }
    else if( !dontRefactor )
    {
      // resids aren't good enough.
      if( freshFactor ) { // We weren't doing iterative refinement,
        // let's do so
        job = 2;
        icntl[8] = 10;
        // Mark this factorization as stale
        freshFactor = false;
        // And grow more pessimistic about the next factorization
        if( kThresholdPivoting >= kThresholdPivotingMax ) {
          // We have already refactored as with a high a pivtol as we
          // are willing to use
          dontRefactor = 1;
        } else {
          // refactor with a higher Threshold Pivoting parameter
          kThresholdPivoting *= kThresholdPivotingFactor;
          if( kThresholdPivoting > kThresholdPivotingMax )
            kThresholdPivoting = kThresholdPivotingMax;
          this->setThresholdPivoting();
          cout << "Setting ThresholdPivoting parameter to "
            << kThresholdPivoting << " for future factorizations" << endl;
        }
      } else if ( dontRefactor ) {
        // We might have tried a refactor, but the pivtol is above our
        // limit.
        done = 1;
      } else {
        // Otherwise, we have already tried iterative refinement, and
        // have already increased the ThresholdPivoting parameter
        std::cout << "Refactoring with Threshold Pivoting parameter "
          << kThresholdPivoting << std::endl;
        this->matrixChanged();
        refactorizations++;
        // be optimistic about the next factorization
        job = 0;
        icntl[8] = 1;
      } // end else we hava already tried iterative refinement
    } // end else resids aren't good enough
  } // end while not done
rhs.copyFrom( *x );
}

Ma57Solver::~Ma57Solver()
{
   freeWorkingArrays();
}

void Ma57Solver::freeWorkingArrays()
{
   if( jcolM )
      delete[] jcolM;
   if( irowM )
      delete[] irowM;
   if( fact )
      delete[] fact;
   if( ifact )
      delete[] ifact;
   if( keep )
      delete[] keep;
   if( iworkn )
      delete[] iworkn;
   if( dworkn )
      delete[] dworkn;

   jcolM = irowM = ifact = keep = iworkn = nullptr;
   fact = dworkn = nullptr;
}

void Ma57Solver::solve(int solveType, OoqpVector& rhs_in)
{
  if( solveType < 1 || solveType > 4 )
    assert("Unknown JOB assigned for use in MA57CD!" && 0);
  else if( solveType == 1 )
  {
    // we prefer iterative refinement.
    solve(rhs_in);
    return;
  } /*else */

  int job = solveType; // Solve using A
  int one = 1;

  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);

  double * drhs   = rhs.elements();
  dworkn  = new_dworkn(n);
  iworkn = new_iworkn(n);

#ifdef HAVE_GETRUSAGE
  rusage before;
  getrusage( RUSAGE_SELF, &before );
#endif

  FNAME(ma57cd)( &job,       &n,
	   fact,       &lfact,    ifact,  &lifact,
	   &one,       drhs,      &n,
	   dworkn,      &n,        iworkn,
	   icntl,      info );
}

void Ma57Solver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N,NRHS;

  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize(NRHS,N);
  assert( n == N );

  // we need checks on the residuals, can't do that with multiple RHS
  for (int i = 0; i < NRHS; i++) {
    SimpleVector v(rhs[i],N);
    solve(v);
  }
}

void Ma57Solver::solve( int nrhss, double* rhss, int* colSparsity )
{
   for (int i = 0; i < nrhss; i++) {
     SimpleVector v(rhss + i * n, n);
     solve(v);
   }
}

int* Ma57Solver::new_iworkn(int dim)
{
   if( niworkn != dim )
   {
      if( iworkn )
         delete[] iworkn;
      iworkn = new int[dim];
      niworkn = dim;
   }
   else
   {
      if( nullptr == iworkn )
         iworkn = new int[dim];
   }
   return iworkn;
}

double* Ma57Solver::new_dworkn(int dim)
{
   if( ndworkn != dim )
   {
      if( dworkn )
         delete[] dworkn;
      dworkn = new double[dim];
      ndworkn = dim;
   }
   else
   {
      if( nullptr == dworkn )
         dworkn = new double[dim];
   }
   return dworkn;
}


/*
void Ma57Solver::Refine( OoqpVector& x_in, OoqpVector& rhs_in )
{
  int job=2; //calculate r=b-Ax, solve A(dx)=r, update solution and exit.

  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
  SimpleVector & x   = dynamic_cast<SimpleVector &>(x_in);

  double * drhs   = rhs.elements();
  double * dx     = x.elements();
  double * dresid = new double[n];
  dworkn  = new double[5*n];

  icntl[8]=2;//steps of iterative refinement

  int * iwork = new_iworkn(n);
  FNAME(ma57dd)( &job,       &n,        &nnz,   M,        irowM,   jcolM,
      fact,       &lfact,    ifact,  &lifact,  drhs,    dx,
      dresid,      dwork,    iwork,  icntl,    cntl,    info, rinfo );

  if(info[0]!=0)
     std::cout << "ma57dd: info[0]=: " << info[0] << std::endl;

  delete[] dwork; delete[] dresid;

}

void Ma57Solver::solve(GenMatrix& rhs_in)
{
   DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
   int N, NRHS;

   // rhs vectors are on the "rows", for continuous memory
   rhs.getSize(NRHS ,N);
   assert( n == N );

   int job = 1;
   const int BLOCKSIZE = 20;

   dworkn  = new_dworkn(n*BLOCKSIZE);
   int dworksize = n*BLOCKSIZE;
   iworkn = new_iworkn(n);

   for (int startcol = 0; startcol < NRHS; startcol += BLOCKSIZE)
   {
      double *drhs = rhs[startcol];
      int endcol = std::min( startcol + BLOCKSIZE, NRHS);
      int numcols = endcol - startcol;

      FNAME(ma57cd)( &job, &n, fact, &lfact, ifact, &lifact,
         &numcols, drhs, &n, dworkn, &dworksize, iworkn,
         icntl, info );
      assert(info[0] >= 0);
      if (info[0] > 0)
         printf("warning from ma57cd, info[0]=%d\n",info[0]);
   }
}

void Ma57Solver::Lsolve( OoqpVector& x )
{
  solve(2,x);
}

void Ma57Solver::Dsolve( OoqpVector& x )
{
  solve(3,x);
}

void Ma57Solver::Ltsolve( OoqpVector& x )
{
  solve(4,x);
}
*/
