/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP
 * Modified by Cosmin Petra to perform solves with the factors.
 */

#include "Ma57Solver.h"

#include <algorithm>
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
  for(int i=0; i<nnz; i++)  printf("%6d %6d %10.2f\n", irow[i], jcol[i], M[i]);
  printf("\n");
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}

Ma57Solver::Ma57Solver( SparseSymMatrix * sgm ) :
      M(mStorage->M),
      n(mStorage->n),
      nnz(mStorage->krowM[n]),
      mStorage(sgm->getStorageHandle())
{
   /* set default parameters */
   FNAME(ma57id)( cntl, icntl );

   cntl[0] = threshold_pivoting; // default is 1e-2 - only use if icntl[6] == 1 (default)
   cntl[1] = treat_pivot_as_zero; // for Ma57BD

//   cntl[2] convergence for Ma57DD iterative refinement (norm improvement by factor..) default 0.5
//   cntl[3] Ma57BD (reuse of old factorization)
//   cntl[4] static pivoting

//icntl[1] = -1; // don't print warning messages
   // 1 -> use pivot order in KEEP, 0 -> AMD ordering using MC47 (no dense rows), 2 -> AMD using MC47,
   // 3 -> min degree MA27, 4 -> METIS ordering, 5 -> auto choice (default)
//   icntl[5] = 5;

   // numerical pivoting: 1 (default) use threshold in cntl[0], 2 -> no pivoting exit on sign change/zero,
   // 3 -> no pivoting exit on pivot < cntl[1], 4 -> no pivoting but alter matrix to have pivots of smae sign only...
//   icntl[6]
//   icntl[7]
   icntl[8] = n_iterative_refinement;
//   icntl[9]
//   icntl[10]
//   icntl[11]
//   icntl[12]
//   icntl[13]

   // default 1 -> scale using symmetrized MC64 version, if != 1 no sclaing
//   icntl[14] = 1;
   // if 0 nothing (default) if 1 remove small entries cntl[1] and place corresponding pivots at the end of factorization -> for highly rank deficient matrices ..
//   icntl[15] = 1

  // set the largest value of ThresholdPivoting parameter we are
  // willing to tolerate.
  threshold_pivoting_max = 1.e-1;

  // set the increase factor for ThresholdPivoting parameter
  threshold_pivoting_factor = 10.0;
}

void Ma57Solver::firstCall()
{
  irowM = new int[nnz];
  jcolM = new int[nnz];

  /* keep irowM, jcolM in case we solve with iterative refinement */
  int * krowM = mStorage->krowM;

  for( int i = 0; i < n; i++ )
  {
    for( int k = krowM[i]; k < krowM[i+1]; k++ ) {
      irowM[k] = i + 1;
    }
  }

  for( int k = 0; k < nnz; k++ )
    jcolM[k] = mStorage->jcolM[k] + 1;

  lkeep = 5 * n + nnz + 2 * std::max(n, nnz) + 42;
  keep = new int[lkeep];

  int * iwork = new int[5 * n];

  int nn = n;
  FNAME(ma57ad)( &nn, &nnz, irowM, jcolM, &lkeep, keep, iwork, icntl, info, rinfo );

  delete [] iwork;

  lfact = info[8];
  lfact = 2 * static_cast<int>(rpessimism * lfact);
  fact  = new double[lfact];

  lifact = info[9];
  lifact = static_cast<int>(ipessimism * lifact);
  ifact  = new int[lifact];

}
void Ma57Solver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void Ma57Solver::matrixChanged()
{
  if( keep == nullptr )
     this->firstCall();

  int* iwork = new int[n];

  bool done = false;
  int tries = 0;
  do
  {

    int nn = n;
    FNAME(ma57bd)( &nn, &nnz, M, fact, &lfact, ifact,
          &lifact, &lkeep, keep, iwork, icntl, cntl,
          info, rinfo );

    if( info[0] != 0 )
       std::cout << "ma57bd: Factorization: info[0]=: " << info[0] << std::endl;

    switch( info[0] )
    {
         case 0:
         {
            done = true;
            break;
         }
         case -3:
         {
            int ic = 0;
            int lnfact = (int) (info[16] * rpessimism);
            double *newfact = new double[lnfact];
            int nn = n;
            FNAME(ma57ed)(&nn, &ic, keep, fact, &lfact, newfact, &lnfact, ifact,
                  &lifact, ifact, &lifact, info);
            delete[] fact;
            fact = newfact;
            lfact = lnfact;
            rpessimism *= 1.1;
            std::cout << "Resizing real part. pessimism = " << rpessimism << std::endl;
            break;
         }
         case -4:
         {
            int ic = 1;
            int lnifact = (int) (info[17] * ipessimism);
            int * nifact = new int[ lnifact ];
            int nn = n;
            FNAME(ma57ed)( &nn, &ic, keep, fact, &lfact, fact, &lfact,
                  ifact, &lifact, nifact, &lnifact, info );
            delete [] ifact;
            ifact = nifact;
            lifact = lnifact;
            ipessimism *= 1.1;
            std::cout << "Resizing int part. pessimism = " << ipessimism << std::endl;
            break;
         }
         default:
         {
            if( info[0] >= 0 )
               done = true;
            assert( info[0] >= 0 );
         }
    }

    tries++;
  }
  while( !done );

  freshFactor = true;

  delete [] iwork;
}

void Ma57Solver::solve( OoqpVector& rhs_in )
{
	int job = 0;

	if( freshFactor )
	   icntl[8] = 1; // No iterative refinement
	else
		icntl[8] = 10; // Iterative refinement

	double* rhs = dynamic_cast<SimpleVector&>(rhs_in).elements();
	const double rhs_inf = rhs_in.infnorm();

	double* x = new double[n];
	double* resid = new double[n];
	double* work = new double[6 * n];

	int* iwork = new int[n];

	bool done = false;
	int refactorizations = 0;
	bool no_refactor =  (threshold_pivoting > threshold_pivoting_max);

	while( !done && refactorizations < 10 )
	{
	   int nn = n;
      FNAME(ma57dd)(&job, &nn, &nnz, M, irowM, jcolM, fact, &lfact, ifact,
            &lifact, rhs, x, resid, work, iwork, icntl, cntl, info, rinfo);

      const double res = *std::max_element(resid, resid + n, [](const double& a, const double& b){ return std::abs(a) < std::abs(b); } );

      if( res < precision * ( 1 + rhs_inf ) )
         done = true;
      else
      {
         if( freshFactor )
         {
            // Switch to iterative refinement
            job = 2;
            icntl[8] = 10;
            freshFactor = 0;

            if( threshold_pivoting >= threshold_pivoting_max )
               no_refactor = true;
            else
            {
               threshold_pivoting *= threshold_pivoting_factor;
               if( threshold_pivoting > threshold_pivoting_max )
                  threshold_pivoting = threshold_pivoting_max;

               this->setThresholdPivoting();
               std::cout << "Setting ThresholdPivoting parameter to " << threshold_pivoting << " for future factorizations" << std::endl;
            }
         }
         else if ( no_refactor )
         {
            done = true;
         }
         else
         {
            std::cout << "Refactoring with Threshold Pivoting parameter" << threshold_pivoting << std::endl;
            this->matrixChanged();
            refactorizations++;

            // be optimistic about the next factorization
            job = 0;
            icntl[8] = 1;
         }
      }
	}

   std::copy(x, x + n, rhs);

   delete[] iwork;
   delete[] work;
   delete[] resid;
   delete[] x;
}

/*void Ma57Solver::Refine( OoqpVector& x_in, OoqpVector& rhs_in )
{
  int job=2; //calculate r=b-Ax, solve A(dx)=r, update solution and exit.

  SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);
  SimpleVector & x   = dynamic_cast<SimpleVector &>(x_in);

  double * drhs   = rhs.elements();
  double * dx     = x.elements();
  double * dresid = new double[n];
  double * dwork  = new double[5*n];

  icntl[8]=2;//steps of iterative refinement

  int * iwork = new_iworkn(n);
  FNAME(ma57dd)( &job,       &n,        &nnz,   M,        irowM,   jcolM,
	   fact,       &lfact,    ifact,  &lifact,  drhs,    dx,
	   dresid,      dwork,    iwork,  icntl,    cntl,    info, rinfo );

  if(info[0]!=0) cout << "ma57dd: info[0]=: " << info[0] << endl;

  delete[] dwork; delete[] dresid;

  }*/

Ma57Solver::~Ma57Solver()
{
  delete [] irowM;
  delete [] jcolM;
  delete [] fact;
  delete [] ifact;
  delete [] keep;
}

/*void Ma57Solver::Lsolve( OoqpVector& x )
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


void Ma57Solver::solve(int solveType, OoqpVector& rhs_in)
{
   if( solveType < 1 || solveType > 4 )
      assert("Unknown JOB assigned for use in MA57CD!" && 0);
   else if( solveType == 1)
   {
      solve(rhs_in);
   }
   else
   {
      int job = solveType;
      int one_rhs = 1;

      double* rhs = dynamic_cast<SimpleVector &>(rhs_in).elements();
      double* dwork = new double[n];
      int* iwork = new int[n];

      int nn = n;
      FNAME(ma57cd)(&job, &nn, fact, &lfact, ifact, &lifact, &one_rhs, rhs, &nn, dwork, &nn, iwork, icntl, info);

      delete[] iwork;
      delete[] dwork;
   }
}

void Ma57Solver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N, NRHS;

  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize( NRHS, N );
  assert(n == N);


  // we need checks on the residuals, can't do that with multiple RHS
  for (int i = 0; i < NRHS; i++)
  {
     SimpleVector v(rhs[i],N);
     solve(v);
  }


//   int job = 1;

//   const int BLOCKSIZE = 20;

//   double * dwork  = new_dworkn(n*BLOCKSIZE);
//   int dworksize = n*BLOCKSIZE;
//   int * iwork     = new_iworkn(n);

//   for (int startcol = 0; startcol < NRHS; startcol += BLOCKSIZE) {
//     double *drhs = rhs[startcol];
//     int endcol = MIN(startcol+BLOCKSIZE,NRHS);
//     int numcols = endcol-startcol;
//     //cout << "MA57 multiple RHS" << endl;
//     FNAME(ma57cd)( &job,       &n,
// 		   fact,       &lfact,    ifact,  &lifact,
// 		   &numcols,       drhs,      &n,
// 		   dwork,      &dworksize,        iwork,
// 		   icntl,      info );
//     assert(info[0] >= 0);
//     if (info[0] > 0) {
//       printf("warning from ma57cd, info[0]=%d\n",info[0]);
//     }
//   }

}
