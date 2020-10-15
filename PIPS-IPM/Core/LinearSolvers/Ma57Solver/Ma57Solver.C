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

extern int gOoqpPrintLevel;

#include <mpi.h>

void dumpdata(int* irow, int* jcol, double*M, int n, int nnz)
{
  printf("======================================================\n");
  for(int i = 0; i < nnz; i++)
     printf("%6d %6d %10.2f\n", irow[i], jcol[i], M[i]);
  printf("\n");
  printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}

Ma57Solver::Ma57Solver( SparseSymMatrix * sgm ) :
      n( sgm->getStorageRef().n)
{

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
   assert( n == mStorage->n );

   M = mStorage->M;

   nnz = mStorage->numberOfNonZeros();

   FNAME(ma57id)( cntl, icntl );

   icntl[1] = -1; // don't print warning messages
   icntl[8] = n_iterative_refinement;

   // 1 -> use pivot order in KEEP, 0 -> AMD ordering using MC47 (no dense rows), 2 -> AMD using MC47,
   // 3 -> min degree MA27, 4 -> METIS ordering, 5 -> auto choice (default)
   icntl[5] = 4;
   // default 1 -> scale using symmetrized MC64 version, if != 1 no sclaing
   //icntl[14] = 1;
   // if 0 nothing (default) if 1 remove small entries cntl[1] and place corresponding pivots at the end of factorization -> for highly rank deficient matrices ..
   icntl[15] = 0;

   // set initial value of "Treat As Zero" parameter
//   cntl[1] = treat_pivot_as_zero; // for Ma57BD
   //   cntl[2] convergence for Ma57DD iterative refinement (norm improvement by factor..) default 0.5
   //   cntl[3] Ma57BD (reuse of old factorization)
   //   cntl[4] static pivoting

   // set initial value of Threshold parameter
   cntl[0] = threshold_pivoting; // default is 1e-2 - only use if icntl[6] == 1 (default)
}

void Ma57Solver::firstCall()
{
  irowM = new int[nnz];
  jcolM = new int[nnz];

  // TODO : move somewhere .. same as in MA27
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

  int N = n;
  FNAME(ma57ad)( &N, &nnz, irowM, jcolM, &lkeep, keep, iworkn, icntl, info, rinfo );
  assert(info[0] >= 0);

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
   if( !keep )
      this->firstCall();

   assert(mStorage->n == mStorage->m);
   assert( n == mStorage->n);

   iworkn = new_iworkn(n);
   bool done = false;
   int tries = 0;

   do
   {
      int N = n;
      FNAME(ma57bd)( &N, &nnz, M, fact, &lfact, ifact,
            &lifact, &lkeep, keep, iworkn, icntl, cntl, info, rinfo );
      done = checkErrorsAndReact();
      ++tries;
   }
   while( !done && tries < max_tries );

   if ( !done && tries > max_tries )
   {
      std::cerr << "ERROR MA27: could not get factorization of matrix after max " << max_tries << " tries" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   freshFactor = true;
}

void Ma57Solver::solve( OoqpVector& rhs_in )
{
	int job = 0;

	if( freshFactor )
		icntl[8] = 1; // No iterative refinement
	else
		icntl[8] = n_iterative_refinement; // Iterative refinement

	SimpleVectorHandle x( new SimpleVector(n) );
	SimpleVectorHandle resid( new SimpleVector(n) );
	SimpleVector& rhs = dynamic_cast<SimpleVector &>(rhs_in);

	double * drhs = rhs.elements();
	double * dx   = x->elements();
	double * dresid = resid->elements();

	const double rhsnorm = rhs.infnorm();

	dworkn = new_dworkn(n * 5);
	iworkn = new_iworkn(n);

	bool done = false;
	int tries = 0;

	while( !done && tries < max_tries )
	{
	   int N = n;
      FNAME(ma57dd)(&job, &N, &nnz, M, irowM, jcolM, fact, &lfact, ifact,
            &lifact, drhs, dx, dresid, dworkn, iworkn, icntl, cntl, info,
            rinfo);

      done = checkErrorsAndReact();

      if( resid->infnorm() < precision * ( 1 + rhsnorm ) )
         done = true;
      else
      {
         if( freshFactor )
         {
            job = 2; // set to iterative refinement
            icntl[8] = 10;

            freshFactor = false;

            if( threshold_pivoting < threshold_pivoting_max )
            {
               // refactor with a higher Threshold Pivoting parameter
               threshold_pivoting = std::min( threshold_pivoting * threshold_pivoting_factor, threshold_pivoting_max);
               cntl[0] = threshold_pivoting;

               if( gOoqpPrintLevel >= ooqp_print_level_warnings )
                  std::cout << "Setting ThresholdPivoting parameter to " << threshold_pivoting << " for future factorizations" << endl;
            }
         }
         else if ( threshold_pivoting == threshold_pivoting_max )
            done = true;
         else
         {
            if( gOoqpPrintLevel >= ooqp_print_level_warnings )
               std::cout << "Refactoring with Threshold Pivoting parameter " << threshold_pivoting << std::endl;

            this->matrixChanged();

            job = 0;
            icntl[8] = 1;
            ++tries;
         }
      }
	}
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
   else if( solveType == 1)
   {
      solve(rhs_in);
   }
   else
   {
      int job = solveType;
      int one_rhs = 1;

      SimpleVector & rhs = dynamic_cast<SimpleVector &>(rhs_in);

      double * drhs = rhs.elements();
      dworkn = new_dworkn(n);
      iworkn = new_iworkn(n);

      int N = n;
      FNAME(ma57cd)(&job, &N, fact, &lfact, ifact, &lifact, &one_rhs, drhs, &N, dworkn, &N, iworkn, icntl, info);
   }
}

void Ma57Solver::solve(GenMatrix& rhs_in)
{
  DenseGenMatrix &rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);
  int N, NRHS;

  // rhs vectors are on the "rows", for continuous memory
  rhs.getSize( NRHS, N );
  assert( n == N );

  // we need checks on the residuals, can't do that with multiple RHS
  for (int i = 0; i < NRHS; i++)
  {
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

bool Ma57Solver::checkErrorsAndReact()
{
   bool error = false;
   const int error_flag = info[0];
   const int error_info = info[1];

   switch ( error_flag )
   {
      case 0 :
         break;
      case -1 :
      {
         std::cerr << "ERROR MA57: N out of range or < -1: " << error_info << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -2 :
      {
         std::cerr << "ERROR MA57: NNZ out of range or < -1 : " << error_info << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -3 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: insufficient space in fact: " << info[1] << " suggest reset to " << info[16] << std::endl;
         rpessimism *= 1.1;

         assert( fact );
         assert( keep );
         int lnew = std::max( static_cast<int>(1.1 * lfact), info[16] );
         double * newfac = new double[lnew];
         assert( lnew > lfact );
         int copy_real_array = 0;
         int dummy = lifact + 1;

         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << lnew << std::endl;

         int N = n;
         FNAME(ma57ed)( &N, &copy_real_array, keep, fact, &lfact, newfac, &lnew, ifact, &lifact,
               nullptr, &dummy, info );

         delete[] fact;
         fact = newfac;
         lfact = lnew;

         error = true;
      }; break;
      case -4 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: insufficient factorization space in ifact: " << info[1] << " suggest reset to " << info[17] << std::endl;;
         ipessimism *= 1.1;

         assert( ifact );
         assert( keep );
         int linew = std::max( static_cast<int>(1.1 * lifact), info[17] );
         int* newifac = new int[linew];
         assert( linew > lifact );
         int copy_int_array = 1;
         int dummy = lfact + 1;

         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << linew << std::endl;

         int N = n;
         FNAME(ma57ed)( &N, &copy_int_array, keep, fact, &lfact, nullptr, &dummy, ifact, &lifact,
               newifac, &linew, info );

         delete[] ifact;
         ifact = newifac;
         lifact = linew;

         error = true;
      } break;
      case -5:
      {
         std::cerr << "WARNING MA57: Small pivot found when pivoting disabled" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -6:
      {
         std::cerr << "ERROR M527: change of sign of pivots detected at stage even though matrix supposedly definiet" << error_info << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -7:
      {
         std::cerr << "ERROR MA57: value new int/real array not bigger than old value in MA57ED " << error_info << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -8:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: iterative refinement failed to converge" << std::endl;
         error = true;
      }; break;
      case -9:
      {
         std::cerr << "ERROR MA57: error in permutation array pos " << error_info << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -10:
      {
         std::cerr << "ERROR MA57: icntl[6] " << error_info << " is out of range" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -11:
      {
         std::cerr << "ERROR MA57: LRHS " << error_info << " smaller than n " << n << " on solve call" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -12:
      {
         std::cerr << "ERROR MA57: invalid value " << error_info << " for job" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -13:
      {
         std::cerr << "ERROR MA57: invalid value " << error_info << " for icntl[8]" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -14:
      {
         std::cerr << "ERROR MA57: MC71AD (icntl > 10) failed inside of MA57DD" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -15:
      {
         std::cerr << "ERROR MA57: lkeep " << error_info << " less than required" << error_info << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -16:
      {
         std::cerr << "ERROR MA57: NRHS " << error_info << " less than one"<< std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -17:
      {
         std::cerr << "ERROR MA57: lwork too small: " << error_info << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -18:
      {
         std::cerr << "ERROR MA57: METIS ordering requested but METIS was not linked" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case 1 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: detected " << info[2] << " entries out of range in irowM and jcolM; ignored" << std::endl;
      }; break;
      case 2 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout<< "WARNING M57: detected " << info[3] << " duplicate entries in user supplied matrix detected - summing them up" << std::endl;
      } break;
      case 3:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: out of range etries and duplicates detected .. ignoring/summing them up" << std::endl;
      }; break;
      case 4:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: rank deficient matrix detected; apparent rank is " << info[24] << std::endl;
      }; break;
      case 5:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: pivots have different sign when factorizing supposedly definite matrix; " << info[25] << " sign changes detected" << std::endl;
      }; break;
      case 8:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: inf norm of computed solution was zero" << std::endl;
      }; break;
      case 10:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: insufficient space in fact: " << lfact << std::endl;
         rpessimism *= 1.5;

         assert( fact );
         assert( keep );
         int lnew = static_cast<int>(1.5 * lfact);
         double * newfac = new double[lnew];
         assert( lnew > lfact );
         int copy_real_array = 0;
         int dummy = lifact + 1;

         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << lnew << std::endl;

         int N = n;
         FNAME(ma57ed)( &N, &copy_real_array, keep, fact, &lfact, newfac, &lnew, ifact, &lifact,
               nullptr, &dummy, info );

         delete[] fact;
         fact = newfac;
         lfact = lnew;

         error = true;
      }; break;
      case 11:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA57: insufficient factorization space in ifact: " << lifact << std::endl;
         ipessimism *= 1.5;

         assert( ifact );
         assert( keep );
         int linew = static_cast<int>(1.5 * lifact);
         int* newifac = new int[linew];
         assert( linew > lifact );
         int copy_int_array = 1;
         int dummy = lfact + 1;

         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << linew << std::endl;

         int N = n;
         FNAME(ma57ed)( &N, &copy_int_array, keep, fact, &lfact, nullptr, &dummy, ifact, &lifact,
               newifac, &linew, info );

         delete[] ifact;
         ifact = newifac;
         lifact = linew;

         error = true;
      }; break;
      default :
      {
         assert( error_flag == 0);
      }; break;
   }

   return error;
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
