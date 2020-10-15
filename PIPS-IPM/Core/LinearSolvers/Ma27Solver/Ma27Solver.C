/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "Ma27Solver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"

#include <fstream>

extern int gOoqpPrintLevel;

Ma27Solver::Ma27Solver(const SparseSymMatrix* sgm) :
       irowM(nullptr), jcolM(nullptr), fact(nullptr), mat(sgm), mat_storage(sgm->getStorageHandle())
{
   init();
}

void Ma27Solver::init()
{
   const double default_small_pivot = 1.0e-10;
   const double default_threshold_pivoting = 0.01;
   /* detecting dense rows during the factorization to preserve sparsity */
   const double default_fratio = 0.5;
   assert( mat_storage->n == mat_storage->m );
   n = mat_storage->n;
   nnz = mat_storage->numberOfNonZeros();

   FNAME(ma27id)(icntl, cntl);

   this->setThresholdPivoting( default_threshold_pivoting );
   cntl[1] = default_fratio;
   this->setSmallPivot( default_small_pivot );

   icntl[0] = 0;
   icntl[1] = 0;
}

void Ma27Solver::firstCall()
{
  irowM = new int[nnz];
  jcolM = new int[nnz];

  this->getIndices( irowM, jcolM );

  liw = static_cast<int>(ipessimism * (2 * nnz + 3 * n + 1));
  iw = new int[liw];
  iw1 = new int[2 * n];
  ikeep = new int[3 * n];

  int iflag = 0; // set to 1 if ikeep contains pivot order
  double ops;

  bool done = false;
  int tries = 0;
  do
  {
     FNAME(ma27ad)( &n, &nnz, irowM, jcolM, iw, &liw, ikeep, iw1, &nsteps, &iflag, icntl, cntl, info, &ops);
     done = !checkErrorsAndReact();
     ++tries;
  }
  while( !done && tries < max_tries );

  if ( !done && tries > max_tries )
  {
     std::cerr << "ERROR MA27: could not get ordering of matrix after max " << max_tries << " tries" << std::endl;
     MPI_Abort(MPI_COMM_WORLD, -1);
  }

  delete [] iw;
  delete [] iw1;

  la = 2 *  this->minimumRealWorkspace();
  fact = new double[la];

  // set iw and iw1 in prep for calls to ma27bd and ma27cd
  liw = 2 *  this->minimumIntWorkspace();
  iw = new int[liw];
  iw1 = new int[n];
}  

void Ma27Solver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

void Ma27Solver::matrixChanged()
{
   if( !fact )
      this->firstCall();

   bool done = false;
   int tries = 0;
   do
   {
      // copy M to fact
      this->copyMatrixElements(fact, la);

      FNAME(ma27bd)(&n, &nnz, irowM, jcolM, fact, &la, iw, &liw, ikeep, &nsteps,
            &maxfrt, iw1, icntl, cntl, info);

      done = !checkErrorsAndReact();
      tries++;
   }
   while( !done && tries < max_tries );

   if ( !done && tries > max_tries )
   {
      std::cerr << "ERROR MA27: could not get factorization of matrix after max " << max_tries << " tries" << std::endl;
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   iw2 = new int[nsteps];
   w = new double[maxfrt];
}

void Ma27Solver::solve( int nrhss, double* rhss, int* colSparsity )
{
   for (int i = 0; i < nrhss; i++) {
     SimpleVector v(rhss + i * n, n);
     solve(v);
   }
}

void Ma27Solver::solve( OoqpVector& rhs_in )
{
   SimpleVector &rhs = dynamic_cast<SimpleVector&>(rhs_in);

   // define structures to save rhs and store residuals
   SimpleVectorHandle iter(new SimpleVector(n));
   iter->setToZero();

   SimpleVectorHandle best_iter(new SimpleVector(n));
   double best_resid = std::numeric_limits<double>::infinity();

   SimpleVectorHandle residual(new SimpleVector(n));
   residual->copyFrom( rhs );

   const double rhsnorm = rhs.twonorm();

   bool done = false;
   int n_iter_ref = 0;

   double rnorm = -1.0;
   double res_last = -1.0;
   /* iterative refinement loop */
   while( !done && n_iter_ref < max_n_iter_refinement )
   {
      /* solve Ax = residual */
      FNAME(ma27cd)(&n, fact, &la, iw, &liw, w, &maxfrt, residual->elements(), iw1,
            &nsteps, icntl, info);
      iter->axpy(1.0, *residual);

      residual->copyFrom(rhs);
      /* calculate residual and possibly new rhs */
      mat->mult( 1.0, *residual, -1.0, *iter);

      /* res = res - A * drhs where A * drhs_out = drhs_in */
      rnorm = residual->twonorm();

      if( PIPSisEQ(rnorm, res_last) || rnorm > 100 * best_resid )
         done = true;
      else
         res_last = rnorm;

      if( rnorm < best_resid )
      {
         best_resid = rnorm;
         best_iter->copyFrom(*iter);
      }

      if( rnorm < precision * ( 1.0 + rhsnorm ) )
         done = true;

      ++n_iter_ref;

      /* refactorize */
      if( done && rnorm >= precision * (1.0 + rhsnorm ) )
      {
         if ( thresholdPivoting() >= threshold_pivoting_max )
         {
            if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            {
               std::cout << "WARNING MA27: threshold_pivoting parameter is already at its max and iterative refinement steps are exceeded with unsifficient precision" << std::endl;
               std::cout << " did not converge but still keeping the iterate" << std::endl;
            }
         }
         else
         {
            setThresholdPivoting( std::min( thresholdPivoting() * threshold_pivoting_factor, threshold_pivoting_max) );

            if( gOoqpPrintLevel >= ooqp_print_level_warnings )
               std::cout << "STATUS Ma27: Setting ThresholdPivoting parameter to " << thresholdPivoting() << " and refactorizing" << std::endl;

            done = false;
            n_iter_ref = 0;

            residual->copyFrom(rhs);
            iter->setToZero();
         }
      }

   }

   rnorm = best_resid;
   if( rnorm >= precision * (1.0 + rhsnorm ) )
   {
      std::cout << "WARNING MA27: big residual after solve : " << rnorm / (1.0 + rhsnorm ) << " > " << precision << std::endl;
      std::fstream file("mat.out");

      mat->writeToStreamDense(file);
      file << " rhs ";
      rhs.writeToStreamAll(file);
      file << " x ";
      best_iter->writeToStreamAll(file);
      residual->copyFrom(rhs);
      file << " resid ";
      mat->mult( 1.0, *residual, -1.0, *best_iter);
      residual->writeToStreamAll(file);
      assert(false);
   }

   rhs.copyFrom(*best_iter);
}

void Ma27Solver::copyMatrixElements( double afact[], int lafact ) const
{
   assert( lafact >= nnz );
   const double * M = mat_storage->M;
   std::copy( M, M + nnz, afact );

   if( lafact > nnz )
      std::fill( afact + nnz, afact + (lafact - nnz), 0.0 );
}

// TODO same as the one in MA57 - move somewhere else, some common MA_Solver thing maybe..
void Ma27Solver::getIndices( int irow[], int jcol[] ) const
{
   const int *krowM = mat_storage->krowM;
   for( int i = 0; i < mat_storage->n; i++ )
   {
      if( mat_storage->fortranIndexed() )
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
      if( !mat_storage->fortranIndexed() )
         jcolM[k] = mat_storage->jcolM[k] + 1;
      else
         jcolM[k] = mat_storage->jcolM[k];
   }
}

Ma27Solver::~Ma27Solver()
{
   freeWorkingArrays();
}

void Ma27Solver::freeWorkingArrays()
{
   if( irowM )
      delete[] irowM;
   if( jcolM )
      delete[] jcolM;
   if( fact )
      delete[] fact;
   if( ikeep )
      delete[] ikeep;
   if( iw )
      delete[] iw;
   if( iw1 )
      delete[] iw1;
   if( iw2 )
      delete[] iw2;
   if( w )
      delete[] w;

   irowM = jcolM = ikeep = iw = iw1 = iw2 = nullptr;
   fact = w = nullptr;
}

bool Ma27Solver::checkErrorsAndReact()
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
         std::cerr << "ERROR MA27: N out of range or < -1: " << n << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -2 :
      {
         std::cerr << "ERROR MA27: NNZ out of range or < -1 : " << nnz << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -3 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27: insufficient space in iw: " << liw << " suggest reset to " << error_info << std::endl;
         ipessimism *= 1.1;

         assert( iw );
         delete[] iw;

         liw = std::max( error_info, static_cast<int>(1.1 * liw) );
         iw = new int[ liw ];
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << liw << std::endl;

         error = true;
      }; break;
      case -4 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27: insufficient factorization space: " << la << std::endl;;
         rpessimism *= 1.1;

         assert( fact );
         delete[] fact;

         la = std::max( error_info, static_cast<int>(1.1 * la) );
         fact = new double[la];

         this->copyMatrixElements(fact, la);
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << la << std::endl;

         error = true;
      } break;
      case -5:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27: matrix apparently numerically singular, detected at stage " << error_info << std::endl;

         if( getSmallPivot() <= threshold_pivtol )
         {
            std::cout << " cannot decrease pivtol anymore -- accepting factorization anyway" << std::endl;
            assert( getSmallPivot() == threshold_pivtol );
         }
         else
         {
            const double curr_pivtol = getSmallPivot();
            const double new_pivtol = std::max( threshold_pivtol, curr_pivtol * threshold_pivtol_factor );
            std::cout << " decreasing pivtol from " << curr_pivtol << " to " << new_pivtol << std::endl;

            setSmallPivot( new_pivtol );
         }
      }; break;
      case -6:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27: change of sign of pivots detected at stage " << error_info << std::endl;
      }; break;
      case -7:
      {
         std::cerr << "ERROR MA27: value of NSTEPS out of range " << nsteps << " (should not happen..) " << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case 1 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27: detected " << error_info << " entries out of range in irowM and jcolM; ignored" << std::endl;
      }; break;
      case 2 :
      {
         std::cerr << "ERROR MA27: change of sign in pivots detected when matrix is supposedly definite" << std::endl;
         MPI_Abort(MPI_COMM_WORLD, -1);
      } break;
      case 3:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27: rank deficient matrix detected; apparent rank is " << error_info << std::endl;
      }; break;
      default :
      {
         assert( 0 == error_flag );
      }; break;
   }

   return error;
}
