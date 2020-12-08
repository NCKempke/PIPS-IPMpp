/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "Ma27Solver.h"
#include "Mc30Scaler.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"

#include <fstream>
#include <algorithm>

extern int gOoqpPrintLevel;

Ma27Solver::Ma27Solver(const SparseSymMatrix* sgm, const std::string& name_) :
       mat(sgm), mat_storage(sgm->getStorageHandle()), name( name_ ), scaler( new Mc30Scaler() )
{
   init();
}

void Ma27Solver::init()
{
   const double default_small_pivot = 1e-10;
   const double default_threshold_pivoting = 0.5;
   /* detecting dense rows during the factorization to preserve sparsity */
   const double default_fratio = 0.5;
   assert( mat_storage->n == mat_storage->m );
   n = mat_storage->n;
   assert( n > 0 );
   nnz = mat_storage->numberOfNonZeros();

   FNAME(ma27id)(icntl.data(), cntl.data());

   this->setThresholdPivoting( default_threshold_pivoting );
   cntl[1] = default_fratio;
   this->setSmallPivot( default_small_pivot );

   icntl[0] = 0;
   icntl[1] = 0;
}

void Ma27Solver::firstCall()
{
  /* convert to MA27 format */
  irowM.resize(nnz);
  jcolM.resize(nnz);
  this->getIndices( irowM, jcolM );

  /* set working arrays for iterative refinement */
  w_ma60.resize(3 * n);
  iw_ma60.resize(2 * n);

  /* set working stuff for factorization routine */
  liw = static_cast<int>(ipessimism * (2 * nnz + 3 * n + 1));
  iw = new int[liw];
  iw1 = new int[2 * n];
  ikeep = new int[3 * n];

  int iflag = 0; // set to 1 if ikeep contains pivot order
  double ops;

  /* do ordering */
  bool done = false;
  int tries = 0;
  do
  {
     FNAME(ma27ad)( &n, &nnz, irowM.data(), jcolM.data(), iw, &liw, ikeep, iw1, &nsteps, &iflag, icntl.data(), cntl.data(), info.data(), &ops);
     done = !checkErrorsAndReact();
     ++tries;
  }
  while( !done && tries < max_tries );

  if ( !done && tries > max_tries )
  {
     std::cerr << "ERROR MA27: could not get ordering of matrix after max " << max_tries << " tries" << "\n";
     MPI_Abort(MPI_COMM_WORLD, -1);
  }

  delete [] iw;
  delete [] iw1;

  la = rpessimism * this->minimumRealWorkspace();
  fact.resize(la);

  // set iw and iw1 in prep for calls to ma27bd and ma27cd
  liw = ipessimism * this->minimumIntWorkspace();
  iw = new int[liw];
  iw1 = new int[n];
}  

void Ma27Solver::diagonalChanged( int /* idiag */, int /* extent */ )
{
  this->matrixChanged();
}

bool new_factor;
void Ma27Solver::matrixChanged()
{
   if( fact.size() == 0 )
      this->firstCall();

   bool done = false;
   int tries = 0;
   new_factor = true;
   do
   {
      // copy M to fact
      this->copyMatrixElements(fact, la);
      scaler->scaleMatrixTripletFormat( n, nnz, fact.data(), irowM.data(), jcolM.data(), true);

      FNAME(ma27bd)(&n, &nnz, irowM.data(), jcolM.data(), fact.data(), &la, iw, &liw, ikeep, &nsteps,
            &maxfrt, iw1, icntl.data(), cntl.data(), info.data());

      done = !checkErrorsAndReact();
      tries++;
   }
   while( !done && tries < max_tries );

   if ( !done && tries > max_tries )
   {
      std::cout << "ERROR MA27: could not get factorization of matrix after max " << max_tries << " tries" << "\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   assert( 0 < nsteps && nsteps < n );
   assert( 0 < maxfrt && nsteps < n );
   iw2 = new int[nsteps];
   w = new double[maxfrt];

   assert( iw2 );
   assert( w );
}

void Ma27Solver::solve( int nrhss, double* rhss, int* colSparsity )
{
   for (int i = 0; i < nrhss; i++) {
     SimpleVector v(rhss + i * n, n);

     solve(v);
     //solveIterRef(v); slower for some reason
   }
}

void Ma27Solver::solveIterRef( OoqpVector& rhs_in )
{
   SimpleVector &rhs = dynamic_cast<SimpleVector&>(rhs_in);

   SimpleVectorHandle x(new SimpleVector(n));
   SimpleVectorHandle y(new SimpleVector(n));
   x->copyFrom( rhs );

   /* init iterative refinement process */
   FNAME(ma60id)(icntl_ma60, keep_ma60, rkeep_ma60);

   /* obtain initial solution */
   scaler->scaleVector(*x);
   FNAME(ma27cd)(&n, fact.data(), &la, iw, &liw, w, &maxfrt, x->elements(), iw1,
         &nsteps, icntl.data(), info.data());
   scaler->scaleVector(*x);

   /* start of iterative refinement */

   /* job[0] <= 0 -> only do backward error estimate
    * job[1] = 1 -> matrix is symmetric
    */
   int job[2] = { 0, 1 };
   int kase = 0;
   double omega[2];
   double error_x;
   double cond[2];
   int n_iters_needed;

   // icntl[0] -> output -> 0 to suppress
   icntl[1] = max_n_iter_refinement; // default 16

   bool done = false;
   while( !done )
   {
      FNAME(ma60ad)( &n, &nnz, mat_storage->M, irowM.data(), jcolM.data(), rhs.elements(), x->elements(), y->elements(),
            scaler->getScaling(), w_ma60.data(), iw_ma60.data(),
            &kase, omega, &error_x, job, cond, &n_iters_needed, icntl_ma60, keep_ma60, rkeep_ma60 );
      assert( kase <= 2 );
      assert( kase >= 0 || kase == -3 );
      if( kase == -3 )
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27 " << name << ": large residual " << error_x <<
               " even after (" << icntl[1] << "/" << max_n_iter_refinement << ") iterative refinement steps" << "\n";
         done = true;
      }
      else if( kase == 0 )
      {
         done = true;
      }
      else if( kase == 1 || kase == 2 )
      {
         FNAME(ma27cd)(&n, fact.data(), &la, iw, &liw, w, &maxfrt, y->elements(), iw1,
               &nsteps, icntl.data(), info.data());
      }
   }

   rhs.copyFrom( *x );
}

void Ma27Solver::solve( OoqpVector& rhs_in )
{
   SimpleVector &rhs = dynamic_cast<SimpleVector&>(rhs_in);

#ifndef NDEBUG
   for( int i = 0; i < rhs.length(); ++i )
      if( std::fabs(rhs[i]) > 1e50 )
      {
         std::cout << "Big entry in right hand side vector..." << "\n";
         break;
      }
#endif
//   SimpleVector* rhs_cpy = dynamic_cast<SimpleVector*>(rhs_in.cloneFull());

//   /* sparsify rhs */
//   for( int i = 0; i < rhs.length(); ++i )
//      if( std::fabs(rhs[i]) < precision )
//         rhs[i] = 0.0;

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
      assert( maxfrt > 0 );
      assert( nsteps > 0 );
      assert( w );
      assert( residual->elements() );
      assert( iw1 );

      /* solve DAD y = D*residual */
      scaler->scaleVector(*residual);
      FNAME(ma27cd)(&n, fact.data(), &la, iw, &liw, w, &maxfrt, residual->elements(), iw1,
            &nsteps, icntl.data(), info.data());
      /* Dy = x */
      scaler->scaleVector(*residual);

      iter->axpy(1.0, *residual);

      residual->copyFrom(rhs);
      /* calculate residual and possibly new rhs */
      mat->mult( 1.0, *residual, -1.0, *iter);

      /* res = res - A * drhs where A * drhs_out = drhs_in */
      rnorm = residual->twonorm();

      if( PIPSisEQ(rnorm, res_last, 1e-7) || rnorm > 100 * best_resid )
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
               std::cout << "WARNING MA27 " << name << ": threshold_pivoting parameter is already at its max and iterative refinement steps are exceeded with unsifficient precision" << "\n";
               std::cout << " did not converge but still keeping the iterate" << "\n";

               std::cout << "Error is " << rnorm << " vs " << precision * (1.0 + rhsnorm) << " required.. \n";
            }
         }
         else
         {
            setThresholdPivoting( std::min( thresholdPivoting() * threshold_pivoting_factor, threshold_pivoting_max) );

            if( gOoqpPrintLevel >= ooqp_print_level_warnings )
               std::cout << "STATUS MA27 " << name << ": Setting ThresholdPivoting parameter to " << thresholdPivoting() << " and refactorizing" << "\n";

            done = false;
            n_iter_ref = 0;

            residual->copyFrom(rhs);
            iter->setToZero();
            matrixChanged();
         }

      }
   }

   rnorm = best_resid;

   if( rnorm >= precision * (1.0 + rhsnorm) )
   {
//      std::cout << "ERROR " << rnorm/(1.0 + rhsnorm) << " > " << precision << " (after " << n_iter_ref << " iter refs)\n";
//      rhs_cpy->writeToStreamAll(std::cout);
//      std::cout << "ERROR " << rnorm << " vs " << precision * (1.0 + rhsnorm) << " required " << "\n";
//      best_iter->writeToStreamAll(std::cout);
//      mat->writeToStreamDense(std::cout);
//      assert( false );
//      std::cout << "Writing K of local schur complement computation..." << "\n";
//      std::ofstream myfile("../test.prb");
//
//      myfile << "n: " << n << "\n";
//      myfile << "nnz: " << nnz << "\n";
//
//      myfile << "ia: ";
//      for( int i = 0; i <= n; i++ )
//         myfile << mat->getStorageRef().krowM[i] << ", ";
//      myfile << "\n";
//
//      myfile << "ja: ";
//      for( int i = 0; i < nnz; i++ )
//         myfile << mat->getStorageRef().jcolM[i] << ", ";
//      myfile << "\n";
//
//      myfile << "a: ";
//      for( int i = 0; i < nnz; i++ )
//         myfile << mat->getStorageRef().M[i] << ", ";
//      myfile << "\n";
//
//      myfile.close();
//
//      std::cout << "Writing rhs from local schur complement computation..." << "\n";
//      myfile.open("../test.rhs");
//
//      std::cout << "sizerhs " << size_t(1) * size_t(n) <<  "\n";
//
//      myfile << "nrhs: " << 1 << "\n";
//
//      myfile << "rhs: ";
//      for( int i = 0; i < n; i++ )
//         myfile << (*rhs_cpy)[i] << ", ";
//      myfile << "\n";
//
//      myfile.close();
//
//      assert(false);
   }

   rhs.copyFrom(*best_iter);

#ifndef NDEBUG
   auto pos = std::find_if( rhs.elements(), rhs.elements() + rhs.length(), []( double el ) { return std::fabs(el) > 1e50; });
   if( pos != rhs.elements() + rhs.length() )
   {
      std::cout << *pos << "\n";
      assert( false && "Big entry in solution vector... " );

   }
#endif

   /* sparsify rhs */
   std::transform( rhs.elements(), rhs.elements() + n, rhs.elements(), []( double el ){
      return std::fabs(el) < 1e-16 ? 0 : el;
   });
}

void Ma27Solver::copyMatrixElements( std::vector<double>& afact, int lafact ) const
{
   assert( lafact >= nnz );
   const double * M = mat_storage->M;
   std::copy( M, M + nnz, afact.begin() );

   if( lafact > nnz )
      std::fill( afact.begin() + nnz, afact.begin() + (lafact - nnz), 0.0 );
}

// TODO same as the one in MA57 - move somewhere else, some common MA_Solver thing maybe..
void Ma27Solver::getIndices( std::vector<int>& irowM, std::vector<int>& jcolM ) const
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
   if( scaler )
      delete scaler;

   freeWorkingArrays();
}

void Ma27Solver::freeWorkingArrays()
{
   irowM.clear();
   jcolM.clear();
   fact.clear();

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

   ikeep = iw = iw1 = iw2 = nullptr;
   w = nullptr;
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
         std::cout << "ERROR MA27 " << name << ": N out of range or < -1: " << n << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -2 :
      {
         std::cout << "ERROR MA27 " << name << ": NNZ out of range or < -1 : " << nnz << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case -3 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27 " << name << ": insufficient space in iw: " << liw << " suggest reset to " << error_info << "\n";
         ipessimism *= 1.1;

         assert( iw );
         delete[] iw;

         liw = std::max( static_cast<int>(ipessimism * error_info), static_cast<int>(ipessimism * liw) );
         iw = new int[ liw ];
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << liw << "\n";

         error = true;
      }; break;
      case -4 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27 " << name << ": insufficient factorization space: " << la << "\n";;
         rpessimism *= 1.1;

         la = std::max( static_cast<int>(rpessimism * error_info), static_cast<int>(rpessimism * la) );
         fact.resize(la);

         this->copyMatrixElements(fact, la);
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << " resetting to " << la << "\n";

         error = true;
      } break;
      case -5:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27 " << name << ": matrix apparently numerically singular, detected at stage " << error_info << "\n";

         if( getSmallPivot() <= threshold_pivtol )
         {
            std::cout << " cannot decrease pivtol anymore -- accepting factorization anyway" << "\n";
            assert( getSmallPivot() == threshold_pivtol );
         }
         else
         {
            const double curr_pivtol = getSmallPivot();
            const double new_pivtol = std::max( threshold_pivtol, curr_pivtol * threshold_pivtol_factor );
            std::cout << " decreasing pivtol from " << curr_pivtol << " to " << new_pivtol << "\n";

            setSmallPivot( new_pivtol );
         }
      }; break;
      case -6:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27 " << name << ": change of sign of pivots detected at stage " << error_info << "\n";
      }; break;
      case -7:
      {
         std::cerr << "ERROR MA27 " << name << ": value of NSTEPS out of range " << nsteps << " (should not happen..) " << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      }; break;
      case 1 :
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27 " << name << ": detected " << error_info << " entries out of range in irowM and jcolM; ignored" << "\n";
      }; break;
      case 2 :
      {
         std::cerr << "ERROR MA27 " << name << ": change of sign in pivots detected when matrix is supposedly definite" << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      } break;
      case 3:
      {
         if( gOoqpPrintLevel >= ooqp_print_level_warnings )
            std::cout << "WARNING MA27 " << name << ": rank deficient matrix detected; apparent rank is " << error_info << " != n : " << this->n << "\n";

         static double last_pert = 1e-14;

         double pert = 0;
         if( new_factor )
            pert = last_pert;
         else
            pert *= 10;

         int n,m; mat_storage->getSize(m,n);
         assert( m == n );
         std::cout << "PERTURBATION! " << pert << "\n";
         for ( int i = 0; i < m; i++ )
         {
            for( int k = mat_storage->krowM[i]; k < mat_storage->krowM[i+1]; k++ )
            {
               int j = mat_storage->jcolM[k];
               if ( i == j )
                  mat_storage->M[k] += pert;
            }
         }
         error = true;

      }; break;
      default :
      {
         assert( 0 == error_flag );
      }; break;
   }

   return error;
}
