/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "Ma27Solver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "SimpleVectorHandle.h"

extern int gOoqpPrintLevel;

Ma27Solver::Ma27Solver( SparseSymMatrix * sgm ) :
      maxMa27Iter(18), precision(1e-7), threshold_pivoting_max(1.e-2), threshold_pivoting_factor(10.0),
      irowM(nullptr), jcolM(nullptr), fact(nullptr),
      ipessimism(2.0), rpessimism(2.0)
{
   mStorage = sgm->getStorageHandle();
   init();
}


void Ma27Solver::init()
{
   assert( mStorage->n == mStorage->m );
   n = mStorage->n;
   nnz = mStorage->numberOfNonZeros();

   FNAME(ma27id)(icntl, cntl);

   this->setTreatAsZero( 1.0e-12 );
   this->setThresholdPivoting( 1.0e-8 );

   icntl[0] = 0;
   icntl[1] = 0;
}

void Ma27Solver::firstCall()
{
  irowM = new int[nnz];
  jcolM = new int[nnz];

  this->getIndices( irowM, jcolM );

  // set array lengths as recommended in ma27 docs
  liw = 2 * (nnz + 3 * n + 1);
  iw = new int[liw];
  iw1 = new int[2 * n];
  ikeep = new int[3 * n];

  int iflag = 0; // set to 1 if ikeep contains pivot order
  double ops;

  FNAME(ma27ad)( &n, &nnz, irowM, jcolM, iw, &liw, ikeep, iw1, &nsteps, &iflag, icntl, cntl, info, &ops);

  delete [] iw;
  delete [] iw1;

  switch ( this->ma27ErrFlg() )
  {
     case -1 :
     {
        std::cerr << "n out of range: " << n << std::endl;
        assert( 0 );
     }; break;
     case -2 :
     {
        std::cerr << "nnz out of range: " << nnz << std::endl;
        assert(0);
     }; break;
     case -3 :
     {
        if( gOoqpPrintLevel >= 100 )
           std::cout << "insufficient space in iw: " << liw << " suggest reset to " << this->ierror() << std::endl;
     }; break;
     case 1 :
     {
        std::cerr << "detected " << this->ierror() << " entries out of range in irowM and jcolM; ignored" << std::endl;
     }; break;
  }

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

   int done = 0, tries = 0;
   do
   {
      // copy M to fact
      this->copyMatrixElements(fact, la);

      FNAME(ma27bd)(&n, &nnz, irowM, jcolM, fact, &la, iw, &liw, ikeep, &nsteps,
            &maxfrt, iw1, icntl, cntl, info);

      switch( this->ma27ErrFlg() )
         {
         case 0:
            done = 1;
            break;
         case -1:
            {
               std::cerr << "n out of range: " << n << std::endl;
               assert(0);
            };
            break;
         case -2:
            {
               std::cerr << "nnz out of range: " << nnz << std::endl;
               assert(0);
            };
            break;
         case -3:
            {
               if( gOoqpPrintLevel >= 100 )
                  std::cout << "insufficient space in iw: " << liw;
               delete[] iw;
               liw = (this->ierror() > ipessimism * liw) ?
                     this->ierror() : static_cast<int>(ipessimism * liw);
               iw = new int[liw];
               if( gOoqpPrintLevel >= 100 )
                  std::cout << " resetting to " << liw << std::endl;

               ipessimism *= 1.1;
            };
            break;
         case -4:
            {
               if( gOoqpPrintLevel >= 100 )
                  std::cout << "insufficient factorization space: " << la;
               delete[] fact;
               la = (this->ierror() > rpessimism * la) ?
                     this->ierror() : static_cast<int>(rpessimism * la);
               fact = new double[la];
               this->copyMatrixElements(fact, la);
               if( gOoqpPrintLevel >= 100 )
                  std::cout << " resetting to " << la << std::endl;
               rpessimism *= 1.1;
            };
            break;
         case -5:
            {
               if( gOoqpPrintLevel >= 100 )
               {
                  std::cout << "matrix apparently numerically singular, detected at stage " << this->ierror() << std::endl;
                  std::cout << "accept this factorization and hope for the best.." << std::endl;
               }
               done = 1;
            };
            break;
         case -6:
            {
               if( gOoqpPrintLevel >= 100 )
               {
                  std::cout << "change of sign of pivots detected at stage " << this->ierror() << std::endl;
                  std::cout << "but who cares " << std::endl;
               }
               done = 1;
            };
            break;
         case -7:
            {
               std::cerr << "value of NSTEPS out of range " << nsteps << std::endl;
               assert(0);
            };
            break;
         case 1:
            {
               if( gOoqpPrintLevel >= 100 )
               {
                  std::cout << "detected " << this->ierror() << " entries out of range in irowM and jcolM; ignored"
                        << std::endl;
               }
               done = 1;
            };
            break;
         case 3:
            {
               if( gOoqpPrintLevel >= 100 )
               {
                  std::cout << "rank deficient matrix detected; apparent rank is "
                        << this->ierror() << std::endl;
               }
               done = 1;
            };
            break;
         default:
            break;
         }
      tries++;
   }
   while( !done && tries < maxMa27Iter );

   if ( !done && maxMa27Iter >= 10)
   {
      if( gOoqpPrintLevel >= 100 )
         std::cout << "we are screwed; did not get a factorization after 10 tries " << std::endl;
   }

  iw2 = new int[nsteps];
  w = new double[maxfrt];
}

void Ma27Solver::solve( OoqpVector& rhs_in )
{
   SimpleVector &rhs = dynamic_cast<SimpleVector&>(rhs_in);

   // define structures to save rhs and store residuals
   SimpleVectorHandle resid(new SimpleVector(n));
   SimpleVectorHandle rhsSave(new SimpleVector(n));

   double *drhs = rhs.elements();
   double *dresid = resid->elements();

   double rnorm = 0.0;

   rhsSave->copyFrom(rhs);
   resid->copyFrom(rhs);

   const double rhsnorm = rhs.infnorm();

   bool done = false;
   int refactorizations = 0;

   /* iterative refinement loop */
   while( !done && refactorizations < 10 )
   {
      FNAME(ma27cd)(&n, fact, &la, iw, &liw, w, &maxfrt, drhs, iw1,
            &nsteps, icntl, info);

      // res = res - A * drhs where A * drhs_out = drhs_in
      mStorage->mult(-1.0, dresid, 1, 1.0, drhs, 1);
      rnorm = resid->infnorm();
    
      if( rnorm < precision * ( 1.0 + rhsnorm ) )
         done = true;
      else if ( this->thresholdPivoting() >= threshold_pivoting_max || refactorizations > 10)
      {
         if( gOoqpPrintLevel >= 10 )
            std::cout << "ThresholdPivoting parameter is already too high" << std::endl;
         done = true;
      }
      else
      {
         double tp = this->thresholdPivoting();
         tp *= threshold_pivoting_factor;
         if( tp > threshold_pivoting_max )
            tp = threshold_pivoting_max;
         this->setThresholdPivoting(tp);

         if( gOoqpPrintLevel >= 10 ) {
            std::cout << "Ma27: Setting ThresholdPivoting parameter to " << this->thresholdPivoting() << " for future factorizations" << std::endl;
      }

      this->matrixChanged();
      refactorizations++;

      resid->copyFrom(*rhsSave);
      rhs.copyFrom(*rhsSave);
    }
  }
}

void Ma27Solver::copyMatrixElements( double afact[], int lafact ) const
{
   assert( lafact >= nnz );
   const double * M = mStorage->M;
   std::copy( M, M + nnz, afact );

   if( lafact > nnz )
      std::fill_n( afact + nnz, afact + (lafact - nnz), 0.0 );
}

// TODO same as the one in MA57 - move somewhere else, some common MA_Solver thing maybe..
void Ma27Solver::getIndices( int irow[], int jcol[] ) const
{
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




