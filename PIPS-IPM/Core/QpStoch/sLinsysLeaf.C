/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeaf.h"

sLinsysLeaf::~sLinsysLeaf()
{
   if( computeBlockwiseSC )
      freeBlockedSolvers();
}

void sLinsysLeaf::freeBlockedSolvers()
{
   assert( solvers_blocked != nullptr );

   #pragma omp parallel num_threads(n_solvers)
   {
      const int id = omp_get_thread_num();

      delete solvers_blocked[id];
      delete problems_blocked[id];
   }

   delete[] solvers_blocked;
   delete[] problems_blocked;
}

void sLinsysLeaf::factor2(sData *prob, Variables *vars)
{
   // Diagonals were already updated, so
   // just trigger a local refactorization (if needed, depends on the type of lin solver).
   stochNode->resMon.recFactTmLocal_start();

   if( computeBlockwiseSC )
   {
      #pragma omp parallel num_threads(n_solvers)
      {
         const SparseStorage& kkt_mod = dynamic_cast<SparseSymMatrix&>(*kkt).getStorageRef();
         const int id = omp_get_thread_num();

         SparseSymMatrix& my_kkt = *problems_blocked[id];
         kkt_mod.copyFrom( my_kkt.krowM(), my_kkt.jcolM(), my_kkt.M() );
         solvers_blocked[id]->matrixChanged();
      }
   }
   else
      solver->matrixChanged();

   stochNode->resMon.recFactTmLocal_stop();
}

void sLinsysLeaf::putXDiagonal( OoqpVector& xdiag_ )
{
  StochVector& xdiag = dynamic_cast<StochVector&>(xdiag_);
  kkt->atPutDiagonal( 0, *xdiag.vec );
}

void sLinsysLeaf::putZDiagonal( OoqpVector& zdiag_)
{
  StochVector& zdiag = dynamic_cast<StochVector&>(zdiag_);
  kkt->atPutDiagonal( locnx + locmy, *zdiag.vec );
}

void sLinsysLeaf::Lsolve (  sData *prob, OoqpVector& x_in )
{
   return;
}

void sLinsysLeaf::Dsolve( sData *prob, OoqpVector& x_in )
{
   StochVector& x = dynamic_cast<StochVector&>(x_in);
   assert(x.children.size()==0);
   stochNode->resMon.recDsolveTmChildren_start();
   solver->Dsolve(*x.vec);
   stochNode->resMon.recDsolveTmChildren_stop();
}

void sLinsysLeaf::Ltsolve (  sData *prob, OoqpVector& x_in )
{
   return;
}

void sLinsysLeaf::Ltsolve2( sData *prob, StochVector& x, SimpleVector& xp)
{
   StochVector& b = dynamic_cast<StochVector&>(x);
   SimpleVector& bi = dynamic_cast<SimpleVector&>(*b.vec);
   assert( 0 == b.children.size() );

#ifdef TIMING
   stochNode->resMon.eLtsolve.clear();
   stochNode->resMon.recLtsolveTmLocal_start();
#endif

   //b_i -= Lni^T x0
   this->LniTransMult(prob, bi, -1.0, xp);
   //  solver->Ltsolve(bi); -> empty
#ifdef TIMING
   stochNode->resMon.recLtsolveTmChildren_stop();
#endif
}

void sLinsysLeaf::sync()
{ assert(false); }


void sLinsysLeaf::deleteChildren()
{ }

void sLinsysLeaf::addTermToSchurComplBlocked(sData *prob, bool sparseSC, SymMatrix& SC)
{
   const bool sc_is_sym = true;

   BorderBiBlock border_right( prob->getLocalCrossHessian(), prob->getLocalA(), prob->getLocalC(), prob->getLocalF().getTranspose(), prob->getLocalG().getTranspose() );
   BorderBiBlock border_left_transp( prob->getLocalCrossHessian().getTranspose(), prob->getLocalA().getTranspose(), prob->getLocalC().getTranspose(),
         prob->getLocalF(), prob->getLocalG() );

   addBiTLeftKiBiRightToResBlockedParallelSolvers( sparseSC, sc_is_sym, border_left_transp, border_right, SC );
}

void sLinsysLeaf::mySymAtPutSubmatrix(SymMatrix& kkt_, 
					     GenMatrix& B_, GenMatrix& D_, 
					     int locnx, int locmy, int locmz)
{
  SparseSymMatrix& kkt = reinterpret_cast<SparseSymMatrix&>(kkt_);
  SparseGenMatrix& B   = reinterpret_cast<SparseGenMatrix&>(B_);
  //SparseGenMatrix& D   = reinterpret_cast<SparseGenMatrix&>(D_);

  int* jcolK = kkt.jcolM(); int* jcolB = B.jcolM(); //int* jcolD = D.jcolM(); 
  int* krowK = kkt.krowM(); int* krowB = B.krowM(); //int* krowD =  D.krowM();
  double* MK = kkt.M();     double* MB = B.M();

  for(int i=0; i<locmy; i++) {
    int itK = krowK[i+locnx];
    int j = krowB[i];

    for(; j<krowB[i+1]; j++) { 

      if(jcolB[j]<i+locnx) {
	jcolK[itK]=jcolB[j]; 
	MK[itK]=MB[j]; 
	itK++;
      }
    }
    jcolK[itK]=i+locnx; MK[itK] = 0.0; itK++;

    assert(j==krowB[i+1]);

    krowK[i+locnx+1]=itK;
  }
  assert(locmz==0);
}

void sLinsysLeaf::addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border )
{
   assert( border.A.children.size() == 0 );
   assert( rhs.children.size() == 0 );

   assert( border.R.mat );
   assert( border.A.mat );
   assert( border.C.mat );

   SparseGenMatrix& Ri_border = *border.R.mat;
   int mRi, nRi; Ri_border.getSize(mRi, nRi);

   SparseGenMatrix& Ai_border = *border.A.mat;
   int mAi, nAi; Ai_border.getSize(mAi, nAi);

   SparseGenMatrix& Ci_border = *border.C.mat;
   int mCi, nCi; Ci_border.getSize(mCi, nCi);

   assert( border.F.mat );
   assert( border.G.mat );
   SparseGenMatrix& Fi_border = *border.F.mat;
   int mFi, nFi; Fi_border.getSize(mFi, nFi);
   SparseGenMatrix& Gi_border = *border.G.mat;
   int mGi, nGi; Gi_border.getSize(mGi, nGi);

   assert( rhs.vec );
   assert( rhs.vec->length() == mRi + mAi + mCi );
   assert( b0.length() == nRi + mFi + mGi );

   SimpleVector& zi = dynamic_cast<SimpleVector&>(*rhs.vec);

   SimpleVector zi1 (&zi[0], mRi);
   SimpleVector zi2 (&zi[mRi], mAi );
   SimpleVector zi3 (&zi[mRi + mAi], mCi);

   SimpleVector b1( &b0[0], nRi );
   SimpleVector b2( &b0[nRi], mFi );
   SimpleVector b3( &b0[nRi + mFi], mGi );

   Ri_border.transMult(1.0, b1, -1.0, zi1);
   Ai_border.transMult(1.0, b1, -1.0, zi2);
   Ci_border.transMult(1.0, b1, -1.0, zi3);

   Fi_border.mult(1.0, b2, -1.0, zi1);
   Gi_border.mult(1.0, b3, -1.0, zi1);
}

void sLinsysLeaf::addBorderX0ToRhs( StochVector& rhs, SimpleVector& x0, BorderLinsys& border )
{
   assert( border.A.children.size() == 0 );
   assert( rhs.children.size() == 0 );

   assert( border.R.mat );
   assert( border.A.mat );
   assert( border.C.mat );

   SparseGenMatrix& Ri_border = *border.R.mat;
   int mRi, nRi; Ri_border.getSize(mRi, nRi);

   SparseGenMatrix& Ai_border = *border.A.mat;
   int mAi, nAi; Ai_border.getSize(mAi, nAi);

   SparseGenMatrix& Ci_border = *border.C.mat;
   int mCi, nCi; Ci_border.getSize(mCi, nCi);

   assert( border.F.mat );
   assert( border.G.mat );
   SparseGenMatrix& Fi_border = *border.F.mat;
   int mFi, nFi; Fi_border.getSize(mFi, nFi);
   SparseGenMatrix& Gi_border = *border.G.mat;
   int mGi, nGi; Gi_border.getSize(mGi, nGi);

   assert( rhs.vec );
   assert( rhs.vec->length() == mRi + mAi + mCi );
   assert( x0.length() == nRi + mFi + mGi );

   SimpleVector& zi = dynamic_cast<SimpleVector&>(*rhs.vec);

   SimpleVector zi1 (&zi[0], mRi);
   SimpleVector zi2 (&zi[mRi], mAi );
   SimpleVector zi3 (&zi[mRi + mAi], mCi);

   SimpleVector x1( &x0[0], nRi );
   SimpleVector x2( &x0[nRi], mFi );
   SimpleVector x3( &x0[nRi + mFi], mGi );

   Ri_border.mult(1.0, zi1, -1.0, x1);
   Ai_border.mult(1.0, zi2, -1.0, x1);
   Ci_border.mult(1.0, zi3, -1.0, x1);

   Fi_border.transMult(1.0, zi1, -1.0, x2);
   Gi_border.transMult(1.0, zi1, -1.0, x3);
}
