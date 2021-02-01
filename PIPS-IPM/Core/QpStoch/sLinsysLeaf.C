/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeaf.h"

void sLinsysLeaf::factor2(sData*, Variables*)
{
   // Diagonals were already updated, so
   // just trigger a local refactorization (if needed, depends on the type of lin solver).
   stochNode->resMon.recFactTmLocal_start();

   if( computeBlockwiseSC )
   {
      #pragma omp parallel num_threads(n_solvers)
      {
         omp_set_num_threads(n_threads_solvers);

         const SparseStorage& kkt_mod = dynamic_cast<SparseSymMatrix&>(*kkt).getStorageRef();
         const int id = omp_get_thread_num();

         SparseSymMatrix& my_kkt = dynamic_cast<SparseSymMatrix&>(*problems_blocked[id].get());
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

void sLinsysLeaf::Dsolve( sData*, OoqpVector& x_in )
{
   StochVector& x = dynamic_cast<StochVector&>(x_in);
   assert(x.children.size()==0);
   stochNode->resMon.recDsolveTmChildren_start();
   solver->Dsolve(*x.vec);
   stochNode->resMon.recDsolveTmChildren_stop();
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

void sLinsysLeaf::deleteChildren()
{ }

void sLinsysLeaf::addTermToSchurComplBlocked( sData *prob, bool sparseSC, SymMatrix& SC, bool use_RAC_inner_border )
{
   assert( use_RAC_inner_border );

   const bool sc_is_sym = true;

   BorderBiBlock border_right( prob->getLocalCrossHessian(), prob->getLocalA(), prob->getLocalC(),
         prob->getLocalF().getTranspose(), prob->getLocalG().getTranspose() );
   BorderBiBlock border_left_transp( prob->getLocalCrossHessian().getTranspose(), prob->getLocalA().getTranspose(), prob->getLocalC().getTranspose(),
         prob->getLocalF(), prob->getLocalG() );

   addBiTLeftKiBiRightToResBlockedParallelSolvers( sparseSC, sc_is_sym, border_left_transp, border_right, SC );
}


/* compute result += B_inner K^-1 Br */
void sLinsysLeaf::addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, bool use_RAC_inner_border)
{
   assert( Br.A.children.size() == 0 );

   const bool result_sparse = false;
   const bool result_sym = false;

   const bool has_RAC = !( Br.R.isEmpty() && Br.A.isEmpty() && Br.C.isEmpty() );
   if( use_RAC_inner_border )
      assert( !has_RAC );

   std::unique_ptr<SparseGenMatrix> dummy( new SparseGenMatrix(0,0,0) );

   SparseGenMatrix& R = has_RAC ? dynamic_cast<SparseGenMatrix&>(*Br.R.mat) : ( use_RAC_inner_border ? data->getLocalCrossHessian() : *dummy );
   SparseGenMatrix& A = has_RAC ? dynamic_cast<SparseGenMatrix&>(*Br.A.mat) : ( use_RAC_inner_border ? data->getLocalA() : *dummy );
   SparseGenMatrix& C = has_RAC ? dynamic_cast<SparseGenMatrix&>(*Br.C.mat) : ( use_RAC_inner_border ? data->getLocalC() : *dummy );

   BorderBiBlock border_right( R, A, C,
         dynamic_cast<SparseGenMatrix&>(*Br.F.mat).getTranspose(),
         dynamic_cast<SparseGenMatrix&>(*Br.G.mat).getTranspose() );

   BorderBiBlock border_left_transp( *dummy, *dummy, *dummy, data->getLocalF(), data->getLocalG() );

   addBiTLeftKiBiRightToResBlockedParallelSolvers( result_sparse, result_sym, border_left_transp, border_right, result);
}

void sLinsysLeaf::mySymAtPutSubmatrix(SymMatrix& kkt_, 
					     GenMatrix& B_, GenMatrix&,
					     int locnx, int locmy, int)
{
  SparseSymMatrix& kkt = dynamic_cast<SparseSymMatrix&>(kkt_);
  SparseGenMatrix& B   = dynamic_cast<SparseGenMatrix&>(B_);
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
}

void sLinsysLeaf::addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border )
{
   assert( border.A.children.size() == 0 );
   assert( rhs.children.size() == 0 );

   assert( border.R.mat );
   assert( border.A.mat );
   assert( border.C.mat );

   SparseGenMatrix& Ri_border = dynamic_cast<SparseGenMatrix&>(*border.R.mat);
   int mRi, nRi; Ri_border.getSize(mRi, nRi);

   SparseGenMatrix& Ai_border = dynamic_cast<SparseGenMatrix&>(*border.A.mat);
   int mAi, nAi; Ai_border.getSize(mAi, nAi);

   SparseGenMatrix& Ci_border = dynamic_cast<SparseGenMatrix&>(*border.C.mat);
   int mCi, nCi; Ci_border.getSize(mCi, nCi);

   assert( border.F.mat );
   assert( border.G.mat );
   SparseGenMatrix& Fi_border = dynamic_cast<SparseGenMatrix&>(*border.F.mat);
   int mFi, nFi; Fi_border.getSize(mFi, nFi);
   SparseGenMatrix& Gi_border = dynamic_cast<SparseGenMatrix&>(*border.G.mat);
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

void sLinsysLeaf::addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border )
{
   assert( border.A.children.size() == 0 );
   assert( rhs.children.size() == 0 );

   assert( border.R.mat );
   assert( border.A.mat );
   assert( border.C.mat );

   SparseGenMatrix& Ri_border = dynamic_cast<SparseGenMatrix&>(*border.R.mat);
   int mRi, nRi; Ri_border.getSize(mRi, nRi);

   SparseGenMatrix& Ai_border = dynamic_cast<SparseGenMatrix&>(*border.A.mat);
   int mAi, nAi; Ai_border.getSize(mAi, nAi);

   SparseGenMatrix& Ci_border = dynamic_cast<SparseGenMatrix&>(*border.C.mat);
   int mCi, nCi; Ci_border.getSize(mCi, nCi);

   assert( border.F.mat );
   assert( border.G.mat );
   SparseGenMatrix& Fi_border = dynamic_cast<SparseGenMatrix&>(*border.F.mat);
   int mFi, nFi; Fi_border.getSize(mFi, nFi);
   SparseGenMatrix& Gi_border = dynamic_cast<SparseGenMatrix&>(*border.G.mat);
   int mGi, nGi; Gi_border.getSize(mGi, nGi);

   assert( rhs.vec );
   assert( rhs.vec->length() == mRi + mAi + mCi );
   assert( x0.length() == nRi + mFi + mGi );

   SimpleVector& rhsi = dynamic_cast<SimpleVector&>(*rhs.vec);

   double* rhsi1 = &rhsi[0];
   double* rhsi2 = &rhsi[mRi];
   double* rhsi3 = &rhsi[mRi + mAi];

   const double* x1 = &x0[0];
   const double* x2 = &x0[nRi];
   const double* x3 = &x0[nRi + mFi];

   //   Ri_border.mult(1.0, rhsi1, -1.0, x1);
   Ri_border.mult( 1.0, rhsi1, 1, -1.0, x1, 1 );
   //   Ai_border.mult(1.0, rhsi2, -1.0, x1);
   Ai_border.mult( 1.0, rhsi2, 1, -1.0, x1, 1 );
   //   Ci_border.mult(1.0, rhsi3, -1.0, x1);
   Ci_border.mult( 1.0, rhsi3, 1, -1.0, x1, 1 );


   //   Fi_border.transMult(1.0, rhsi1, -1.0, x2);
   Fi_border.transMult( 1.0, rhsi1, 1, -1.0, x2, 1 );
   //   Gi_border.transMult(1.0, rhsi1, -1.0, x3);
   Gi_border.transMult( 1.0, rhsi1, 1, -1.0, x3, 1 );
}
