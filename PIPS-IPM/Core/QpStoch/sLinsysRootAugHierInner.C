/*
 * sLinsysRootAugHierInner.C
 *
 *  Created on: 27.01.2021
 *      Author: bzfkempk
 */

#include "sLinsysRootAugHierInner.h"

sLinsysRootAugHierInner::sLinsysRootAugHierInner(sFactory *factory,
      sData *prob_, OoqpVector *dd_, OoqpVector *dq_, OoqpVector *nomegaInv_, OoqpVector *rhs_) :
      sLinsysRootAug(factory, prob_, dynamic_cast<StochVector*>(dd_)->vec,
            dynamic_cast<StochVector*>(dq_)->vec,
            dynamic_cast<StochVector*>(nomegaInv_)->vec, rhs_)
{
   assert( locnx == 0 );
   assert( locmy == 0 );
}


void sLinsysRootAugHierInner::assembleLocalKKT( sData* prob )
{
   for( size_t c = 0; c < children.size(); ++c )
   {
#ifdef STOCH_TESTING
      g_scenNum = c;
#endif
      if( children[c]->mpiComm == MPI_COMM_NULL )
         continue;

      children[c]->stochNode->resMon.recFactTmChildren_start();
      //---------------------------------------------
      addTermToSchurCompl(prob, c, false);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();
   }
}


void sLinsysRootAugHierInner::addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border )
{
   assert( dynamic_cast<StochGenMatrix&>(*data->A).Blmat->isKindOf(kStringGenMatrix) );
   assert( dynamic_cast<StochGenMatrix&>(*data->C).Blmat->isKindOf(kStringGenMatrix) );

   BorderLinsys Bl( dynamic_cast<StringGenMatrix&>(*dynamic_cast<StochGenMatrix&>(*data->A).Blmat),
         dynamic_cast<StringGenMatrix&>(*dynamic_cast<StochGenMatrix&>(*data->C).Blmat), false );
   std::vector<BorderMod> border_mod;

   addBTKiInvBToSC( result, Bl, Br, border_mod, false, false );
}

/* buffer is still transposed ..*/
void sLinsysRootAugHierInner::finalizeZ0Hierarchical( DenseGenMatrix& buffer, BorderLinsys&, std::vector<BorderMod>& Br_mod_border )
{
   assert( Br_mod_border.empty() );

   const SparseGenMatrix& F = data->getLocalF().getTranspose();
   const SparseGenMatrix& G = data->getLocalG().getTranspose();

   int mF, nF;
   F.getSize(mF, nF);
   int mG, nG;
   G.getSize(mG, nG);

   int mBuf, nBuf;
   buffer.getSize(mBuf, nBuf);

   assert( nBuf == nF + nG );
   assert( mF == mG );
   assert( mBuf >= mF );

   if( mF > 0 )
   {
      for( int rowF = 0; rowF < mF; ++rowF )
      {
         const double* MF = F.M();
         const int* krowF = F.krowM();
         const int* jcolF = F.jcolM();

         const int rowF_start = krowF[rowF];
         const int rowF_end = krowF[rowF + 1];

         for( int k = rowF_start; k < rowF_end; ++k )
         {
            const int col = jcolF[k];
            const double val_F = MF[k];

            assert( 0 <= col && col < nBuf - nG );
            buffer[rowF][col] += val_F;
         }
      }
   }

   if( mG > 0 )
   {
      for( int rowG = 0; rowG < mG; ++rowG )
      {
         const double* MG = G.M();
         const int* krowG = G.krowM();
         const int* jcolG = G.jcolM();

         const int rowG_start = krowG[rowG];
         const int rowG_end = krowG[rowG + 1];

         for( int k = rowG_start; k < rowG_end; ++k )
         {
            const int col = jcolG[k] + nF;
            const double val_G = MG[k];

            assert( nF <= col && col < nBuf);
            buffer[rowG][col] += val_G;
         }
      }
   }
}


void sLinsysRootAugHierInner::addTermToSchurComplBlocked(sData* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC )
{
   BorderLinsys Bl( prob->getLocalFBorder(), prob->getLocalGBorder(), use_local_RAC );
   BorderLinsys Br( prob->getLocalFBorder(), prob->getLocalGBorder(), use_local_RAC );
   std::vector<BorderMod> border_mod;
   addBTKiInvBToSC( SC, Bl, Br, border_mod, true, sparseSC );
}

void sLinsysRootAugHierInner::LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0,
      BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, int parent_nx, int parent_my, int parent_mz, bool sparse_res, bool sym_res )
{
   // TODO add own border to vector and pass on to solve
   assert(false);
}

void sLinsysRootAugHierInner::putXDiagonal( OoqpVector& xdiag_ )
{
  assert( dynamic_cast<StochVector&>(xdiag_).vec->isKindOf(kStochVector) );
  StochVector& xdiag = dynamic_cast<StochVector&>(*dynamic_cast<StochVector&>(xdiag_).vec);

  assert(children.size() == xdiag.children.size());

  xDiag = xdiag.vec;

  for(size_t it = 0; it < children.size(); it++)
    children[it]->putXDiagonal(*xdiag.children[it]);
}


void sLinsysRootAugHierInner::putZDiagonal( OoqpVector& zdiag_ )
{
  assert( dynamic_cast<StochVector&>(zdiag_).vec->isKindOf(kStochVector) );
  StochVector& zdiag = dynamic_cast<StochVector&>(*dynamic_cast<StochVector&>(zdiag_).vec);

  assert(children.size() == zdiag.children.size());

  zDiag = zdiag.vec;
  zDiagLinkCons = zdiag.vecl;

  for(size_t it = 0; it < children.size(); it++)
    children[it]->putZDiagonal(*zdiag.children[it]);
}
