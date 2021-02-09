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
         dynamic_cast<StringGenMatrix&>(*dynamic_cast<StochGenMatrix&>(*data->C).Blmat), true );

   GenMatrix* BlFbuf = Bl.F.mat;
   GenMatrix* BlGbuf = Bl.G.mat;

   Bl.F.mat = &data->getLocalF();
   Bl.G.mat = &data->getLocalG();

   addBTKiInvBToSC( result, Bl, Br, Br_mod_border, false, false );

   Bl.F.mat = BlFbuf;
   Bl.G.mat = BlGbuf;
}

void sLinsysRootAugHierInner::addTermToSchurComplBlocked(sData* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC )
{
   assert( data == prob );

   BorderLinsys Bl( prob->getLocalFBorder(), prob->getLocalGBorder(), use_local_RAC );
   BorderLinsys Br( prob->getLocalFBorder(), prob->getLocalGBorder(), use_local_RAC );

   GenMatrix* BlFbuf = Bl.F.mat;
   GenMatrix* BlGbuf = Bl.G.mat;

   Bl.F.mat = &prob->getLocalF();
   Bl.G.mat = &prob->getLocalG();

   std::vector<BorderMod> border_mod;
   addBTKiInvBToSC( SC, Bl, Br, border_mod, true, sparseSC );

   Bl.F.mat = BlFbuf;
   Bl.G.mat = BlGbuf;
}

void sLinsysRootAugHierInner::LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0,
      BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sparse_res, bool sym_res )
{
   BorderLinsys B_inner( data->getLocalFBorder(), data->getLocalGBorder(), true );

   GenMatrix* B_inner_Fbuf = Bl.F.mat;
   GenMatrix* B_inner_Gbuf = Bl.G.mat;

   B_inner.F.mat = &data->getLocalF();
   B_inner.G.mat = &data->getLocalG();

   BorderMod B_inner_mod(B_inner, X0);

   std::vector<BorderMod> border_mod( Br_mod_border );
   border_mod.push_back( B_inner_mod );

   addBTKiInvBToSC( res, Bl, Br, border_mod, sym_res, sparse_res );

   B_inner.F.mat = B_inner_Fbuf;
   B_inner.G.mat = B_inner_Gbuf;

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
