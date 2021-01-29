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

void sLinsysRootAugHierInner::addTermToSchurComplBlocked(sData* prob, bool sparseSC, SymMatrix& SC )
{
   std::unique_ptr<StringGenMatrix> dummy( new StringGenMatrix() );
   for( unsigned int i = 0; i < children.size(); ++i )
   {
      BorderLinsys Bl( *dummy, *dummy, *dummy, prob->getLocalFBorder(), prob->getLocalGBorder() );
      BorderLinsys Br( *dummy, *dummy, *dummy, prob->getLocalFBorder(), prob->getLocalGBorder() );
      addBTKiInvBToSC(SC, Bl, Br, true, sparseSC );
   }
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
