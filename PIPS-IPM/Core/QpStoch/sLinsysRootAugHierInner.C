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
}
