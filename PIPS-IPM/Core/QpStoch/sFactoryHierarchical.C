/*
 * sFactoryHierarchical.C
 *
 *  Created on: 10.09.2020
 *      Author: bzfkempk
 */

#include "sFactoryHierarchical.h"
#include <cassert>

sLinsysRoot* sFactoryHierarchical::newLinsysRoot()
{
   // TODO return new root linear system (linear system with inner matrix and border only)
   //return new sLinsysRootBordered(this, data);

   assert("not implemented" && 0 );
   return nullptr;
}

sLinsysRoot* sFactoryHierarchical::newLinsysRoot(sData *prob, OoqpVector *dd, OoqpVector *dq,
      OoqpVector *nomegaInv, OoqpVector *rhs)
{
   // TODO return new root linear system (linear system with inner matrix and border only)
   //return new sLinsysRootBordered(this, prob, dd, dq, nomegaInv, rhs);

   assert("not implemented" && 0 );
   return nullptr;
}

