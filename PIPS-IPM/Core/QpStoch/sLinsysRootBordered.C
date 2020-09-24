/*
 * sLinsysRootBordered.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "sLinsysRootBordered.h"
#include "sLinsysRootAug.h"
#include "sFactory.h"

sLinsysRootBordered::sLinsysRootBordered(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_, true)
{
   assert(locmyl >= 0 && locmzl >= 0);

   kkt = createKKT(prob_);
   solver = createSolver(prob_, kkt);

   //   buffer = NULL;
}

sLinsysRootBordered::~sLinsysRootBordered()
{
   // delete buffer;
   assert( children.size() == 0 );
   delete children[0];
}

void sLinsysRootBordered::finalizeKKT(sData* prob, Variables* vars)
{
   assert( 0 && "TODO: implement..");
}

void sLinsysRootBordered::solveReduced( sData *prob, SimpleVector& b)
{
   assert( 0 && "TODO: implement..");
}

SymMatrix* sLinsysRootBordered::createKKT(sData* prob)
{
   assert( 0 && "TODO: implement..");
   return nullptr;
}

DoubleLinearSolver* sLinsysRootBordered::createSolver(sData* prob, SymMatrix* kktmat)
{
   assert( 0 && "TODO: implement..");
   return nullptr;
}
