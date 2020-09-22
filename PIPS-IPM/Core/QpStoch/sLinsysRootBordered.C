/*
 * sLinsysRootBordered.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "sLinsysRootBordered.h"

sLinsysRootBordered::sLinsysRootBordered(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_)
{
//   buffer = NULL;
}

sLinsysRootBordered::~sLinsysRootBordered()
{
//   delete buffer;
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
