/*
 * sLinsysRootBordered.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "sLinsysRootBordered.h"

#include "DeSymIndefSolver.h"
#include "DeSymIndefSolver2.h"
#include "DeSymPSDSolver.h"
#include "sLinsysRootAug.h"
#include "sFactory.h"

sLinsysRootBordered::sLinsysRootBordered(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_, true)
{
   assert(locmyl >= 0 && locmzl >= 0);

   kkt = createKKT(prob_);
   solver = createSolver(prob_, kkt);
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

/* create kkt used to store Schur Complement of border layer */
SymMatrix* sLinsysRootBordered::createKKT(sData* prob)
{
   const int n = locmy + locmyl + locmzl;

   return new DenseSymMatrix(n);
}

void sLinsysRootBordered::assembleLocalKKT(sData* prob)
{
   assert( is_hierarchy_root );
   assert( !hasSparseKkt );
   assert( children.size() == 1 );

   DenseSymMatrix& kktd = dynamic_cast<DenseSymMatrix&>(*kkt);

   this->children[0]->addInnerToHierarchicalSchurComplement(kktd, prob);

   assert( 0 && "TODO : implement..");
}

DoubleLinearSolver* sLinsysRootBordered::createSolver(sData* prob, SymMatrix* kktmat_)
{
   DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
   return new DeSymIndefSolver(kktmat);
   //return new DeSymIndefSolver2(kktmat, locnx); // saddle point solver
   //return new DeSymPSDSolver(kktmat);
}
