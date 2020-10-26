/*
 * sLinsysRootBordered.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "sLinsysRootBordered.h"

#include "BorderedSymMatrix.h"
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
   assert( children.size() == 0 );
   delete children[0];
}

void sLinsysRootBordered::finalizeKKT(/* const */sData* prob, Variables* vars)
{
   /* Add corner block
    * [ Q0 F0T G0T  ]
    * [ F0  0   0   ]
    * [ G0  0 OmN+1 ]
    */
   assert( prob->isHierarchieRoot() );

   const SparseGenMatrix& F0 = *dynamic_cast<const BorderedGenMatrix&>(*prob->A).bottom_left_block;
   const SparseGenMatrix& G0 = *dynamic_cast<const BorderedGenMatrix&>(*prob->C).bottom_left_block;
   const SparseSymMatrix& Q0 = dynamic_cast<const SparseSymMatrix&>(*dynamic_cast<const BorderedSymMatrix&>(*prob->Q).top_left_block);

   DenseSymMatrix& SC = dynamic_cast<DenseSymMatrix&>(*kkt);
   //   SC.symAtPutDense()

   /////////////////////////////////////////////////////////////
   // update the KKT with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   const int* krowQ0 = Q0.krowM();
   const int* jcolQ0 = Q0.jcolM();
   const double* MQ = Q0.M();

   for( int rowQ = 0; rowQ < locnx; rowQ++)
   {
      for( int k = krowQ0[rowQ]; k < krowQ0[rowQ + 1]; ++k)
      {
         const int colQ = jcolQ0[k];
         if( rowQ == colQ )
            continue;
         double val = MQ[k];
         SC[rowQ][colQ] += val;
         SC[colQ][rowQ] += val;

         assert(0 && "non-empty Q currently not supported");
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q) + X^{-1} S
   /////////////////////////////////////////////////////////////
   if( xDiag )
   {
      SimpleVector& sxDiag = dynamic_cast<SimpleVector&>(*xDiag);
      for( int i = 0; i < locnx; i++)
         SC[i][i] += sxDiag[i];
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with F
   /////////////////////////////////////////////////////////////
   if( locmyl > 0 )
   {
      const double* MF0 = F0.M();
      const int* krowF0 = F0.krowM();
      const int* jcolF0 = F0.jcolM();

      for( int rowF0 = 0; rowF0 < locmyl; ++rowF0)
      {
         for( int k = krowF0[rowF0]; k < krowF0[rowF0 + 1]; ++k )
         {
            const int colF0 = jcolF0[k];
            assert(colF0 < locnx);

            const double valF0 = MF0[k];
            SC[locnx + rowF0][colF0] += valF0;
         }
      }
   }


   /////////////////////////////////////////////////////////////
   // update the KKT with G and put z diagonal
   /////////////////////////////////////////////////////////////
   if( locmzl > 0 )
   {
      assert(zDiagLinkCons);
      SimpleVector& szDiagLinkCons = dynamic_cast<SimpleVector&>(*zDiagLinkCons);

      const double* MG0 = G0.M();
      const int* krowG0 = G0.krowM();
      const int* jcolG0 = G0.jcolM();

      for( int rowG0 = 0; rowG0 < locmzl; ++rowG0 )
      {
         SC[locnx + locmyl + rowG0][locnx + locmyl + rowG0] += szDiagLinkCons[rowG0];

         for( int k = krowG0[rowG0]; k < krowG0[rowG0 + 1]; ++k )
         {
            const int colG0 = jcolG0[k];
            assert( colG0 < locnx);

            const double valG0 = MG0[k];
            SC[locnx + locmyl + rowG0][colG0] += valG0;
         }
      }
   }
}

void sLinsysRootBordered::solveReduced( sData *prob, SimpleVector& b)
{
   assert( 0 && "TODO: implement..");
}

void sLinsysRootBordered::computeSchurCompRightHandSide( StochVector& rhs )
{
   /* solve inner system */
   this->children[0]->solveCompressed( *rhs.children[0] );

   assert( rhs.vec );
   assert( !rhs.vecl );

   if( PIPS_MPIgetRank(mpiComm) != 0 )
      rhs.vec->setToZero();

   BorderLinsys border( *dynamic_cast<BorderedSymMatrix&>(*data->Q).border_vertical,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_bottom,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_bottom);

   this->children[0]->addBorderTimesRhsToB0( *rhs.children[0], dynamic_cast<SimpleVector&>(*rhs.vec), border );

   PIPS_MPIsumArrayInPlace( dynamic_cast<SimpleVector&>(*rhs.vec).elements(), rhs.vec->length(), mpiComm );
}


void sLinsysRootBordered::Lsolve(sData *prob, OoqpVector& x)
{
   assert( is_hierarchy_root );
   assert( children.size() == 1 );
   assert( prob );

   StochVector& xs = dynamic_cast<StochVector&>(x);
   assert( xs.children.size() == 1 );
   assert( data->children.size() == 1 );
   assert( xs.vec );

   computeSchurCompRightHandSide( xs );
   solver->solve( *xs.vec );
}

void sLinsysRootBordered::Dsolve(sData *prob, OoqpVector& x)
{
   assert( false && "TODO: implement" );
}

/* create kkt used to store Schur Complement of border layer */
SymMatrix* sLinsysRootBordered::createKKT(sData* prob)
{
   const int n = locnx + locmyl + locmzl;

   return new DenseSymMatrix(n);
}

void sLinsysRootBordered::assembleLocalKKT(sData* prob)
{
   assert( is_hierarchy_root );
   assert( !hasSparseKkt );
   assert( children.size() == 1 );

   // assemble complete inner KKT from children
   DenseSymMatrix& SC = dynamic_cast<DenseSymMatrix&>(*kkt);

   this->children[0]->addInnerToHierarchicalSchurComplement(SC, prob);
}

/* since we have only one child we will not allreduce anything */
void sLinsysRootBordered::reduceKKT(sData* prob)
{
   return;
}

DoubleLinearSolver* sLinsysRootBordered::createSolver(sData* prob, SymMatrix* kktmat_)
{
   DenseSymMatrix* kktmat = dynamic_cast<DenseSymMatrix*>(kktmat_);
   return new DeSymIndefSolver(kktmat);
   //return new DeSymIndefSolver2(kktmat, locnx); // saddle point solver
   //return new DeSymPSDSolver(kktmat);
}
