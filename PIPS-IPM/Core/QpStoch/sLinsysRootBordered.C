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
#include "StochOptions.h"

#include "sLinsysRootAug.h"
#include "sFactory.h"

sLinsysRootBordered::sLinsysRootBordered(sFactory * factory_, sData * prob_)
  : sLinsysRoot(factory_, prob_, true)
{
   assert(locmyl >= 0 && locmzl >= 0);

   kkt.reset(createKKT(prob_));

   solver.reset( createSolver(prob_, kkt.get()) );
}

void sLinsysRootBordered::finalizeKKT(/* const */sData* prob, Variables*)
{

   /* Add corner block
    * [ Q0 F0T G0T  ]
    * [ F0  0   0   ]
    * [ G0  0 OmN+1 ]
    */
   assert( prob->isHierarchyRoot() );

   const SparseGenMatrix& F0 = *dynamic_cast<const BorderedGenMatrix&>(*prob->A).bottom_left_block;
   const SparseGenMatrix& G0 = *dynamic_cast<const BorderedGenMatrix&>(*prob->C).bottom_left_block;
   const SparseSymMatrix& Q0 = dynamic_cast<const SparseSymMatrix&>(*dynamic_cast<const BorderedSymMatrix&>(*prob->Q).top_left_block);

   DenseSymMatrix& SC = dynamic_cast<DenseSymMatrix&>(*kkt);
   int mSC, nSC; SC.getSize(mSC, nSC);
   assert( mSC == nSC );

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
            SC[colF0][locnx + rowF0] += valF0;
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
            SC[colG0][locnx + locmyl + rowG0] += valG0;
         }
      }
   }
}

void sLinsysRootBordered::computeSchurCompRightHandSide( const StochVector& rhs_inner, SimpleVector& b0 )
{
   if( !sol_inner )
      sol_inner.reset(dynamic_cast<StochVector*>(rhs_inner.cloneFull()));
   else
      sol_inner->copyFrom(rhs_inner);

   /* solve inner system */
   children[0]->solveCompressed( *sol_inner );

   if( PIPS_MPIgetRank(mpiComm) != 0 )
      b0.setToZero();

   BorderLinsys border( *dynamic_cast<BorderedSymMatrix&>(*data->Q).border_vertical,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_bottom,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_bottom);

   children[0]->addBorderTimesRhsToB0( *sol_inner, b0, border );

   PIPS_MPIsumArrayInPlace( b0.elements(), b0.length(), mpiComm );
}

void sLinsysRootBordered::computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& b0, bool )
{
   BorderLinsys border( *dynamic_cast<BorderedSymMatrix&>(*data->Q).border_vertical,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_bottom,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_bottom);

   children[0]->addBorderX0ToRhs( rhs_inner, b0, border );
}


/*  The solve :
 *    [  K  B  ] [  x  ] = [  b  ]
 *    [ B^T K0 ] [ x_0 ] = [ b_0 ]
 */
/* forms right hand side for schur system \tilda{b_0} = b_0 - B^T * K^-1 b and in doing so solves K^-1 b */
void sLinsysRootBordered::Lsolve(sData* , OoqpVector& x)
{
   assert( is_hierarchy_root );
   assert( children.size() == 1 );

   StochVector& xs = dynamic_cast<StochVector&>(x);
   assert( xs.children.size() == 1 );
   assert( data->children.size() == 1 );
   StochVector& b = *dynamic_cast<StochVector&>(x).children[0];

   assert( xs.vec );
   assert( !xs.vecl );
   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*xs.vec);

   computeSchurCompRightHandSide( b, b0 );
}

/* does Schur Complement solve and computes SC x_0 = \tilda{b_0} = ( K0 - B^T K B ) x_0 */
void sLinsysRootBordered::Dsolve(sData*, OoqpVector& x)
{
   assert( is_hierarchy_root );
   assert( children.size() == 1 );

   StochVector& xs = dynamic_cast<StochVector&>(x);
   assert( xs.children.size() == 1 );
   assert( data->children.size() == 1 );
   assert( xs.vec );
   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*xs.vec);

   solver->solve( b0 );
}

/* back substitute x_0 : K x = b - B x_0 and solve for x */
void sLinsysRootBordered::Ltsolve(sData*, OoqpVector& x)
{
   assert( is_hierarchy_root );
   assert( children.size() == 1 );

   StochVector& xs = dynamic_cast<StochVector&>(x);
   assert( xs.children.size() == 1 );
   assert( data->children.size() == 1 );
   StochVector& b = *dynamic_cast<StochVector&>(x).children[0];

   assert( xs.vec );
   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*xs.vec);

   computeInnerSystemRightHandSide( b, b0, false );

   children[0]->solveCompressed( b );
}

/* create kkt used to store Schur Complement of border layer */
SymMatrix* sLinsysRootBordered::createKKT(sData*)
{
   const int n = locnx + locmyl + locmzl;

   if( PIPS_MPIgetRank(mpiComm) == 0 )
      std::cout << "sLinsysRootBordered: getSchurCompMaxNnz " << n * n << "\n";

   return new DenseSymMatrix(n);
}

void sLinsysRootBordered::assembleLocalKKT(sData* prob)
{
   assert(allreduce_kkt);
   assert( is_hierarchy_root );
   assert( !hasSparseKkt );
   assert( children.size() == 1 );

   // assemble complete inner KKT from children
   DenseSymMatrix& SC = dynamic_cast<DenseSymMatrix&>(*kkt);

   assert( prob->children.size() == 1 );

   BorderLinsys B( *dynamic_cast<BorderedSymMatrix&>(*prob->Q).border_vertical,
         *dynamic_cast<BorderedGenMatrix&>(*prob->A).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*prob->C).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*prob->A).border_bottom,
         *dynamic_cast<BorderedGenMatrix&>(*prob->C).border_bottom);
   std::vector<BorderMod> border_mod;

   children[0]->addBTKiInvBToSC(SC, B, B, border_mod, true, false);
}

/* since we have only one child we will not allreduce anything */
void sLinsysRootBordered::reduceKKT(sData*)
{
   if( iAmDistrib )
      allreduceMatrix(*kkt, false, true, mpiComm );
}

DoubleLinearSolver* sLinsysRootBordered::createSolver(sData*, const SymMatrix* kktmat_)
{
   const SolverTypeDense solver = pips_options::getSolverDense();
   const DenseSymMatrix* kktmat = dynamic_cast<const DenseSymMatrix*>(kktmat_);

   static bool printed = false;
   if( !printed && 0 == PIPS_MPIgetRank(mpiComm) )
      std::cout << "sLinsysRootBordered: using " << solver << " for dense border Schur complement\n";
   printed = true;

   if( solver == SolverTypeDense::SOLVER_DENSE_SYM_INDEF )
      return new DeSymIndefSolver(kktmat);
   else if( solver == SolverTypeDense::SOLVER_DENSE_SYM_INDEF_SADDLE_POINT )
      return new DeSymIndefSolver2(kktmat, locnx);
   else
   {
      assert( solver == SolverTypeDense::SOLVER_DENSE_SYM_PSD );
      return new DeSymPSDSolver(kktmat);
   }
}
