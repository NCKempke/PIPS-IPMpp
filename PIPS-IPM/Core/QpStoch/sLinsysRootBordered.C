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
  : sLinsysRootAug(factory_, prob_, true)
{
}

void sLinsysRootBordered::computeSchurCompRightHandSide( const StochVector& rhs_inner, SimpleVector& b0 )
{
   if( !sol_inner )
      sol_inner.reset(dynamic_cast<StochVector*>(rhs_inner.cloneFull()));
   else
      sol_inner->copyFrom(rhs_inner);

   /* solve inner system */
   this->children[0]->solveCompressed( *sol_inner );

   if( PIPS_MPIgetRank(mpiComm) != 0 )
      b0.setToZero();

   BorderLinsys border( *dynamic_cast<BorderedSymMatrix&>(*data->Q).border_vertical,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_bottom,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_bottom);

   this->children[0]->addBorderTimesRhsToB0( *sol_inner, b0, border );

   PIPS_MPIsumArrayInPlace( b0.elements(), b0.length(), mpiComm );
}

void sLinsysRootBordered::computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& b0 )
{
   BorderLinsys border( *dynamic_cast<BorderedSymMatrix&>(*data->Q).border_vertical,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_left,
         *dynamic_cast<BorderedGenMatrix&>(*data->A).border_bottom,
         *dynamic_cast<BorderedGenMatrix&>(*data->C).border_bottom);

   this->children[0]->addBorderX0ToRhs( rhs_inner, b0, border );
}


/*  The solve :
 *    [  K  B  ] [  x  ] = [  b  ]
 *    [ B^T K0 ] [ x_0 ] = [ b_0 ]
 */
/* forms right hand side for schur system \tilda{b_0} = b_0 - B^T * K^-1 b and in doing so solves K^-1 b */
void sLinsysRootBordered::Lsolve(sData *prob, OoqpVector& x)
{
   assert( is_hierarchy_root );
   assert( children.size() == 1 );
   assert( prob );

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
void sLinsysRootBordered::Dsolve(sData *prob, OoqpVector& x)
{
   assert( is_hierarchy_root );
   assert( children.size() == 1 );
   assert( prob );

   StochVector& xs = dynamic_cast<StochVector&>(x);
   assert( xs.children.size() == 1 );
   assert( data->children.size() == 1 );
   assert( xs.vec );
   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*xs.vec);

   solver->solve( b0 );
}

/* back substitute x_0 : K x = b - B x_0 and solve for x */
void sLinsysRootBordered::Ltsolve(sData* prob, OoqpVector& x)
{
   assert( is_hierarchy_root );
   assert( children.size() == 1 );
   assert( prob );

   StochVector& xs = dynamic_cast<StochVector&>(x);
   assert( xs.children.size() == 1 );
   assert( data->children.size() == 1 );
   StochVector& b = *dynamic_cast<StochVector&>(x).children[0];

   assert( xs.vec );
   SimpleVector& b0 = dynamic_cast<SimpleVector&>(*xs.vec);

   computeInnerSystemRightHandSide( b, b0 );

   this->children[0]->solveCompressed( b );
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
