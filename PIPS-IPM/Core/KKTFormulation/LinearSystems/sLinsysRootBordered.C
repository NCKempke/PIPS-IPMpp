/*
 * sLinsysRootBordered.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "sLinsysRootBordered.h"
#include "BorderedSymmetricMatrix.h"
#include "PIPSIPMppOptions.h"
#include "DistributedFactory.hpp"

sLinsysRootBordered::sLinsysRootBordered(const DistributedFactory& factory_, DistributedProblem* prob_) : sLinsysRootAug(factory_, prob_, true) {
   assert(locmz == 0);
   assert(!hasSparseKkt);

   if (PIPS_MPIgetRank(mpiComm) == 0) {
      print_solver_regularization_and_sc_info("sLinsysRootBordered");
   }
}

void sLinsysRootBordered::computeSchurCompRightHandSide(const DistributedVector<double>& rhs_inner, SimpleVector<double>& b0) {
   if (!sol_inner)
      sol_inner.reset(dynamic_cast<DistributedVector<double>*>(rhs_inner.clone_full()));
   else
      sol_inner->copyFrom(rhs_inner);

   /* solve inner system */
   children[0]->solveCompressed(*sol_inner);

   if (PIPS_MPIgetRank(mpiComm) != 0)
      b0.setToZero();

   BorderLinsys border(*dynamic_cast<const BorderedSymmetricMatrix&>(*data->hessian).border_vertical, *dynamic_cast<const BorderedMatrix&>(*data->equality_jacobian).border_left,
         *dynamic_cast<const BorderedMatrix&>(*data->inequality_jacobian).border_left, 0, *dynamic_cast<const BorderedMatrix&>(*data->equality_jacobian).border_bottom,
         *dynamic_cast<const BorderedMatrix&>(*data->inequality_jacobian).border_bottom);

   children[0]->addBorderTimesRhsToB0(*sol_inner, b0, border);

   PIPS_MPIsumArrayInPlace(b0.elements(), b0.length(), mpiComm);
}

void sLinsysRootBordered::computeInnerSystemRightHandSide(DistributedVector<double>& rhs_inner, const SimpleVector<double>& b0, bool) {
   BorderLinsys border(*dynamic_cast<const BorderedSymmetricMatrix&>(*data->hessian).border_vertical, *dynamic_cast<const BorderedMatrix&>(*data->equality_jacobian).border_left,
         *dynamic_cast<const BorderedMatrix&>(*data->inequality_jacobian).border_left, 0, *dynamic_cast<const BorderedMatrix&>(*data->equality_jacobian).border_bottom,
         *dynamic_cast<const BorderedMatrix&>(*data->inequality_jacobian).border_bottom);

   children[0]->addBorderX0ToRhs(rhs_inner, b0, border);
}

/*  The solve :
 *    [  K  B  ] [  x  ] = [  b  ]
 *    [ B^T K0 ] [ x_0 ] = [ b_0 ]
 */
/* forms right hand side for schur system \tilda{b_0} = b_0 - B^T * K^-1 b and in doing so solves K^-1 b */
void sLinsysRootBordered::Lsolve(Vector<double>& x) {
   assert(is_hierarchy_root);
   assert(children.size() == 1);

   auto& xs = dynamic_cast<DistributedVector<double>&>(x);
   assert(xs.children.size() == 1);
   assert(data->children.size() == 1);
   DistributedVector<double>& b = *dynamic_cast<DistributedVector<double>&>(x).children[0];

   assert(xs.first);
   assert(!xs.last);
   auto& b0 = dynamic_cast<SimpleVector<double>&>(*xs.first);

   computeSchurCompRightHandSide(b, b0);
}

/* does Schur Complement solve and computes SC x_0 = \tilda{b_0} = ( K0 - B^T K B ) x_0 */
void sLinsysRootBordered::Dsolve(Vector<double>& x) {
   assert(is_hierarchy_root);
   assert(children.size() == 1);

   auto& xs = dynamic_cast<DistributedVector<double>&>(x);
   assert(xs.children.size() == 1);
   assert(data->children.size() == 1);
   assert(xs.first);
   auto& b0 = dynamic_cast<SimpleVector<double>&>(*xs.first);

   solver->solve(b0);
}

/* back substitute x_0 : K x = b - B x_0 and solve for x */
void sLinsysRootBordered::Ltsolve(Vector<double>& x) {
   assert(is_hierarchy_root);
   assert(children.size() == 1);

   auto& xs = dynamic_cast<DistributedVector<double>&>(x);
   assert(xs.children.size() == 1);
   assert(data->children.size() == 1);
   DistributedVector<double>& b = *dynamic_cast<DistributedVector<double>&>(x).children[0];

   assert(xs.first);
   auto& b0 = dynamic_cast<SimpleVector<double>&>(*xs.first);

   computeInnerSystemRightHandSide(b, b0, false);

   children[0]->solveCompressed(b);
}

void sLinsysRootBordered::assembleLocalKKT() {
   assert(allreduce_kkt);
   assert(is_hierarchy_root);
   assert(!hasSparseKkt);
   assert(children.size() == 1);
   assert(locmy >= 0);
   assert(data->children.size() == 1);

   // assemble complete inner KKT from children
   auto& SC = dynamic_cast<DenseSymmetricMatrix&>(*kkt);

   BorderLinsys B(*dynamic_cast<const BorderedSymmetricMatrix&>(*data->hessian).border_vertical, *dynamic_cast<const BorderedMatrix&>(*data->equality_jacobian).border_left,
         *dynamic_cast<const BorderedMatrix&>(*data->inequality_jacobian).border_left, 0, *dynamic_cast<const BorderedMatrix&>(*data->equality_jacobian).border_bottom,
         *dynamic_cast<const BorderedMatrix&>(*data->inequality_jacobian).border_bottom);
   assert(data->children.size() == 1);

   std::vector<BorderMod> border_mod;

   children[0]->addBlTKiInvBrToRes(SC, B, B, border_mod, true, false);
}
