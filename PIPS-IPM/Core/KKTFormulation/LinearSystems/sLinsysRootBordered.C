/*
 * sLinsysRootBordered.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "sLinsysRootBordered.h"
#include "BorderedSymmetricMatrix.h"
#include "DeSymIndefSolver.h"
#include "DeSymIndefSolver2.h"
#include "DeSymPSDSolver.h"
#include "PIPSIPMppOptions.h"
#include "DistributedFactory.hpp"

sLinsysRootBordered::sLinsysRootBordered(DistributedFactory* factory_, DistributedProblem* prob_) : DistributedRootLinearSystem(factory_, prob_, true) {
   assert(locmyl >= 0 && locmzl >= 0);

   if (apply_regularization) {
      assert(false && "TODO : regularization not implemented for hierarchical approach currently");
   }
   kkt.reset(createKKT());

   solver.reset(createSolver(kkt.get()));
}

void sLinsysRootBordered::add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization,
      double dual_inequality_regularization) {
   assert(dynamic_cast<DistributedVector<double>*>(primal_regularization_diagonal.get()));
   assert(dynamic_cast<DistributedVector<double>*>(dual_equality_regularization_diagonal.get()));
   assert(dynamic_cast<DistributedVector<double>*>(dual_inequality_regularization_diagonal.get()));

   if (locnx > 0) {
      assert(dynamic_cast<DistributedVector<double>&>(*primal_regularization_diagonal).first);
      dynamic_cast<DistributedVector<double>&>(*primal_regularization_diagonal).first->addConstant(primal_regularization);
      kkt->diagonal_add_constant_from(0, locnx, primal_regularization);
   }

   if (locmyl > 0) {
      assert(dynamic_cast<DistributedVector<double>&>(*dual_equality_regularization_diagonal).last);
      dynamic_cast<DistributedVector<double>&>(*dual_equality_regularization_diagonal).last->addConstant(dual_equality_regularization);
      kkt->diagonal_add_constant_from(locnx, locmyl, dual_equality_regularization);
   }

   if (locmzl > 0) {
      assert(dynamic_cast<DistributedVector<double>&>(*dual_inequality_regularization_diagonal).last);
      dynamic_cast<DistributedVector<double>&>(*dual_inequality_regularization_diagonal).last->addConstant(dual_inequality_regularization);
      kkt->diagonal_add_constant_from(locnx + locmyl, locmzl, dual_inequality_regularization);
   }
}

void sLinsysRootBordered::finalizeKKT() {
   /* Add corner block
    * [ Q0   A0^T  F0T   G0T  ]
    * [ A0   xReg   0     0   ]
    * [ F0    0   xReg_l  0   ]
    * [ G0    0     0   OmN+1 ]
    */
   assert(data->isHierarchyRoot());

   const auto& F0 = *dynamic_cast<const BorderedMatrix&>(*data->A).bottom_left_block;
   const auto& G0 = *dynamic_cast<const BorderedMatrix&>(*data->C).bottom_left_block;
   const auto& Q0 = dynamic_cast<const SparseSymmetricMatrix&>(*dynamic_cast<const BorderedSymmetricMatrix&>(*data->Q).top_left_block);

   auto& SC = dynamic_cast<DenseSymmetricMatrix&>(*kkt);

   /////////////////////////////////////////////////////////////
   // update the KKT with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   const int* krowQ0 = Q0.krowM();
   const int* jcolQ0 = Q0.jcolM();
   const double* MQ = Q0.M();

   for (int rowQ = 0; rowQ < locnx; rowQ++) {
      for (int k = krowQ0[rowQ]; k < krowQ0[rowQ + 1]; ++k) {
         const int colQ = jcolQ0[k];
         if (rowQ == colQ)
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
   if (locnx > 0) {
      assert(xDiag);
      SC.atAddDiagonal(0, *xDiag);
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with F
   /////////////////////////////////////////////////////////////
   if (locmyl > 0) {
      SC.add_matrix_at(dynamic_cast<const SparseMatrix&>(F0), locnx, 0);
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with G and put z diagonal
   /////////////////////////////////////////////////////////////
   if (locmzl > 0) {
      assert(zDiagLinkCons);

      SC.add_matrix_at(dynamic_cast<const SparseMatrix&>(G0), locnx + locmyl, 0);
      SC.atAddDiagonal(locnx + locmyl, *zDiagLinkCons);
   }
}

void sLinsysRootBordered::computeSchurCompRightHandSide(const DistributedVector<double>& rhs_inner, SimpleVector<double>& b0) {
   if (!sol_inner)
      sol_inner.reset(dynamic_cast<DistributedVector<double>*>(rhs_inner.cloneFull()));
   else
      sol_inner->copyFrom(rhs_inner);

   /* solve inner system */
   children[0]->solveCompressed(*sol_inner);

   if (PIPS_MPIgetRank(mpiComm) != 0)
      b0.setToZero();

   BorderLinsys border(*dynamic_cast<const BorderedSymmetricMatrix&>(*data->Q).border_vertical, *dynamic_cast<const BorderedMatrix&>(*data->A).border_left,
         *dynamic_cast<const BorderedMatrix&>(*data->C).border_left, 0, *dynamic_cast<const BorderedMatrix&>(*data->A).border_bottom,
         *dynamic_cast<const BorderedMatrix&>(*data->C).border_bottom);

   children[0]->addBorderTimesRhsToB0(*sol_inner, b0, border);

   PIPS_MPIsumArrayInPlace(b0.elements(), b0.length(), mpiComm);
}

void sLinsysRootBordered::computeInnerSystemRightHandSide(DistributedVector<double>& rhs_inner, const SimpleVector<double>& b0, bool) {
   BorderLinsys border(*dynamic_cast<const BorderedSymmetricMatrix&>(*data->Q).border_vertical, *dynamic_cast<const BorderedMatrix&>(*data->A).border_left,
         *dynamic_cast<const BorderedMatrix&>(*data->C).border_left, 0, *dynamic_cast<const BorderedMatrix&>(*data->A).border_bottom,
         *dynamic_cast<const BorderedMatrix&>(*data->C).border_bottom);

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

/* create kkt used to store Schur Complement of border layer */
SymmetricMatrix* sLinsysRootBordered::createKKT() {
   assert(stochNode->mz() == -1);
   //assert(locmy >= 0);
   const int n = locnx + locmyl + locmzl;

   if (PIPS_MPIgetRank(mpiComm) == 0)
      std::cout << "sLinsysRootBordered: getSchurCompMaxNnz " << n * n << "\n";

   return new DenseSymmetricMatrix(n);
}

void sLinsysRootBordered::assembleLocalKKT() {
   assert(allreduce_kkt);
   assert(is_hierarchy_root);
   assert(!hasSparseKkt);
   assert(children.size() == 1);

   // assemble complete inner KKT from children
   auto& SC = dynamic_cast<DenseSymmetricMatrix&>(*kkt);

   assert(data->children.size() == 1);

   BorderLinsys B(*dynamic_cast<const BorderedSymmetricMatrix&>(*data->Q).border_vertical, *dynamic_cast<const BorderedMatrix&>(*data->A).border_left,
         *dynamic_cast<const BorderedMatrix&>(*data->C).border_left, 0, *dynamic_cast<const BorderedMatrix&>(*data->A).border_bottom,
         *dynamic_cast<const BorderedMatrix&>(*data->C).border_bottom);
   std::vector<BorderMod> border_mod;

   children[0]->addBlTKiInvBrToRes(SC, B, B, border_mod, true, false);
}

/* since we have only one child we will not allreduce anything */
void sLinsysRootBordered::reduceKKT() {
   if (iAmDistrib)
      allreduceMatrix(*kkt, false, true, mpiComm);
}

DoubleLinearSolver* sLinsysRootBordered::createSolver(const SymmetricMatrix* kktmat_) {
   const SolverTypeDense solver = pipsipmpp_options::get_solver_dense();
   const auto& kktmat = dynamic_cast<const DenseSymmetricMatrix&>(*kktmat_);

   static bool printed = false;
   if (!printed && 0 == PIPS_MPIgetRank(mpiComm))
      std::cout << "sLinsysRootBordered: using " << solver << " for dense border Schur complement\n";
   printed = true;

   if (solver == SolverTypeDense::SOLVER_DENSE_SYM_INDEF)
      return new DeSymIndefSolver(kktmat);
   else if (solver == SolverTypeDense::SOLVER_DENSE_SYM_INDEF_SADDLE_POINT)
      return new DeSymIndefSolver2(kktmat, locnx);
   else {
      assert(solver == SolverTypeDense::SOLVER_DENSE_SYM_PSD);
      return new DeSymPSDSolver(kktmat);
   }
}
