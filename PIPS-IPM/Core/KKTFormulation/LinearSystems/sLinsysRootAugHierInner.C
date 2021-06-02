/*
 * sLinsysRootAugHierInner.C
 *
 *  Created on: 27.01.2021
 *      Author: bzfkempk
 */

#include "sLinsysRootAugHierInner.h"

sLinsysRootAugHierInner::sLinsysRootAugHierInner(DistributedFactory* factory, DistributedQP* prob_,
   std::shared_ptr<Vector<double>> dd_, std::shared_ptr<Vector<double>> dq_,
   std::shared_ptr<Vector<double>> nomegaInv_, std::shared_ptr<Vector<double>> regP_,
   std::shared_ptr<Vector<double>> regDy_, std::shared_ptr<Vector<double>> regDz_, std::shared_ptr<Vector<double>> rhs_)
   : sLinsysRootAug(
   factory, prob_, dynamic_cast<DistributedVector<double>&>(*dd_).first,
   dynamic_cast<DistributedVector<double>&>(*dq_).first,
   dynamic_cast<DistributedVector<double>&>(*nomegaInv_).first, dynamic_cast<DistributedVector<double>&>(*regP_).first,
   dynamic_cast<DistributedVector<double>&>(*regDy_).first, dynamic_cast<DistributedVector<double>&>(*regDz_).first,
   rhs_, true) {
   assert(locnx == 0 && locmy == 0 && locmz == 0);
   assert(hasSparseKkt);

   static bool printed = false;
   if (!printed && PIPS_MPIgetRank() == 0) {
      print_solver_regularization_and_sc_info("sLinsysRootAugHierInner");
      printed = true;
   }
}

void sLinsysRootAugHierInner::assembleLocalKKT() {
   for (size_t c = 0; c < children.size(); ++c) {
#ifdef STOCH_TESTING
      g_scenNum = c;
#endif
      if (children[c]->mpiComm == MPI_COMM_NULL)
         continue;

      children[c]->stochNode->resMon.recFactTmChildren_start();
      //---------------------------------------------
      addTermToSchurCompl(c, false);
      //---------------------------------------------
      children[c]->stochNode->resMon.recFactTmChildren_stop();
   }
}

void sLinsysRootAugHierInner::Ltsolve(Vector<double>& x) {
   auto& b = dynamic_cast<DistributedVector<double>&>(x);
   auto& b0 = dynamic_cast<SimpleVector<double>&>(*b.first);

   //dumpRhs(0, "sol",  b0);
   SimpleVector<double>& z0 = b0; //just another name, for clarity

   for (size_t it = 0; it < children.size(); it++)
      children[it]->Ltsolve2(*b.children[it], z0, false);
}

void sLinsysRootAugHierInner::Ltsolve2(DistributedVector<double>& x, SimpleVector<double>& x0, bool use_local_RAC) {
   assert(pipsipmpp_options::get_bool_parameter("HIERARCHICAL"));

   auto& b = dynamic_cast<DistributedVector<double>&>(x);

   computeInnerSystemRightHandSide(b, x0, use_local_RAC);
   solveCompressed(x);
}

void
sLinsysRootAugHierInner::LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool two_link_border,
   int begin_cols, int end_cols) {
   LsolveHierarchyBorder(result, Br, Br_mod_border, false, two_link_border, begin_cols, end_cols);
}

void sLinsysRootAugHierInner::LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl,
   BorderLinsys& Br,
   std::vector<BorderMod>& br_mod_border, bool sym_res, bool sparse_res, int begin_cols, int end_cols) {
   if (Bl.isEmpty() || (Br.isEmpty() && br_mod_border.empty()))
      return;

   LtsolveHierarchyBorder(res, X0, Bl, Br, br_mod_border, sym_res, sparse_res, false, begin_cols, end_cols);
}

void
sLinsysRootAugHierInner::computeInnerSystemRightHandSide(DistributedVector<double>& rhs_inner,
   const SimpleVector<double>& b0, bool use_local_RAC) {
   BorderLinsys Border(0, dynamic_cast<StripMatrix&>(*dynamic_cast<const DistributedMatrix&>(*data->A).Blmat),
      dynamic_cast<StripMatrix&>(*dynamic_cast<const DistributedMatrix&>(*data->C).Blmat), use_local_RAC);

   if (Border.isEmpty())
      return;

   addBorderX0ToRhs(rhs_inner, b0, Border);
}

/* compute Schur rhs b0 - sum Bi^T Ki^-1 bi for all children */
void sLinsysRootAugHierInner::Lsolve(Vector<double>& x) {
   assert(!is_hierarchy_root);

   auto& b = dynamic_cast<DistributedVector<double>&>(x);
   assert(children.size() == b.children.size());

   auto& b0 = dynamic_cast<SimpleVector<double>&>(*b.first);
   assert(!b.last);

   if (iAmDistrib && PIPS_MPIgetRank(mpiComm) > 0)
      b0.setToZero();

   // compute Bi^T Ki^-1 rhs_i and sum it up
   for (size_t it = 0; it < children.size(); it++)
      children[it]->addLniziLinkCons(b0, *b.children[it], false);

   if (iAmDistrib)
      PIPS_MPIsumArrayInPlace(b0.elements(), b0.length(), mpiComm);
}

void sLinsysRootAugHierInner::addLniziLinkCons(Vector<double>& z0_, Vector<double>& zi, bool use_local_RAC) {
   assert(zi.isKindOf(kStochVector));
   auto& z0 = dynamic_cast<SimpleVector<double>&>(z0_);

   if (!sol_inner)
      sol_inner.reset(dynamic_cast<DistributedVector<double>*>(zi.cloneFull()));
   else
      sol_inner->copyFrom(zi);

   /* solve system */
   solveCompressed(*sol_inner);

   BorderLinsys Bl(0, dynamic_cast<StripMatrix&>(*dynamic_cast<const DistributedMatrix&>(*data->A).Blmat),
      dynamic_cast<StripMatrix&>(*dynamic_cast<const DistributedMatrix&>(*data->C).Blmat), use_local_RAC);

   addBorderTimesRhsToB0(*sol_inner, z0, Bl);
}

void sLinsysRootAugHierInner::addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0,
   BorderLinsys& border) {
   assert(rhs.children.size() == children.size());
   assert(border.F.children.size() == children.size());
   for (size_t i = 0; i < children.size(); ++i) {
      BorderLinsys child_border = getChild(border, i);
      if (child_border.isEmpty())
         continue;

      children[i]->addBorderTimesRhsToB0(*rhs.children[i], b0, child_border);
   }

   /* add schur complement part */
   if ((border.has_RAC || border.use_local_RAC) && (PIPS_MPIgetSize(mpiComm) == 0 || PIPS_MPIgetRank(mpiComm) == 0)) {
      if (border.has_RAC) {
         assert(border.A.last);
         assert(border.C.last);
      }

      const SparseMatrix& F_border = border.has_RAC ? dynamic_cast<SparseMatrix&>(*border.A.last).getTranspose()
         : data->getLocalF().getTranspose();
      const SparseMatrix& G_border = border.has_RAC ? dynamic_cast<SparseMatrix&>(*border.C.last).getTranspose()
         : data->getLocalG().getTranspose();

      const auto[mFb, nFb] = F_border.n_rows_columns();
      int nGb = G_border.n_columns();

      assert(mFb == G_border.n_rows());
      assert(rhs.first);
      assert(rhs.first->length() == nFb + nGb);

      assert(b0.length() >= mFb);

      auto& zi = dynamic_cast<SimpleVector<double>&>(*rhs.first);

      SimpleVector<double> zi1(&zi[0], nFb);
      SimpleVector<double> zi2(&zi[nFb], nGb);

      SimpleVector<double> b1(&b0[0], mFb);

      F_border.mult(1.0, b1, -1.0, zi1);
      G_border.mult(1.0, b1, -1.0, zi2);
   }
}

void sLinsysRootAugHierInner::addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0,
   BorderLinsys& border) {
   assert(rhs.children.size() == children.size());
   assert(border.F.children.size() == children.size());

   for (size_t i = 0; i < children.size(); ++i) {
      BorderLinsys child_border = getChild(border, i);
      if (child_border.isEmpty())
         continue;
      children[i]->addBorderX0ToRhs(*rhs.children[i], x0, child_border);
   }

   if (border.has_RAC) {
      assert(border.A.last);
      assert(border.C.last);
   }

   /* add schur complement part */
   if (border.has_RAC || border.use_local_RAC) {
      const SparseMatrix& F_border = border.has_RAC ? dynamic_cast<SparseMatrix&>(*border.A.last) : data->getLocalF();
      const SparseMatrix& G_border = border.has_RAC ? dynamic_cast<SparseMatrix&>(*border.C.last) : data->getLocalG();
      const auto[mFb, nFb] = F_border.n_rows_columns();
      const auto mGb = G_border.n_rows();

      assert(rhs.first);
      assert(rhs.first->length() == mFb + mGb);
      assert(nFb == G_border.n_columns());
      assert(x0.length() >= nFb);

      auto& rhs0 = dynamic_cast<SimpleVector<double>&>(*rhs.first);

      SimpleVector<double> rhs01(&rhs0[0], mFb);
      SimpleVector<double> rhs02(&rhs0[mFb], mGb);

      SimpleVector<double> x1((double*) &x0[0], nFb);

      F_border.mult(1.0, rhs01, -1.0, x1);
      G_border.mult(1.0, rhs02, -1.0, x1);
   }
}

/* res += [ B_inner^T K_i^{-1} (Br - SUM_j Brmodj Xj) ]^T = (Br^T - SUM_j Xj^T Brmodj^T) K_i^{-1} B_inner from begin_cols to end_cols in */
void
sLinsysRootAugHierInner::addInnerBorderKiInvBrToRes(AbstractMatrix& result, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool use_local_RAC,
   bool sparse_res, bool sym_res, int begin_cols, int end_cols, int n_empty_rows_inner_border) {
   const auto mres = result.n_rows();
   assert(dynamic_cast<const DistributedMatrix&>(*data->A).Blmat->is_a(kStripMatrix));
   assert(dynamic_cast<const DistributedMatrix&>(*data->C).Blmat->is_a(kStripMatrix));

   BorderLinsys Bl(n_empty_rows_inner_border, data->getLocalFBorder(), data->getLocalGBorder(), use_local_RAC);
   if (Bl.isEmpty() || (Br.isEmpty() && Br_mod_border.empty()))
      return;

   const int n_buffer = locnx + locmy + locmz + locmyl + locmzl;
   // buffer b0 for blockwise computation of Br0 - SUM_i  Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ), stored in transposed form (for quick access of cols in solve)
   // dense since we have no clue about any structure in the system and Xij are dense
#ifndef NDEBUG
   const int m_buffer = allocateAndZeroBlockedComputationsBuffer(mres, n_buffer);
   assert(end_cols - begin_cols <= m_buffer);
#else
   allocateAndZeroBlockedComputationsBuffer(mres, n_buffer);
#endif

   addBlTKiInvBrToResBlockwise(result, Bl, Br, Br_mod_border, sym_res, sparse_res, *buffer_blocked_hierarchical,
      begin_cols, end_cols);
}

void sLinsysRootAugHierInner::addTermToSchurComplBlocked(bool sparseSC, SymmetricMatrix& SC, bool use_local_RAC,
   int n_empty_rows_inner_border) {
   assert(data->getLocalFBorder().n_columns() == data->getLocalGBorder().n_columns());
   assert(dynamic_cast<const DistributedMatrix&>(*data->A).Blmat->is_a(kStripMatrix));
   assert(dynamic_cast<const DistributedMatrix&>(*data->C).Blmat->is_a(kStripMatrix));

   BorderLinsys Bl(n_empty_rows_inner_border, data->getLocalFBorder(), data->getLocalGBorder(), use_local_RAC);
   BorderLinsys Br(n_empty_rows_inner_border, data->getLocalFBorder(), data->getLocalGBorder(), use_local_RAC);
   if (Bl.isEmpty())
      return;

   std::vector<BorderMod> border_mod;

   const bool sc_is_symmetric = true;
   addBlTKiInvBrToRes(SC, Bl, Br, border_mod, sc_is_symmetric, sparseSC);
}

/* compute (Bli^T X_i = Bli^T Ki^-1 (Bri - Bi_{inner} X0))^T and add it to SC */
void sLinsysRootAugHierInner::LniTransMultHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl,
   BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool sparse_res, bool sym_res, bool use_local_RAC, int begin_cols,
   int end_cols,
   int n_empty_rows_inner_border) {

   if (Bl.isEmpty() || (Br.isEmpty() && Br_mod_border.empty()))
      return;

   BorderLinsys B_inner(n_empty_rows_inner_border, data->getLocalFBorder(), data->getLocalGBorder(), use_local_RAC);

   std::vector<BorderMod> border_mod(Br_mod_border);
   if (!B_inner.isEmpty()) {
      BorderMod B_inner_mod(B_inner, X0);
      border_mod.push_back(B_inner_mod);
   }

   const int n_buffer = locnx + locmy + locmz + locmyl + locmzl;

   // buffer b0 for blockwise computation of Br0 - SUM_i  Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ), stored in transposed form (for quick access of cols in solve)
   // dense since we have no clue about any structure in the system and Xij are dense
   const int m_buffer = allocateAndZeroBlockedComputationsBuffer(res.n_rows(), n_buffer);
   (void) m_buffer;
   assert(end_cols - begin_cols <= m_buffer);

   addBlTKiInvBrToResBlockwise(res, Bl, Br, border_mod, sym_res, sparse_res, *buffer_blocked_hierarchical, begin_cols,
      end_cols);
}

void sLinsysRootAugHierInner::put_primal_diagonal() {
   assert(primal_diagonal);
   assert(primal_diagonal->isKindOf(kStochVector));

   auto& primal_diagonal_loc = dynamic_cast<DistributedVector<double>&>(*primal_diagonal);
   xDiag = primal_diagonal_loc.first.get();

   assert(children.size() == primal_diagonal_loc.children.size());
   for (auto& it : children)
      it->put_primal_diagonal();
}


void sLinsysRootAugHierInner::put_dual_inequalites_diagonal() {
   assert(nomegaInv);
   assert(nomegaInv->isKindOf(kStochVector));

   auto& dual_inequality_diagonal_loc = dynamic_cast<DistributedVector<double>&>(*nomegaInv);

   zDiag = dual_inequality_diagonal_loc.first.get();
   zDiagLinkCons = dual_inequality_diagonal_loc.last.get();

   assert(children.size() == dual_inequality_diagonal_loc.children.size());
   for (auto& it : children)
      it->put_dual_inequalites_diagonal();
}
