/*
 * StochPresolverBase.cpp
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

//#define PIPS_DEBUG
#include "StochPresolverBase.h"

#include "StochOptions.h"
#include "pipsdef.h"
#include "StochVectorUtilities.h"
#include <cassert>

StochPresolverBase::StochPresolverBase(PresolveData& presolve_data, const DistributedQP& origProb) : my_rank(PIPS_MPIgetRank(MPI_COMM_WORLD)),
      distributed(PIPS_MPIgetDistributed(MPI_COMM_WORLD)), verbosity(pips_options::getIntParameter("PRESOLVE_VERBOSITY")),
      INF_NEG(-pips_options::getDoubleParameter("PRESOLVE_INFINITY")), INF_POS(pips_options::getDoubleParameter("PRESOLVE_INFINITY")),
      n_linking_vars(dynamic_cast<const DistributedVector<double>&>(*origProb.g).first->length()), n_linking_rows_eq(
            dynamic_cast<const DistributedVector<double>&>(*origProb.bA).last
            ? dynamic_cast<const DistributedVector<double>&>(*origProb.bA).last->length() : 0), n_linking_rows_ineq(
            dynamic_cast<const DistributedVector<double>&>(*origProb.iclow).last
            ? dynamic_cast<const DistributedVector<double>&>(*origProb.iclow).last->length() : 0), presolve_data(presolve_data), origProb(origProb) {
   localNelims = 0;
   nChildren = presolve_data.getNChildren();

   setPointersToNull();
}

StochPresolverBase::~StochPresolverBase() {
}

void StochPresolverBase::setPointersToNull() {
   currAmat = currAmatTrans = currBmat = currBmatTrans = currBlmat = currBlmatTrans = nullptr;

   currxlowParent = currIxlowParent = currxuppParent = currIxuppParent = currxlowChild = currIxlowChild = currxuppChild = currIxuppChild = currEqRhs = currIneqLhs = currIclow = currIneqRhs = currIcupp = currEqRhsLink = currIneqLhsLink = currIclowLink = currIneqRhsLink = currIcuppLink = currgParent = currgChild = nullptr;

   currNnzRow = currNnzRowLink = currNnzColParent = currNnzColChild = nullptr;
}

void StochPresolverBase::countRowsCols()// method is const but changes pointers
{
   if (verbosity <= 1)
      return;

   std::vector<int> count(17, 0);

   int zero_dummy = 0;

   int& n_rows_eq = count[0];
   int& n_rows_ineq = count[1];
   int& n_rows_empty_eq = count[2];
   int& n_rows_empty_ineq = count[3];

   int n_rows_fixed_eq = 0;
   int& n_rows_fixed_ineq = count[4];
   int& n_rows_boxed_ineq = count[5];
   int& n_rows_onsided_ineq = count[6];

   int& n_rows_singleton_eq = count[7];
   int& n_rows_singleton_ineq = count[8];

   int& n_cols = count[9];
   int& n_cols_empty = count[10];
   int& n_cols_free = count[11];
   int& n_cols_onesided = count[12];
   int& n_cols_boxed = count[13];
   int& n_cols_singleton = count[14];
   int& n_cols_orig_free = count[15];
   int& n_cols_orig_free_removed = count[16];

   /* root nodes of equality and inequality system - linking varaiables and linking constraints */
   if (my_rank == 0) {

      int n_rows_linking_eq = 0;
      int n_rows_linking_ineq = 0;
      int n_rows_empty_linking_eq = 0;
      int n_rows_empty_linking_ineq = 0;
      int n_rows_singleton_linking_eq = 0;
      int n_rows_singleton_linking_ineq = 0;

      updatePointersForCurrentNode(-1, EQUALITY_SYSTEM);

      countRowsBlock(n_rows_eq, n_rows_empty_eq, zero_dummy, zero_dummy, n_rows_fixed_eq, n_rows_singleton_eq, EQUALITY_SYSTEM, B_MAT);
      assert(n_rows_eq - n_rows_empty_eq == n_rows_fixed_eq);

      countRowsBlock(n_rows_linking_eq, n_rows_empty_linking_eq, zero_dummy, zero_dummy, n_rows_fixed_eq, n_rows_singleton_linking_eq,
            EQUALITY_SYSTEM, BL_MAT);

      assert(zero_dummy == 0);

      updatePointersForCurrentNode(-1, INEQUALITY_SYSTEM);
      countRowsBlock(n_rows_ineq, n_rows_empty_ineq, n_rows_onsided_ineq, n_rows_boxed_ineq, n_rows_fixed_ineq, n_rows_singleton_ineq,
            INEQUALITY_SYSTEM, B_MAT);
      countRowsBlock(n_rows_linking_ineq, n_rows_empty_linking_ineq, n_rows_onsided_ineq, n_rows_boxed_ineq, n_rows_fixed_ineq,
            n_rows_singleton_linking_ineq, INEQUALITY_SYSTEM, BL_MAT);

      const SimpleVector<double>& ixlow_orig = dynamic_cast<const SimpleVector<double>&>(*dynamic_cast<const DistributedVector<double>& >(*origProb.ixlow).first);
      const SimpleVector<double>& ixupp_orig = dynamic_cast<const SimpleVector<double>&>(*dynamic_cast<const DistributedVector<double>& >(*origProb.ixupp).first);

      countBoxedColumns(n_cols, n_cols_empty, n_cols_free, n_cols_onesided, n_cols_boxed, n_cols_singleton, n_cols_orig_free,
            n_cols_orig_free_removed, ixlow_orig, ixupp_orig, true);

      assert(n_cols - n_cols_empty == n_cols_free + n_cols_onesided + n_cols_boxed);

      std::cout << "#linking_vars:\t\t" << n_cols << "\t(#empty: n_cols_empty " << n_cols_empty << ", #singleton: " << n_cols_singleton << ", #free: "
                << n_cols_free << ", #onesided: " << n_cols_onesided << ", #boxed: " << n_cols_boxed << " #orig_free_non_empty: "
                << n_cols_orig_free - n_cols_orig_free_removed << ")" << "\n";

      std::cout << "#rows B0:\t\t" << n_rows_eq << "\t(#empty: n_rows_empty " << n_rows_empty_eq << ", #singleton: " << n_rows_singleton_eq << ")"
                << "\n";
      std::cout << "#rows Bl_0:\t\t" << n_rows_linking_eq << "\t(#empty: n_rows_empty " << n_rows_empty_linking_eq << ", #singleton: "
                << n_rows_singleton_linking_eq << ")" << "\n";
      std::cout << "#rows D0:\t\t" << n_rows_ineq << "\t(#empty: n_rows_empty " << n_rows_empty_ineq << ", #singleton: " << n_rows_singleton_ineq
                << ")" << "\n";
      std::cout << "#rows Dl_0:\t\t" << n_rows_linking_ineq << "\t(#empty: n_rows_empty " << n_rows_empty_linking_ineq << ", #singleton: "
                << n_rows_singleton_linking_ineq << ")" << "\n";

      n_rows_eq += n_rows_linking_eq;
      n_rows_empty_eq += n_rows_empty_linking_eq;
      n_rows_singleton_eq += n_rows_singleton_linking_eq;

      n_rows_ineq += n_rows_linking_ineq;
      n_rows_empty_ineq += n_rows_empty_linking_ineq;
      n_rows_singleton_ineq += n_rows_singleton_linking_ineq;

      assert(n_rows_eq - n_rows_empty_eq == n_rows_fixed_eq);
      assert(n_rows_ineq - n_rows_empty_ineq == n_rows_onsided_ineq + n_rows_boxed_ineq + n_rows_fixed_ineq);
   }

   /* child nodes in both systems */
   for (int node = 0; node < nChildren; node++) {
      if (presolve_data.nodeIsDummy(node))
         continue;

      /* equality system */
      updatePointersForCurrentNode(node, EQUALITY_SYSTEM);
      countRowsBlock(n_rows_eq, n_rows_empty_eq, zero_dummy, zero_dummy, n_rows_fixed_eq, n_rows_singleton_eq, EQUALITY_SYSTEM, A_MAT);
      assert(zero_dummy == 0);
      assert(n_rows_eq - n_rows_empty_eq == n_rows_fixed_eq);

      /* inequality system */
      updatePointersForCurrentNode(node, INEQUALITY_SYSTEM);
      countRowsBlock(n_rows_ineq, n_rows_empty_ineq, n_rows_onsided_ineq, n_rows_boxed_ineq, n_rows_fixed_ineq, n_rows_singleton_ineq,
            INEQUALITY_SYSTEM, A_MAT);
      assert(n_rows_ineq - n_rows_empty_ineq == n_rows_onsided_ineq + n_rows_boxed_ineq + n_rows_fixed_ineq);

      const SimpleVector<double>& ixlow_orig = dynamic_cast<const SimpleVector<double>&>(*dynamic_cast<const DistributedVector<double>& >(*origProb.ixlow).children[node]->first);
      const SimpleVector<double>& ixupp_orig = dynamic_cast<const SimpleVector<double>&>(*dynamic_cast<const DistributedVector<double>& >(*origProb.ixupp).children[node]->first);

      countBoxedColumns(n_cols, n_cols_empty, n_cols_free, n_cols_onesided, n_cols_boxed, n_cols_singleton, n_cols_orig_free,
            n_cols_orig_free_removed, ixlow_orig, ixupp_orig, false);
      assert(n_cols - n_cols_empty == n_cols_free + n_cols_onesided + n_cols_boxed);
   }

#if 0//TIMING // TODO
   // count how many linking rows do not really link two blocks:
#endif

   /* sync data */
   if (distributed)
      PIPS_MPIsumArrayInPlace(count, MPI_COMM_WORLD);

   n_rows_fixed_eq = n_rows_eq - n_rows_empty_eq;

   if (my_rank == 0) {
      std::cout << "#rows_total:\t\t" << n_rows_eq + n_rows_ineq << "\t(#empty: " << n_rows_empty_eq + n_rows_empty_ineq << ", #fixed: "
                << n_rows_fixed_eq + n_rows_fixed_ineq << ", #boxed: " << n_rows_boxed_ineq << ", #onesided: " << n_rows_onsided_ineq
                << ", #singleton: " << n_rows_singleton_eq + n_rows_singleton_ineq << ")" << "\n";
      std::cout << "#rows non-empty A:\t" << n_rows_eq - n_rows_empty_eq << "\t(#singleton: " << n_rows_singleton_eq << ")" << "\n";
      std::cout << "#rows non-empty C:\t" << n_rows_ineq - n_rows_empty_ineq << "\t(#singleton: " << n_rows_singleton_ineq << ")" << "\n";

      std::cout << "#vars_total:\t\t" << n_cols << "\t(#empty: " << n_cols_empty << ", #free: " << n_cols_free << ", #onesided: " << n_cols_onesided
                << ", #boxed: " << n_cols_boxed << ", #singleton: " << n_cols_singleton << ", #orig_free_non_empty: "
                << n_cols_orig_free - n_cols_orig_free_removed << ")" << "\n";
      std::cout << "#vars non_empty:\t" << n_cols - n_cols_empty << "\n";
   }
}

void StochPresolverBase::countRowsBlock(int& n_rows_total, int& n_rows_empty, int& n_rows_onesided, int& n_rows_boxed, int& n_rows_fixed,
      int& n_rows_singleton, SystemType system_type, BlockType block_type) const {
   if (block_type == BL_MAT)
      if (!presolve_data.hasLinking(system_type))
         return;

   const SimpleVector<int>* nnz_row = (block_type != BL_MAT) ? currNnzRow : currNnzRowLink;
   const SimpleVector<double>* iclow = (block_type != BL_MAT) ? currIclow : currIclowLink;
   const SimpleVector<double>* lhs = (block_type != BL_MAT) ? currIneqLhs : currIneqLhsLink;
   const SimpleVector<double>* icupp = (block_type != BL_MAT) ? currIcupp : currIcuppLink;
   const SimpleVector<double>* rhs = (block_type != BL_MAT) ? currIneqRhs : currIneqRhsLink;
   if (system_type == EQUALITY_SYSTEM)
      rhs = (block_type != BL_MAT) ? currEqRhs : currEqRhsLink;

#ifndef NDEBUG
   if (system_type == EQUALITY_SYSTEM) {
      assert(rhs);
      assert(lhs == nullptr);
      assert(iclow == nullptr);
      assert(icupp == nullptr);
      assert(nnz_row->length() == rhs->length());
   }
   else {
      assert(lhs);
      assert(rhs);
      assert(iclow);
      assert(icupp);
      assert(lhs->length() == rhs->length());
      assert(iclow->length() == icupp->length());
      assert(lhs->length() == nnz_row->length());
      assert(iclow->length() == nnz_row->length());
   }
#endif

   n_rows_total += rhs->length();

   for (int i = 0; i < rhs->length(); ++i) {
      if ((*nnz_row)[i] != 0) {
         if ((*nnz_row)[i] == 1)
            ++n_rows_singleton;

         if (system_type == EQUALITY_SYSTEM) {
            /* row with rhs = lhs */
            ++n_rows_fixed;
         }
         else {
            if (!PIPSisZero((*iclow)[i]) && !PIPSisZero((*icupp)[i])) {
               if (PIPSisEQ((*lhs)[i], (*rhs)[i]))
                  ++n_rows_fixed;
               else
                  ++n_rows_boxed;
            }
            else {
               assert(!PIPSisZero((*iclow)[i]) || !PIPSisZero((*icupp)[i]));
               ++n_rows_onesided;
            }
         }
      }
      else {
         ++n_rows_empty;
      }
   }
}

void StochPresolverBase::countBoxedColumns(int& n_cols_total, int& n_cols_empty, int& n_cols_free, int& n_cols_onesided, int& n_cols_boxed,
      int& n_cols_singleton, int& n_cols_orig_free, int& n_cols_orig_free_removed, const SimpleVector<double>& ixlow_orig,
      const SimpleVector<double>& ixupp_orig, bool at_root_node) const {
   const SimpleVector<double>& ixlow = (at_root_node) ? *currIxlowParent : *currIxlowChild;
   const SimpleVector<double>& ixupp = (at_root_node) ? *currIxuppParent : *currIxuppChild;
   const SimpleVector<int>& curr_nnz = (at_root_node) ? *currNnzColParent : *currNnzColChild;

#ifndef NDEBUG
   const SimpleVector<double>& xupp = (at_root_node) ? *currxuppParent : *currxuppChild;
   const SimpleVector<double>& xlow = (at_root_node) ? *currxlowParent : *currxlowChild;
   assert(ixlow.length() == ixupp.length());
   assert(ixlow.length() == xlow.length());
   assert(xlow.length() == xupp.length());
   assert(ixlow.length() == curr_nnz.length());
#endif

   n_cols_total += ixlow.length();

   for (int i = 0; i < ixlow.length(); i++) {
      if (PIPSisZero(ixupp_orig[i]) && PIPSisZero(ixlow_orig[i]))
         ++n_cols_orig_free;

      if (curr_nnz[i] != 0) {
         if (curr_nnz[i] == 1) {
            ++n_cols_singleton;
         }

         if (!PIPSisZero(ixlow[i]) && !PIPSisZero(ixupp[i]))
            ++n_cols_boxed;
         else if (PIPSisZero(ixlow[i]) && PIPSisZero(ixupp[i]))
            ++n_cols_free;
         else
            ++n_cols_onesided;
      }
      else {
         if (PIPSisZero(ixupp_orig[i]) && PIPSisZero(ixlow_orig[i]))
            ++n_cols_orig_free_removed;
         ++n_cols_empty;
      }
   }
}

/**
 * set all pointers to the currently necessary data
 * If node == -1 we are in the root node
 */
void StochPresolverBase::updatePointersForCurrentNode(int node, SystemType system_type) {
   assert(!presolve_data.nodeIsDummy(node));
   assert(-1 <= node && node < nChildren);
   assert(system_type == EQUALITY_SYSTEM || system_type == INEQUALITY_SYSTEM);

   const SmartPointer<GenMatrix> matrix = (system_type == EQUALITY_SYSTEM) ? presolve_data.getPresProb().A : presolve_data.getPresProb().C;

   /* set matrix pointers for A B and Bl */
   setPointersMatrices(matrix, node);

   /* set lhs rhs for equations */
   setPointersMatrixBoundsActivities(system_type, node);

   /* set x lower upper bounds */
   setPointersVarBounds(node);

   /* set objective function pointers */
   setPointersObjective(node);

   /* set reduction pointers columns and rows */
   setReductionPointers(system_type, node);
}

// todo : set pointers nullptr if no linking constraints?
void StochPresolverBase::setPointersMatrices(const SmartPointer<GenMatrix> mat, int node) {
   assert(-1 <= node && node < nChildren);
   const StochGenMatrix& smat = dynamic_cast<const StochGenMatrix&>(*mat);

   /* in root node only B0 and Bl0 are present */
   if (node == -1) {
      currAmat = nullptr;
      currAmatTrans = nullptr;

      currBmat = dynamic_cast<const SparseGenMatrix*>(smat.Bmat)->getStorageDynamic();
      currBmatTrans = dynamic_cast<const SparseGenMatrix*>(smat.Bmat)->getStorageDynamicTransposed();

      currBlmat = dynamic_cast<const SparseGenMatrix*>(smat.Blmat)->getStorageDynamic();
      currBlmatTrans = dynamic_cast<const SparseGenMatrix*>(smat.Blmat)->getStorageDynamicTransposed();
   }
   else {
      currAmat = dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Amat)->getStorageDynamic();
      currAmatTrans = dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Amat)->getStorageDynamicTransposed();
      currBmat = dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Bmat)->getStorageDynamic();
      currBmatTrans = dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Bmat)->getStorageDynamicTransposed();
      currBlmat = dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Blmat)->getStorageDynamic();
      currBlmatTrans = dynamic_cast<const SparseGenMatrix*>(smat.children[node]->Blmat)->getStorageDynamicTransposed();
   }
}

// todo make one for bounds one for activities
void StochPresolverBase::setPointersMatrixBoundsActivities(SystemType system_type, int node) {
   assert(-1 <= node && node < nChildren);

   if (system_type == EQUALITY_SYSTEM) {
      currEqRhs = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().bA, node, false);

      currIneqLhs = currIclow = currIneqRhs = currIcupp = currIneqLhsLink = currIclowLink = currIneqRhsLink = currIcuppLink = nullptr;

      if (presolve_data.hasLinking(system_type))
         currEqRhsLink = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().bA, node, true);
      else
         currEqRhsLink = nullptr;

   }
   else {
      currIneqLhs = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().bl, node, false);
      currIclow = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().iclow, node, false);
      currIneqRhs = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().bu, node, false);
      currIcupp = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().icupp, node, false);

      currIneqLhsLink = currIclowLink = currIneqRhsLink = currIcuppLink = nullptr;

      if (presolve_data.hasLinking(system_type)) {
         currIneqLhsLink = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().bl, node, true);
         currIclowLink = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().iclow, node, true);
         currIneqRhsLink = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().bu, node, true);
         currIcuppLink = &getSimpleVecFromRowStochVec(*presolve_data.getPresProb().icupp, node, true);
      }
      else
         currEqRhs = currEqRhsLink = nullptr;
   }
}

void StochPresolverBase::setPointersVarBounds(int node) {
   assert(-1 <= node && node < nChildren);

   currxlowParent = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().blx, -1);
   currIxlowParent = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().ixlow, -1);
   currxuppParent = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().bux, -1);
   currIxuppParent = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().ixupp, -1);

   if (node != -1) {
      currxlowChild = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().blx, node);
      currIxlowChild = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().ixlow, node);
      currxuppChild = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().bux, node);
      currIxuppChild = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().ixupp, node);
   }
   else
      currxlowChild = currxuppChild = currIxlowChild = currIxuppChild = nullptr;
}

void StochPresolverBase::setPointersObjective(int node) {
   currgParent = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().g, -1);
   if (node != -1)
      currgChild = &getSimpleVecFromColStochVec(*presolve_data.getPresProb().g, node);
   else
      currgChild = nullptr;
}

void StochPresolverBase::setReductionPointers(SystemType system_type, int node) {
   assert(-1 <= node && node < nChildren);

   const DistributedVector<int>& row_nnz = (system_type == EQUALITY_SYSTEM) ? presolve_data.getNnzsRowA() : presolve_data.getNnzsRowC();

   /* rows */
   currNnzRow = &getSimpleVecFromRowStochVec(row_nnz, node, false);

   if (presolve_data.hasLinking(system_type))
      currNnzRowLink = &getSimpleVecFromRowStochVec(row_nnz, node, true);
   else
      currNnzRowLink = nullptr;

   /* columns */
   currNnzColParent = &getSimpleVecFromColStochVec(presolve_data.getNnzsCol(), -1);

   if (node != -1)
      currNnzColChild = &getSimpleVecFromColStochVec(presolve_data.getNnzsCol(), node);
   else
      currNnzColChild = nullptr;
}
