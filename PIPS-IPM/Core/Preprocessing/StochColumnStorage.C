/*
 * StochColumnStorage.C
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 *
 *
 *
 *  stored_cols_eq/ineq:
 *      | 0   |
 *      | Bl1 | B1 |
 *      | Bl2 |    | B2 |
 *         .            | B3 |
 *         .                   ...
 *         .
 *      | BlN |                    | BN |
 *      | Bl0 | A1 | A2 | A3 | ... | AN |
 *
 *  b0_block_linking_cols_eq/ineq
 *      | B0  |
 *
 *
 */
// todo improve description

#include <memory>

#include "StochColumnStorage.h"
#include "DoubleMatrixTypes.h"
#include "SystemType.h"
#include "DistributedMatrixUtilities.h"
#include "DistributedVectorUtilities.h"

StochColumnStorage::StochColumnStorage(const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part) : nChildren(
      matrix_eq_part.children.size()) {
   assert(matrix_ineq_part.children.size() == nChildren);

   createStorageMatrix(EQUALITY_SYSTEM, matrix_eq_part);
   createStorageMatrix(INEQUALITY_SYSTEM, matrix_ineq_part);
}

void StochColumnStorage::createStorageMatrix(SystemType system_type, const DistributedMatrix& sys_matrix) {
   std::unique_ptr<SparseMatrix>& b0_block_storage = (system_type == EQUALITY_SYSTEM) ? B0_eq : B0_ineq;
   std::unique_ptr<DistributedMatrix>& col_storage = (system_type == EQUALITY_SYSTEM) ? stored_cols_eq : stored_cols_ineq;

   /* extra storage for b0mat entries in linking variable column we want to store */
   assert(sys_matrix.Bmat);
   b0_block_storage = dynamic_cast<const SparseMatrix&>(*sys_matrix.Bmat).cloneEmptyRowsTransposed(true);

   // todo : n, m are wrong here but the counters are broken anyways?


   /* the rest of the matrix is going to get transposed - note that this would normally reverse the dummy non-dummy structure of the matrix
    * we will not reverse it here but rather permute the child array such that the dummy non dummy structure stays the same as before
    * this will make the implementation of multiplication etc easier
    * e.g. A1 will not be swapped with Bln but Bl1
    * Bi will stay at their position in the matrix and just get transposed
    */

   /* clone submatrices */
   /* Amat is empty in root node */
   auto Amat_clone = dynamic_cast<const SparseMatrix&>(*sys_matrix.Amat).cloneEmptyRows(true);
   /* B0mat will not be used and stay empty */
   auto Bmat_clone = dynamic_cast<const SparseMatrix&>(*sys_matrix.Blmat).cloneEmptyRowsTransposed(true);
   auto Blmat_clone = dynamic_cast<const SparseMatrix&>(*sys_matrix.Blmat).cloneEmptyRowsTransposed(true);

   col_storage = std::make_unique<DistributedMatrix>(std::move(Amat_clone), std::move(Bmat_clone), std::move(Blmat_clone), sys_matrix.mpiComm);

   for (const auto & it : sys_matrix.children) {
      const DistributedMatrix& child = *it;

      /* clone submatrices */
      auto Amat_child_clone = dynamic_cast<const SparseMatrix&>(*child.Blmat).cloneEmptyRowsTransposed(true);
      auto Bmat_child_clone = dynamic_cast<const SparseMatrix&>(*child.Bmat).cloneEmptyRowsTransposed(true);
      auto Blmat_child_clone = dynamic_cast<const SparseMatrix&>(*child.Amat).cloneEmptyRowsTransposed(true);

      /* create child */
      std::shared_ptr<DistributedMatrix> child_clone;
      if (child.is_a(kStochGenDummyMatrix))
         child_clone = std::make_unique<StochGenDummyMatrix>();
      else
         child_clone = std::make_unique<DistributedMatrix>(std::move(Amat_child_clone), std::move(Bmat_child_clone), std::move(Blmat_child_clone), child.mpiComm);

      col_storage->children.push_back(child_clone);
   }
}

int StochColumnStorage::storeCol(const INDEX& col, const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part) {
   assert(col.isCol());
   assert(matrix_eq_part.children.size() == matrix_ineq_part.children.size());

   if (col.getNode() == -1)
      return storeLinkingCol(col.getIndex(), matrix_eq_part, matrix_ineq_part);
   else
      return storeLocalCol(col, matrix_eq_part, matrix_ineq_part);
}

int StochColumnStorage::storeLinkingCol(int col, const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part) {
   assert(matrix_eq_part.children.size() == matrix_ineq_part.children.size());
   assert(matrix_eq_part.children.size() == stored_cols_eq->children.size());

   assert(!matrix_eq_part.is_a(kStochGenDummyMatrix));
   assert(!matrix_eq_part.is_a(kStochGenDummyMatrix));

   const SparseMatrix& B0 = *getSparseGenMatrixFromStochMat(matrix_eq_part, -1, B_MAT);
   const SparseMatrix& D0 = *getSparseGenMatrixFromStochMat(matrix_ineq_part, -1, B_MAT);
   const SparseMatrix& Bl0 = *getSparseGenMatrixFromStochMat(matrix_eq_part, -1, BL_MAT);
   const SparseMatrix& Dl0 = *getSparseGenMatrixFromStochMat(matrix_ineq_part, -1, BL_MAT);

   SparseMatrix& Bl0_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, -1, BL_MAT);
   SparseMatrix& Dl0_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, -1, BL_MAT);

   const int B0_index = B0_eq->appendCol(B0, col);
#ifndef NDEBUG
   const int D0_index = B0_ineq->appendCol(D0, col);
   const int Bl0_index = Bl0_storage.appendCol(Bl0, col);
   const int Dl0_index = Dl0_storage.appendCol(Dl0, col);
   assert(B0_index == D0_index);
   assert(Bl0_index == Dl0_index);
   assert(B0_index == Bl0_index);
#else
   B0_ineq->appendCol(D0, col);
   Bl0_storage.appendCol(Bl0, col);
   Dl0_storage.appendCol(Dl0, col);
#endif

   for (unsigned int i = 0; i < nChildren; ++i) {
      if (matrix_eq_part.children[i]->is_a(kStochGenDummyMatrix)) {
         assert(matrix_ineq_part.children[i]->is_a(kStochGenDummyMatrix));
         continue;
      }
      else {
         assert(!matrix_ineq_part.children[i]->is_a(kStochGenDummyMatrix));
      }

      const SparseMatrix& Ai = *getSparseGenMatrixFromStochMat(matrix_eq_part, i, A_MAT);
      const SparseMatrix& Ci = *getSparseGenMatrixFromStochMat(matrix_ineq_part, i, A_MAT);

      SparseMatrix& Ai_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, i, BL_MAT);
      SparseMatrix& Ci_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, i, BL_MAT);

#ifndef NDEBUG
      const int Ai_index = Ai_storage.appendCol(Ai, col);
      const int Ci_index = Ci_storage.appendCol(Ci, col);
      assert(Ai_index == Ci_index);
      assert(Ai_index == B0_index);
#else
      Ai_storage.appendCol(Ai, col);
      Ci_storage.appendCol(Ci, col);
#endif
   }

   return B0_index;
}

int StochColumnStorage::storeLocalCol(const INDEX& col, const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part) {
   const int node = col.getNode();
   const int col_index = col.getIndex();

   assert(0 <= node && node < static_cast<int>(matrix_eq_part.children.size()));
   assert(matrix_eq_part.children.size() == matrix_ineq_part.children.size());
   assert(matrix_eq_part.children.size() == stored_cols_eq->children.size());

   assert(!matrix_eq_part.children[node]->is_a(kStochGenDummyMatrix));
   assert(!matrix_ineq_part.children[node]->is_a(kStochGenDummyMatrix));

   const SparseMatrix& Bi = *getSparseGenMatrixFromStochMat(matrix_eq_part, node, B_MAT);
   const SparseMatrix& Di = *getSparseGenMatrixFromStochMat(matrix_ineq_part, node, B_MAT);
   const SparseMatrix& Bli = *getSparseGenMatrixFromStochMat(matrix_eq_part, node, BL_MAT);
   const SparseMatrix& Dli = *getSparseGenMatrixFromStochMat(matrix_ineq_part, node, BL_MAT);

   SparseMatrix& Bi_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, B_MAT);
   SparseMatrix& Di_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, B_MAT);
   SparseMatrix& Bli_storage = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, A_MAT);
   SparseMatrix& Dli_storage = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, A_MAT);

   const int Bi_index = Bi_storage.appendCol(Bi, col_index);
#ifndef NDEBUG
   const int Di_index = Di_storage.appendCol(Di, col_index);
   const int Bli_index = Bli_storage.appendCol(Bli, col_index);
   const int Dli_index = Dli_storage.appendCol(Dli, col_index);

   assert(Bi_index == Di_index);
   assert(Bli_index == Dli_index);
   assert(Bi_index == Bli_index);
#else
   Di_storage.appendCol(Di, col_index);
   Bli_storage.appendCol(Bli, col_index);
   Dli_storage.appendCol(Dli, col_index);
#endif

   return Bi_index;
}

void StochColumnStorage::axpyAtCol(double beta, DistributedVector<double>* eq_vec, DistributedVector<double>* ineq_vec, SimpleVector<double>* eq_link,
      SimpleVector<double>* ineq_link, double alpha, const INDEX& col) const {
   assert(col.isCol());

   if (!PIPSisEQ(1.0, beta)) {
      if (eq_vec)
         eq_vec->scale(beta);
      if (ineq_vec)
         ineq_vec->scale(beta);
   }

   if (eq_vec)
      axpyAtCol(*eq_vec, eq_link, alpha, col, EQUALITY_SYSTEM);

   if (ineq_vec)
      axpyAtCol(*ineq_vec, ineq_link, alpha, col, INEQUALITY_SYSTEM);
}

void StochColumnStorage::axpyAtCol(DistributedVector<double>& vec, SimpleVector<double>* vec_link, double alpha, const INDEX& col,
      SystemType system_type) const {
   assert(col.isCol());
   assert(!vec.isKindOf(kStochDummy));

   const SparseMatrix& B0_mat = (system_type == EQUALITY_SYSTEM) ? *B0_eq : *B0_ineq;
   const DistributedMatrix& mat = (system_type == EQUALITY_SYSTEM) ? *stored_cols_eq : *stored_cols_ineq;

   if (col.isLinkingCol()) {
      /* B0 */
      SimpleVector<double>& b0_vec = getSimpleVecFromRowStochVec(vec, -1, false);
      B0_mat.axpyWithRowAt(alpha, b0_vec, col.getIndex());

      /* Bl0 */
      const SparseMatrix& Bl0_mat = *getSparseGenMatrixFromStochMat(mat, -1, BL_MAT);
      SimpleVector<double>& bl0_vec = (vec_link == nullptr) ? getSimpleVecFromRowStochVec(vec, -1, true) : *vec_link;
      Bl0_mat.axpyWithRowAt(alpha, bl0_vec, col.getIndex());

      /* Amats */
      for (unsigned int node = 0; node < nChildren; ++node) {
         if (!vec.children[node]->isKindOf(kStochDummy)) {
            assert(!mat.children[node]->is_a(kStochGenDummyMatrix));

            const SparseMatrix& A_mat = *getSparseGenMatrixFromStochMat(mat, node, BL_MAT);
            SimpleVector<double>& a_vec = getSimpleVecFromRowStochVec(vec, node, false);

            A_mat.axpyWithRowAt(alpha, a_vec, col.getIndex());
         }
      }
   }
   else {
      assert(!mat.children[col.getNode()]->is_a(kStochGenDummyMatrix));

      const SparseMatrix& Bi_mat = *getSparseGenMatrixFromStochMat(mat, col.getNode(), B_MAT);
      SimpleVector<double>& bi_vec = getSimpleVecFromRowStochVec(vec, col.getNode(), false);

      const SparseMatrix& Bli_mat = *getSparseGenMatrixFromStochMat(mat, col.getNode(), A_MAT);
      SimpleVector<double>& bli_vec = (vec_link == nullptr) ? getSimpleVecFromRowStochVec(vec, -1, true) : *vec_link;

      Bi_mat.axpyWithRowAt(alpha, bi_vec, col.getIndex());
      Bli_mat.axpyWithRowAt(alpha, bli_vec, col.getIndex());
   }
}


double
StochColumnStorage::multColTimesVec(const INDEX& col, const DistributedVector<double>& vec_eq, const DistributedVector<double>& vec_ineq) const {
   assert(col.isCol());
   assert(col.hasValidNode(nChildren));
   assert(nChildren == vec_eq.children.size());
   assert(nChildren == vec_ineq.children.size());

   double res{0.0};

   if (col.isLinkingCol()) {
      assert(PIPS_MPIisValueEqual(col.getNode()));
      /* we need to synchronize the column times duals in this case */
      if (PIPS_MPIgetRank() == 0)
         res = multiplyLinkingColTimesVec(col.getIndex(), vec_eq, vec_ineq);
      else
         res = multiplyLinkingColTimesVecWithoutRootNode(col.getIndex(), vec_eq, vec_ineq);

      PIPS_MPIgetSumInPlace(res);
   }
   else
      res = multiplyLocalColTimesVec(col, vec_eq, vec_ineq);

   return res;
}

double
StochColumnStorage::multiplyLinkingColTimesVec(int col, const DistributedVector<double>& vec_eq, const DistributedVector<double>& vec_ineq) const {
   assert(!vec_eq.isKindOf(kStochDummy) && !vec_ineq.isKindOf(kStochDummy));

   double res = 0.0;
   /* equality system */
   const SparseMatrix& Bl0_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, -1, BL_MAT);

   res += B0_eq->localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, -1, false), col);
   res += Bl0_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, -1, true), col);

   /* inequality system */
   const SparseMatrix& Dl0_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, -1, BL_MAT);

   res += B0_ineq->localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, -1, false), col);
   res += Dl0_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, -1, true), col);

   res += multiplyLinkingColTimesVecWithoutRootNode(col, vec_eq, vec_ineq);

   return res;
}

double StochColumnStorage::multiplyLinkingColTimesVecWithoutRootNode(int col, const DistributedVector<double>& vec_eq,
      const DistributedVector<double>& vec_ineq) const {
   double res = 0.0;

   /* Amat equality and inequality */
   for (unsigned int node = 0; node < nChildren; ++node) {
      if (!vec_eq.children[node]->isKindOf(kStochDummy)) {
         assert(!vec_ineq.children[node]->isKindOf(kStochDummy));
         assert(!stored_cols_eq->children[node]->is_a(kStochGenDummyMatrix));
         assert(!stored_cols_ineq->children[node]->is_a(kStochGenDummyMatrix));

         const SparseMatrix& A_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, BL_MAT);
         assert(vec_eq.children[node]->first);

         const SimpleVector<double>& a_vec = getSimpleVecFromRowStochVec(vec_eq, node, false);
         res += A_mat.localRowTimesVec(a_vec, col);

         const SparseMatrix& C_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, BL_MAT);
         assert(vec_ineq.children[node]->first);

         const SimpleVector<double>& c_vec = getSimpleVecFromRowStochVec(vec_ineq, node, false);
         res += C_mat.localRowTimesVec(c_vec, col);
      }
   }

   return res;
}

double StochColumnStorage::multiplyLocalColTimesVec(const INDEX& col, const DistributedVector<double>& vec_eq,
      const DistributedVector<double>& vec_ineq) const {
   const int node = col.getNode();
   const int col_index = col.getIndex();
   double res = 0.0;

   assert(0 <= node && node < static_cast<int>(nChildren));

   assert(!vec_eq.isKindOf(kStochDummy) && !vec_ineq.isKindOf(kStochDummy));
   assert(!stored_cols_eq->children[node]->is_a(kStochGenDummyMatrix));
   assert(!stored_cols_ineq->children[node]->is_a(kStochGenDummyMatrix));

   /* equality system */
   const SparseMatrix& Bi_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, B_MAT);
   const SparseMatrix& Bli_mat = *getSparseGenMatrixFromStochMat(*stored_cols_eq, node, A_MAT);

   res += Bi_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, node, false), col_index);
   res += Bli_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_eq, -1, true), col_index);

   /* inequality system */
   const SparseMatrix& Di_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, B_MAT);
   const SparseMatrix& Dli_mat = *getSparseGenMatrixFromStochMat(*stored_cols_ineq, node, A_MAT);

   res += Di_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, node, false), col_index);
   res += Dli_mat.localRowTimesVec(getSimpleVecFromRowStochVec(vec_ineq, -1, true), col_index);

   return res;
}

