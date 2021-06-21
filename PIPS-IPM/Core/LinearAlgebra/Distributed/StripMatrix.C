/*
 * StripMatrix.C
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */

#include "StripMatrix.h"

#include "DistributedVector.h"
#include "SimpleVector.hpp"
#include "DoubleMatrixTypes.h"
#include "DistributedTreeCallbacks.h"

#include "pipsdef.h"
#include <algorithm>

StripMatrix::StripMatrix(bool is_vertical, std::unique_ptr<GeneralMatrix> first_in, std::unique_ptr<GeneralMatrix> last_in,
   MPI_Comm mpi_comm_, bool is_view) : first(std::move(first_in)),
   last(std::move(last_in)), is_vertical(is_vertical), mpi_comm(mpi_comm_), distributed(PIPS_MPIgetDistributed(mpi_comm)),
   rank(PIPS_MPIgetRank(mpi_comm)), is_view{is_view} {
   assert(first);

   m = first->n_rows();
   n = first->n_columns();

   nonzeros += first->numberOfNonZeros();

   if (last) {
      const auto[ml, nl] = last->n_rows_columns();
      nonzeros += last->numberOfNonZeros();

      if (is_vertical) {
         assert(n == nl);
         m += ml;
      } else {
         assert(m == ml);
         n += nl;
      }
   }
}

StripMatrix::~StripMatrix() {
   if(is_view)
   {
      first.release();
      last.release();

      for( auto& child : children )
         child.release();
   }
}

void StripMatrix::addChild(std::unique_ptr<StripMatrix> child_in) {
   assert(child_in);

   StripMatrix* child_ptr = child_in.get();
   children.push_back(std::move(child_in));

   nonzeros += child_ptr->numberOfNonZeros();

   const auto[m_, n_] = child_ptr->n_rows_columns();

   assert(child_ptr->is_vertical == this->is_vertical || child_ptr->is_a(kStringGenDummyMatrix));

   if (!child_ptr->is_a(kStringGenDummyMatrix)) {
      if (is_vertical) {
         assert(n == n_);
         m += m_;
      } else {
         assert(m == m_);
         n += n_;
      }
   }
}

std::pair<long long, long long> StripMatrix::n_rows_columns() const {
   return {m, n};
}

long long StripMatrix::n_rows() const {
   return m;
}

long long StripMatrix::n_columns() const {
   return n;
}

bool StripMatrix::isEmpty() const {
   return !first && !last && children.empty() && m == 0 && n == 0;
}

int StripMatrix::is_a(int type) const {
   return (type == kStripMatrix || type == kStringMatrix || type == kGenMatrix);
}

double StripMatrix::inf_norm() const {
   double norm = 0.0;

   for (const auto& it : children)
      norm = std::max(norm, it->inf_norm());

   if (distributed)
      norm = PIPS_MPIgetMax(norm, mpi_comm);

   norm = std::max(norm, first->inf_norm());

   if (last)
      norm = std::max(norm, last->inf_norm());

   return norm;
}

void StripMatrix::scalarMult(double num) {
   first->scalarMult(num);

   if (last)
      last->scalarMult(num);

   for (auto& it : children)
      it->scalarMult(num);
}

void StripMatrix::mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   if (is_vertical)
      multVertical(beta, y, alpha, x, &AbstractMatrix::mult);
   else
      multHorizontal(beta, y, alpha, x, true, &AbstractMatrix::mult);
}

void StripMatrix::mult_transform(double beta, Vector<double>& y, double alpha, const Vector<double>& x,
   const std::function<double(const double&)>& transform) const {

   auto mult = [&capture0 = std::as_const(transform)](const GeneralMatrix* mat, auto&& PH1, auto&& PH2, auto&& PH3, auto&& PH4) {
      mat->mult_transform(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
         std::forward<decltype(PH3)>(PH3), std::forward<decltype(PH4)>(PH4), capture0);
   };

   if (is_vertical) {
      multVertical(beta, y, alpha, x, mult);
   } else {
      multHorizontal(beta, y, alpha, x, true, mult);
   }
}

void StripMatrix::transpose_mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   if (is_vertical)
      transMultVertical(beta, y, alpha, x, true, &AbstractMatrix::transpose_mult);
   else
      transMultHorizontal(beta, y, alpha, x, &AbstractMatrix::transpose_mult);
}

void StripMatrix::transpose_mult_transform(double beta, Vector<double>& y, double alpha, const Vector<double>& x,
   const std::function<double(const double&)>& transform) const {

   auto transpose_mult = [&capture0 = std::as_const(transform)](const GeneralMatrix* mat, auto&& PH1, auto&& PH2, auto&& PH3, auto&& PH4) {
      mat->transpose_mult_transform(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
         std::forward<decltype(PH3)>(PH3), std::forward<decltype(PH4)>(PH4), capture0);
   };

   if (is_vertical)
      transMultVertical(beta, y, alpha, x, true, transpose_mult);
   else
      transMultHorizontal(beta, y, alpha, x, transpose_mult);
}

int StripMatrix::numberOfNonZeros() const {
#ifndef NDEBUG
   const int nonzeros_now = static_cast<int>(nonzeros);
   const_cast<StripMatrix*>(this)->recomputeNonzeros();
   assert(nonzeros_now == nonzeros);
#endif
   return static_cast<int>(nonzeros);
}

void StripMatrix::recomputeNonzeros() {
   nonzeros = 0;

   for (auto& child : children) {
      child->recomputeNonzeros();
      if (PIPS_MPIgetRank(child->mpi_comm) == 0)
         nonzeros += child->numberOfNonZeros();
#ifndef NDEBUG
         // !in DEBUG the nonzeros get recomputed so all processes have to join in - a terrible hack but for safety..
      else
         static_cast<void>(child->numberOfNonZeros());
#endif
   }

   if (dynamic_cast<const StripMatrix*>(first.get())) {
      auto& matstr = dynamic_cast<StripMatrix&>(*first);
      matstr.recomputeNonzeros();
      if (PIPS_MPIgetRank(matstr.mpi_comm) == 0)
         nonzeros += matstr.numberOfNonZeros();
#ifndef NDEBUG
         // !in DEBUG the nonzeros get recomputed so all processes have to join in - a terrible hack but for safety..
      else
         static_cast<void>(matstr.numberOfNonZeros());
#endif
   } else if (PIPS_MPIiAmSpecial(distributed, mpi_comm)) {
      if (first)
         nonzeros += first->numberOfNonZeros();
      if (last)
         nonzeros += last->numberOfNonZeros();
   }

   PIPS_MPIgetSumInPlace(nonzeros, mpi_comm);
}

void StripMatrix::multVertical(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   const std::function<void(const GeneralMatrix*, double, Vector<double>&, double, const Vector<double>&)>& mult) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   assert(is_vertical);
   assert(y.children.size() == children.size());

   mult(first.get(), beta, *y.first, alpha, x);

   if (last) {
      assert(y.last);
      mult(last.get(), beta, *y.last, alpha, x);
   }

   for (size_t i = 0; i < children.size(); ++i) {
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(y.children[i]->isKindOf(kStochDummy));

      children[i]->multVertical(beta, *y.children[i], alpha, x, mult);
   }
}

void StripMatrix::multHorizontal(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   bool root, const std::function<void(const GeneralMatrix*, double, Vector<double>&, double, const Vector<double>&)>& mult) const {
   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(!is_vertical);
   assert(x.children.size() == children.size());
   assert((x.last && last) || (x.last == nullptr && last == nullptr));

   if (first->is_a(kStripMatrix)) {
      assert(!last);
      assert(children.empty());
      assert(!root);
      dynamic_cast<const StripMatrix&>(*first).multHorizontal(1.0, y, alpha, *x.first, false, mult);
   } else {
      if (PIPS_MPIiAmSpecial(distributed, mpi_comm)) {
         mult(first.get(), root ? beta : 1.0, y, alpha, *x.first);
         if (last)
            mult(last.get(), 1.0, y, alpha, *x.last);
      } else if (root)
         y.setToZero();

      for (size_t i = 0; i < children.size(); ++i)
         children[i]->multHorizontal(1.0, y, alpha, *x.children[i], false, mult);

      if (distributed && root)
         PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);
   }
}

void StripMatrix::transMultVertical(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   bool root, const std::function<void(const GeneralMatrix*, double, Vector<double>&, double, const Vector<double>&)>& transpose_mult) const {
   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(is_vertical);
   assert(x.children.size() == children.size());
   assert(x.first && first);
   assert((x.last && last) || (x.last == nullptr && last == nullptr));

   if (first->is_a(kStripMatrix)) {
      assert(!last);
      assert(children.empty());
      assert(!root);
      dynamic_cast<StripMatrix&>(*first).transMultVertical(1.0, y, alpha, *x.first, false, transpose_mult);
   } else {
      if (rank == 0) {
         transpose_mult(first.get(), root ? beta : 1.0, y, alpha, *x.first);
         if (last)
            transpose_mult(last.get(), 1.0, y, alpha, *x.last);
      } else if (root)
         y.setToZero();

      for (size_t i = 0; i < children.size(); i++)
         children[i]->transMultVertical(1.0, y, alpha, *x.children[i], false, transpose_mult);

      if (distributed && root)
         PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);
   }
}

void
StripMatrix::transMultHorizontal(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   const std::function<void(const GeneralMatrix*, double, Vector<double>&, double, const Vector<double>&)>& transpose_mult) const {

   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   assert(!is_vertical);
   assert(y.children.size() == children.size());
   if (y.last)
      assert(last);
   else
      assert(last == nullptr);

   transpose_mult(first.get(), beta, *y.first, alpha, x);

   if (last)
      transpose_mult(last.get(), beta, *y.last, alpha, x);

   for (size_t i = 0; i < children.size(); ++i) {
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(y.children[i]->isKindOf(kStochDummy));

      children[i]->transMultHorizontal(beta, *y.children[i], alpha, x, transpose_mult);
   }
}

void StripMatrix::getColMinMaxVecHorizontal(bool get_min, bool initialize_vec, const Vector<double>* row_scale,
   Vector<double>& minmax_in) const {
   assert(!is_vertical);
   auto& minmax = dynamic_cast<DistributedVector<double>&>(minmax_in);

   assert(minmax.first && first);
   assert(minmax.children.size() == children.size());
   assert((minmax.last && last) || (minmax.last == nullptr && last == nullptr));


   for (size_t i = 0; i < children.size(); ++i) {
      assert(minmax.children[i]);
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(minmax.children[i]->isKindOf(kStochDummy));

      children[i]->getColMinMaxVecHorizontal(get_min, initialize_vec, row_scale, *minmax.children[i]);
   }

   first->getColMinMaxVec(get_min, initialize_vec, row_scale, *minmax.first);

   if (last)
      last->getColMinMaxVec(get_min, initialize_vec, row_scale, *minmax.last);
}

void StripMatrix::getColMinMaxVecVertical(bool get_min, bool initialize_vec, const Vector<double>* row_scale_in,
   Vector<double>& minmax_) const {
   assert(is_vertical);
   const bool has_rowscale = (row_scale_in != nullptr);

   const auto* row_scale = dynamic_cast<const DistributedVector<double>*>(row_scale_in);
   auto& minmax = dynamic_cast<SimpleVector<double>&>(minmax_);

   assert(!has_rowscale || row_scale->children.size() == children.size());
   if (has_rowscale)
      assert((row_scale->last && last) || (row_scale->last == nullptr && last == nullptr));


   for (size_t i = 0; i < children.size(); i++) {
      if (has_rowscale)
         if (children[i]->is_a(kStringGenDummyMatrix))
            assert(row_scale->children[i]->isKindOf(kStochDummy));

      children[i]->getColMinMaxVecVertical(get_min, false, has_rowscale ? row_scale->children[i].get() : nullptr, minmax);
   }

   if (distributed) {
      if (get_min)
         PIPS_MPIminArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
      else
         PIPS_MPImaxArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
   }

   first->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->first.get() : nullptr, minmax);

   if (last)
      last->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->last.get() : nullptr, minmax);
}

/** DistributedVector<double> colScaleVec, SimpleVector<double> minmaxVec */
void
StripMatrix::getRowMinMaxVecHorizontal(bool get_min, bool initialize_vec, const Vector<double>* col_scale_in,
   Vector<double>& minmax_) const {
   assert(!is_vertical);
   const bool has_colscale = (col_scale_in != nullptr);

   const auto* col_scale = dynamic_cast<const DistributedVector<double>*>(col_scale_in);
   auto& minmax = dynamic_cast<SimpleVector<double>&>(minmax_);
   assert(!has_colscale || col_scale->children.size() == children.size());
   if (has_colscale)
      assert((col_scale->last && last) || (col_scale->last == nullptr && last == nullptr));


   for (size_t i = 0; i < children.size(); i++) {
      if (has_colscale)
         if (children[i]->is_a(kStringGenDummyMatrix))
            assert(col_scale->children[i]->isKindOf(kStochDummy));

      children[i]->getRowMinMaxVecHorizontal(get_min, false, has_colscale ? col_scale->children[i].get() : nullptr, minmax);
   }

   if (distributed) {
      if (get_min)
         PIPS_MPIminArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
      else
         PIPS_MPImaxArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
   }

   first->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->first.get() : nullptr, minmax);

   if (last)
      last->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->last.get() : nullptr, minmax);
}

/** DistributedVector<double> minmaxVec, SimpleVector<double> colScaleVec */
void StripMatrix::getRowMinMaxVecVertical(bool get_min, bool initialize_vec, const Vector<double>* col_scale,
   Vector<double>& minmax_in) const {
   assert(is_vertical);

   auto& minmax = dynamic_cast<DistributedVector<double>&>(minmax_in);

   assert(minmax.first && first);
   assert(minmax.children.size() == children.size());
   assert((minmax.last && last) || (minmax.last == nullptr && last == nullptr));

   first->getRowMinMaxVec(get_min, initialize_vec, col_scale, *minmax.first);

   for (size_t i = 0; i < children.size(); ++i) {
      assert(minmax.children[i]);
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(minmax.children[i]->isKindOf(kStochDummy));

      children[i]->getRowMinMaxVecVertical(get_min, initialize_vec, col_scale, *minmax.children[i]);
   }

   if (last)
      last->getRowMinMaxVec(get_min, initialize_vec, col_scale, *minmax.last);
}

void StripMatrix::getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec,
   Vector<double>& minmaxVec) const {
   if (is_vertical)
      getRowMinMaxVecVertical(getMin, initializeVec, colScaleVec, minmaxVec);
   else
      getRowMinMaxVecHorizontal(getMin, initializeVec, colScaleVec, minmaxVec);
}


void StripMatrix::getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec,
   Vector<double>& minmaxVec) const {
   if (is_vertical)
      getColMinMaxVecVertical(getMin, initializeVec, rowScaleVec, minmaxVec);
   else
      getColMinMaxVecHorizontal(getMin, initializeVec, rowScaleVec, minmaxVec);
}

void StripMatrix::columnScaleVertical(const Vector<double>& vec) {
   assert(is_vertical);

   first->columnScale(vec);

   for (auto& i : children)
      i->columnScaleVertical(vec);

   if (last)
      last->columnScale(vec);
}

void StripMatrix::columnScaleHorizontal(const Vector<double>& vec_in) {
   assert(!is_vertical);

   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && last) || (vec.last == nullptr && last == nullptr));

   first->columnScale(*vec.first);

   for (size_t i = 0; i < children.size(); ++i) {
      assert(vec.children[i]);
      if (children[i]->is_a(kStochGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->columnScaleHorizontal(*vec.children[i]);
   }

   if (last)
      last->columnScale(*vec.last);
}

void StripMatrix::rowScaleVertical(const Vector<double>& vec_in) {
   assert(is_vertical);

   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && last) || (vec.last == nullptr && last == nullptr));

   first->rowScale(*vec.first);

   for (size_t i = 0; i < children.size(); i++) {
      assert(vec.children[i]);
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->rowScaleVertical(*vec.children[i]);
   }

   if (last)
      last->rowScale(*vec.last);
}

void StripMatrix::rowScaleHorizontal(const Vector<double>& vec) {
   assert(!is_vertical);
   first->rowScale(vec);

   if (last)
      last->rowScale(vec);

   for (auto& i : children)
      i->rowScaleHorizontal(vec);
}

void StripMatrix::columnScale(const Vector<double>& vec) {
   if (is_vertical)
      columnScaleVertical(vec);
   else
      columnScaleHorizontal(vec);
}

void StripMatrix::rowScale(const Vector<double>& vec) {
   if (is_vertical)
      rowScaleVertical(vec);
   else
      rowScaleHorizontal(vec);
}

void StripMatrix::addRowSumsVertical(Vector<double>& vec_in) const {
   assert(is_vertical);

   auto& vec = dynamic_cast<DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && last) || (vec.last == nullptr && last == nullptr));

   first->addRowSums(*vec.first);

   for (size_t i = 0; i < children.size(); i++) {
      assert(vec.children[i]);
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->addRowSumsVertical(*vec.children[i]);
   }

   if (last)
      last->addRowSums(*vec.last);
}

void StripMatrix::addRowSumsHorizontal(Vector<double>& vec_in) const {
   assert(!is_vertical);
   auto& vec = dynamic_cast<SimpleVector<double>&>(vec_in);

   for (const auto& i : children)
      i->addRowSumsHorizontal(vec);

   if (distributed)
      PIPS_MPIsumArrayInPlace(vec.elements(), vec.length(), mpi_comm);

   first->addRowSums(vec);

   if (last)
      last->addRowSums(vec);
}

void StripMatrix::addColSumsVertical(Vector<double>& vec_in) const {
   assert(is_vertical);
   auto& vec = dynamic_cast<SimpleVector<double>&>(vec_in);

   for (const auto& i : children)
      i->addColSumsVertical(vec);

   if (distributed)
      PIPS_MPIsumArrayInPlace(vec.elements(), vec.length(), mpi_comm);

   first->addColSums(vec);

   if (last)
      last->addColSums(vec);
}

void StripMatrix::addColSumsHorizontal(Vector<double>& vec_in) const {
   assert(!is_vertical);

   auto& vec = dynamic_cast<DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && last) || (vec.last == nullptr && last == nullptr));

   first->addColSums(*vec.first);

   for (size_t i = 0; i < children.size(); i++) {
      assert(vec.children[i]);
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->addColSumsHorizontal(*vec.children[i]);
   }

   if (last)
      last->addColSums(*vec.last);
}

void StripMatrix::sum_transform_rows_vertical(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   assert(is_vertical);

   auto& result = dynamic_cast<DistributedVector<double>&>(result_);

   assert(result.first);
   assert(result.children.size() == children.size());
   assert((result.last && last) || (!result.last && !last));

   first->sum_transform_rows(*result.first, transform);

   for (size_t i = 0; i < children.size(); i++) {
      assert(result.children[i]);
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(result.children[i]->isKindOf(kStochDummy));

      children[i]->sum_transform_rows(*result.children[i], transform);
   }

   if (last)
      last->sum_transform_rows(*result.last, transform);
}

void StripMatrix::sum_transform_rows_horizontal(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   assert(!is_vertical);
   auto& result = dynamic_cast<SimpleVector<double>&>(result_);

   for (const auto& i : children)
      i->sum_transform_rows(result, transform);

   if (distributed)
      PIPS_MPIsumArrayInPlace(result.elements(), result.length(), mpi_comm);

   first->sum_transform_rows(result, transform);

   if (last)
      last->sum_transform_rows(result, transform);
}

void StripMatrix::sum_transform_rows(Vector<double>& result, const std::function<double(const double&)>& transform) const {
   if (is_vertical) {
      sum_transform_rows_vertical(result, transform);
   } else {
      sum_transform_rows_horizontal(result, transform);
   }
}

void StripMatrix::sum_transform_columns(Vector<double>& result, const std::function<double(const double&)>& transform) const {
   if (is_vertical) {
      sum_transform_columns_vertical(result, transform);
   } else {
      sum_transform_columns_horizontal(result, transform);
   }
}

void StripMatrix::sum_transform_columns_vertical(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   assert(is_vertical);
   auto& result = dynamic_cast<SimpleVector<double>&>(result_);

   for (const auto& i : children)
      i->sum_transform_columns_vertical(result, transform);

   if (distributed)
      PIPS_MPIsumArrayInPlace(result.elements(), result.length(), mpi_comm);

   first->sum_transform_columns(result, transform);

   if (last)
      last->sum_transform_columns(result, transform);
}

void StripMatrix::sum_transform_columns_horizontal(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   assert(!is_vertical);

   auto& result = dynamic_cast<DistributedVector<double>&>(result_);

   assert(result.first);
   assert(result.children.size() == children.size());
   assert((result.last && last) || (!result.last && !last));

   first->sum_transform_columns(*result.first, transform);

   for (size_t i = 0; i < children.size(); i++) {
      assert(result.children[i]);
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(result.children[i]->isKindOf(kStochDummy));

      children[i]->sum_transform_columns_horizontal(*result.children[i], transform);
   }

   if (last)
      last->sum_transform_columns(*result.last, transform);
}

void StripMatrix::addRowSums(Vector<double>& vec) const {
   if (is_vertical)
      addRowSumsVertical(vec);
   else
      addRowSumsHorizontal(vec);
}

void StripMatrix::addColSums(Vector<double>& vec) const {
   if (is_vertical)
      addColSumsVertical(vec);
   else
      addColSumsHorizontal(vec);
}

void StripMatrix::combineChildrenInNewChildren(const std::vector<unsigned int>& map_child_subchild,
   const std::vector<MPI_Comm>& child_comms) {

#ifndef NDEBUG
   const unsigned int n_new_children = getNDistinctValues(map_child_subchild);
   assert(child_comms.size() == n_new_children);
   assert(children.size() == map_child_subchild.size());
#endif

   unsigned int n_children{0};
   for (unsigned int i = 0; i < map_child_subchild.size(); ++i) {
      if (child_comms[n_children] == MPI_COMM_NULL) {
         addChild(std::make_unique<StringGenDummyMatrix>());

         while (i + 1 != map_child_subchild.size() && map_child_subchild[i] == map_child_subchild[i + 1]) {
            ++i;
         }
      } else {
         std::unique_ptr<SparseMatrix> empty_filler = is_vertical ? std::make_unique<SparseMatrix>(0, n, 0) : std::make_unique<SparseMatrix>(m, 0, 0);
         std::unique_ptr<StripMatrix> new_child = std::make_unique<StripMatrix>(is_vertical, std::move(empty_filler), nullptr, child_comms[n_children]);

         StripMatrix* new_child_ptr = new_child.get();
         /* will not change size of StringGenMat since new_child is of size zero */
         addChild(std::move(new_child));

         new_child_ptr->addChild(std::move(children[i]));

         while (i + 1 != map_child_subchild.size() && map_child_subchild[i] == map_child_subchild[i + 1]) {
            ++i;
            new_child_ptr->addChild(std::move(children[i]));
         }
      }

      ++n_children;
   }

   assert(n_children == n_new_children);
   assert(children.size() == n_new_children + map_child_subchild.size());

   children.erase(children.begin(), children.begin() + map_child_subchild.size());
   assert(children.size() == n_new_children);

   recomputeNonzeros();
}

std::unique_ptr<GeneralMatrix> StripMatrix::shaveBottom(int n_rows) {
   assert(!is_vertical);
   assert(first);

   std::unique_ptr<GeneralMatrix> mat_border{first->shaveBottom(n_rows)};
   std::unique_ptr<GeneralMatrix> matlink_border{last ? last->shaveBottom(n_rows) : nullptr};

   std::unique_ptr<StripMatrix> border = std::make_unique<StripMatrix>(false, std::move(mat_border), std::move(matlink_border), mpi_comm);

#ifndef NDEBUG
   const auto[mB, nB] = border->n_rows_columns();
   (void) nB;
   assert(mB == n_rows);
#endif

   for (auto& child : children)
      border->addChild(std::unique_ptr<StripMatrix>(dynamic_cast<StripMatrix*>(child->shaveBottom(n_rows).release())));

   m -= n_rows;

   recomputeNonzeros();
   border->recomputeNonzeros();

   return border;
}

void StripMatrix::write_to_streamDense(std::ostream& out) const {
   assert(!is_vertical);
   for (int i = 0; i < m; ++i) {
      write_to_streamDenseRow(out, i);
      out << "\n";
   }
}

void StripMatrix::write_to_streamDenseRow(std::ostream& out, int row) const {
   assert(!is_vertical);

   const int my_rank = PIPS_MPIgetRank(mpi_comm);

   std::ostringstream row_stream{};
   first->write_to_streamDenseRow(row_stream, row);

   if (my_rank != 0) {
      row_stream.str("");
      row_stream.clear();
   }

   for (auto& child : children)
      child->write_to_streamDenseRow(row_stream, row);

   const std::string my_row_part = row_stream.str();
   const std::string full_row = PIPS_MPIallgatherString(my_row_part, mpi_comm);

   if (my_rank == 0)
      out << full_row;

   if (last && my_rank == 0)
      last->write_to_streamDenseRow(out, row);
}

void StripMatrix::writeDashedLineToStream(std::ostream& out) const {
   assert(!is_vertical);

   std::stringstream row_stream{};

   first->writeDashedLineToStream(row_stream);

   if (PIPS_MPIgetRank(mpi_comm) != 0) {
      row_stream.str("");
      row_stream.clear();
   }

   for (auto& child : children)
      child->writeDashedLineToStream(row_stream);

   const std::string my_row_part = row_stream.str();
   const std::string full_row = PIPS_MPIallgatherString(my_row_part, mpi_comm);

   if (PIPS_MPIgetRank(mpi_comm) == 0)
      out << full_row;

   if (last && PIPS_MPIgetRank(mpi_comm) == 0)
      last->writeDashedLineToStream(out);
}

void StripMatrix::splitAlongTree(const DistributedTreeCallbacks& tree) {
   if (tree.getMapBlockSubTrees().empty())
      return;
   combineChildrenInNewChildren(tree.getMapBlockSubTrees(), tree.getChildComms());

   assert(tree.getChildComms().size() == children.size());

   for (size_t i = 0; i < children.size(); ++i) {
      children[i] = std::make_unique<StripMatrix>(is_vertical, std::move(children[i]), nullptr, tree.getChildComms()[i]);
   }

   assert(children.size() == tree.getChildren().size());
   const auto& tree_children = tree.getChildren();
   for (size_t i = 0; i < tree_children.size(); ++i) {
      const auto& tree_child = tree_children[i];
      if (tree_child->getCommWorkers() == MPI_COMM_NULL) {
         children[i] = std::make_unique<StringGenDummyMatrix>();
      } else if (tree_child->getSubRoot()) {
         assert(children[i]->first->is_a(kStripMatrix));
         dynamic_cast<StripMatrix&>(*children[i]->first).splitAlongTree(
            dynamic_cast<const DistributedTreeCallbacks&>(*tree_child->getSubRoot()));
      }
   }

   recomputeNonzeros();
}
