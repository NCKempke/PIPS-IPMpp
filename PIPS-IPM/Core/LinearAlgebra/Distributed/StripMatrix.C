/*
 * StripMatrix.C
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */

#include "StripMatrix.h"

#include "DistributedVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"
#include "DistributedTreeCallbacks.h"

#include "pipsdef.h"
#include <algorithm>

StripMatrix::StripMatrix(bool is_vertical, GeneralMatrix* first, GeneralMatrix* last, MPI_Comm mpi_comm_, bool is_view) : first(first),
      last(last), is_vertical(is_vertical), mpi_comm(mpi_comm_), distributed(PIPS_MPIgetDistributed(mpi_comm)),
      rank(PIPS_MPIgetRank(mpi_comm)), is_view{is_view} {
   assert(first);

   first->getSize(m, n);

   nonzeros += first->numberOfNonZeros();

   if (last) {
      long long ml, nl;
      last->getSize(ml, nl);
      nonzeros += last->numberOfNonZeros();

      if (is_vertical) {
         assert(n == nl);
         m += ml;
      }
      else {
         assert(m == ml);
         n += nl;
      }
   }
}

StripMatrix::~StripMatrix() {
   if (is_view)
      return;

   for (StripMatrix* child : children)
      delete child;

   delete first;
   delete last;
}

void StripMatrix::addChild(StripMatrix* child) {
   children.push_back(child);

   nonzeros += child->numberOfNonZeros();

   long long m_, n_;
   child->getSize(m_, n_);

   assert(child->is_vertical == this->is_vertical || child->is_a(kStringGenDummyMatrix));

   if (!child->is_a(kStringGenDummyMatrix)) {
      if (is_vertical) {
         assert(n == n_);
         m += m_;
      }
      else {
         assert(m == m_);
         n += n_;
      }
   }
}

bool StripMatrix::isEmpty() const {
   return !first && !last && children.empty() && m == 0 && n == 0;
}

int StripMatrix::is_a(int type) const {
   return (type == kStripMatrix || type == kStringMatrix || type == kGenMatrix);
}

double StripMatrix::inf_norm() const {
   double norm = 0.0;

   for (auto it : children)
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

   for (auto & it : children)
      it->scalarMult(num);
}


void StripMatrix::mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   if (is_vertical)
      multVertical(beta, y, alpha, x);
   else
      multHorizontal(beta, y, alpha, x, true);
}

void StripMatrix::transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   if (is_vertical)
      transMultVertical(beta, y, alpha, x, true);
   else
      transMultHorizontal(beta, y, alpha, x);
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
      else
         child->numberOfNonZeros();
   }

   if (dynamic_cast<const StripMatrix*>(first)) {
      auto& matstr = dynamic_cast<StripMatrix&>(*first);
      matstr.recomputeNonzeros();
      if (PIPS_MPIgetRank(matstr.mpi_comm) == 0)
         nonzeros += matstr.numberOfNonZeros();
      else
         matstr.numberOfNonZeros();
   }
   else if (PIPS_MPIiAmSpecial(distributed, mpi_comm)) {
      if (first)
         nonzeros += first->numberOfNonZeros();
      if (last)
         nonzeros += last->numberOfNonZeros();
   }

   PIPS_MPIgetSumInPlace(nonzeros, mpi_comm);
}

void StripMatrix::multVertical(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   assert(is_vertical);
   assert(y.children.size() == children.size());

   first->mult(beta, *y.first, alpha, x);

   if (last) {
      assert(y.last);
      last->mult(beta, *y.last, alpha, x);
   }

   for (size_t i = 0; i < children.size(); ++i) {
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(y.children[i]->isKindOf(kStochDummy));

      children[i]->multVertical(beta, *y.children[i], alpha, x);
   }
}

void StripMatrix::multHorizontal(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in, bool root) const {
   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(!is_vertical);
   assert(x.children.size() == children.size());
   assert((x.last && last) || (x.last == nullptr && last == nullptr));

   if (first->is_a(kStripMatrix)) {
      assert(!last);
      assert(children.empty());
      assert(!root);
      dynamic_cast<const StripMatrix*>(first)->multHorizontal(1.0, y, alpha, *x.first, false);
   }
   else {
      if (PIPS_MPIiAmSpecial(distributed, mpi_comm)) {
         first->mult(root ? beta : 1.0, y, alpha, *x.first);
         if (last)
            last->mult(1.0, y, alpha, *x.last);
      }
      else if (root)
         y.setToZero();

      for (size_t i = 0; i < children.size(); ++i)
         children[i]->multHorizontal(1.0, y, alpha, *x.children[i], false);

      if (distributed && root)
         PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);
   }
}

void StripMatrix::transMultVertical(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in, bool root) const {
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
      dynamic_cast<StripMatrix*>(first)->transMultVertical(1.0, y, alpha, *x.first, false);
   }
   else {
      if (rank == 0) {
         first->transMult(root ? beta : 1.0, y, alpha, *x.first);
         if (last)
            last->transMult(1.0, y, alpha, *x.last);
      }
      else if (root)
         y.setToZero();

      for (size_t i = 0; i < children.size(); i++)
         children[i]->transMultVertical(1.0, y, alpha, *x.children[i], false);

      if (distributed && root)
         PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);
   }
}

void StripMatrix::transMultHorizontal(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {

   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   assert(!is_vertical);
   assert(y.children.size() == children.size());
   if (y.last)
      assert(last);
   else
      assert(last == nullptr);

   first->transMult(beta, *y.first, alpha, x);

   if (last)
      last->transMult(beta, *y.last, alpha, x);

   for (size_t i = 0; i < children.size(); ++i) {
      if (children[i]->is_a(kStringGenDummyMatrix))
         assert(y.children[i]->isKindOf(kStochDummy));

      children[i]->transMultHorizontal(beta, *y.children[i], alpha, x);
   }
}

void StripMatrix::getColMinMaxVecHorizontal(bool get_min, bool initialize_vec, const Vector<double>* row_scale, Vector<double>& minmax_in) const {
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

void StripMatrix::getColMinMaxVecVertical(bool get_min, bool initialize_vec, const Vector<double>* row_scale_in, Vector<double>& minmax_) const {
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

      children[i]->getColMinMaxVecVertical(get_min, false, has_rowscale ? row_scale->children[i] : nullptr, minmax);
   }

   if (distributed) {
      if (get_min)
         PIPS_MPIminArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
      else
         PIPS_MPImaxArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
   }

   first->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->first : nullptr, minmax);

   if (last)
      last->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->last : nullptr, minmax);
}

/** DistributedVector<double> colScaleVec, SimpleVector<double> minmaxVec */
void
StripMatrix::getRowMinMaxVecHorizontal(bool get_min, bool initialize_vec, const Vector<double>* col_scale_in, Vector<double>& minmax_) const {
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

      children[i]->getRowMinMaxVecHorizontal(get_min, false, has_colscale ? col_scale->children[i] : nullptr, minmax);
   }

   if (distributed) {
      if (get_min)
         PIPS_MPIminArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
      else
         PIPS_MPImaxArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
   }

   first->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->first : nullptr, minmax);

   if (last)
      last->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->last : nullptr, minmax);
}

/** DistributedVector<double> minmaxVec, SimpleVector<double> colScaleVec */
void StripMatrix::getRowMinMaxVecVertical(bool get_min, bool initialize_vec, const Vector<double>* col_scale, Vector<double>& minmax_in) const {
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

void StripMatrix::getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) const {
   if (is_vertical)
      getRowMinMaxVecVertical(getMin, initializeVec, colScaleVec, minmaxVec);
   else
      getRowMinMaxVecHorizontal(getMin, initializeVec, colScaleVec, minmaxVec);
}


void StripMatrix::getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) const {
   if (is_vertical)
      getColMinMaxVecVertical(getMin, initializeVec, rowScaleVec, minmaxVec);
   else
      getColMinMaxVecHorizontal(getMin, initializeVec, rowScaleVec, minmaxVec);
}

void StripMatrix::columnScaleVertical(const Vector<double>& vec) {
   assert(is_vertical);

   first->columnScale(vec);

   for (auto & i : children)
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

   for (auto & i : children)
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

   for (auto i : children)
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

   for (auto i : children)
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

void StripMatrix::combineChildrenInNewChildren(const std::vector<unsigned int>& map_child_subchild, const std::vector<MPI_Comm>& child_comms) {

#ifndef NDEBUG
   const unsigned int n_new_children = getNDistinctValues(map_child_subchild);
   assert(child_comms.size() == n_new_children);
   assert(children.size() == map_child_subchild.size());
#endif

   unsigned int n_children{0};
   for (unsigned int i = 0; i < map_child_subchild.size(); ++i) {
      if (child_comms[n_children] == MPI_COMM_NULL) {
         addChild(new StringGenDummyMatrix());
         delete children[i];
         while (i + 1 != map_child_subchild.size() && map_child_subchild[i] == map_child_subchild[i + 1]) {
            ++i;
            delete children[i];
         }
      }
      else {
         SparseMatrix* empty_filler = is_vertical ? new SparseMatrix(0, n, 0) : new SparseMatrix(m, 0, 0);
         auto* new_child = new StripMatrix(is_vertical, empty_filler, nullptr, child_comms[n_children]);

         /* will not change size of StringGenMat since new_child is of size zero */
         addChild(new_child);

         new_child->addChild(children[i]);

         while (i + 1 != map_child_subchild.size() && map_child_subchild[i] == map_child_subchild[i + 1]) {
            ++i;
            new_child->addChild(children[i]);
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

GeneralMatrix* StripMatrix::shaveBottom(int n_rows) {
   assert(!is_vertical);
   assert(first);

   GeneralMatrix* mat_border = first->shaveBottom(n_rows);
   GeneralMatrix* matlink_border = last ? last->shaveBottom(n_rows) : nullptr;

   auto* border = new StripMatrix(false, mat_border, matlink_border, mpi_comm);

#ifndef NDEBUG
   int mB, nB;
   border->getSize(mB, nB);
   assert(mB == n_rows);
#endif

   for (auto& child : children)
      border->addChild(dynamic_cast<StripMatrix*>(child->shaveBottom(n_rows)));

   m -= n_rows;

   recomputeNonzeros();
   border->recomputeNonzeros();

   return border;
}

void StripMatrix::writeToStreamDense(std::ostream& out) const {
   assert(!is_vertical);
   for (int i = 0; i < m; ++i) {
      writeToStreamDenseRow(out, i);
      out << "\n";
   }
}

void StripMatrix::writeToStreamDenseRow(std::ostream& out, int row) const {
   assert(!is_vertical);

   const int my_rank = PIPS_MPIgetRank(mpi_comm);

   std::ostringstream row_stream{};
   first->writeToStreamDenseRow(row_stream, row);

   if (my_rank != 0) {
      row_stream.str("");
      row_stream.clear();
   }

   for (auto& child : children)
      child->writeToStreamDenseRow(row_stream, row);

   const std::string my_row_part = row_stream.str();
   const std::string full_row = PIPS_MPIallgatherString(my_row_part, mpi_comm);

   if (my_rank == 0)
      out << full_row;

   if (last && my_rank == 0)
      last->writeToStreamDenseRow(out, row);
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
      auto& child = children[i];
      auto* new_child = new StripMatrix(is_vertical, child, nullptr, tree.getChildComms()[i]);
      child = new_child;
   }

   assert(children.size() == tree.getChildren().size());
   const auto& tree_children = tree.getChildren();
   for (size_t i = 0; i < tree_children.size(); ++i) {
      const auto& tree_child = tree_children[i];
      if (tree_child->getCommWorkers() == MPI_COMM_NULL) {
         delete children[i];
         children[i] = new StringGenDummyMatrix();
      }
      else if (tree_child->getSubRoot()) {
         assert(children[i]->first->is_a(kStripMatrix));
         dynamic_cast<StripMatrix*>(children[i]->first)->splitAlongTree(dynamic_cast<const DistributedTreeCallbacks&>(*tree_child->getSubRoot()));
      }
   }

   recomputeNonzeros();
}
