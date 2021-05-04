/*
 * StringGenMatrix.C
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */

#include "StringGenMatrix.h"

#include "DistributedVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"
#include "DistributedTreeCallbacks.h"

#include "pipsdef.h"
#include <algorithm>

StringGenMatrix::StringGenMatrix(bool is_vertical, GenMatrix* mat, GenMatrix* mat_link, MPI_Comm mpi_comm_, bool is_view) : mat(mat),
      mat_link(mat_link), is_vertical(is_vertical), mpi_comm(mpi_comm_), distributed(PIPS_MPIgetDistributed(mpi_comm)),
      rank(PIPS_MPIgetRank(mpi_comm)), is_view{is_view} {
   assert(mat);

   mat->getSize(m, n);

   nonzeros += mat->numberOfNonZeros();

   if (mat_link) {
      long long ml, nl;
      mat_link->getSize(ml, nl);
      nonzeros += mat_link->numberOfNonZeros();

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

StringGenMatrix::~StringGenMatrix() {
   if (is_view)
      return;

   for (StringGenMatrix* child : children)
      delete child;

   delete mat;
   delete mat_link;
}

void StringGenMatrix::addChild(StringGenMatrix* child) {
   children.push_back(child);

   nonzeros += child->numberOfNonZeros();

   long long m_, n_;
   child->getSize(m_, n_);

   assert(child->is_vertical == this->is_vertical || child->isKindOf(kStringGenDummyMatrix));

   if (!child->isKindOf(kStringGenDummyMatrix)) {
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

bool StringGenMatrix::isEmpty() const {
   return !mat && !mat_link && children.empty() && m == 0 && n == 0;
}

int StringGenMatrix::isKindOf(int type) const {
   return (type == kStringGenMatrix || type == kStringMatrix || type == kGenMatrix);
}

double StringGenMatrix::abmaxnorm() const {
   double norm = 0.0;

   for (size_t it = 0; it < children.size(); it++)
      norm = std::max(norm, children[it]->abmaxnorm());

   if (distributed)
      norm = PIPS_MPIgetMax(norm, mpi_comm);

   norm = std::max(norm, mat->abmaxnorm());

   if (mat_link)
      norm = std::max(norm, mat_link->abmaxnorm());

   return norm;
}

void StringGenMatrix::scalarMult(double num) {
   mat->scalarMult(num);

   if (mat_link)
      mat_link->scalarMult(num);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->scalarMult(num);
}


void StringGenMatrix::mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   if (is_vertical)
      multVertical(beta, y, alpha, x);
   else
      multHorizontal(beta, y, alpha, x, true);
}

void StringGenMatrix::transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   if (is_vertical)
      transMultVertical(beta, y, alpha, x, true);
   else
      transMultHorizontal(beta, y, alpha, x);
}

int StringGenMatrix::numberOfNonZeros() const {
#ifndef NDEBUG
   const int nonzeros_now = nonzeros;
   const_cast<StringGenMatrix*>(this)->recomputeNonzeros();
   assert(nonzeros_now == nonzeros);
#endif
   return nonzeros;
}

void StringGenMatrix::recomputeNonzeros() {
   nonzeros = 0;

   for (auto& child : children) {
      child->recomputeNonzeros();
      if (PIPS_MPIgetRank(child->mpi_comm) == 0)
         nonzeros += child->numberOfNonZeros();
      else
         child->numberOfNonZeros();
   }

   if (dynamic_cast<const StringGenMatrix*>(mat)) {
      StringGenMatrix& matstr = dynamic_cast<StringGenMatrix&>(*mat);
      matstr.recomputeNonzeros();
      if (PIPS_MPIgetRank(matstr.mpi_comm) == 0)
         nonzeros += matstr.numberOfNonZeros();
      else
         matstr.numberOfNonZeros();
   }
   else if (PIPS_MPIiAmSpecial(distributed, mpi_comm)) {
      if (mat)
         nonzeros += mat->numberOfNonZeros();
      if (mat_link)
         nonzeros += mat_link->numberOfNonZeros();
   }

   PIPS_MPIgetSumInPlace(nonzeros, mpi_comm);
}

void StringGenMatrix::multVertical(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const SimpleVector<double>& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   DistributedVector<double>& y = dynamic_cast<DistributedVector<double>&>(y_in);

   assert(is_vertical);
   assert(y.children.size() == children.size());

   mat->mult(beta, *y.first, alpha, x);

   if (mat_link) {
      assert(y.last);
      mat_link->mult(beta, *y.last, alpha, x);
   }

   for (size_t i = 0; i < children.size(); ++i) {
      if (children[i]->isKindOf(kStringGenDummyMatrix))
         assert(y.children[i]->isKindOf(kStochDummy));

      children[i]->multVertical(beta, *y.children[i], alpha, x);
   }
}

void StringGenMatrix::multHorizontal(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in, bool root) const {
   const DistributedVector<double>& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   SimpleVector<double>& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(!is_vertical);
   assert(x.children.size() == children.size());
   assert((x.last && mat_link) || (x.last == nullptr && mat_link == nullptr));

   if (mat->isKindOf(kStringGenMatrix)) {
      assert(!mat_link);
      assert(children.empty());
      assert(!root);
      dynamic_cast<const StringGenMatrix*>(mat)->multHorizontal(1.0, y, alpha, *x.first, false);
   }
   else {
      if (PIPS_MPIiAmSpecial(distributed, mpi_comm)) {
         mat->mult(root ? beta : 1.0, y, alpha, *x.first);
         if (mat_link)
            mat_link->mult(1.0, y, alpha, *x.last);
      }
      else if (root)
         y.setToZero();

      for (size_t i = 0; i < children.size(); ++i)
         children[i]->multHorizontal(1.0, y, alpha, *x.children[i], false);

      if (distributed && root)
         PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);
   }
}

void StringGenMatrix::transMultVertical(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in, bool root) const {
   const DistributedVector<double>& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   SimpleVector<double>& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(is_vertical);
   assert(x.children.size() == children.size());
   assert(x.first && mat);
   assert((x.last && mat_link) || (x.last == nullptr && mat_link == nullptr));

   if (mat->isKindOf(kStringGenMatrix)) {
      assert(!mat_link);
      assert(children.empty());
      assert(!root);
      dynamic_cast<StringGenMatrix*>(mat)->transMultVertical(1.0, y, alpha, *x.first, false);
   }
   else {
      if (rank == 0) {
         mat->transMult(root ? beta : 1.0, y, alpha, *x.first);
         if (mat_link)
            mat_link->transMult(1.0, y, alpha, *x.last);
      }
      else if (root)
         y.setToZero();

      for (size_t i = 0; i < children.size(); i++)
         children[i]->transMultVertical(1.0, y, alpha, *x.children[i], false);

      if (distributed && root)
         PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);
   }
}

void StringGenMatrix::transMultHorizontal(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {

   const SimpleVector<double>& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   DistributedVector<double>& y = dynamic_cast<DistributedVector<double>&>(y_in);

   assert(!is_vertical);
   assert(y.children.size() == children.size());
   if (y.last)
      assert(mat_link);
   else
      assert(mat_link == nullptr);

   mat->transMult(beta, *y.first, alpha, x);

   if (mat_link)
      mat_link->transMult(beta, *y.last, alpha, x);

   for (size_t i = 0; i < children.size(); ++i) {
      if (children[i]->isKindOf(kStringGenDummyMatrix))
         assert(y.children[i]->isKindOf(kStochDummy));

      children[i]->transMultHorizontal(beta, *y.children[i], alpha, x);
   }
}

void StringGenMatrix::getColMinMaxVecHorizontal(bool get_min, bool initialize_vec, const Vector<double>* row_scale, Vector<double>& minmax_in) const {
   assert(!is_vertical);
   DistributedVector<double>& minmax = dynamic_cast<DistributedVector<double>&>(minmax_in);

   assert(minmax.first && mat);
   assert(minmax.children.size() == children.size());
   assert((minmax.last && mat_link) || (minmax.last == nullptr && mat_link == nullptr));


   for (size_t i = 0; i < children.size(); ++i) {
      assert(minmax.children[i]);
      if (children[i]->isKindOf(kStringGenDummyMatrix))
         assert(minmax.children[i]->isKindOf(kStochDummy));

      children[i]->getColMinMaxVecHorizontal(get_min, initialize_vec, row_scale, *minmax.children[i]);
   }

   mat->getColMinMaxVec(get_min, initialize_vec, row_scale, *minmax.first);

   if (mat_link)
      mat_link->getColMinMaxVec(get_min, initialize_vec, row_scale, *minmax.last);
}

void StringGenMatrix::getColMinMaxVecVertical(bool get_min, bool initialize_vec, const Vector<double>* row_scale_in, Vector<double>& minmax_) const {
   assert(is_vertical);
   const bool has_rowscale = (row_scale_in != nullptr);

   const DistributedVector<double>* row_scale = dynamic_cast<const DistributedVector<double>*>(row_scale_in);
   SimpleVector<double>& minmax = dynamic_cast<SimpleVector<double>&>(minmax_);

   assert(!has_rowscale || row_scale->children.size() == children.size());
   if (has_rowscale)
      assert((row_scale->last && mat_link) || (row_scale->last == nullptr && mat_link == nullptr));


   for (size_t i = 0; i < children.size(); i++) {
      if (has_rowscale)
         if (children[i]->isKindOf(kStringGenDummyMatrix))
            assert(row_scale->children[i]->isKindOf(kStochDummy));

      children[i]->getColMinMaxVecVertical(get_min, false, has_rowscale ? row_scale->children[i] : nullptr, minmax);
   }

   if (distributed) {
      if (get_min)
         PIPS_MPIminArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
      else
         PIPS_MPImaxArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
   }

   mat->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->first : nullptr, minmax);

   if (mat_link)
      mat_link->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->last : nullptr, minmax);
}

/** DistributedVector<double> colScaleVec, SimpleVector<double> minmaxVec */
void
StringGenMatrix::getRowMinMaxVecHorizontal(bool get_min, bool initialize_vec, const Vector<double>* col_scale_in, Vector<double>& minmax_) const {
   assert(!is_vertical);
   const bool has_colscale = (col_scale_in != nullptr);

   const DistributedVector<double>* col_scale = dynamic_cast<const DistributedVector<double>*>(col_scale_in);
   SimpleVector<double>& minmax = dynamic_cast<SimpleVector<double>&>(minmax_);
   assert(!has_colscale || col_scale->children.size() == children.size());
   if (has_colscale)
      assert((col_scale->last && mat_link) || (col_scale->last == nullptr && mat_link == nullptr));


   for (size_t i = 0; i < children.size(); i++) {
      if (has_colscale)
         if (children[i]->isKindOf(kStringGenDummyMatrix))
            assert(col_scale->children[i]->isKindOf(kStochDummy));

      children[i]->getRowMinMaxVecHorizontal(get_min, false, has_colscale ? col_scale->children[i] : nullptr, minmax);
   }

   if (distributed) {
      if (get_min)
         PIPS_MPIminArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
      else
         PIPS_MPImaxArrayInPlace(minmax.elements(), minmax.length(), mpi_comm);
   }

   mat->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->first : nullptr, minmax);

   if (mat_link)
      mat_link->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->last : nullptr, minmax);
}

/** DistributedVector<double> minmaxVec, SimpleVector<double> colScaleVec */
void StringGenMatrix::getRowMinMaxVecVertical(bool get_min, bool initialize_vec, const Vector<double>* col_scale, Vector<double>& minmax_in) const {
   assert(is_vertical);

   DistributedVector<double>& minmax = dynamic_cast<DistributedVector<double>&>(minmax_in);

   assert(minmax.first && mat);
   assert(minmax.children.size() == children.size());
   assert((minmax.last && mat_link) || (minmax.last == nullptr && mat_link == nullptr));

   mat->getRowMinMaxVec(get_min, initialize_vec, col_scale, *minmax.first);

   for (size_t i = 0; i < children.size(); ++i) {
      assert(minmax.children[i]);
      if (children[i]->isKindOf(kStringGenDummyMatrix))
         assert(minmax.children[i]->isKindOf(kStochDummy));

      children[i]->getRowMinMaxVecVertical(get_min, initialize_vec, col_scale, *minmax.children[i]);
   }

   if (mat_link)
      mat_link->getRowMinMaxVec(get_min, initialize_vec, col_scale, *minmax.last);
}

void StringGenMatrix::getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) {
   if (is_vertical)
      getRowMinMaxVecVertical(getMin, initializeVec, colScaleVec, minmaxVec);
   else
      getRowMinMaxVecHorizontal(getMin, initializeVec, colScaleVec, minmaxVec);
}


void StringGenMatrix::getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) {
   if (is_vertical)
      getColMinMaxVecVertical(getMin, initializeVec, rowScaleVec, minmaxVec);
   else
      getColMinMaxVecHorizontal(getMin, initializeVec, rowScaleVec, minmaxVec);
}

void StringGenMatrix::columnScaleVertical(const Vector<double>& vec) {
   assert(is_vertical);

   mat->columnScale(vec);

   for (size_t i = 0; i < children.size(); i++)
      children[i]->columnScaleVertical(vec);

   if (mat_link)
      mat_link->columnScale(vec);
}

void StringGenMatrix::columnScaleHorizontal(const Vector<double>& vec_in) {
   assert(!is_vertical);

   const DistributedVector<double>& vec = dynamic_cast<const DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && mat_link) || (vec.last == nullptr && mat_link == nullptr));

   mat->columnScale(*vec.first);

   for (size_t i = 0; i < children.size(); ++i) {
      assert(vec.children[i]);
      if (children[i]->isKindOf(kStochGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->columnScaleHorizontal(*vec.children[i]);
   }

   if (mat_link)
      mat_link->columnScale(*vec.last);
}

void StringGenMatrix::rowScaleVertical(const Vector<double>& vec_in) {
   assert(is_vertical);

   const DistributedVector<double>& vec = dynamic_cast<const DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && mat_link) || (vec.last == nullptr && mat_link == nullptr));

   mat->rowScale(*vec.first);

   for (size_t i = 0; i < children.size(); i++) {
      assert(vec.children[i]);
      if (children[i]->isKindOf(kStringGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->rowScaleVertical(*vec.children[i]);
   }

   if (mat_link)
      mat_link->rowScale(*vec.last);
}

void StringGenMatrix::rowScaleHorizontal(const Vector<double>& vec) {
   assert(!is_vertical);
   mat->rowScale(vec);

   if (mat_link)
      mat_link->rowScale(vec);

   for (size_t i = 0; i < children.size(); ++i)
      children[i]->rowScaleHorizontal(vec);
}

void StringGenMatrix::columnScale(const Vector<double>& vec) {
   if (is_vertical)
      columnScaleVertical(vec);
   else
      columnScaleHorizontal(vec);
}

void StringGenMatrix::rowScale(const Vector<double>& vec) {
   if (is_vertical)
      rowScaleVertical(vec);
   else
      rowScaleHorizontal(vec);
}

void StringGenMatrix::addRowSumsVertical(Vector<double>& vec_in) const {
   assert(is_vertical);

   DistributedVector<double>& vec = dynamic_cast<DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && mat_link) || (vec.last == nullptr && mat_link == nullptr));

   mat->addRowSums(*vec.first);

   for (size_t i = 0; i < children.size(); i++) {
      assert(vec.children[i]);
      if (children[i]->isKindOf(kStringGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->addRowSumsVertical(*vec.children[i]);
   }

   if (mat_link)
      mat_link->addRowSums(*vec.last);
}

void StringGenMatrix::addRowSumsHorizontal(Vector<double>& vec_in) const {
   assert(!is_vertical);
   SimpleVector<double>& vec = dynamic_cast<SimpleVector<double>&>(vec_in);

   for (size_t i = 0; i < children.size(); i++)
      children[i]->addRowSumsHorizontal(vec);

   if (distributed)
      PIPS_MPIsumArrayInPlace(vec.elements(), vec.length(), mpi_comm);

   mat->addRowSums(vec);

   if (mat_link)
      mat_link->addRowSums(vec);
}

void StringGenMatrix::addColSumsVertical(Vector<double>& vec_in) const {
   assert(is_vertical);
   SimpleVector<double>& vec = dynamic_cast<SimpleVector<double>&>(vec_in);

   for (size_t i = 0; i < children.size(); i++)
      children[i]->addColSumsVertical(vec);

   if (distributed)
      PIPS_MPIsumArrayInPlace(vec.elements(), vec.length(), mpi_comm);

   mat->addColSums(vec);

   if (mat_link)
      mat_link->addColSums(vec);
}

void StringGenMatrix::addColSumsHorizontal(Vector<double>& vec_in) const {
   assert(!is_vertical);

   DistributedVector<double>& vec = dynamic_cast<DistributedVector<double>&>(vec_in);

   assert(vec.first);
   assert(vec.children.size() == children.size());
   assert((vec.last && mat_link) || (vec.last == nullptr && mat_link == nullptr));

   mat->addColSums(*vec.first);

   for (size_t i = 0; i < children.size(); i++) {
      assert(vec.children[i]);
      if (children[i]->isKindOf(kStringGenDummyMatrix))
         assert(vec.children[i]->isKindOf(kStochDummy));

      children[i]->addColSumsHorizontal(*vec.children[i]);
   }

   if (mat_link)
      mat_link->addColSums(*vec.last);
}

void StringGenMatrix::addRowSums(Vector<double>& vec) const {
   if (is_vertical)
      addRowSumsVertical(vec);
   else
      addRowSumsHorizontal(vec);
}

void StringGenMatrix::addColSums(Vector<double>& vec) const {
   if (is_vertical)
      addColSumsVertical(vec);
   else
      addColSumsHorizontal(vec);
}

void StringGenMatrix::combineChildrenInNewChildren(const std::vector<unsigned int>& map_child_subchild, const std::vector<MPI_Comm>& child_comms) {

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
         SparseGenMatrix* empty_filler = is_vertical ? new SparseGenMatrix(0, n, 0) : new SparseGenMatrix(m, 0, 0);
         StringGenMatrix* new_child = new StringGenMatrix(is_vertical, empty_filler, nullptr, child_comms[n_children]);

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

GenMatrix* StringGenMatrix::shaveBottom(int n_rows) {
   assert(!is_vertical);
   assert(mat);

   GenMatrix* mat_border = mat->shaveBottom(n_rows);
   GenMatrix* matlink_border = mat_link ? mat_link->shaveBottom(n_rows) : nullptr;

   StringGenMatrix* border = new StringGenMatrix(false, mat_border, matlink_border, mpi_comm);

#ifndef NDEBUG
   int mB, nB;
   border->getSize(mB, nB);
   assert(mB == n_rows);
#endif

   for (auto& child : children)
      border->addChild(dynamic_cast<StringGenMatrix*>(child->shaveBottom(n_rows)));

   m -= n_rows;

   recomputeNonzeros();
   border->recomputeNonzeros();

   return border;
}

void StringGenMatrix::writeToStreamDense(std::ostream& out) const {
   assert(!is_vertical);
   for (int i = 0; i < m; ++i) {
      writeToStreamDenseRow(out, i);
      out << "\n";
   }
}

void StringGenMatrix::writeToStreamDenseRow(std::ostream& out, int row) const {
   assert(!is_vertical);

   const int my_rank = PIPS_MPIgetRank(mpi_comm);

   std::ostringstream row_stream{};
   mat->writeToStreamDenseRow(row_stream, row);

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

   if (mat_link && my_rank == 0)
      mat_link->writeToStreamDenseRow(out, row);
}

void StringGenMatrix::writeDashedLineToStream(std::ostream& out) const {
   assert(!is_vertical);

   std::stringstream row_stream{};

   mat->writeDashedLineToStream(row_stream);

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

   if (mat_link && PIPS_MPIgetRank(mpi_comm) == 0)
      mat_link->writeDashedLineToStream(out);
}

void StringGenMatrix::splitAlongTree(const DistributedTreeCallbacks& tree) {
   if (tree.getMapBlockSubTrees().empty())
      return;
   combineChildrenInNewChildren(tree.getMapBlockSubTrees(), tree.getChildComms());

   assert(tree.getChildComms().size() == children.size());

   for (size_t i = 0; i < children.size(); ++i) {
      auto& child = children[i];
      StringGenMatrix* new_child = new StringGenMatrix(is_vertical, child, nullptr, tree.getChildComms()[i]);
      child = new_child;
   }

   assert(children.size() == tree.getChildren().size());
   const auto tree_children = tree.getChildren();
   for (size_t i = 0; i < tree_children.size(); ++i) {
      const auto& tree_child = tree_children[i];
      if (tree_child->getCommWorkers() == MPI_COMM_NULL) {
         delete children[i];
         children[i] = new StringGenDummyMatrix();
      }
      else if (tree_child->getSubRoot()) {
         assert(children[i]->mat->isKindOf(kStringGenMatrix));
         dynamic_cast<StringGenMatrix*>(children[i]->mat)->splitAlongTree(dynamic_cast<const DistributedTreeCallbacks&>(*tree_child->getSubRoot()));
      }
   }

   recomputeNonzeros();
}
