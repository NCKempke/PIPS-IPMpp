/*
 * BorderedGenMatrix.C
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */


#include "BorderedMatrix.h"
#include "Vector.hpp"
#include "DistributedMatrix.h"
#include "SparseMatrix.h"
#include "DistributedVector.h"
#include "DoubleMatrixTypes.h"

#include "pipsdef.h"

#include <vector>
#include <cassert>
#include "StripMatrix.h"

BorderedMatrix::BorderedMatrix(std::shared_ptr<DistributedMatrix> inner_matrix_, std::unique_ptr<StripMatrix> border_left_, std::unique_ptr<StripMatrix> border_bottom_,
      std::unique_ptr<GeneralMatrix> bottom_left_block_, MPI_Comm mpi_comm_) : inner_matrix{std::move(inner_matrix_)}, border_left{std::move(border_left_)},
      border_bottom{std::move(border_bottom_)}, bottom_left_block{std::move(bottom_left_block_)}, mpi_comm(mpi_comm_), distributed(mpi_comm == MPI_COMM_NULL),
      rank(PIPS_MPIgetRank(mpi_comm)) {
   assert(inner_matrix);
   assert(border_left);
   assert(border_bottom);
   assert(bottom_left_block);

   m = bottom_left_block->n_rows();
   n = bottom_left_block->n_columns();

#ifndef NDEBUG
   const auto [m_bottom, n_bottom] = border_bottom->n_rows_columns();
   const auto [m_left, n_left] = border_left->n_rows_columns();
   const auto [m_inner, n_inner] = inner_matrix->n_rows_columns();

   assert(n == n_left);
   assert(m == m_bottom);

   assert(m_inner == m_left);
   assert(n_inner == n_bottom);
#endif

   m += border_left->n_rows();
   n += border_bottom->n_columns();
}

int BorderedMatrix::is_a(int type) const {
   return type == kBorderedGenMatrix || type == kBorderedMatrix || type == kGenMatrix;
}

void BorderedMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   this->mult(beta, y_in, alpha, x_in, &AbstractMatrix::mult);
}

void BorderedMatrix::mult_transform(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   const std::function<double(const double&)>& transform) const {

   auto mult = [&capture0 = std::as_const(transform)](const GeneralMatrix* mat, auto&& PH1, auto&& PH2, auto&& PH3, auto&& PH4) {
      mat->mult_transform(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
         std::forward<decltype(PH3)>(PH3), std::forward<decltype(PH4)>(PH4), capture0);
   };

   this->mult(beta, y_in, alpha, x_in, mult);
}

void BorderedMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   const std::function<void(const GeneralMatrix*, double, Vector<double>&, double, const Vector<double>&)>& mult) const {
   /* x row, y column shaped */
   assert(hasVecStructureForBorderedMat(x_in, true));
   assert(hasVecStructureForBorderedMat(y_in, false));

   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   mult(border_left.get(), beta, *y.children[0], alpha, *x.first);
   mult(inner_matrix.get(), 1.0, *y.children[0], alpha, *x.children[0]);

   mult(bottom_left_block.get(), beta, *y.last, alpha, *x.first);
   mult(border_bottom.get(), 1.0, *y.last, alpha, *x.children[0]);

}

/** y = beta * y + alpha * this^T * x */
void BorderedMatrix::transpose_mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   this->transpose_mult(beta, y_in, alpha, x_in, &AbstractMatrix::transpose_mult);
}

/** y = beta * y + alpha * this^T * x */
void BorderedMatrix::transpose_mult_transform(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   const std::function<double(const double&)>& transform) const {

   auto transpose_mult = [&capture0 = std::as_const(transform)](const GeneralMatrix* mat, auto&& PH1, auto&& PH2, auto&& PH3, auto&& PH4) {
      mat->transpose_mult_transform(std::forward<decltype(PH1)>(PH1), std::forward<decltype(PH2)>(PH2),
         std::forward<decltype(PH3)>(PH3), std::forward<decltype(PH4)>(PH4), capture0);
   };

   this->transpose_mult(beta, y_in, alpha, x_in, transpose_mult);
}

void BorderedMatrix::transpose_mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in,
   const std::function<void(const GeneralMatrix*, double, Vector<double>&, double, const Vector<double>&)>& transpose_mult) const {

   /* x column, y row shaped */
   assert(hasVecStructureForBorderedMat(x_in, false));
   assert(hasVecStructureForBorderedMat(y_in, true));

   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   transpose_mult(border_left.get(), beta, *y.first, alpha, *x.children[0]);
   transpose_mult(bottom_left_block.get(), 1.0, *y.first, alpha, *x.last);

   transpose_mult(inner_matrix.get(), beta, *y.children[0], alpha, *x.children[0]);
   transpose_mult(border_bottom.get(), 1.0, *y.children[0], alpha, *x.last);
}


double BorderedMatrix::inf_norm() const {
   double norm = -std::numeric_limits<double>::infinity();

   norm = std::max(norm, inner_matrix->inf_norm());
   norm = std::max(norm, border_left->inf_norm());
   norm = std::max(norm, border_bottom->inf_norm());
   norm = std::max(norm, bottom_left_block->inf_norm());

   return norm;
}

void BorderedMatrix::columnScale(const Vector<double>& vec) {
   assert(hasVecStructureForBorderedMat(vec, true));

   const auto& svec = dynamic_cast<const DistributedVector<double>&>(vec);

   border_left->columnScale(*svec.first);
   bottom_left_block->columnScale(*svec.first);

   inner_matrix->columnScale(*svec.children[0]);
   border_bottom->columnScale(*svec.children[0]);
}

void BorderedMatrix::rowScale(const Vector<double>& vec) {
   assert(hasVecStructureForBorderedMat(vec, false));

   const auto& svec = dynamic_cast<const DistributedVector<double>&>(vec);

   border_left->rowScale(*svec.children[0]);
   inner_matrix->rowScale(*svec.children[0]);

   bottom_left_block->rowScale(*svec.last);
   border_bottom->rowScale(*svec.last);
}

void BorderedMatrix::scalarMult(double num) {
   inner_matrix->scalarMult(num);
   border_left->scalarMult(num);
   border_bottom->scalarMult(num);
   bottom_left_block->scalarMult(num);
}

std::pair<long long, long long> BorderedMatrix::n_rows_columns() const {
   return {m, n};
}

long long BorderedMatrix::n_rows() const {
   return m;
}

long long BorderedMatrix::n_columns() const {
   return n;
}

void BorderedMatrix::getRowMinMaxVec(bool get_min, bool initialize_vec, const Vector<double>* col_scale_in, Vector<double>& minmax_in) const {
   assert(hasVecStructureForBorderedMat(minmax_in, false));

   const bool has_colscale = (col_scale_in != nullptr);
   if (has_colscale)
      assert(hasVecStructureForBorderedMat(*col_scale_in, true));

   auto& minmax = dynamic_cast<DistributedVector<double>&>(minmax_in);

   const DistributedVector<double>* col_scale = has_colscale ? dynamic_cast<const DistributedVector<double>*>(col_scale_in) : nullptr;

   border_left->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->first.get() : nullptr, *minmax.children[0]);
   inner_matrix->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->children[0].get() : nullptr, *minmax.children[0]);

   bottom_left_block->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->first.get() : nullptr, *minmax.last);
   border_bottom->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->children[0].get() : nullptr, *minmax.last);
}

void BorderedMatrix::getColMinMaxVec(bool get_min, bool initialize_vec, const Vector<double>* row_scale_in, Vector<double>& minmax_in) const {
   assert(hasVecStructureForBorderedMat(minmax_in, true));

   const bool has_rowscale = (row_scale_in != nullptr);
   if (has_rowscale)
      assert(hasVecStructureForBorderedMat(*row_scale_in, false));

   auto& minmax = dynamic_cast<DistributedVector<double>&>(minmax_in);
   const DistributedVector<double>* row_scale = has_rowscale ? dynamic_cast<const DistributedVector<double>*>(row_scale_in) : nullptr;

   border_left->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->children[0].get() : nullptr, *minmax.first);
   bottom_left_block->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->last.get() : nullptr, *minmax.first);

   inner_matrix->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->children[0].get() : nullptr, *minmax.children[0]);
   border_bottom->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->last.get() : nullptr, *minmax.children[0]);
}

void BorderedMatrix::sum_transform_rows(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   assert(hasVecStructureForBorderedMat(result_, false));

   auto& result = dynamic_cast<DistributedVector<double>&>(result_);

   border_left->sum_transform_rows(*result.children[0], transform);
   inner_matrix->sum_transform_rows(*result.children[0], transform);

   bottom_left_block->sum_transform_rows(*result.last, transform);
   border_bottom->sum_transform_rows(*result.last, transform);
}

void BorderedMatrix::sum_transform_columns(Vector<double>& result_,
   const std::function<double(const double&)>& transform) const {
   assert(hasVecStructureForBorderedMat(result_, true));

   auto& result = dynamic_cast<DistributedVector<double>&>(result_);

   border_left->sum_transform_columns(*result.first, transform);
   bottom_left_block->sum_transform_columns(*result.first, transform);

   inner_matrix->sum_transform_columns(*result.children[0], transform);
   border_bottom->sum_transform_columns(*result.children[0], transform);
}

void BorderedMatrix::addRowSums(Vector<double>& vec_) const {
   assert(hasVecStructureForBorderedMat(vec_, false));

   auto& vec = dynamic_cast<DistributedVector<double>&>(vec_);

   border_left->addRowSums(*vec.children[0]);
   inner_matrix->addRowSums(*vec.children[0]);

   bottom_left_block->addRowSums(*vec.last);
   border_bottom->addRowSums(*vec.last);
}

void BorderedMatrix::addColSums(Vector<double>& vec_) const {
   assert(hasVecStructureForBorderedMat(vec_, true));

   auto& vec = dynamic_cast<DistributedVector<double>&>(vec_);

   border_left->addColSums(*vec.first);
   bottom_left_block->addColSums(*vec.first);

   inner_matrix->addColSums(*vec.children[0]);
   border_bottom->addColSums(*vec.children[0]);
}

template<typename T>
bool BorderedMatrix::hasVecStructureForBorderedMat(const Vector<T>& vec, bool row_vec) const {
   const auto& vecs = dynamic_cast<const DistributedVector<T>&>(vec);

   if (vecs.children.size() != 1) {
      std::cout << "children\n";
      return false;
   }

   if (!vecs.children[0]) {
      std::cout << "child[0]\n";
      return false;
   }

   if (row_vec) {
      if (vecs.last) {
         std::cout << "row-first but root.last\n";
         return false;
      }
      if (!vecs.first) {
         std::cout << "row-first but NO root.first\n";
         return false;
      }
   }
   else {
      if (vecs.first) {
         std::cout << "col-first but root.first\n";
         return false;
      }
      if (!vecs.last) {
         std::cout << "col-first but NO root.last\n";
         return false;
      }
   }

   if (row_vec && vecs.length() != n) {
      std::cout << "ROW: root.length = " << vecs.length() << " != " << n << " = border.n \n";
      return false;
   }
   if (!row_vec && vecs.length() != m) {
      std::cout << "COL: root.length = " << vecs.length() << " != " << m << " = border.m \n";
      return false;
   }

   const auto [m_border, n_border] = border_left->n_rows_columns();

   if (row_vec && vecs.first->length() != n_border) {
      std::cout << "ROW: root.first.length = " << vecs.first->length() << " != " << n_border << " = border.n \n";
      return false;
   }
   if (!row_vec && vecs.length() != m) {
      std::cout << "COL: root.last.length = " << vecs.last->length() << " != " << m_border << " = border.m \n";
      return false;
   }

   return true;
}

void BorderedMatrix::write_to_streamDense(std::ostream& out) const {
   const int my_rank = PIPS_MPIgetRank(mpi_comm);
   const int size = PIPS_MPIgetSize(mpi_comm);
   const bool iAmDistrib = (size != 0);

   inner_matrix->write_to_streamDenseBordered(*border_left, out,0);

   const auto mL = this->bottom_left_block->n_rows();
   if (mL > 0) {
      if (iAmDistrib)
         MPI_Barrier(mpi_comm);

      // for each row r do:
      for (int r = 0; r < mL; r++) {
         MPI_Barrier(mpi_comm);

         // process Zero collects all the information and then prints it.
         if (my_rank == 0) {
            bottom_left_block->write_to_streamDenseRow(out, r);
            out << "|\t";
         }
         border_bottom->write_to_streamDenseRow(out, r);

         if (my_rank == 0)
            out << "\n";

         MPI_Barrier(mpi_comm);
      }
   }
   if (iAmDistrib)
      MPI_Barrier(mpi_comm);
}