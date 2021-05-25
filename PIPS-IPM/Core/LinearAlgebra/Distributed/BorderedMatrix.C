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

BorderedMatrix::BorderedMatrix(DistributedMatrix* inner_matrix_, StripMatrix* border_left_, StripMatrix* border_bottom_,
      SparseMatrix* bottom_left_block_, MPI_Comm mpi_comm_) : inner_matrix{inner_matrix_}, border_left{border_left_},
      border_bottom{border_bottom_}, bottom_left_block{bottom_left_block_}, mpi_comm(mpi_comm_), distributed(mpi_comm == MPI_COMM_NULL),
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

BorderedMatrix::~BorderedMatrix() {
   delete bottom_left_block;
   delete border_bottom;
   delete border_left;
   delete inner_matrix;
}

int BorderedMatrix::is_a(int type) const {
   return type == kBorderedGenMatrix || type == kBorderedMatrix || type == kGenMatrix;
}

void BorderedMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   /* x row, y column shaped */
   assert(hasVecStructureForBorderedMat(x_in, true));
   assert(hasVecStructureForBorderedMat(y_in, false));

   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   border_left->mult(beta, *y.children[0], alpha, *x.first);
   inner_matrix->mult(1.0, *y.children[0], alpha, *x.children[0]);

   bottom_left_block->mult(beta, *y.last, alpha, *x.first);
   border_bottom->mult(1.0, *y.last, alpha, *x.children[0]);
}

/** y = beta * y + alpha * this^T * x */
void BorderedMatrix::transMult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   /* x column, y row shaped */
   assert(hasVecStructureForBorderedMat(x_in, false));
   assert(hasVecStructureForBorderedMat(y_in, true));

   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_in);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_in);

   border_left->transMult(beta, *y.first, alpha, *x.children[0]);
   bottom_left_block->transMult(1.0, *y.first, alpha, *x.last);

   inner_matrix->transMult(beta, *y.children[0], alpha, *x.children[0]);
   border_bottom->transMult(1.0, *y.children[0], alpha, *x.last);
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

   border_left->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->first : nullptr, *minmax.children[0]);
   inner_matrix->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->children[0] : nullptr, *minmax.children[0]);

   bottom_left_block->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->first : nullptr, *minmax.last);
   border_bottom->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->children[0] : nullptr, *minmax.last);
}

void BorderedMatrix::getColMinMaxVec(bool get_min, bool initialize_vec, const Vector<double>* row_scale_in, Vector<double>& minmax_in) const {
   assert(hasVecStructureForBorderedMat(minmax_in, true));

   const bool has_rowscale = (row_scale_in != nullptr);
   if (has_rowscale)
      assert(hasVecStructureForBorderedMat(*row_scale_in, false));

   auto& minmax = dynamic_cast<DistributedVector<double>&>(minmax_in);
   const DistributedVector<double>* row_scale = has_rowscale ? dynamic_cast<const DistributedVector<double>*>(row_scale_in) : nullptr;

   border_left->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->children[0] : nullptr, *minmax.first);
   bottom_left_block->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->last : nullptr, *minmax.first);

   inner_matrix->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->children[0] : nullptr, *minmax.children[0]);
   border_bottom->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->last : nullptr, *minmax.children[0]);
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
      std::cout << "children" << std::endl;
      return false;
   }

   if (vecs.children[0] == nullptr) {
      std::cout << "child[0]" << std::endl;
      return false;
   }

   if (row_vec) {
      if (vecs.last != nullptr) {
         std::cout << "row-first but root.last" << std::endl;
         return false;
      }
      if (vecs.first == nullptr) {
         std::cout << "row-first but NO root.first" << std::endl;
         return false;
      }
   }
   else {
      if (vecs.first != nullptr) {
         std::cout << "col-first but root.first" << std::endl;
         return false;
      }
      if (vecs.last == nullptr) {
         std::cout << "col-first but NO root.last" << std::endl;
         return false;
      }

   }

   if (row_vec && vecs.length() != n) {
      std::cout << "ROW: root.length = " << vecs.length() << " != " << n << " = border.n " << std::endl;
      return false;
   }
   if (!row_vec && vecs.length() != m) {
      std::cout << "COL: root.length = " << vecs.length() << " != " << m << " = border.m " << std::endl;
      return false;
   }

   const auto [m_border, n_border] = border_left->n_rows_columns();

   if (row_vec && vecs.first->length() != n_border) {
      std::cout << "ROW: root.first.length = " << vecs.first->length() << " != " << n_border << " = border.n " << std::endl;
      return false;
   }
   if (!row_vec && vecs.length() != m) {
      std::cout << "COL: root.last.length = " << vecs.last->length() << " != " << m_border << " = border.m " << std::endl;
      return false;
   }

   return true;
}

void BorderedMatrix::writeToStreamDense(std::ostream& out) const {
   const int my_rank = PIPS_MPIgetRank(mpi_comm);
   const int size = PIPS_MPIgetSize(mpi_comm);
   const bool iAmDistrib = (size != 0);

   inner_matrix->writeToStreamDenseBordered(*border_left, out,0);

   const auto mL = this->bottom_left_block->n_rows();
   if (mL > 0) {
      if (iAmDistrib)
         MPI_Barrier(mpi_comm);

      // for each row r do:
      for (int r = 0; r < mL; r++) {
         MPI_Barrier(mpi_comm);

         // process Zero collects all the information and then prints it.
         if (my_rank == 0) {
            bottom_left_block->writeToStreamDenseRow(out, r);
            out << "|\t";
         }
         border_bottom->writeToStreamDenseRow(out, r);

         if (my_rank == 0)
            out << "\n";

         MPI_Barrier(mpi_comm);
      }
   }
   if (iAmDistrib)
      MPI_Barrier(mpi_comm);
}