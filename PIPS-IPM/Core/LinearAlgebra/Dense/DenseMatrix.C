/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP 
 */

#include <cassert>
#include "DenseMatrix.h"
#include "OoqpBlas.h"
#include "SimpleVector.hpp"
#include "DoubleMatrixTypes.h"
#include "SparseMatrix.h"

int DenseMatrix::is_a(int type) const {
   return type == kDenseGenMatrix || type == kGenMatrix;
}

DenseMatrix::DenseMatrix(int size) {
   mStorage = std::make_shared<DenseStorage>(size, size);
}

DenseMatrix::DenseMatrix(double A[], int m, int n) {
   mStorage = std::make_shared<DenseStorage>(A, m, n);
}

DenseMatrix::DenseMatrix(int m, int n) {
   mStorage = std::make_shared<DenseStorage>(m, n);
}

void DenseMatrix::atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) {
   mStorage->atPutDense(row, col, A, lda, rowExtent, colExtent);
}

void DenseMatrix::atPutZeros(int row, int col, int rowExtent, int colExtent) {
   mStorage->atPutZeros(row, col, rowExtent, colExtent);
}

void DenseMatrix::putZeros() {
   mStorage->putZeros();
}

void DenseMatrix::sum_transform_rows(Vector<double>& result, const std::function<double(const double&)>& transform) const {
   assert(result.length() == mStorage->n_rows());
   mStorage->sum_transform_rows(result, transform);
}


void DenseMatrix::putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) {
   mStorage->putSparseTriple(irow, len, jcol, A, info);
}

void DenseMatrix::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const {
   mStorage->fromGetSpRow(row, col, A, lenA, jcolA, nnz, colExtent, info);
}

void DenseMatrix::atPutSpRow(int row, const double* A, int lenA, const int* jcolA, int& info) {
   mStorage->atPutSpRow(row, A, lenA, jcolA, info);
}

void DenseMatrix::getDiagonal(Vector<double>& vec) const {
   mStorage->getDiagonal(vec);
}

void DenseMatrix::setToDiagonal(const Vector<double>& vec) {
   mStorage->setToDiagonal(vec);
}

std::pair<long long, long long> DenseMatrix::n_rows_columns() const {
   return mStorage->n_rows_columns();
}

long long DenseMatrix::n_rows() const {
   return mStorage->n_rows();
}

long long DenseMatrix::n_columns() const {
   return mStorage->n_columns();
}

void DenseMatrix::atPutSubmatrix(int destRow, int destCol, const AbstractMatrix& Mat, int srcRow, int srcCol, int rowExtent, int colExtent) {
   int m = mStorage->m, n = mStorage->n;
   double** M = mStorage->M;

   assert(destRow >= 0 && destRow + rowExtent <= m);
   assert(destCol >= 0 && destCol + colExtent <= n);

   // If assertions are turned off, clip to the actual size of this matrix
   destRow = (destRow >= 0) ? destRow : 0;
   destCol = (destCol >= 0) ? destCol : 0;
   rowExtent = (destRow + rowExtent <= m) ? rowExtent : m - destRow;
   colExtent = (destCol + colExtent <= n) ? colExtent : n - destCol;

   Mat.fromGetDense(srcRow, srcCol, &M[destRow][destCol], n, rowExtent, colExtent);
}

void DenseMatrix::mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   char fortranTrans = 'T';
   int n = mStorage->n, m = mStorage->m;

   dgemv(&fortranTrans, &n, &m, &alpha, &mStorage->M[0][0], &n, x, &incx, &beta, y, &incy);
}

void DenseMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   char fortranTrans = 'T';
   int n = mStorage->n, m = mStorage->m;
   double** M = mStorage->M;
   int incx = 1, incy = 1;

   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   if (n != 0 && m != 0) {
      dgemv(&fortranTrans, &n, &m, &alpha, &M[0][0], &n, &x[0], &incx, &beta, &y[0], &incy);
   }
   else {
      if (m != 0)
         y.scale(beta);
   }
}


void DenseMatrix::transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   char fortranTrans = 'N';
   int n = mStorage->n, m = mStorage->m;
   double** M = mStorage->M;

   dgemv(&fortranTrans, &n, &m, &alpha, &M[0][0], &n, x, &incx, &beta, y, &incy);
}

void DenseMatrix::transMult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   char fortranTrans = 'N';
   int n = mStorage->n, m = mStorage->m;
   double** M = mStorage->M;
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);
   int incx = 1, incy = 1;

   if (m != 0 && n != 0) {
      dgemv(&fortranTrans, &n, &m, &alpha, &M[0][0], &n, &x[0], &incx, &beta, &y[0], &incy);
   }
   else {
      if (n != 0) {
         y.scale(beta);
      }
   }
}

double DenseMatrix::inf_norm() const {
   assert(mStorage != nullptr);
   return mStorage->inf_norm();
}
double DenseMatrix::abminnormNonZero(double tol) const {
   assert(mStorage != nullptr);
   return mStorage->abminnormNonZero(tol);
}

void DenseMatrix::write_to_stream(std::ostream& out) const {
   for (int i = 0; i < mStorage->m; i++) {
      for (int j = 0; j < mStorage->n; j++)
         out << mStorage->M[i][j] << "\t";

      out << "\n";
   }
}

void DenseMatrix::write_to_streamDense(std::ostream& out) const {
   write_to_stream(out);
}

void DenseMatrix::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const {
   const int m = mStorage->m;
   const int n = mStorage->n;

   assert(row >= 0 && row + rowExtent <= m);
   assert(col >= 0 && col + colExtent <= n);

   // If assertions are turned off, clip to the actual size of this matrix
   row = (row >= 0) ? row : 0;
   col = (col >= 0) ? col : 0;
   int lrow = (row + rowExtent <= m) ? rowExtent : m - row;
   int lcol = (col + colExtent <= n) ? colExtent : n - col;

   mStorage->fromGetDense(row, col, A, lda, lrow, lcol);
}

void DenseMatrix::atPutDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atPutDiagonal(idiag, v);
}

void DenseMatrix::atAddDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atAddDiagonal(idiag, v);
}

void DenseMatrix::fromGetDiagonal(int idiag, Vector<double>& v) const {
   mStorage->fromGetDiagonal(idiag, v);
}

void DenseMatrix::getRow(int rowIndex, Vector<double>& v_in) {
   assert (rowIndex >= 0 && rowIndex <= mStorage->m);
   auto& v = dynamic_cast<SimpleVector<double>&>(v_in);

   mStorage->fromGetDense(rowIndex, 0, &v[0], 1, 1, mStorage->n);
}

void DenseMatrix::columnScale(const Vector<double>& vec) {
   mStorage->columnScale(vec);
}

void DenseMatrix::symmetricScale(const Vector<double>& vec) {
   mStorage->symmetricScale(vec);
}

void DenseMatrix::rowScale(const Vector<double>& vec) {
   mStorage->columnScale(vec);
}

void DenseMatrix::scalarMult(double num) {
   mStorage->scalarMult(num);
}

// TODO : probably move to some utility class..
/* compute beta * res += alpha * this * mat where mat gets multiplied to the submatrix
 * starting at mul_start and the results gets added starting at res_start */
void
DenseMatrix::multMatAt(int row_start, int row_end, int col_offset_this, double beta, int row_start_res, int col_offset_result, DenseMatrix& res,
      double alpha, const SparseMatrix& mat) const {
   assert(0 <= col_offset_this);
   assert(0 <= row_start);
   assert(row_start <= row_end);
   assert(row_end <= mStorage->m);

   const auto mat_n = mat.n_columns();

   assert(col_offset_this <= mStorage->n && col_offset_this + mat.n_rows() <= mStorage->n);

   const int n_rows = row_end - row_start;
   assert(col_offset_result >= 0);
   assert(col_offset_result <= res.mStorage->n && col_offset_result + mat_n <= res.mStorage->n);
   assert(n_rows + row_start_res <= res.mStorage->m);

   const SparseStorage& mat_tp = mat.getTranspose().getStorage();
   for (int j = 0; j < n_rows; ++j) {
      for (int i = 0; i < mat_n; ++i) {
         if (beta != 1.0)
            res[row_start + j][col_offset_result + i] *= beta;

         const int col_start = mat_tp.krowM[i];
         const int col_end = mat_tp.krowM[i + 1];

         for (int k = col_start; k < col_end; ++k) {
            const int row = mat_tp.jcolM[k];
            const double val = mat_tp.M[k];

            assert(col_offset_result + i < res.mStorage->n);
            assert(row < mat.n_rows());
            assert(row + col_offset_this < mStorage->n);
            res[row_start_res + j][col_offset_result + i] += mStorage->M[row_start + j][row + col_offset_this] * val * alpha;
         }
      }
   }
}

void DenseMatrix::addMatAt(const SparseMatrix& mat, int mat_row_start, int mat_row_end, int this_row_0, int this_col_0) {
   const auto [mmat, nmat] = mat.n_rows_columns();
   if (mmat <= 0 || nmat <= 0)
      return;

   assert(0 <= mat_row_start && mat_row_start <= mat_row_end && mat_row_end - 1 < mmat);

   const int n_rows = mat_row_end - mat_row_start;
   assert(0 <= this_row_0 && this_row_0 + n_rows - 1 < mStorage->m);
   assert(0 <= this_col_0 && this_col_0 + nmat - 1 < mStorage->n);

   for (int row = 0; row < n_rows; ++row) {
      const int row_in_mat = mat_row_start + row;
      const int row_start = mat.krowM()[row_in_mat];
      const int row_end = mat.krowM()[row_in_mat + 1];

      for (int j = row_start; j < row_end; ++j) {
         const int col = mat.jcolM()[j];
         const double val = mat.M()[j];

         assert(col < nmat);
         assert(this_col_0 + col < mStorage->n);
         assert(this_row_0 + row < mStorage->m);
         this->operator[](this_row_0 + row)[this_col_0 + col] += val;
      }
   }
}
