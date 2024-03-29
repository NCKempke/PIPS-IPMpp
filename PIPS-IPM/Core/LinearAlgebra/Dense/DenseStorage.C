/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cassert>
#include <numeric>

#include "OoqpBlas.h"
#include "DenseStorage.h"
#include "SparseStorage.h"
#include "Vector.hpp"
#include "DenseVector.hpp"

int DenseStorageInstances = 0;

void DenseStorage::fromGetDiagonal(int idiag, Vector<double>& vec) const {
   const int extent = vec.length();

   assert(idiag + extent <= n);
   assert(idiag + extent <= m);

   auto& sv = (DenseVector<double>&) vec;

   for (int k = idiag; k < idiag + extent; k++) {
      sv[k] = M[k][k];
   }
}

void DenseStorage::getDiagonal(Vector<double>& vec) const {
   this->fromGetDiagonal(0, vec);
}


void DenseStorage::setToDiagonal(const Vector<double>& vec) {
   const int extent = vec.length();

   assert(extent <= n);
   assert(extent <= m);

   auto& sv = (DenseVector<double>&) vec;
   for (int i = 0; i < m; i++) {
      for (int k = 0; k < n; k++) {
         M[i][k] = 0.0;
      }
   }

   for (int k = 0; k < extent; k++) {
      M[k][k] = sv[k];
   }
}

std::pair<int,int> DenseStorage::n_rows_columns() const {
   return {m, n};
}

int DenseStorage::n_rows() const {
   return m;
}

int DenseStorage::n_columns() const {
   return n;
}

DenseStorage::DenseStorage(int min, int nin) {
   DenseStorageInstances++;
   m = min;
   n = nin;

   int mbar = (m > 0) ? m : 1; // We always allocate one row.
   try {
      neverDeleteElts = 0;

      M = new double* [mbar];
      if (m > 0) {
         M[0] = new double[m * n];
      }
      else {
         M[0] = nullptr;
      }
      int i;
      for (i = 1; i < m; i++)
         M[i] = M[0] + i * n;
   }
   catch (...) {
      std::cerr << "Out of memory in DenseStorage::DenseStorage(" << m << ", " << n << ")\n";
      throw;
   }
}

DenseStorage::DenseStorage(double A[], int min, int nin) {
   DenseStorageInstances++;
   m = min;
   n = nin;

   M = new double* [m];
   int i;
   for (i = 0; i < m; i++) {
      M[i] = A + i * n;
   }

   neverDeleteElts = 1;
}

void DenseStorage::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const {
   assert(col >= 0 && col + colExtent <= n);
   assert(row >= 0 && row < m);

   int k = 0;
   info = 0;

   for (int j = col; j < col + colExtent; j++) {
      if (M[row][j] != 0.0) {
         if (k < lenA) {
            // Add the element to A
            A[k] = M[row][j];
            jcolA[k] = j;
            k++;
         }
         else {
            // Count the number of additional elements needed in A
            info++;
         }
      }
   }
   nnz = k;
}

void DenseStorage::atPutSpRow(int row, const double A[], int lenA, const int jcolA[], int& info) {
   info = 0;
   int k;
   for (k = 0; k < lenA; k++) {
      M[row][jcolA[k]] = A[k];
   }
}


void DenseStorage::putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) {

   for (int i = 0; i < m; i++) {
      for (int k = 0; k < n; k++) {
         M[i][k] = 0.0;
      }
   }

   for (int k = 0; k < len; k++) {
      assert(irow[k] >= 0 && irow[k] < m);
      assert(jcol[k] >= 0 && jcol[k] < n);
      M[irow[k]][jcol[k]] = A[k];
   }
   info = 0;
}


void DenseStorage::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const {
   assert(row >= 0 && row + rowExtent <= m);
   assert(col >= 0 && col + colExtent <= n);

   for (int i = 0; i < rowExtent; i++) {
      //printf("\tcopying at %d in number of %d\n", i*lda, colExtent);
      memcpy(&A[i * lda], &M[i + row][col], colExtent * sizeof(double));
   }
}


DenseStorage::~DenseStorage() {
   DenseStorageInstances--;
   if (!neverDeleteElts) {
      delete[] M[0];
   }
   delete[] M;
}

void DenseStorage::atPutZeros(int row, int col, int rowExtent, int colExtent) {
   assert(row >= 0 && row + rowExtent <= m);
   assert(col >= 0 && col + colExtent <= n);

   // If assertions are turned off, clip to the actual size of this matrix
   row = (row >= 0) ? row : 0;
   col = (col >= 0) ? col : 0;
   int mrow = (row + rowExtent <= m) ? row + rowExtent : m;
   int ncol = (col + colExtent <= n) ? col + colExtent : n;

   for (int i = row; i < mrow; i++) {
      for (int j = col; j < ncol; j++) {
         this->M[i][j] = 0.0;
      }
   }
}

void DenseStorage::putZeros() {
   std::fill(M[0], M[0] + n * m, 0.0);
}

void DenseStorage::sum_transform_rows(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   assert(result_.length() == this->n_rows());

   auto& result = dynamic_cast<DenseVector<double>&>(result_);

   auto accumulate = [&transform] (const double& sum, const double& other) {
      return sum + transform(other);
   };

   for (int i = 0; i < m; ++i) {
      result[i] = std::accumulate(M[i], M[i] + n, result[i], accumulate);
   }
}

void DenseStorage::sum_transform_columns(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   assert(result_.length() == this->n_columns());

   auto& result = dynamic_cast<DenseVector<double>&>(result_);

   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         result[j] += transform(M[i][j]);
      }
   }
}

void DenseStorage::atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) {
   assert(row >= 0 && row + rowExtent <= m);
   assert(col >= 0 && col + colExtent <= n);

   for (int i = 0; i < rowExtent; i++) {
      memcpy(&M[i + row][col], &A[i * lda], colExtent * sizeof(double));
   }
}

void DenseStorage::atAddOuterProductOf(int row, int col, double alpha, double* x, int incx, int nx) {
   assert(row >= 0 && row + nx <= m);
   assert(col >= 0 && col + nx <= n);

   // If assertions are turned off, clip to the actual size of this matrix
   row = (row >= 0) ? row : 0;
   col = (col >= 0) ? col : 0;

   nx = (row + nx <= m) ? nx : m - row;
   nx = (col + nx <= n) ? nx : n - col;

   char fortranUplo = 'U';

   dsyr_(&fortranUplo, &nx, &alpha, x, &incx, &M[row][col], &n);
}


void DenseStorage::addToDiagonalAt(double alpha, double x[], int incx, int idiag, int extent) {
   assert(idiag + extent <= n);
   assert(idiag + extent <= m);

   // If assertions are off, clip to the actual size of this matrix
   if (idiag + extent < n)
      extent = n - idiag;
   if (idiag + extent < m)
      extent = m - idiag;

   int incy = n + 1;
   daxpy_(&extent, &alpha, x, &incx, &M[idiag][idiag], &incy);

}


void DenseStorage::atPutDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const DenseVector<double>&>(vvec);
   atPutDiagonal(idiag, &v[0], 1, v.length());
}

void DenseStorage::atAddDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const DenseVector<double>&>(vvec);
   atAddDiagonal(idiag, &v[0], 1, v.length());
}

void DenseStorage::atPutDiagonal(int idiag, const double x[], int incx, int extent) {
   for (int i = 0; i < extent; i++) {
      M[i + idiag][i + idiag] = x[i * incx];
   }
}

void DenseStorage::atAddDiagonal(int idiag, const double x[], int incx, int extent) {
   for (int i = 0; i < extent; i++) {
      M[i + idiag][i + idiag] += x[i * incx];
   }
}

void DenseStorage::diagonal_set_to_constant_from(int from, int length, double value) {
   assert(from + length < this->m);
   assert(from + length < this->n);

   for(int i = 0; i < from + length; ++i)
   {
      M[i][i] = value;
   }
}

void DenseStorage::fill_from_sparse(const SparseStorage& other)
{
   assert(this->n_rows_columns() == other.n_rows_columns());

   this->putZeros();
   for (int row = 0; row < other.n_rows(); ++row) {
      for(int j = other.krowM[row]; j < other.krowM[row + 1]; ++j)
      {
         const int col = other.jcolM[j];
         this->M[row][col] = other.M[j];
      }
   }
}

void DenseStorage::fill_from_dense(const DenseStorage& other)
{
   assert(this->n_rows_columns() == other.n_rows_columns());
   std::copy(other.M[0], other.M[0] + n * m, this->M[0]);
}

void DenseStorage::diagonal_add_constant_from(int from, int length, double value) {
   assert(from + length <= this->m);
   assert(from + length <= this->n);

   for (int i = from; i < from + length; ++i) {
      M[i][i] += value;
   }
}

double DenseStorage::inf_norm() const {
   double max = 0.0;

   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         if (std::fabs(M[i][j]) > max)
            max = std::fabs(M[i][j]);
      }
   }
   return max;
}

double DenseStorage::abminnormNonZero(double tol) const {
   double min = std::numeric_limits<double>::infinity();

   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         if (std::fabs(M[i][j]) < min && tol < std::fabs(M[i][j]))
            min = std::fabs(M[i][j]);
      }
   }
   return min;
}

int DenseStorage::non_zeros() const {
   int nnzs{0};
   for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
         if (M[i][j] != 0)
            ++nnzs;
      }
   }
   return nnzs;
}

void DenseStorage::columnScale(const Vector<double>& scale_in) {
   const auto& scale = dynamic_cast<const DenseVector<double>&>(scale_in);

   assert(scale.length() == n);

   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++)
         M[i][j] = M[i][j] * scale[j];
   }
}

void DenseStorage::scalarMult(double num) {
   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++)
         M[i][j] = M[i][j] * num;
   }
}

void DenseStorage::rowScale(const Vector<double>& scale_in) {
   const auto& scale = dynamic_cast<const DenseVector<double>&>(scale_in);

   assert(scale.length() == m);

   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++)
         M[i][j] = M[i][j] * scale[i];
   }
}

void DenseStorage::symmetricScale(const Vector<double>& scale_in) {
   const auto& scale = dynamic_cast<const DenseVector<double>&>(scale_in);

   assert(scale.length() == n);
   assert(scale.length() == m);

   for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++)
         M[i][j] = M[i][j] * scale[i] * scale[j];
   }
}