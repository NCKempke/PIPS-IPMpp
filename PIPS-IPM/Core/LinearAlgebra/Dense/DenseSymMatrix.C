/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include <cassert>
#include "DenseSymMatrix.h"
#include "OoqpBlas.h"
#include "SimpleVector.h"

#include "DenseGenMatrix.h"
#include "SparseGenMatrix.h"
#include "SparseSymMatrix.h"

#include "DoubleMatrixTypes.h"

// TODO : move to Blas header..
extern "C" void dsyrk_(char* UPLO, char* TRANS, int* N, int* K, double* alpha, double* A, int* lda, double* beta, double* C, int* ldc);

int DenseSymMatrix::isKindOf(int matrixType) const {
   return matrixType == kDenseSymMatrix || matrixType == kSymMatrix;
}

DenseSymMatrix::DenseSymMatrix(int size) {
   mStorage = std::make_shared<DenseStorage>(size, size);
}

DenseSymMatrix::DenseSymMatrix(double Q[], int size) {
   mStorage = std::make_shared<DenseStorage>(Q, size, size);
}

void DenseSymMatrix::symAtPutZeros(int row, int col, int rowExtent, int colExtent) {
   mStorage->atPutZeros(row, col, rowExtent, colExtent);
}

void DenseSymMatrix::putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) {
   mStorage->putSparseTriple(irow, len, jcol, A, info);
}

void DenseSymMatrix::atAddOuterProductOf(int row, int col, double alpha, double* x, int incx, int nx) {
   mStorage->atAddOuterProductOf(row, col, alpha, x, incx, nx);
}

void DenseSymMatrix::getDiagonal(Vector<double>& vec) const {
   mStorage->getDiagonal(vec);
}

void DenseSymMatrix::setToDiagonal(const Vector<double>& vec) {
   mStorage->setToDiagonal(vec);
}

void DenseSymMatrix::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const {
   if (col + colExtent < row + 1) {
      mStorage->fromGetSpRow(row, col, A, lenA, jcolA, nnz, colExtent, info);
   }
   else {
      if (col <= row) {
         mStorage->fromGetSpRow(row, col, A, lenA, jcolA, nnz, row - col + 1, info);
      }
   }
}

void DenseSymMatrix::getSize(long long& m, long long& n) const {
   m = mStorage->m;
   n = mStorage->n;
}

void DenseSymMatrix::getSize(int& m, int& n) const {
   m = mStorage->m;
   n = mStorage->n;
}

long long DenseSymMatrix::size() const {
   return mStorage->m;
}

void DenseSymMatrix::symAtPutSubmatrix(int destRow, int destCol, const DoubleMatrix& Mat, int srcRow, int srcCol, int rowExtent, int colExtent) {
   const int m = mStorage->m;
   const int n = mStorage->n;
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

void DenseSymMatrix::mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   char fortranUplo = 'U';
   int n = mStorage->n;

   dsymv_(&fortranUplo, &n, &alpha, &mStorage->M[0][0], &n, x, &incx, &beta, y, &incy);
}

void DenseSymMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   char fortranUplo = 'U';
   int n = mStorage->n;
   auto& y = (SimpleVector<double>&) y_in;
   auto& x = (SimpleVector<double>&) x_in;
   int incx = 1, incy = 1;

   if (n != 0) {
      dsymv_(&fortranUplo, &n, &alpha, &mStorage->M[0][0], &n, &x[0], &incx, &beta, &y[0], &incy);
   }
}

void DenseSymMatrix::transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const {
   this->mult(beta, y, alpha, x);
}

void DenseSymMatrix::transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   this->mult(beta, y, incy, alpha, x, incx);
}

double DenseSymMatrix::abmaxnorm() const {
   return mStorage->abmaxnorm();
}

double DenseSymMatrix::abminnormNonZero(double tol) const {
   return mStorage->abminnormNonZero(tol);
}

void DenseSymMatrix::writeToStream(std::ostream& out) const {
   for (int i = 0; i < mStorage->m; i++) {
      for (int j = 0; j < mStorage->n; j++)
         out << mStorage->M[i][j] << "\t";

      out << "\n";
   }
}

void DenseSymMatrix::writeToStreamDense(std::ostream& out) const {
   writeToStream(out);
}

void DenseSymMatrix::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const {
   const int m = mStorage->m;
   const int n = mStorage->n;
   double** M = mStorage->M;

   assert(row >= 0 && row + rowExtent <= m);
   assert(col >= 0 && col + colExtent <= n);

   // If assertions are turned off, clip to the actual size of this matrix
   row = (row >= 0) ? row : 0;
   col = (col >= 0) ? col : 0;
   rowExtent = (row + rowExtent <= m) ? rowExtent : m - row;
   colExtent = (col + colExtent <= n) ? colExtent : n - col;

   for (int i = row; i < row + rowExtent; i++) {
      int j;
      for (j = col; j <= i && j < col + colExtent; j++) {
         A[(i - row) * lda + j - col] = M[i][j];
      }
      for (; j < col + colExtent; j++) {
         A[(i - row) * lda + j - col] = M[j][i];
      }
   }
}


void DenseSymMatrix::atPutDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atPutDiagonal(idiag, v);
}

void DenseSymMatrix::atAddDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atAddDiagonal(idiag, v);
}

void DenseSymMatrix::fromGetDiagonal(int idiag, Vector<double>& v) const {
   mStorage->fromGetDiagonal(idiag, v);
}

void DenseSymMatrix::diagonal_add_constant_from(int from, int length, double value) {
   assert(0 <= from);
   assert(0 <= length);
   assert(from + length < this->size());
   mStorage->diagonal_add_constant_from(from, length, value);
}

void DenseSymMatrix::diagonal_set_to_constant_from(int from, int length, double value) {
   assert(0 <= from);
   assert(0 <= length);
   assert(from + length < this->size());
   mStorage->diagonal_set_to_constant_from(from, length, value);
}

void DenseSymMatrix::symAtPutSpRow(int row, const double A[], int lenA, const int jcolA[], int& info) {
   // Lower triangular put
   int lA = lenA;
   while (lA > 0 && jcolA[lA - 1] > row)
      lA--;
   if (lA > 0) {
      mStorage->atPutSpRow(row, A, lA, jcolA, info);
   }
   else {
      info = 0;
   }
}

void DenseSymMatrix::symAtPutDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) {
   mStorage->atPutDense(row, col, A, lda, rowExtent, colExtent);
}

void DenseSymMatrix::symmetricScale(const Vector<double>& vec) {
   mStorage->symmetricScale(vec);
}

void DenseSymMatrix::columnScale(const Vector<double>& vec) {
   mStorage->columnScale(vec);
}

void DenseSymMatrix::rowScale(const Vector<double>& vec) {
   mStorage->rowScale(vec);
}

void DenseSymMatrix::scalarMult(double num) {
   mStorage->scalarMult(num);
}

/* updates the upper left block only  if sizes does not matches */
void DenseSymMatrix::matMult(double alpha, GenMatrix& A_, int transA, GenMatrix& B_, int transB, double beta) {

   auto& A = dynamic_cast<DenseGenMatrix&>(A_);
   auto& B = dynamic_cast<DenseGenMatrix&>(B_);

   // the other way around since fortran stores column-wise and we store row-wise
   char forTransA = (transA == 0 ? 'T' : 'N');
   char forTransB = (transB == 0 ? 'T' : 'N');

   DenseSymMatrix& C = *this;

   int m, n, k, kB;
   int ldc = mStorage->m;

   if (!transA) {
      A.getSize(m, k);
   }
   else {
      A.getSize(k, m);
   }

   if (!transB) {
      B.getSize(kB, n);
   }
   else {
      B.getSize(n, kB);
   }

   assert(k == kB);

   assert(mStorage->m >= m);
   assert(mStorage->n >= n);

   double** AA = A.mStorage->M;
   double** BB = B.mStorage->M;
   double** CC = C.mStorage->M;

   dgemm_(&forTransA, &forTransB, &m, &n, &k, &alpha, &AA[0][0], &m, &BB[0][0], &n, &beta, &CC[0][0], &ldc);
}


void DenseSymMatrix::symAtPutSubmatrix(int destRow, int destCol, const DoubleMatrix& Mat, int srcRow, int srcCol, int rowExtent, int colExtent,
      int forceSymUpdate) {

   if (forceSymUpdate == 0) {
      symAtPutSubmatrix(destRow, destCol, Mat, srcRow, srcCol, rowExtent, colExtent);

      return;
   }

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

   for (int i = destRow; i < destRow + rowExtent; i++) {
      for (int j = destCol; j < destCol + colExtent; j++) {
         M[j][i] = M[i][j];
      }
   }
}

void DenseSymMatrix::atRankkUpdate(double alpha, double beta, DenseGenMatrix& U, int trans) {
   //-----------------------------------------------
   // setup if the U is stored in column-major form
   // (FORTRAN Style)
   //   char UPLO  = 'U'; char TRANS = trans==0?'N':'T';

   //   U.getSize(n,k); ldu=n;
   //   if(trans) k=n;

   //   n = mStorage->n;
   //   lda=n;
   //----------------------------------------------

   // U and 'this' are stored in row-major form -> a little change in passing params to FORTRAN is needed
   char UPLO = 'U'; //update LOWER triagular part for symmetric matrix 'this'
   //trans=1 -> this += U'*U -> tell BLAS to do U*U'
   //trans=0 -> this += U*U' -> tell BLAS to do U'*U
   char TRANS = trans == 0 ? 'T' : 'N';

   int m, k;
   U.getSize(m, k);
   int ldu = k; // change leading dim so that U in row-major  in col-major
   if (trans)
      k = m;

   int n = mStorage->n;
   int lda = n;

#ifdef DEBUG
   //TRANS = 'N', k specifies the number of columns of the matrix U
   //we pass U' instead of U, so k should be the number of rows of U
   int r,c; U.getSize(r,c);
   if(TRANS=='N') assert(k==r);
   else if(TRANS=='T') assert(k==c);
   else assert(false);
#endif


   dsyrk_(&UPLO, &TRANS, &n, &k, &beta, &U.getStorageRef().M[0][0], &ldu, &alpha, &mStorage->M[0][0], &lda);
}

int DenseSymMatrix::getNumberOfNonZeros() const {
   assert(mStorage->m == mStorage->n);
   int nnz = 0;
   for (int i = 0; i < mStorage->m; i++) {
      for (int j = i + 1; j < mStorage->n; j++)
         if (mStorage->M[i][j] != 0.0)
            nnz++;
      nnz++; //always have diags
   }
   return nnz;
}
