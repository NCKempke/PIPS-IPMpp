/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "SparseStorage.h"
#include <cassert>
#include <cmath>
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"

int SparseSymMatrix::isKindOf(int type) const {
   return type == kSparseSymMatrix || type == kSymMatrix;
}

SparseSymMatrix::SparseSymMatrix() : isLower(true) {
   mStorage = nullptr;
}

SparseSymMatrix::SparseSymMatrix(const SparseSymMatrix& mat) : isLower(mat.isLower) {
   mStorage = SparseStorageHandle(new SparseStorage(mat.mStorage->m, mat.mStorage->n, mat.mStorage->len));
   mat.getStorageRef().copyFrom(mStorage->krowM, mStorage->jcolM, mStorage->M);
}

SparseSymMatrix::SparseSymMatrix(int size, int nnz, bool isLower) : isLower(isLower) {
   mStorage = SparseStorageHandle(new SparseStorage(size, size, nnz));
}

SparseSymMatrix::SparseSymMatrix(SparseStorage* m_storage, bool is_lower_) : isLower(is_lower_) {
   mStorage = SparseStorageHandle(m_storage);
}

SparseSymMatrix::SparseSymMatrix(int size, int nnz, int krowM[], int jcolM[], double M[], int deleteElts, bool isLower) : isLower(isLower) {
   mStorage = SparseStorageHandle(new SparseStorage(size, size, nnz, krowM, jcolM, M, deleteElts));
}


void SparseSymMatrix::putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) {
   mStorage->putSparseTriple(irow, len, jcol, A, info);
}

SymMatrix* SparseSymMatrix::clone() const {
   return new SparseSymMatrix(*this);
}

void SparseSymMatrix::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const {
   mStorage->fromGetDense(row, col, A, lda, rowExtent, colExtent);
}


void SparseSymMatrix::getDiagonal(Vector<double>& vec) const {
   mStorage->getDiagonal(vec);
}

void SparseSymMatrix::setToDiagonal(const Vector<double>& vec) {
   mStorage->setToDiagonal(vec);
}

void SparseSymMatrix::diagonal_add_constant_from(int from, int length, double value) {
   assert(0 <= from);
   assert(from + length <= this->size());
   mStorage->diagonal_add_constant_from(from, length, value);
}


void SparseSymMatrix::symAtPutSpRow(int row, const double A[], int lenA, const int jcolA[], int& info) {
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

void SparseSymMatrix::symPutZeroes() {
   assert(mStorage);
   mStorage->clear();
}

void SparseSymMatrix::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const {
   mStorage->fromGetSpRow(row, col, A, lenA, jcolA, nnz, colExtent, info);
}

void SparseSymMatrix::symAtPutSubmatrix(int destRow, int destCol, const DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) {
   int i, k;
   int info, nnz;

   auto* ja = new int[colExtent];
   auto* a = new double[colExtent];

   nnz = 0;
   for (i = 0; i < rowExtent; i++) {
      M.fromGetSpRow(srcRow + i, srcCol, a, colExtent, ja, nnz, colExtent, info);
      for (k = 0; k < nnz; k++) {
         ja[k] += (destCol - srcCol);
      }
      this->symAtPutSpRow(destRow + i, a, nnz, ja, info);
      assert(info == 0);
   }

   delete[] a;
   delete[] ja;
}

// Pass these to storage
void SparseSymMatrix::getSize(long long& m, long long& n) const {
   int mint, nint;
   mStorage->getSize(mint, nint);
   m = mint;
   n = nint;
}

void SparseSymMatrix::getSize(int& m, int& n) const {
   mStorage->getSize(m, n);
}


long long SparseSymMatrix::size() const {
   return mStorage->rows();
}

void SparseSymMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() == mStorage->n && y.length() == mStorage->m);

   const double* xv = nullptr;
   double* yv = nullptr;
   if (x.length() > 0)
      xv = &x[0];
   if (y.length() > 0)
      yv = &y[0];

   this->mult(beta, yv, 1, alpha, xv, 1);
}

void SparseSymMatrix::transMult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() == mStorage->n && y.length() == mStorage->m);

   const double* xv = nullptr;
   double* yv = nullptr;
   if (x.length() > 0)
      xv = &x[0];
   if (y.length() > 0)
      yv = &y[0];

   this->mult(beta, yv, 1, alpha, xv, 1);
}

void SparseSymMatrix::transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   this->mult(beta, y, incy, alpha, x, incx);
}

double SparseSymMatrix::abmaxnorm() const {
   return mStorage->abmaxnorm();
}

double SparseSymMatrix::abminnormNonZero(double tol) const {
   if (mStorage.notNil())
      return mStorage->abminnormNonZero(tol);
   else
      return std::numeric_limits<double>::infinity();
}

void SparseSymMatrix::writeToStream(std::ostream& out) const {
   mStorage->writeToStream(out);
}

void SparseSymMatrix::writeNNZpatternToStreamDense(std::ostream& out) const {
   mStorage->writeNNZpatternToStreamDense(out);
}

void SparseSymMatrix::writeToStreamDense(std::ostream& out) const {
   mStorage->writeToStreamDense(out);
}


void SparseSymMatrix::writeToStreamDenseRow(std::ostream& out, int row) const {
   if (mStorage->n > 0) {
      assert(row < mStorage->m);
      mStorage->writeToStreamDenseRow(out, row);
   }
}

void SparseSymMatrix::mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   mStorage->multSym(beta, y, incy, alpha, x, incx);
}

void SparseSymMatrix::atPutDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atPutDiagonal(idiag, v);
}

void SparseSymMatrix::atAddDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atAddDiagonal(idiag, v);
}

void SparseSymMatrix::fromGetDiagonal(int idiag, Vector<double>& v) const {
   mStorage->fromGetDiagonal(idiag, v);
}

void SparseSymMatrix::symmetricScale(const Vector<double>& vec) {
   mStorage->symmetricScale(vec);
}

void SparseSymMatrix::columnScale(const Vector<double>& vec) {
   mStorage->columnScale(vec);
}

void SparseSymMatrix::rowScale(const Vector<double>& vec) {
   mStorage->rowScale(vec);
}

void SparseSymMatrix::scalarMult(double num) {
   mStorage->scalarMult(num);
}

void SparseSymMatrix::reduceToLower() {
   mStorage->reduceToLower();
}

void SparseSymMatrix::deleteEmptyRowsCols(const Vector<int>& nnzVec) {
   const auto& vec = dynamic_cast<const SimpleVector<int>&>(nnzVec);
#ifndef NDEBUG
   int m, n;
   mStorage->getSize(m, n);
   assert(nnzVec.length() == m);
   assert(nnzVec.length() == n);
#endif
   mStorage->deleteEmptyRowsCols(vec.elements(), vec.elements());
}

void SparseSymMatrix::getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const {
   mStorage->getSparseTriplet_c2fortran(irn, jcn, val);
}

void SparseSymMatrix::getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const {
   mStorage->getSparseTriplet_fortran2fortran(irn, jcn, val);
}


void SparseSymMatrix::deleteZeroRowsCols(int*& new2orgIdx) {
   mStorage->deleteZeroRowsColsSym(new2orgIdx);
}

SparseGenMatrix* SparseSymMatrix::shaveSymLeftBottom(int n_vars) {
   assert(n_vars <= mStorage->n);
   assert(n_vars >= 0);
   //   SparseStorage* border = mStorage->shaveSymLeftBottom( n_vars );
   // TODO : not implemented properly ..
   assert(mStorage->len == 0);
   assert(mStorage->m == mStorage->n);

   mStorage->m -= n_vars;
   mStorage->n -= n_vars;

   auto* m_border = new SparseStorage(mStorage->m, n_vars, 0);
   auto* border = new SparseGenMatrix(m_border);

   return border;
}
