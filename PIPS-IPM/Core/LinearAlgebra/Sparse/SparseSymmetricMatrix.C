/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseSymmetricMatrix.h"
#include "SparseMatrix.h"
#include "SparseStorage.h"
#include <cassert>
#include <cmath>
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"

int SparseSymmetricMatrix::is_a(int type) const {
   return type == kSparseSymMatrix || type == kSymMatrix;
}

SparseSymmetricMatrix::SparseSymmetricMatrix() : isLower(true) {
   mStorage = nullptr;
}

// TODO implement mStorage copy ctor..
SparseSymmetricMatrix::SparseSymmetricMatrix(const SparseSymmetricMatrix& mat) : mStorage{
   std::make_unique<SparseStorage>(mat.mStorage->m, mat.mStorage->n, mat.mStorage->len)}, isLower(mat.isLower) {
   mat.getStorage().copyFrom(mStorage->krowM, mStorage->jcolM, mStorage->M);
}

SparseSymmetricMatrix::SparseSymmetricMatrix(int size, int nnz, bool isLower) : mStorage{
   std::make_unique<SparseStorage>(size, size, nnz)}, isLower(isLower) {}

SparseSymmetricMatrix::SparseSymmetricMatrix(std::unique_ptr<SparseStorage> m_storage, bool is_lower_) :
   mStorage{std::move(m_storage)}, isLower(is_lower_) {
}

SparseSymmetricMatrix::SparseSymmetricMatrix(int size, int nnz, int krowM[], int jcolM[], double M[], int deleteElts,
   bool isLower) : mStorage{std::make_unique<SparseStorage>(size, size, nnz, krowM, jcolM, M, deleteElts)}, isLower(isLower) {}


void SparseSymmetricMatrix::putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) {
   mStorage->putSparseTriple(irow, len, jcol, A, info);
}

std::unique_ptr<SymmetricMatrix> SparseSymmetricMatrix::clone() const {
   return std::make_unique<SparseSymmetricMatrix>(*this);
}

void SparseSymmetricMatrix::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const {
   mStorage->fromGetDense(row, col, A, lda, rowExtent, colExtent);
}


void SparseSymmetricMatrix::getDiagonal(Vector<double>& vec) const {
   mStorage->getDiagonal(vec);
}

void SparseSymmetricMatrix::setToDiagonal(const Vector<double>& vec) {
   mStorage->setToDiagonal(vec);
}

void SparseSymmetricMatrix::diagonal_add_constant_from(int from, int length, double value) {
   assert(0 <= from);
   assert(from + length <= this->size());
   mStorage->diagonal_add_constant_from(from, length, value);
}

void SparseSymmetricMatrix::diagonal_set_to_constant_from(int from, int length, double value) {
   assert(0 <= from);
   assert(from + length <= this->size());
   mStorage->diagonal_set_to_constant_from(from, length, value);
}

void SparseSymmetricMatrix::symAtPutSpRow(int row, const double A[], int lenA, const int jcolA[], int& info) {
   // Lower triangular put
   int lA = lenA;
   while (lA > 0 && jcolA[lA - 1] > row)
      lA--;
   if (lA > 0) {
      mStorage->atPutSpRow(row, A, lA, jcolA, info);
   } else {
      info = 0;
   }
}

void SparseSymmetricMatrix::symPutZeroes() {
   assert(mStorage);
   mStorage->clear();
}

void SparseSymmetricMatrix::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent,
   int& info) const {
   mStorage->fromGetSpRow(row, col, A, lenA, jcolA, nnz, colExtent, info);
}

void SparseSymmetricMatrix::symAtPutSubmatrix(int destRow, int destCol, const AbstractMatrix& M, int srcRow, int srcCol,
   int rowExtent, int colExtent) {
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

std::pair<long long, long long> SparseSymmetricMatrix::n_rows_columns() const {
   return mStorage->n_rows_columns();
}

long long SparseSymmetricMatrix::size() const {
   return mStorage->n_rows();
}

long long SparseSymmetricMatrix::n_rows() const {
   return mStorage->n_rows();
}

long long SparseSymmetricMatrix::n_columns() const {
   return mStorage->n_columns();
}

void SparseSymmetricMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
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

void
SparseSymmetricMatrix::transMult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
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

void
SparseSymmetricMatrix::transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   this->mult(beta, y, incy, alpha, x, incx);
}

double SparseSymmetricMatrix::inf_norm() const {
   return mStorage->inf_norm();
}

double SparseSymmetricMatrix::abminnormNonZero(double tol) const {
   if (mStorage)
      return mStorage->abminnormNonZero(tol);
   else
      return std::numeric_limits<double>::infinity();
}

void SparseSymmetricMatrix::writeToStream(std::ostream& out) const {
   mStorage->writeToStream(out);
}

void SparseSymmetricMatrix::writeNNZpatternToStreamDense(std::ostream& out) const {
   mStorage->writeNNZpatternToStreamDense(out);
}

void SparseSymmetricMatrix::writeToStreamDense(std::ostream& out) const {
   mStorage->writeToStreamDense(out);
}


void SparseSymmetricMatrix::writeToStreamDenseRow(std::ostream& out, int row) const {
   if (mStorage->n > 0) {
      assert(row < mStorage->m);
      mStorage->writeToStreamDenseRow(out, row);
   }
}

void SparseSymmetricMatrix::mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   mStorage->multSym(beta, y, incy, alpha, x, incx);
}

void SparseSymmetricMatrix::atPutDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atPutDiagonal(idiag, v);
}

void SparseSymmetricMatrix::atAddDiagonal(int idiag, const Vector<double>& v) {
   mStorage->atAddDiagonal(idiag, v);
}

void SparseSymmetricMatrix::fromGetDiagonal(int idiag, Vector<double>& v) const {
   mStorage->fromGetDiagonal(idiag, v);
}

void SparseSymmetricMatrix::symmetricScale(const Vector<double>& vec) {
   mStorage->symmetricScale(vec);
}

void SparseSymmetricMatrix::columnScale(const Vector<double>& vec) {
   mStorage->columnScale(vec);
}

void SparseSymmetricMatrix::rowScale(const Vector<double>& vec) {
   mStorage->rowScale(vec);
}

void SparseSymmetricMatrix::scalarMult(double num) {
   mStorage->scalarMult(num);
}

void SparseSymmetricMatrix::reduceToLower() {
   mStorage->reduceToLower();
}

void SparseSymmetricMatrix::deleteEmptyRowsCols(const Vector<int>& nnzVec) {
   const auto& vec = dynamic_cast<const SimpleVector<int>&>(nnzVec);
#ifndef NDEBUG
   const auto[m, n] = mStorage->n_rows_columns();
   assert(nnzVec.length() == m);
   assert(nnzVec.length() == n);
#endif
   mStorage->deleteEmptyRowsCols(vec.elements(), vec.elements());
}

void SparseSymmetricMatrix::getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const {
   mStorage->getSparseTriplet_c2fortran(irn, jcn, val);
}

void SparseSymmetricMatrix::getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const {
   mStorage->getSparseTriplet_fortran2fortran(irn, jcn, val);
}


void SparseSymmetricMatrix::deleteZeroRowsCols(int*& new2orgIdx) {
   mStorage->deleteZeroRowsColsSym(new2orgIdx);
}

std::unique_ptr<SparseMatrix> SparseSymmetricMatrix::shaveSymLeftBottom(int n_vars) {
   assert(n_vars <= mStorage->n);
   assert(n_vars >= 0);
   //   SparseStorage* border = mStorage->shaveSymLeftBottom( n_vars );
   // TODO : not implemented properly ..
   assert(mStorage->len == 0);
   assert(mStorage->m == mStorage->n);

   mStorage->m -= n_vars;
   mStorage->n -= n_vars;

   auto m_border = std::make_unique<SparseStorage>(mStorage->m, n_vars, 0);
   auto border = std::make_unique<SparseMatrix>(std::move(m_border));

   return border;
}

void SparseSymmetricMatrix::append_empty_diagonal(int n_values) {
   const int old_m = mStorage->m;

   mStorage->m += n_values;
   mStorage->n += n_values;

   int* new_krowM = new int[mStorage->m + 1];
   std::copy(mStorage->krowM, mStorage->krowM + old_m + 1, new_krowM);
   std::fill(new_krowM + old_m + 1, new_krowM + mStorage->m + 2, mStorage->krowM[mStorage->m]);
   std::swap(new_krowM, mStorage->krowM);
   delete[] new_krowM;
}

