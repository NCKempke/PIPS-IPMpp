/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseGenMatrix.h"
#include "SparseStorage.h"
#include <cassert>
#include "SimpleVector.h"

#include "DoubleMatrixTypes.h"
#include <limits>
#include "SparseSymMatrix.h"

int SparseGenMatrix::isKindOf(int type) const {
   return type == kSparseGenMatrix || type == kGenMatrix;
}

SparseGenMatrix::SparseGenMatrix(int rows, int cols, int nnz) {
   mStorage = SparseStorageHandle(new SparseStorage(rows, cols, nnz));
}


SparseGenMatrix::SparseGenMatrix(int rows, int cols, int nnz, int krowM[], int jcolM[], double M[], int deleteElts) {
   mStorage = SparseStorageHandle(new SparseStorage(rows, cols, nnz, krowM, jcolM, M, deleteElts));
}

SparseGenMatrix::SparseGenMatrix(SparseStorage* m_storage) {
   mStorage = SparseStorageHandle(m_storage);
}

SparseGenMatrix::~SparseGenMatrix() {
   delete m_Mt;
   delete mStorageDynamic;
}

/* create a matrix with the same amount of columns but no rows in it */
GenMatrix* SparseGenMatrix::cloneEmptyRows(bool switchToDynamicStorage) const {
   SparseGenMatrix* clone;

   if (switchToDynamicStorage) {
      clone = new SparseGenMatrix();
      clone->mStorageDynamic = new SparseStorageDynamic(0, mStorage->n, 0);
      assert(!m_Mt);
   }
   else
      clone = new SparseGenMatrix(0, mStorage->n, 0);

   if (m_Mt) {
      assert(clone->m_Mt == NULL);
      clone->m_Mt = new SparseGenMatrix(0, m_Mt->getStorageRef().length(), 0);
   }

   return clone;
}

/* same as clone empty rows but transposes the matrix first */
SparseGenMatrix* SparseGenMatrix::cloneEmptyRowsTransposed(bool switchToDynamicStorage) const {
   SparseGenMatrix* clone;
   if (switchToDynamicStorage) {
      clone = new SparseGenMatrix();
      clone->mStorageDynamic = new SparseStorageDynamic(0, mStorage->m, 0);
      assert(!m_Mt);
   }
   else
      clone = new SparseGenMatrix(0, mStorage->m, 0);

   if (m_Mt) {
      assert(clone->m_Mt == nullptr);
      clone->m_Mt = new SparseGenMatrix(0, m_Mt->getStorageRef().m, 0);
   }

   return clone;
}

GenMatrix* SparseGenMatrix::cloneFull(bool switchToDynamicStorage) const {
   SparseGenMatrix* clone;

   if (switchToDynamicStorage) {
      clone = new SparseGenMatrix();
      clone->mStorageDynamic = new SparseStorageDynamic(*mStorage);
      assert(!m_Mt);
   }
   else {
      clone = new SparseGenMatrix(mStorage->m, mStorage->n, mStorage->len);
      mStorage->copyFrom(clone->krowM(), clone->jcolM(), clone->M());
   }

   if (m_Mt) {
      assert(clone->m_Mt == nullptr);

      SparseStorage& storage_t = m_Mt->getStorageRef();
      clone->m_Mt = new SparseGenMatrix(storage_t.m, storage_t.length(), storage_t.len);

      SparseGenMatrix* clone_t = clone->m_Mt;

      storage_t.copyFrom(clone_t->krowM(), clone_t->jcolM(), clone_t->M());
   }

   return clone;
}


void SparseGenMatrix::atPutDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) {
   mStorage->atPutDense(row, col, A, lda, rowExtent, colExtent);
   assert(m_Mt == nullptr);
}


void SparseGenMatrix::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) {
   mStorage->fromGetDense(row, col, A, lda, rowExtent, colExtent);
}


void SparseGenMatrix::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) {
   mStorage->fromGetSpRow(row, col, A, lenA, jcolA, nnz, colExtent, info);
}


void SparseGenMatrix::putSparseTriple(int irow[], int len, int jcol[], double A[], int& info) {
   mStorage->putSparseTriple(irow, len, jcol, A, info);
   assert(m_Mt == nullptr);
}


void SparseGenMatrix::writeToStream(std::ostream& out) const {
   mStorage->writeToStream(out);
}

void SparseGenMatrix::writeToStreamDense(std::ostream& out) const {
   assert(mStorage);
   if (mStorageDynamic != nullptr)
      mStorageDynamic->writeToStreamDense(out);
   else
      mStorage->writeToStreamDense(out);

}

void SparseGenMatrix::writeToStreamDenseRow(std::ostream& out, int rowidx) const {
   if (mStorageDynamic != nullptr) {
      if (mStorageDynamic->getN() > 0) {
         assert(rowidx < mStorageDynamic->getM());
         mStorageDynamic->writeToStreamDenseRow(out, rowidx);
      }
   }
   else {
      if (mStorage->n > 0) {
         assert(rowidx < mStorage->m);
         mStorage->writeToStreamDenseRow(out, rowidx);
      }
   }
}

void SparseGenMatrix::writeDashedLineToStream(std::ostream& out) const {
   if (mStorageDynamic != nullptr) {
      for (int i = 0; i < mStorageDynamic->getN(); ++i)
         out << "-\t";
   }
   else {
      for (int i = 0; i < mStorage->n; ++i)
         out << "-\t";
   }
}

void SparseGenMatrix::getDiagonal(Vector<double>& vec) {
   mStorage->getDiagonal(vec);
}


void SparseGenMatrix::setToDiagonal(const Vector<double>& vec) {
   mStorage->setToDiagonal(vec);
   assert(m_Mt == nullptr);
}


void SparseGenMatrix::atPutSpRow(int row, double A[], int lenA, int jcolA[], int& info) {
   mStorage->atPutSpRow(row, A, lenA, jcolA, info);
   assert(m_Mt == nullptr);
}


int SparseGenMatrix::numberOfNonZeros() const {
   return mStorage->numberOfNonZeros();
}


void SparseGenMatrix::symmetrize(int& info) {
   mStorage->symmetrize(info);
   assert(m_Mt == nullptr);
}


void SparseGenMatrix::getSize(long long& m, long long& n) const {
   if (mStorageDynamic != nullptr) {
      m = mStorageDynamic->getM();
      n = mStorageDynamic->getN();
   }
   else {
      if (mStorage) {
         m = mStorage->m;
         n = mStorage->n;
      }
      else {
         m = -1;
         n = -1;
      }
   }
}
void SparseGenMatrix::getSize(int& m, int& n) const {
   if (mStorageDynamic != nullptr) {
      m = mStorageDynamic->getM();
      n = mStorageDynamic->getN();
   }
   else {
      if (mStorage) {
         m = mStorage->m;
         n = mStorage->n;
      }
      else {
         m = -1;
         n = -1;
      }
   }
}


void SparseGenMatrix::atPutSubmatrix(int destRow, int destCol, DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) {
   int i, k;
   int info, nnz;

   int* ja = new int[colExtent];
   double* a = new double[colExtent];

   nnz = 0;
   for (i = 0; i < rowExtent; i++) {
      M.fromGetSpRow(srcRow + i, srcCol, a, colExtent, ja, nnz, colExtent, info);
      for (k = 0; k < nnz; k++) {
         ja[k] += (destCol - srcCol);
      }
      mStorage->atPutSpRow(destRow + i, a, nnz, ja, info);
      assert(m_Mt == nullptr);
   }

   delete[] ja;
   delete[] a;
}


void SparseGenMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const SimpleVector<double>& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   SimpleVector<double>& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() == mStorage->n && y.length() == mStorage->m);

   const double* xv = 0;
   double* yv = 0;

   if (x.length() > 0)
      xv = &x[0];
   if (y.length() > 0)
      yv = &y[0];

   mStorage->mult(beta, yv, 1, alpha, xv, 1);
}

void SparseGenMatrix::mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   mStorage->mult(beta, y, incy, alpha, x, incx);
}

void SparseGenMatrix::multMatSymUpper(double beta, SymMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const {
   SparseSymMatrix& y_sparse = dynamic_cast<SparseSymMatrix&>(y);
   assert(!y_sparse.isLower);

   mStorage->multMatSymUpper(beta, y_sparse.getStorageRef(), alpha, x, yrowstart, ycolstart);
}

void SparseGenMatrix::transmultMatSymUpper(double beta, SymMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const {
   assert(m_Mt);

   SparseSymMatrix& y_sparse = dynamic_cast<SparseSymMatrix&>(y);
   assert(!y_sparse.isLower);

   m_Mt->getStorageRef().multMatSymUpper(beta, y_sparse.getStorageRef(), alpha, x, yrowstart, ycolstart);
}

void SparseGenMatrix::transMult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const SimpleVector<double>& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   SimpleVector<double>& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() == mStorage->m && y.length() == mStorage->n);

   mStorage->transMult(beta, y.elements(), 1, alpha, x.elements(), 1);
}

void SparseGenMatrix::transMult(double beta, Vector<double>& y_in, int incy, double alpha, const Vector<double>& x_in, int incx) const {
   const SimpleVector<double>& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   SimpleVector<double>& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() > 0 && y.length() > 0);
   assert(x.length() >= incx * mStorage->m);
   assert(y.length() >= incy * mStorage->n);

   mStorage->transMult(beta, y.elements(), incy, alpha, x.elements(), incx);
}

void SparseGenMatrix::transMult(double beta, double yv[], int incy, double alpha, const double xv[], int incx) const {
   mStorage->transMult(beta, yv, incy, alpha, xv, incx);
}

void SparseGenMatrix::transMultD(double beta, Vector<double>& y_, double alpha, const Vector<double>& x_, const Vector<double>& d_) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_);
   const auto& d = dynamic_cast<const SimpleVector<double>&>(d_);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_);

   assert(x.length() == d.length());
   if (x.length() == 0 || y.length() == 0)
      return;
   assert(x.length() > 0 && y.length() > 0);
   assert(x.length() == mStorage->m);
   assert(y.length() == mStorage->n);

   mStorage->transMultD(beta, y.elements(), 1, alpha, x.elements(), d.elements(), 1);
}

double SparseGenMatrix::abmaxnorm() const {
   if (mStorage.notNil())
      return mStorage->abmaxnorm();
   else
      return 0.0;
}

double SparseGenMatrix::abminnormNonZero(double tol) const {
   if (mStorage.notNil())
      return mStorage->abminnormNonZero(tol);
   else
      return std::numeric_limits<double>::infinity();
}

void SparseGenMatrix::atPutDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const SimpleVector<double>&>(vvec);

   mStorage->atPutDiagonal(idiag, &v[0], 1, v.length());

   assert(m_Mt == nullptr);
}

void SparseGenMatrix::atAddDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const SimpleVector<double>&>(vvec);
   mStorage->atAddDiagonal(idiag, &v[0], 1, v.length());

   assert(m_Mt == nullptr);
}

void SparseGenMatrix::fromGetDiagonal(int idiag, Vector<double>& vvec) {
   mStorage->fromGetDiagonal(idiag, vvec);
}

void SparseGenMatrix::columnScale(const Vector<double>& vec) {
   mStorage->columnScale(vec);

   if (m_Mt != nullptr)
      m_Mt->rowScale(vec);
}

void SparseGenMatrix::symmetricScale(const Vector<double>& vec) {
   mStorage->symmetricScale(vec);

   if (m_Mt != nullptr)
      m_Mt->symmetricScale(vec);
}

void SparseGenMatrix::rowScale(const Vector<double>& vec) {
   mStorage->rowScale(vec);

   if (m_Mt != nullptr)
      m_Mt->columnScale(vec);
}

void SparseGenMatrix::scalarMult(double num) {
   mStorage->scalarMult(num);

   if (m_Mt != nullptr)
      m_Mt->scalarMult(num);
}

void SparseGenMatrix::matTransDMultMat(Vector<double>& d_, SymMatrix** res) {
   SimpleVector<double>& d = dynamic_cast<SimpleVector<double>&>(d_);

   int m = mStorage->m;
   int n = mStorage->n;
   int nnz = mStorage->numberOfNonZeros();

   if (*res == nullptr) {
      assert(m_Mt == nullptr);
      //we need to form the transpose
      m_Mt = new SparseGenMatrix(n, m, nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());

      //find the sparsity pattern of the product -> the buffers for result will be allocated
      int* krowMtM = nullptr;
      int* jcolMtM = nullptr;
      double* dMtM = nullptr;
      mStorage->matTransDSymbMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), &krowMtM, &jcolMtM, &dMtM);

      *res = new SparseSymMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
   }

   assert(res);
   assert(m_Mt);

   SparseSymMatrix* MtDM = dynamic_cast<SparseSymMatrix*>(*res);

   mStorage->matTransDMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}

void SparseGenMatrix::initTransposed(bool dynamic) {
   assert(m_Mt == nullptr);

   if (dynamic) {
      assert(mStorageDynamic != nullptr);
      m_Mt = new SparseGenMatrix();

      m_Mt->mStorageDynamic = mStorageDynamic->getTranspose();
   }
   else {
      const int m = mStorage->m;
      const int n = mStorage->n;
      const int nnz = mStorage->numberOfNonZeros();

      m_Mt = new SparseGenMatrix(n, m, nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
   }
}

void SparseGenMatrix::matTransDinvMultMat(Vector<double>& d_, SymMatrix** res) {
   SimpleVector<double>& d = dynamic_cast<SimpleVector<double>&>(d_);
   const int m = mStorage->m;
   const int n = mStorage->n;
   const int nnz = mStorage->numberOfNonZeros();

   if (!*res) {

      // we need to form the transpose
      if (!m_Mt) {
         m_Mt = new SparseGenMatrix(n, m, nnz);
         mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
      }

      //find the sparsity pattern of the product -> the buffers for result will be allocated
      int* krowMtM = nullptr;
      int* jcolMtM = nullptr;
      double* dMtM = nullptr;
      mStorage->matTransDSymbMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), &krowMtM, &jcolMtM, &dMtM);

      *res = new SparseSymMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
   }

   assert(res);
   assert(m_Mt);

   SparseSymMatrix* MtDM = dynamic_cast<SparseSymMatrix*>(*res);

   mStorage->matTransDinvMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}
void SparseGenMatrix::matMultTrans(SymMatrix** res) {
   int m = mStorage->m;
   int n = mStorage->n;
   int nnz = mStorage->numberOfNonZeros();

   SimpleVector<double> d(n);
   d.setToConstant(1.0);

   if (*res == nullptr) {
      //we need to form the transpose
      if (!m_Mt) {
         m_Mt = new SparseGenMatrix(n, m, nnz);
         mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
      }
      //find the sparsity pattern of the product -> the buffers for result will be allocated
      int* krowMtM = nullptr;
      int* jcolMtM = nullptr;
      double* dMtM = nullptr;

      m_Mt->mStorage->matTransDSymbMultMat(&d[0], krowM(), jcolM(), M(), &krowMtM, &jcolMtM, &dMtM);
      *res = new SparseSymMatrix(m, krowMtM[m], krowMtM, jcolMtM, dMtM, 1);
   }
   SparseSymMatrix* MMt = dynamic_cast<SparseSymMatrix*>(*res);
   m_Mt->mStorage->matTransDMultMat(&d[0], krowM(), jcolM(), M(), MMt->krowM(), MMt->jcolM(), MMt->M());
}


void SparseGenMatrix::addNnzPerRow(Vector<int>& nnzVec) const {
   SimpleVector<int>& vec = dynamic_cast<SimpleVector<int>&>(nnzVec);

   if (mStorageDynamic != nullptr) {
      assert(vec.length() == mStorageDynamic->getM());
      mStorageDynamic->addNnzPerRow(vec.elements());
   }
   else {
      assert(vec.length() == mStorage->m);
      mStorage->addNnzPerRow(vec.elements());
   }
}

void SparseGenMatrix::addNnzPerCol(Vector<int>& nnzVec) {
   SimpleVector<int>& vec = dynamic_cast<SimpleVector<int>&>(nnzVec);

   if (!m_Mt)
      initTransposed();

   if (m_Mt->mStorageDynamic != nullptr) {
      assert(vec.length() == m_Mt->mStorageDynamic->getM());
      m_Mt->mStorageDynamic->addNnzPerRow(vec.elements());
   }
   else {
      assert(vec.length() == m_Mt->mStorage->m);
      m_Mt->mStorage->addNnzPerRow(vec.elements());
   }
}

void SparseGenMatrix::addNnzPerCol(Vector<int>& nnzVec, int begin_cols, int end_cols) {
   assert(0 <= begin_cols && begin_cols <= end_cols);
   assert(end_cols - begin_cols <= nnzVec.length());

   SimpleVector<int>& vec = dynamic_cast<SimpleVector<int>&>(nnzVec);

   if (!m_Mt)
      initTransposed();

   if (m_Mt->mStorageDynamic != nullptr) {
      assert(end_cols <= m_Mt->mStorageDynamic->getM());
      m_Mt->mStorageDynamic->addNnzPerRow(vec.elements(), begin_cols, end_cols);
   }
   else {
      assert(end_cols <= m_Mt->mStorage->m);
      m_Mt->mStorage->addNnzPerRow(vec.elements(), begin_cols, end_cols);
   }
}

void SparseGenMatrix::addRowSums(Vector<double>& sumVec) const {
   SimpleVector<double>& vec = dynamic_cast<SimpleVector<double>&>(sumVec);

   assert(mStorageDynamic == nullptr && mStorage != nullptr);
   assert(vec.length() == mStorage->m);

   mStorage->addRowSums(vec.elements());
}

void SparseGenMatrix::addColSums(Vector<double>& sumVec) const {
   assert(m_Mt);
   assert(m_Mt->mStorageDynamic == nullptr && m_Mt->mStorage != nullptr);

   SimpleVector<double>& vec = dynamic_cast<SimpleVector<double>&>(sumVec);
   assert(vec.length() == m_Mt->mStorage->m);

   m_Mt->mStorage->addRowSums(vec.elements());
}

void SparseGenMatrix::getMinMaxVec(bool getMin, bool initializeVec, const SparseStorageDynamic* storage_dynamic, const Vector<double>* coScaleVec,
      Vector<double>& minmaxVec) {
   SimpleVector<double>& mvec = dynamic_cast<SimpleVector<double>&>(minmaxVec);

   assert(mvec.length() == storage_dynamic->getM());

   if (initializeVec) {
      if (getMin)
         mvec.setToConstant(std::numeric_limits<double>::max());
      else
         mvec.setToConstant(0.0);
   }

   if (coScaleVec) {
      const SimpleVector<double>* covec = dynamic_cast<const SimpleVector<double>*>(coScaleVec);

      storage_dynamic->getRowMinMaxVec(getMin, covec->elements(), mvec.elements());
   }
   else {
      storage_dynamic->getRowMinMaxVec(getMin, nullptr, mvec.elements());
   }
}

void SparseGenMatrix::getMinMaxVec(bool getMin, bool initializeVec, const SparseStorage* storage, const Vector<double>* coScaleVec,
      Vector<double>& minmaxVec) {
   SimpleVector<double>& mvec = dynamic_cast<SimpleVector<double>&>(minmaxVec);

   assert(mvec.length() == storage->m);
   if (initializeVec) {
      if (getMin)
         mvec.setToConstant(std::numeric_limits<double>::max());
      else
         mvec.setToZero();
   }

   if (coScaleVec) {
      const SimpleVector<double>* covec = dynamic_cast<const SimpleVector<double>*>(coScaleVec);

      storage->getRowMinMaxVec(getMin, covec->elements(), mvec.elements());
   }
   else {
      storage->getRowMinMaxVec(getMin, nullptr, mvec.elements());
   }
}

void SparseGenMatrix::getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) {
   if (hasDynamicStorage())
      getMinMaxVec(getMin, initializeVec, mStorageDynamic, colScaleVec, minmaxVec);
   else
      getMinMaxVec(getMin, initializeVec, mStorage, colScaleVec, minmaxVec);
}

void SparseGenMatrix::getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) {
   if (hasDynamicStorage()) {
      assert(getStorageDynamic());
      assert(getStorageDynamicTransposed());

      getMinMaxVec(getMin, initializeVec, getStorageDynamicTransposed(), rowScaleVec, minmaxVec);
   }
   else {
      // we may need to form the transpose
      if (!m_Mt)
         initTransposed();

      getMinMaxVec(getMin, initializeVec, m_Mt->mStorage, rowScaleVec, minmaxVec);
   }
}

void SparseGenMatrix::initStaticStorageFromDynamic(const Vector<int>& rowNnzVec, const Vector<int>* colNnzVec) {
   assert(mStorageDynamic != nullptr);

   const SimpleVector<int>& rowNnzVecSimple = dynamic_cast<const SimpleVector<int>&>(rowNnzVec);
   const SimpleVector<int>* colNnzVecSimple = dynamic_cast<const SimpleVector<int>*>(colNnzVec);

   mStorageDynamic->restoreOrder();
   SparseStorageHandle staticStorage(
         mStorageDynamic->getStaticStorage(rowNnzVecSimple.elements(), (colNnzVecSimple == nullptr) ? nullptr : colNnzVecSimple->elements()));

   mStorage = staticStorage;

   assert(mStorage->refs() == 2);
}

void SparseGenMatrix::deleteEmptyRowsCols(const Vector<int>& rowNnzVec, const Vector<int>& colNnzVec) {
   const SimpleVector<int>& rowNnzVecSimple = dynamic_cast<const SimpleVector<int>&>(rowNnzVec);
   const SimpleVector<int>& colNnzVecSimple = dynamic_cast<const SimpleVector<int>&>(colNnzVec);

   mStorage->deleteEmptyRowsCols(rowNnzVecSimple.elements(), colNnzVecSimple.elements());
}

void SparseGenMatrix::fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset, double* rowsArrayDense,
      int* rowSparsity) const {

   mStorage->fromGetRowsBlock(rowIndices, nRows, arrayLineSize, arrayLineOffset, rowsArrayDense, rowSparsity);
}

void SparseGenMatrix::fromGetRowsBlock(int row_start, int n_rows, int array_line_size, int array_line_offset, double* rows_array_dense,
      int* row_sparsity) const {
   assert(0 <= row_start && 0 <= n_rows && 0 <= array_line_size && 0 <= array_line_offset);
   mStorage->fromGetRowsBlock(rows_array_dense, row_start, n_rows, array_line_size, array_line_offset, row_sparsity);
}

void SparseGenMatrix::deleteEmptyRows(int*& orgIndex) {
   mStorage->deleteEmptyRows(orgIndex);

   if (m_Mt)
      this->updateTransposed();
}

void SparseGenMatrix::fromGetColsBlock(const int* colIndices, int nCols, int arrayLineSize, int arrayLineOffset, double* colsArrayDense,
      int* rowSparsity) {
   if (!m_Mt)
      updateTransposed();

   m_Mt->getStorageRef().fromGetRowsBlock(colIndices, nCols, arrayLineSize, arrayLineOffset, colsArrayDense, rowSparsity);
}

void SparseGenMatrix::fromGetColsBlock(int col_start, int n_cols, int array_line_size, int array_line_offset, double* cols_array_dense,
      int* row_sparsity) {
   if (!m_Mt)
      updateTransposed();

   m_Mt->getStorageRef().fromGetRowsBlock(cols_array_dense, col_start, n_cols, array_line_size, array_line_offset, row_sparsity);
}

bool SparseGenMatrix::hasTransposed() const {
   return (m_Mt != nullptr);
}

void SparseGenMatrix::freeDynamicStorage() {
   delete mStorageDynamic;
   mStorageDynamic = nullptr;
}

void SparseGenMatrix::updateTransposed() {
   if (m_Mt) {
      delete m_Mt;
      m_Mt = nullptr;
   }
   const int nnz = mStorage->numberOfNonZeros();

   const int m = mStorage->m;
   const int n = mStorage->n;
   m_Mt = new SparseGenMatrix(n, m, nnz);

   mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
}

void SparseGenMatrix::deleteTransposed() {
   if (m_Mt) {
      delete m_Mt;
      m_Mt = nullptr;
   }
}

void SparseGenMatrix::getLinkVarsNnz(std::vector<int>& vec) const {
   mStorage->getLinkVarsNnz(vec);
}

void SparseGenMatrix::updateNonEmptyRowsCount(std::vector<int>& rowcount) const {
   const int m = mStorage->m;
   const int* const rowStart = mStorage->krowM;

   assert(unsigned(m) == rowcount.size());

   for (int i = 0; i < m; i++)
      if (rowStart[i] != rowStart[i + 1])
         rowcount[i]++;
}

void SparseGenMatrix::updateNonEmptyRowsCountNew(int block_id, std::vector<int>& n_blocks_per_row, std::vector<int>& row_start_block,
      std::vector<int>& row_end_block) const {
   const int m = mStorage->m;
   const int* const rowStart = mStorage->krowM;

   assert(block_id >= 0 && m >= 0);
   assert(unsigned(m) == n_blocks_per_row.size());
   assert(n_blocks_per_row.size() == row_start_block.size() && n_blocks_per_row.size() == row_end_block.size());

   for (int i = 0; i < m; i++)
      if (rowStart[i] != rowStart[i + 1]) {
         row_start_block[i] = std::min(row_start_block[i], block_id);
         row_end_block[i] = std::max(row_end_block[i], block_id);

         n_blocks_per_row[i]++;
      }
}

void SparseGenMatrix::updateNonEmptyRowsCount(int block_id, std::vector<int>& n_blocks_per_row, std::vector<int>& row_start_block,
      std::vector<int>& row_end_block) const {
   const int m = mStorage->m;
   const int* const rowStart = mStorage->krowM;

   assert(block_id >= 0 && m >= 0);
   assert(unsigned(m) == n_blocks_per_row.size());
   assert(n_blocks_per_row.size() == row_start_block.size() && n_blocks_per_row.size() == row_end_block.size());

   for (int i = 0; i < m; i++)
      if (rowStart[i] != rowStart[i + 1]) {
         if (row_start_block[i] < 0)
            row_start_block[i] = block_id;
         else
            row_end_block[i] = block_id;

         n_blocks_per_row[i]++;
      }
}

SparseGenMatrix& SparseGenMatrix::getTranspose() {
   if (m_Mt)
      return *m_Mt;

   updateTransposed();

   assert(m_Mt);

   return *m_Mt;
}

void SparseGenMatrix::permuteRows(const std::vector<unsigned int>& permvec) {
   mStorage->permuteRows(permvec);

   if (m_Mt)
      m_Mt->mStorage->permuteCols(permvec);
}

void SparseGenMatrix::permuteCols(const std::vector<unsigned int>& permvec) {
   mStorage->permuteCols(permvec);

   if (m_Mt)
      m_Mt->mStorage->permuteRows(permvec);
}

int SparseGenMatrix::appendRow(const SparseGenMatrix& matrix_row, int row) {
   assert(hasDynamicStorage());
   assert(!hasTransposed());
   assert(matrix_row.hasDynamicStorage());

   mStorageDynamic->appendRow(matrix_row.getStorageDynamicRef(), row);

   return mStorageDynamic->getM() - 1;
}

int SparseGenMatrix::appendCol(const SparseGenMatrix& matrix_col, int col) {
   assert(matrix_col.hasTransposed());
   assert(matrix_col.hasDynamicStorage());
   assert(hasDynamicStorage());
   assert(!hasTransposed());

   mStorageDynamic->appendRow(matrix_col.getStorageDynamicTransposedRef(), col);

   return mStorageDynamic->getM() - 1;
}

void SparseGenMatrix::axpyWithRowAt(double alpha, SimpleVector<double>& y, int row) const {
   assert(hasDynamicStorage());

   mStorageDynamic->axpyWithRowAt(alpha, y.elements(), y.length(), row);
}

void SparseGenMatrix::axpyWithRowAtPosNeg(double alpha, SimpleVector<double>& y_pos, SimpleVector<double>& y_neg, int row) const {
   assert(hasDynamicStorage());
   assert(y_pos.length() == y_neg.length());

   mStorageDynamic->axpyWithRowAtPosNeg(alpha, y_pos.elements(), y_neg.elements(), y_pos.length(), row);
}

double SparseGenMatrix::localRowTimesVec(const SimpleVector<double>& vec, int row) const {
   assert(hasDynamicStorage());

   return mStorageDynamic->rowTimesVec(vec.elements(), vec.length(), row);
}

void SparseGenMatrix::removeRow(int row) {
   assert(hasDynamicStorage());
   if (hasTransposed()) {
      assert(m_Mt->getStorageDynamic());
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
      removeRowUsingTransposed(row, m_Mt->getStorageDynamicRef());
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
   }
   else
      mStorageDynamic->clearRow(row);
}

void SparseGenMatrix::removeRowUsingTransposed(int row, SparseStorageDynamic& mat_trans) {
   assert(hasDynamicStorage());

   SparseStorageDynamic& mat = *mStorageDynamic;

   const int start = mat.getRowPtr(row).start;
   const int end = mat.getRowPtr(row).end;

   for (int i = start; i < end; ++i) {
      const int col = mat.getJcolM(i);
      mat_trans.removeEntryAtRowCol(col, row);
   }

   mat.clearRow(row);
}

void SparseGenMatrix::removeCol(int col) {
   assert(hasDynamicStorage());
   if (hasTransposed()) {
      assert(m_Mt->getStorageDynamic());
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
      m_Mt->removeRowUsingTransposed(col, *mStorageDynamic);

      m_Mt->removeRow(col);
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
   }
   else
      mStorageDynamic->clearCol(col);
}

void SparseGenMatrix::removeEntryAtRowCol(int row, int col) {
   assert(hasDynamicStorage());

   if (hasTransposed())
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());

   mStorageDynamic->removeEntryAtRowCol(row, col);

   if (hasTransposed()) {
      m_Mt->removeEntryAtRowCol(col, row);
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
   }
}

void SparseGenMatrix::removeEntryAtRowColIndex(int row, int col_index) {
   assert(hasDynamicStorage());

   if (hasTransposed())
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());

   const int col = mStorageDynamic->getJcolM(col_index);
   mStorageDynamic->removeEntryAtIndex(row, col_index);

   if (hasTransposed()) {
      m_Mt->removeEntryAtRowCol(col, row);
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
   }
}

void SparseGenMatrix::addColToRow(double coeff, int col, int row) {
   assert(hasDynamicStorage());

   if (hasTransposed())
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());

   mStorageDynamic->addColToRow(coeff, col, row);

   if (hasTransposed()) {
      m_Mt->addColToRow(coeff, row, col);
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
   }
}

SparseGenMatrix* SparseGenMatrix::shaveLeft(int n_cols) {
   assert(!hasDynamicStorage());
   assert(n_cols <= mStorage->n);

   SparseStorage* border = mStorage->shaveLeft(n_cols);

   if (n_cols != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }

   return new SparseGenMatrix(border);
}

GenMatrix* SparseGenMatrix::shaveBottom(int n_rows) {
   assert(!hasDynamicStorage());
   assert(n_rows <= mStorage->m);

   SparseStorage* border = mStorage->shaveBottom(n_rows);

   if (n_rows != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }

   return new SparseGenMatrix(border);
}

void SparseGenMatrix::dropNEmptyRowsBottom(int n_rows) {
   assert(!hasDynamicStorage());
   assert(n_rows <= mStorage->m);

   mStorage->dropNEmptyRowsBottom(n_rows);

   if (n_rows != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }
}

void SparseGenMatrix::dropNEmptyRowsTop(int n_rows) {
   assert(!hasDynamicStorage());
   assert(n_rows <= mStorage->m);

   mStorage->dropNEmptyRowsTop(n_rows);

   if (n_rows != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }
}
