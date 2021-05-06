/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "SparseMatrix.h"
#include "SparseStorage.h"
#include <cassert>
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"
#include <limits>
#include "SparseSymmetricMatrix.h"

int SparseMatrix::is_a(int type) const {
   return type == kSparseGenMatrix || type == kGenMatrix;
}

SparseMatrix::SparseMatrix(int rows, int cols, int nnz) {
   mStorage = SmartPointer<SparseStorage>(new SparseStorage(rows, cols, nnz));
}


SparseMatrix::SparseMatrix(int rows, int cols, int nnz, int krowM[], int jcolM[], double M[], int deleteElts) {
   mStorage = SmartPointer<SparseStorage>(new SparseStorage(rows, cols, nnz, krowM, jcolM, M, deleteElts));
}

SparseMatrix::SparseMatrix(SparseStorage* m_storage) {
   mStorage = SmartPointer<SparseStorage>(m_storage);
}

SparseMatrix::~SparseMatrix() {
   delete m_Mt;
   delete mStorageDynamic;
}

/* create a matrix with the same amount of columns but no rows in it */
GeneralMatrix* SparseMatrix::cloneEmptyRows(bool switchToDynamicStorage) const {
   SparseMatrix* clone;

   if (switchToDynamicStorage) {
      clone = new SparseMatrix();
      clone->mStorageDynamic = new SparseStorageDynamic(0, mStorage->n, 0);
      assert(!m_Mt);
   }
   else
      clone = new SparseMatrix(0, mStorage->n, 0);

   if (m_Mt) {
      assert(!clone->m_Mt);
      clone->m_Mt = new SparseMatrix(0, m_Mt->getStorageRef().length(), 0);
   }

   return clone;
}

/* same as clone empty rows but transposes the matrix first */
SparseMatrix* SparseMatrix::cloneEmptyRowsTransposed(bool switchToDynamicStorage) const {
   SparseMatrix* clone;
   if (switchToDynamicStorage) {
      clone = new SparseMatrix();
      clone->mStorageDynamic = new SparseStorageDynamic(0, mStorage->m, 0);
      assert(!m_Mt);
   }
   else
      clone = new SparseMatrix(0, mStorage->m, 0);

   if (m_Mt) {
      assert(clone->m_Mt == nullptr);
      clone->m_Mt = new SparseMatrix(0, m_Mt->getStorageRef().m, 0);
   }

   return clone;
}

GeneralMatrix* SparseMatrix::cloneFull(bool switchToDynamicStorage) const {
   SparseMatrix* clone;

   if (switchToDynamicStorage) {
      clone = new SparseMatrix();
      clone->mStorageDynamic = new SparseStorageDynamic(*mStorage);
      assert(!m_Mt);
   }
   else {
      clone = new SparseMatrix(mStorage->m, mStorage->n, mStorage->len);
      mStorage->copyFrom(clone->krowM(), clone->jcolM(), clone->M());
   }

   if (m_Mt) {
      assert(clone->m_Mt == nullptr);

      SparseStorage& storage_t = m_Mt->getStorageRef();
      clone->m_Mt = new SparseMatrix(storage_t.m, storage_t.length(), storage_t.len);

      SparseMatrix* clone_t = clone->m_Mt;

      storage_t.copyFrom(clone_t->krowM(), clone_t->jcolM(), clone_t->M());
   }

   return clone;
}


void SparseMatrix::atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) {
   mStorage->atPutDense(row, col, A, lda, rowExtent, colExtent);
   assert(m_Mt == nullptr);
}


void SparseMatrix::fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const {
   mStorage->fromGetDense(row, col, A, lda, rowExtent, colExtent);
}


void SparseMatrix::fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const {
   mStorage->fromGetSpRow(row, col, A, lenA, jcolA, nnz, colExtent, info);
}


void SparseMatrix::putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) {
   mStorage->putSparseTriple(irow, len, jcol, A, info);
   assert(m_Mt == nullptr);
}


void SparseMatrix::writeToStream(std::ostream& out) const {
   mStorage->writeToStream(out);
}

void SparseMatrix::writeToStreamDense(std::ostream& out) const {
   assert(mStorage);
   if (mStorageDynamic != nullptr)
      mStorageDynamic->writeToStreamDense(out);
   else
      mStorage->writeToStreamDense(out);

}

void SparseMatrix::writeToStreamDenseRow(std::ostream& out, int rowidx) const {
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

void SparseMatrix::writeDashedLineToStream(std::ostream& out) const {
   if (mStorageDynamic != nullptr) {
      for (int i = 0; i < mStorageDynamic->getN(); ++i)
         out << "-\t";
   }
   else {
      for (int i = 0; i < mStorage->n; ++i)
         out << "-\t";
   }
}

void SparseMatrix::getDiagonal(Vector<double>& vec) const {
   mStorage->getDiagonal(vec);
}


void SparseMatrix::setToDiagonal(const Vector<double>& vec) {
   mStorage->setToDiagonal(vec);
   assert(m_Mt == nullptr);
}


void SparseMatrix::atPutSpRow(int row, const double A[], int lenA, const int jcolA[], int& info) {
   mStorage->atPutSpRow(row, A, lenA, jcolA, info);
   assert(m_Mt == nullptr);
}


int SparseMatrix::numberOfNonZeros() const {
   return mStorage->numberOfNonZeros();
}


void SparseMatrix::symmetrize(int& info) {
   mStorage->symmetrize(info);
   assert(m_Mt == nullptr);
}


void SparseMatrix::getSize(long long& m, long long& n) const {
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
void SparseMatrix::getSize(int& m, int& n) const {
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


void SparseMatrix::atPutSubmatrix(int destRow, int destCol, const AbstractMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) {
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
      mStorage->atPutSpRow(destRow + i, a, nnz, ja, info);
      assert(m_Mt == nullptr);
   }

   delete[] ja;
   delete[] a;
}


void SparseMatrix::mult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() == mStorage->n && y.length() == mStorage->m);

   const double* xv = 0;
   double* yv = 0;

   if (x.length() > 0)
      xv = &x[0];
   if (y.length() > 0)
      yv = &y[0];

   mStorage->mult(beta, yv, 1, alpha, xv, 1);
}

void SparseMatrix::mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const {
   mStorage->mult(beta, y, incy, alpha, x, incx);
}

void SparseMatrix::multMatSymUpper(double beta, SymmetricMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const {
   auto& y_sparse = dynamic_cast<SparseSymmetricMatrix&>(y);
   assert(!y_sparse.isLower);

   mStorage->multMatSymUpper(beta, y_sparse.getStorageRef(), alpha, x, yrowstart, ycolstart);
}

void SparseMatrix::transmultMatSymUpper(double beta, SymmetricMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const {
   if(!m_Mt)
      initTransposed();
   assert(m_Mt);

   auto& y_sparse = dynamic_cast<SparseSymmetricMatrix&>(y);
   assert(!y_sparse.isLower);

   m_Mt->getStorageRef().multMatSymUpper(beta, y_sparse.getStorageRef(), alpha, x, yrowstart, ycolstart);
}

void SparseMatrix::transMult(double beta, Vector<double>& y_in, double alpha, const Vector<double>& x_in) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() == mStorage->m && y.length() == mStorage->n);

   mStorage->transMult(beta, y.elements(), 1, alpha, x.elements(), 1);
}

void SparseMatrix::transMult(double beta, Vector<double>& y_in, int incy, double alpha, const Vector<double>& x_in, int incx) const {
   const auto& x = dynamic_cast<const SimpleVector<double>&>(x_in);
   auto& y = dynamic_cast<SimpleVector<double>&>(y_in);

   assert(x.length() > 0 && y.length() > 0);
   assert(x.length() >= incx * mStorage->m);
   assert(y.length() >= incy * mStorage->n);

   mStorage->transMult(beta, y.elements(), incy, alpha, x.elements(), incx);
}

void SparseMatrix::transMult(double beta, double yv[], int incy, double alpha, const double xv[], int incx) const {
   mStorage->transMult(beta, yv, incy, alpha, xv, incx);
}

void SparseMatrix::transMultD(double beta, Vector<double>& y_, double alpha, const Vector<double>& x_, const Vector<double>& d_) const {
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

double SparseMatrix::inf_norm() const {
   if (mStorage)
      return mStorage->inf_norm();
   else
      return 0.0;
}

double SparseMatrix::abminnormNonZero(double tol) const {
   if (mStorage)
      return mStorage->abminnormNonZero(tol);
   else
      return std::numeric_limits<double>::infinity();
}

void SparseMatrix::atPutDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const SimpleVector<double>&>(vvec);

   mStorage->atPutDiagonal(idiag, &v[0], 1, v.length());

   if(m_Mt)
      m_Mt->atPutDiagonal(idiag, vvec);
}

void SparseMatrix::atAddDiagonal(int idiag, const Vector<double>& vvec) {
   const auto& v = dynamic_cast<const SimpleVector<double>&>(vvec);
   mStorage->atAddDiagonal(idiag, &v[0], 1, v.length());

   if(m_Mt)
      m_Mt->atAddDiagonal(idiag, vvec);
}

void SparseMatrix::fromGetDiagonal(int idiag, Vector<double>& vvec) const {
   mStorage->fromGetDiagonal(idiag, vvec);
}

void SparseMatrix::columnScale(const Vector<double>& vec) {
   mStorage->columnScale(vec);

   if (m_Mt)
      m_Mt->rowScale(vec);
}

void SparseMatrix::symmetricScale(const Vector<double>& vec) {
   mStorage->symmetricScale(vec);

   if (m_Mt)
      m_Mt->symmetricScale(vec);
}

void SparseMatrix::rowScale(const Vector<double>& vec) {
   mStorage->rowScale(vec);

   if (m_Mt)
      m_Mt->columnScale(vec);
}

void SparseMatrix::scalarMult(double num) {
   mStorage->scalarMult(num);

   if (m_Mt != nullptr)
      m_Mt->scalarMult(num);
}

void SparseMatrix::matTransDMultMat(const Vector<double>& d_, SymmetricMatrix** res) const {
   const auto& d = dynamic_cast<const SimpleVector<double>&>(d_);

   int m = mStorage->m;
   int n = mStorage->n;
   int nnz = mStorage->numberOfNonZeros();

   if (*res == nullptr) {
      assert(m_Mt == nullptr);
      //we need to form the transpose
      m_Mt = new SparseMatrix(n, m, nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());

      //find the sparsity pattern of the product -> the buffers for result will be allocated
      int* krowMtM = nullptr;
      int* jcolMtM = nullptr;
      double* dMtM = nullptr;
      mStorage->matTransDSymbMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), &krowMtM, &jcolMtM, &dMtM);

      *res = new SparseSymmetricMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
   }

   assert(res);
   assert(m_Mt);

   auto* MtDM = dynamic_cast<SparseSymmetricMatrix*>(*res);

   mStorage->matTransDMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}

void SparseMatrix::initTransposed(bool dynamic) const {
   if(m_Mt)
      return;

   if (dynamic) {
      assert(mStorageDynamic != nullptr);
      m_Mt = new SparseMatrix();

      m_Mt->mStorageDynamic = mStorageDynamic->getTranspose();
   }
   else {
      const int m = mStorage->m;
      const int n = mStorage->n;
      const int nnz = mStorage->numberOfNonZeros();

      m_Mt = new SparseMatrix(n, m, nnz);
      mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
   }
}

void SparseMatrix::matTransDinvMultMat(const Vector<double>& d_, SymmetricMatrix** res) const {
   const auto& d = dynamic_cast<const SimpleVector<double>&>(d_);
   const int m = mStorage->m;
   const int n = mStorage->n;
   const int nnz = mStorage->numberOfNonZeros();

   if (!*res) {

      // we need to form the transpose
      if (!m_Mt) {
         m_Mt = new SparseMatrix(n, m, nnz);
         mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
      }

      //find the sparsity pattern of the product -> the buffers for result will be allocated
      int* krowMtM = nullptr;
      int* jcolMtM = nullptr;
      double* dMtM = nullptr;
      mStorage->matTransDSymbMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), &krowMtM, &jcolMtM, &dMtM);

      *res = new SparseSymmetricMatrix(n, krowMtM[n], krowMtM, jcolMtM, dMtM, 1);
   }

   assert(res);
   assert(m_Mt);

   auto* MtDM = dynamic_cast<SparseSymmetricMatrix*>(*res);

   mStorage->matTransDinvMultMat(&d[0], m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M(), MtDM->krowM(), MtDM->jcolM(), MtDM->M());
}
void SparseMatrix::matMultTrans(SymmetricMatrix** res) const {
   int m = mStorage->m;
   int n = mStorage->n;
   int nnz = mStorage->numberOfNonZeros();

   SimpleVector<double> d(n);
   d.setToConstant(1.0);

   if (*res == nullptr) {
      //we need to form the transpose
      if (!m_Mt) {
         m_Mt = new SparseMatrix(n, m, nnz);
         mStorage->transpose(m_Mt->krowM(), m_Mt->jcolM(), m_Mt->M());
      }
      //find the sparsity pattern of the product -> the buffers for result will be allocated
      int* krowMtM = nullptr;
      int* jcolMtM = nullptr;
      double* dMtM = nullptr;

      m_Mt->mStorage->matTransDSymbMultMat(&d[0], krowM(), jcolM(), M(), &krowMtM, &jcolMtM, &dMtM);
      *res = new SparseSymmetricMatrix(m, krowMtM[m], krowMtM, jcolMtM, dMtM, 1);
   }
   auto* MMt = dynamic_cast<SparseSymmetricMatrix*>(*res);
   m_Mt->mStorage->matTransDMultMat(&d[0], krowM(), jcolM(), M(), MMt->krowM(), MMt->jcolM(), MMt->M());
}


void SparseMatrix::addNnzPerRow(Vector<int>& nnzVec) const {
   auto& vec = dynamic_cast<SimpleVector<int>&>(nnzVec);

   if (mStorageDynamic != nullptr) {
      assert(vec.length() == mStorageDynamic->getM());
      mStorageDynamic->addNnzPerRow(vec.elements());
   }
   else {
      assert(vec.length() == mStorage->m);
      mStorage->addNnzPerRow(vec.elements());
   }
}

void SparseMatrix::addNnzPerCol(Vector<int>& nnzVec) const {
   auto& vec = dynamic_cast<SimpleVector<int>&>(nnzVec);

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

void SparseMatrix::addNnzPerCol(Vector<int>& nnzVec, int begin_cols, int end_cols) const {
   assert(0 <= begin_cols && begin_cols <= end_cols);
   assert(end_cols - begin_cols <= nnzVec.length());

   auto& vec = dynamic_cast<SimpleVector<int>&>(nnzVec);

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

void SparseMatrix::addRowSums(Vector<double>& sumVec) const {
   auto& vec = dynamic_cast<SimpleVector<double>&>(sumVec);

   assert(mStorageDynamic == nullptr && mStorage != nullptr);
   assert(vec.length() == mStorage->m);

   mStorage->addRowSums(vec.elements());
}

void SparseMatrix::addColSums(Vector<double>& sumVec) const {
   assert(m_Mt);
   assert(m_Mt->mStorageDynamic == nullptr && m_Mt->mStorage != nullptr);

   auto& vec = dynamic_cast<SimpleVector<double>&>(sumVec);
   assert(vec.length() == m_Mt->mStorage->m);

   m_Mt->mStorage->addRowSums(vec.elements());
}

void SparseMatrix::getMinMaxVec(bool getMin, bool initializeVec, const SparseStorageDynamic* storage_dynamic, const Vector<double>* coScaleVec,
      Vector<double>& minmaxVec) {
   auto& mvec = dynamic_cast<SimpleVector<double>&>(minmaxVec);

   assert(mvec.length() == storage_dynamic->getM());

   if (initializeVec) {
      if (getMin)
         mvec.setToConstant(std::numeric_limits<double>::max());
      else
         mvec.setToConstant(0.0);
   }

   if (coScaleVec) {
      const auto* covec = dynamic_cast<const SimpleVector<double>*>(coScaleVec);

      storage_dynamic->getRowMinMaxVec(getMin, covec->elements(), mvec.elements());
   }
   else {
      storage_dynamic->getRowMinMaxVec(getMin, nullptr, mvec.elements());
   }
}

void SparseMatrix::getMinMaxVec(bool getMin, bool initializeVec, const SparseStorage* storage, const Vector<double>* coScaleVec,
      Vector<double>& minmaxVec) {
   auto& mvec = dynamic_cast<SimpleVector<double>&>(minmaxVec);

   assert(mvec.length() == storage->m);
   if (initializeVec) {
      if (getMin)
         mvec.setToConstant(std::numeric_limits<double>::max());
      else
         mvec.setToZero();
   }

   if (coScaleVec) {
      const auto* covec = dynamic_cast<const SimpleVector<double>*>(coScaleVec);

      storage->getRowMinMaxVec(getMin, covec->elements(), mvec.elements());
   }
   else {
      storage->getRowMinMaxVec(getMin, nullptr, mvec.elements());
   }
}

void SparseMatrix::getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) const {
   if (hasDynamicStorage())
      getMinMaxVec(getMin, initializeVec, mStorageDynamic, colScaleVec, minmaxVec);
   else
      getMinMaxVec(getMin, initializeVec, mStorage, colScaleVec, minmaxVec);
}

void SparseMatrix::getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) const {
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

void SparseMatrix::initStaticStorageFromDynamic(const Vector<int>& rowNnzVec, const Vector<int>* colNnzVec) {
   assert(mStorageDynamic != nullptr);

   const auto& rowNnzVecSimple = dynamic_cast<const SimpleVector<int>&>(rowNnzVec);
   const auto* colNnzVecSimple = dynamic_cast<const SimpleVector<int>*>(colNnzVec);

   mStorageDynamic->restoreOrder();
   SmartPointer<SparseStorage> staticStorage(
         mStorageDynamic->getStaticStorage(rowNnzVecSimple.elements(), (colNnzVecSimple == nullptr) ? nullptr : colNnzVecSimple->elements()));

   mStorage = staticStorage;

   assert(mStorage->refs() == 2);
}

void SparseMatrix::deleteEmptyRowsCols(const Vector<int>& rowNnzVec, const Vector<int>& colNnzVec) {
   const auto& rowNnzVecSimple = dynamic_cast<const SimpleVector<int>&>(rowNnzVec);
   const auto& colNnzVecSimple = dynamic_cast<const SimpleVector<int>&>(colNnzVec);

   mStorage->deleteEmptyRowsCols(rowNnzVecSimple.elements(), colNnzVecSimple.elements());
}

void SparseMatrix::fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset, double* rowsArrayDense,
      int* rowSparsity) const {

   mStorage->fromGetRowsBlock(rowIndices, nRows, arrayLineSize, arrayLineOffset, rowsArrayDense, rowSparsity);
}

void SparseMatrix::fromGetRowsBlock(int row_start, int n_rows, int array_line_size, int array_line_offset, double* rows_array_dense,
      int* row_sparsity) const {
   assert(0 <= row_start && 0 <= n_rows && 0 <= array_line_size && 0 <= array_line_offset);
   mStorage->fromGetRowsBlock(rows_array_dense, row_start, n_rows, array_line_size, array_line_offset, row_sparsity);
}

void SparseMatrix::deleteEmptyRows(int*& orgIndex) {
   mStorage->deleteEmptyRows(orgIndex);

   if (m_Mt)
      this->updateTransposed();
}

void SparseMatrix::fromGetColsBlock(const int* colIndices, int nCols, int arrayLineSize, int arrayLineOffset, double* colsArrayDense,
      int* rowSparsity) {
   if (!m_Mt)
      initTransposed();

   m_Mt->getStorageRef().fromGetRowsBlock(colIndices, nCols, arrayLineSize, arrayLineOffset, colsArrayDense, rowSparsity);
}

void SparseMatrix::fromGetColsBlock(int col_start, int n_cols, int array_line_size, int array_line_offset, double* cols_array_dense,
      int* row_sparsity) {
   if (!m_Mt)
      initTransposed();

   m_Mt->getStorageRef().fromGetRowsBlock(cols_array_dense, col_start, n_cols, array_line_size, array_line_offset, row_sparsity);
}

bool SparseMatrix::hasTransposed() const {
   return (m_Mt != nullptr);
}

void SparseMatrix::freeDynamicStorage() {
   delete mStorageDynamic;
   mStorageDynamic = nullptr;
}

void SparseMatrix::updateTransposed() const {
   deleteTransposed();
   initTransposed();
}

void SparseMatrix::deleteTransposed() const {
   if (m_Mt) {
      delete m_Mt;
      m_Mt = nullptr;
   }
}

void SparseMatrix::getLinkVarsNnz(std::vector<int>& vec) const {
   mStorage->getLinkVarsNnz(vec);
}

void SparseMatrix::updateNonEmptyRowsCount(std::vector<int>& rowcount) const {
   const int m = mStorage->m;
   const int* const rowStart = mStorage->krowM;

   assert(unsigned(m) == rowcount.size());

   for (int i = 0; i < m; i++)
      if (rowStart[i] != rowStart[i + 1])
         rowcount[i]++;
}

void SparseMatrix::updateNonEmptyRowsCountNew(int block_id, std::vector<int>& n_blocks_per_row, std::vector<int>& row_start_block,
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

void SparseMatrix::updateNonEmptyRowsCount(int block_id, std::vector<int>& n_blocks_per_row, std::vector<int>& row_start_block,
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

SparseMatrix& SparseMatrix::getTranspose() {
   if (!m_Mt) {
      initTransposed();
   }

   return *m_Mt;
}

const SparseMatrix& SparseMatrix::getTranspose() const {
   if (!m_Mt) {
      initTransposed();
   }

   return *m_Mt;
}

void SparseMatrix::permuteRows(const std::vector<unsigned int>& permvec) {
   mStorage->permuteRows(permvec);

   if (m_Mt)
      m_Mt->mStorage->permuteCols(permvec);
}

void SparseMatrix::permuteCols(const std::vector<unsigned int>& permvec) {
   mStorage->permuteCols(permvec);

   if (m_Mt)
      m_Mt->mStorage->permuteRows(permvec);
}

int SparseMatrix::appendRow(const SparseMatrix& matrix_row, int row) {
   assert(hasDynamicStorage());
   assert(!hasTransposed());
   assert(matrix_row.hasDynamicStorage());

   mStorageDynamic->appendRow(matrix_row.getStorageDynamicRef(), row);

   return mStorageDynamic->getM() - 1;
}

int SparseMatrix::appendCol(const SparseMatrix& matrix_col, int col) {
   assert(matrix_col.hasTransposed());
   assert(matrix_col.hasDynamicStorage());
   assert(hasDynamicStorage());
   assert(!hasTransposed());

   mStorageDynamic->appendRow(matrix_col.getStorageDynamicTransposedRef(), col);

   return mStorageDynamic->getM() - 1;
}

void SparseMatrix::axpyWithRowAt(double alpha, SimpleVector<double>& y, int row) const {
   assert(hasDynamicStorage());

   mStorageDynamic->axpyWithRowAt(alpha, y.elements(), y.length(), row);
}

void SparseMatrix::axpyWithRowAtPosNeg(double alpha, SimpleVector<double>& y_pos, SimpleVector<double>& y_neg, int row) const {
   assert(hasDynamicStorage());
   assert(y_pos.length() == y_neg.length());

   mStorageDynamic->axpyWithRowAtPosNeg(alpha, y_pos.elements(), y_neg.elements(), y_pos.length(), row);
}

double SparseMatrix::localRowTimesVec(const SimpleVector<double>& vec, int row) const {
   assert(hasDynamicStorage());

   return mStorageDynamic->rowTimesVec(vec.elements(), vec.length(), row);
}

void SparseMatrix::removeRow(int row) {
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

void SparseMatrix::removeRowUsingTransposed(int row, SparseStorageDynamic& mat_trans) {
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

void SparseMatrix::removeCol(int col) {
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

void SparseMatrix::removeEntryAtRowCol(int row, int col) {
   assert(hasDynamicStorage());

   if (hasTransposed())
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());

   mStorageDynamic->removeEntryAtRowCol(row, col);

   if (hasTransposed()) {
      m_Mt->removeEntryAtRowCol(col, row);
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
   }
}

void SparseMatrix::removeEntryAtRowColIndex(int row, int col_index) {
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

void SparseMatrix::addColToRow(double coeff, int col, int row) {
   assert(hasDynamicStorage());

   if (hasTransposed())
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());

   mStorageDynamic->addColToRow(coeff, col, row);

   if (hasTransposed()) {
      m_Mt->addColToRow(coeff, row, col);
      assert(mStorageDynamic->getNVals() == m_Mt->getStorageDynamic()->getNVals());
   }
}

SparseMatrix* SparseMatrix::shaveLeft(int n_cols) {
   assert(!hasDynamicStorage());
   assert(n_cols <= mStorage->n);

   SparseStorage* border = mStorage->shaveLeft(n_cols);

   if (n_cols != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }

   return new SparseMatrix(border);
}

GeneralMatrix* SparseMatrix::shaveBottom(int n_rows) {
   assert(!hasDynamicStorage());
   assert(n_rows <= mStorage->m);

   SparseStorage* border = mStorage->shaveBottom(n_rows);

   if (n_rows != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }

   return new SparseMatrix(border);
}

void SparseMatrix::dropNEmptyRowsBottom(int n_rows) {
   assert(!hasDynamicStorage());
   assert(n_rows <= mStorage->m);

   mStorage->dropNEmptyRowsBottom(n_rows);

   if (n_rows != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }
}

void SparseMatrix::dropNEmptyRowsTop(int n_rows) {
   assert(!hasDynamicStorage());
   assert(n_rows <= mStorage->m);

   mStorage->dropNEmptyRowsTop(n_rows);

   if (n_rows != 0 && m_Mt) {
      delete m_Mt;
      this->initTransposed(false);
   }
}