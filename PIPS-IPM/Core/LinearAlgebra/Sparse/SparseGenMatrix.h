/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSEGENMATRIX_H
#define SPARSEGENMATRIX_H

#include "Vector.hpp"
#include "SmartPointer.h"
#include "DoubleMatrix.h"
#include "SparseStorage.h"
#include "SparseStorageDynamic.h"
#include "SparseGenMatrixHandle.h"
#include "SimpleVector.h"
#include <vector>
#include "pipsport.h"

/** Represents sparse non-symmetric, possibly non-square matrices stored in
 *  row-major Harwell-Boeing format.
 *  @ingroup SparseLinearAlgebra
 */
class SparseGenMatrix : public GenMatrix {
private:
   static void
   getMinMaxVec(bool getMin, bool initializeVec, const SparseStorage* storage, const Vector<double>* coScaleVec, Vector<double>& minmaxVec);
   static void getMinMaxVec(bool getMin, bool initializeVec, const SparseStorageDynamic* storage_dynamic, const Vector<double>* coScaleVec,
         Vector<double>& minmaxVec);
protected:
   SparseStorageHandle mStorage;
   SparseStorageDynamic* mStorageDynamic{};

   // in the case of A'*A we internally form the transpose only once
   SparseGenMatrix* m_Mt{};

public:

   SparseGenMatrix() = default;

   void updateTransposed();
   void deleteTransposed();

   SparseGenMatrix(int rows, int cols, int nnz);
   SparseGenMatrix(int rows, int cols, int nnz, int krowM[], int jcolM[], double M[], int deleteElts = 0);
   SparseGenMatrix(SparseStorage* m_storage);

   GenMatrix* cloneFull(bool switchToDynamicStorage = false) const override;
   GenMatrix* cloneEmptyRows(bool switchToDynamicStorage = false) const override;
   virtual SparseGenMatrix* cloneEmptyRowsTransposed(bool switchToDynamicStorage = false) const;

   void getSize(long long& m, long long& n) const override;
   void getSize(int& m, int& n) const override;

   /** The actual number of structural non-zero elements in this sparse
    *  matrix. This includes so-called "accidental" zeros, elements that
    *  are treated as non-zero even though their value happens to be zero.
    */
   int numberOfNonZeros() const override;

   int isKindOf(int matType) const override;

   void atPutDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) override;
   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) override;

   virtual bool isEmpty() const { return mStorage->n == 0 && mStorage->m == 0 && mStorage->len == 0; };

   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void symmetricScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   void fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) override;
   void atPutSubmatrix(int destRow, int destCol, DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) override;
   void atPutSpRow(int col, double A[], int lenA, int jcolA[], int& info) override;

   void putSparseTriple(int irow[], int len, int jcol[], double A[], int& info) override;

   void getDiagonal(Vector<double>& vec) override;
   void setToDiagonal(const Vector<double>& vec) override;

   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   virtual void mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;

   virtual void multMatSymUpper(double beta, SymMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const;

   virtual void transmultMatSymUpper(double beta, SymMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const;

   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   virtual void transMult(double beta, Vector<double>& y_in, int incy, double alpha, const Vector<double>& x_in, int incx) const;
   virtual void transMult(double beta, double y_in[], int incy, double alpha, const double x_in[], int incx) const;

   /** y = beta * y + this^T diag(d)^-1 x */
   virtual void transMultD(double beta, Vector<double>& y, double alpha, const Vector<double>& x, const Vector<double>& d) const;

   /** res = this^T * D * this where D=diag(d) is a diagonal matrix. */
   void matTransDMultMat(Vector<double>& d, SymMatrix** res) override;
   /** res = this^T * inv(D) * this where D=diag(d) is a diagonal matrix. */
   void matTransDinvMultMat(Vector<double>& d, SymMatrix** res) override;
   /** C = this^T * D * this where D=diag(d) is a diagonal matrix. */

   /** initialize (dynamic) transposed matrix */
   virtual void initTransposed(bool dynamic = false);

   /** res = this * this^T */
   virtual void matMultTrans(SymMatrix** res);

   double abmaxnorm() const override;
   double abminnormNonZero(double tol = 1e-30) const override;

   void writeToStream(std::ostream& out) const override;
   void writeToStreamDense(std::ostream& out) const override;
   void writeToStreamDenseRow(std::ostream& out, int rowidx) const override;
   void writeDashedLineToStream(std::ostream& out) const override;

   /** Make the elements in this matrix symmetric. The elements of interest
    *  must be in the lower triangle, and the upper triangle must be empty.
    *  @param info zero if the operation succeeded. Otherwise, insufficient
    *  space was allocated to symmetrize the matrix.
    */
   virtual void symmetrize(int& info);

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int idiag, Vector<double>& v) override;

   SparseStorageHandle getStorageHandle() { return mStorage; }
   const SparseStorageHandle getStorageHandle() const { return mStorage; }

   SparseStorage& getStorageRef() { return *mStorage; }
   const SparseStorage& getStorageRef() const { return *mStorage; }

   int* krowM() { return mStorage->krowM; }
   const int* krowM() const { return mStorage->krowM; }

   int* jcolM() { return mStorage->jcolM; }
   const int* jcolM() const { return mStorage->jcolM; }

   double* M() { return mStorage->M; }
   const double* M() const { return mStorage->M; }

   SparseStorageDynamic* getStorageDynamic() {
      assert(mStorageDynamic != nullptr);
      return mStorageDynamic;
   }
   const SparseStorageDynamic* getStorageDynamic() const {
      assert(mStorageDynamic != nullptr);
      return mStorageDynamic;
   }
   SparseStorageDynamic& getStorageDynamicRef() {
      assert(mStorageDynamic != nullptr);
      return *mStorageDynamic;
   }
   const SparseStorageDynamic& getStorageDynamicRef() const {
      assert(mStorageDynamic != nullptr);
      return *mStorageDynamic;
   }
   SparseStorageDynamic* getStorageDynamicTransposed() {
      assert(m_Mt != nullptr && m_Mt->hasDynamicStorage());
      return m_Mt->getStorageDynamic();
   }
   const SparseStorageDynamic* getStorageDynamicTransposed() const {
      assert(m_Mt != nullptr && m_Mt->hasDynamicStorage());
      return m_Mt->getStorageDynamic();
   }
   SparseStorageDynamic& getStorageDynamicTransposedRef() {
      assert(m_Mt != nullptr && m_Mt->hasDynamicStorage());
      return m_Mt->getStorageDynamicRef();
   }
   const SparseStorageDynamic& getStorageDynamicTransposedRef() const {
      assert(m_Mt != nullptr && m_Mt->hasDynamicStorage());
      return m_Mt->getStorageDynamicRef();
   }
   bool hasDynamicStorage() const { return (mStorageDynamic != nullptr); };

   virtual void addNnzPerRow(Vector<int>& nnzVec) const;
   virtual void addNnzPerCol(Vector<int>& nnzVec);
   void addNnzPerCol(Vector<int>& nnzVec, int begin_cols, int end_cols);

   /** fill vector with absolute minimum/maximum value of each row */
   void getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) override;

   /** fill vector with absolute minimum/maximum value of each column */
   void getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) override;

   void addRowSums(Vector<double>& sumVec) const override;
   void addColSums(Vector<double>& sumVec) const override;

   void initStaticStorageFromDynamic(const Vector<int>& rowNnzVec, const Vector<int>* colNnzVec);

   void permuteRows(const std::vector<unsigned int>& permvec);

   void permuteCols(const std::vector<unsigned int>& permvec);

   void getLinkVarsNnz(std::vector<int>& vec) const;

   void updateNonEmptyRowsCount(std::vector<int>& rowcount) const;

   void updateNonEmptyRowsCountNew(int blockPosition, std::vector<int>& n_blocks_per_row, std::vector<int>& row_start_block,
         std::vector<int>& row_end_block) const;
   void updateNonEmptyRowsCount(int blockPosition, std::vector<int>& n_blocks_per_row, std::vector<int>& row_start_block,
         std::vector<int>& row_end_block) const;

   SparseGenMatrix& getTranspose();

   void deleteEmptyRowsCols(const Vector<int>& rowNnzVec, const Vector<int>& colNnzVec);

   void deleteEmptyRows(int*& orgIndex);

   void fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset, double* rowsArrayDense,
         int* rowSparsity = nullptr) const;

   void fromGetRowsBlock(int row_start, int n_rows, int array_line_size, int array_line_offset, double* rows_array_dense,
         int* row_sparsity = nullptr) const;

   void
   fromGetColsBlock(const int* colIndices, int nCols, int arrayLineSize, int arrayLineOffset, double* colsArrayDense, int* rowSparsity = nullptr);

   void
   fromGetColsBlock(int col_start, int n_cols, int array_line_size, int array_line_offset, double* cols_array_dense, int* row_sparsity = nullptr);

   bool hasTransposed() const;

   void freeDynamicStorage();

   virtual int appendRow(const SparseGenMatrix& matrix_row, int row);
   /** appends col - need matrix_col to have a transposed! */
   virtual int appendCol(const SparseGenMatrix& matrix_col, int col);

   virtual double localRowTimesVec(const SimpleVector<double>& vec, int row) const;
   virtual void axpyWithRowAt(double alpha, SimpleVector<double>& y, int row) const;
   virtual void axpyWithRowAtPosNeg(double alpha, SimpleVector<double>& y_pos, SimpleVector<double>& y_neg, int row) const;

   virtual void removeRow(int row);
   virtual void removeCol(int col);

   virtual void removeRowUsingTransposed(int row, SparseStorageDynamic& mat_trans);

   virtual void removeEntryAtRowCol(int row, int col);
   virtual void removeEntryAtRowColIndex(int row, int col_index);

   virtual void addColToRow(double coeff, int col, int row);

   virtual SparseGenMatrix* shaveLeft(int n_cols);
   GenMatrix* shaveBottom(int n_rows) override;
   virtual void dropNEmptyRowsBottom(int n_rows);
   virtual void dropNEmptyRowsTop(int n_rows);

   virtual ~SparseGenMatrix();
};

#endif
