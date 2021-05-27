/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSEMATRIX_H
#define SPARSEMATRIX_H

#include "../Abstract/Vector.hpp"
#include "../Abstract/AbstractMatrix.h"

#include "../Dense/SimpleVector.h"

#include "SparseStorage.h"
#include "SparseStorageDynamic.h"

#include <vector>

/** Represents sparse non-symmetric, possibly non-square matrices stored in
 *  row-major Harwell-Boeing format.
 *  @ingroup SparseLinearAlgebra
 */
class SparseMatrix : public GeneralMatrix {
private:
   static void
   getMinMaxVec(bool getMin, bool initializeVec, const SparseStorage* storage, const Vector<double>* coScaleVec, Vector<double>& minmaxVec);
   static void getMinMaxVec(bool getMin, bool initializeVec, const SparseStorageDynamic* storage_dynamic, const Vector<double>* coScaleVec,
         Vector<double>& minmaxVec);

   std::unique_ptr<SparseStorage> mStorage{};
   std::unique_ptr<SparseStorageDynamic> mStorageDynamic{};

   /* transposed will be initialized when necessary */
   mutable std::unique_ptr<SparseMatrix> m_Mt{};
public:

   void updateTransposed() const;
   void deleteTransposed() const;

   SparseMatrix() : SparseMatrix(0,0,0) {};
   SparseMatrix(int rows, int cols, int nnz);
   SparseMatrix(int rows, int cols, int nnz, int krowM[], int jcolM[], double M[], int deleteElts = 0);
   explicit SparseMatrix(std::unique_ptr<SparseStorage> m_storage);

   using GeneralMatrix::cloneFull;
   std::unique_ptr<GeneralMatrix> cloneFull(bool switchToDynamicStorage) const override;

   using GeneralMatrix::cloneEmptyRows;
   std::unique_ptr<GeneralMatrix> cloneEmptyRows(bool switchToDynamicStorage) const override;

   std::unique_ptr<SparseMatrix> cloneEmptyRowsTransposed(bool switchToDynamicStorage = false) const;

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override;
   [[nodiscard]] long long n_rows() const override;
   [[nodiscard]] long long n_columns() const override;

   /** The actual number of structural non-zero elements in this sparse
    *  matrix. This includes so-called "accidental" zeros, elements that
    *  are treated as non-zero even though their value happens to be zero.
    */
   int numberOfNonZeros() const override;

   int is_a(int matType) const override;

   void atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) override;
   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const override;

   virtual bool isEmpty() const { return mStorage->n == 0 && mStorage->m == 0 && mStorage->len == 0; };

   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void symmetricScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   void fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const override;
   void atPutSubmatrix(int destRow, int destCol, const AbstractMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) override;
   void atPutSpRow(int col, const double A[], int lenA, const int jcolA[], int& info) override;

   void putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) override;

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;

   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   void mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;

   void multMatSymUpper(double beta, SymmetricMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const;

   void transmultMatSymUpper(double beta, SymmetricMatrix& y, double alpha, const double x[], int yrowstart, int ycolstart) const;

   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   void transMult(double beta, Vector<double>& y_in, int incy, double alpha, const Vector<double>& x_in, int incx) const;
   void transMult(double beta, double y_in[], int incy, double alpha, const double x_in[], int incx) const;

   /** y = beta * y + this^T diag(d)^-1 x */
   void transMultD(double beta, Vector<double>& y, double alpha, const Vector<double>& x, const Vector<double>& d) const;

   /** res = this^T * D * this where D=diag(d) is a diagonal matrix. */
   void matTransDMultMat(const Vector<double>& d, SymmetricMatrix** res) const override;
   /** res = this^T * inv(D) * this where D=diag(d) is a diagonal matrix. */
   void matTransDinvMultMat(const Vector<double>& d, SymmetricMatrix** res) const override;
   /** C = this^T * D * this where D=diag(d) is a diagonal matrix. */

   /** initialize (dynamic) transposed matrix */
   void initTransposed(bool dynamic = false) const;

   /** res = this * this^T */
   void matMultTrans(SymmetricMatrix** res) const;

   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;

   void writeToStream(std::ostream& out) const override;
   void writeToStreamDense(std::ostream& out) const override;
   void writeToStreamDenseRow(std::ostream& out, int rowidx) const override;
   void writeDashedLineToStream(std::ostream& out) const override;

   /** Make the elements in this matrix symmetric. The elements of interest
    *  must be in the lower triangle, and the upper triangle must be empty.
    *  @param info zero if the operation succeeded. Otherwise, insufficient
    *  space was allocated to symmetrize the matrix.
    */
   void symmetrize(int& info);

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int idiag, Vector<double>& v) const override;

   [[nodiscard]] int* krowM() { return mStorage->krowM; }
   [[nodiscard]] const int* krowM() const { return mStorage->krowM; }

   [[nodiscard]] int* jcolM() { return mStorage->jcolM; }
   [[nodiscard]] const int* jcolM() const { return mStorage->jcolM; }

   [[nodiscard]] double* M() { return mStorage->M; }
   [[nodiscard]] const double* M() const { return mStorage->M; }

   [[nodiscard]] SparseStorageDynamic& getStorageDynamic() {
      assert(mStorageDynamic);
      return *mStorageDynamic;
   }

   [[nodiscard]] const SparseStorageDynamic& getStorageDynamic() const {
      assert(mStorageDynamic);
      return *mStorageDynamic;
   }

   [[nodiscard]] const SparseStorageDynamic* getStorageDynamicPtr() const {
      return mStorageDynamic.get();
   }

   [[nodiscard]] SparseStorage& getStorage() {
      assert(mStorage);
      return *mStorage;
   }

   [[nodiscard]] const SparseStorage& getStorage() const {
      assert(mStorage);
      return *mStorage;
   }

   [[nodiscard]] SparseStorageDynamic& getStorageDynamicTransposed() {
      assert(m_Mt && m_Mt->hasDynamicStorage());
      return m_Mt->getStorageDynamic();
   }

   [[nodiscard]] const SparseStorageDynamic* getStorageDynamicTransposedPtr() const {
      assert(m_Mt && m_Mt->hasDynamicStorage());
      return m_Mt->getStorageDynamicPtr();
   }

   [[nodiscard]] const SparseStorageDynamic& getStorageDynamicTransposed() const {
      assert(m_Mt && m_Mt->hasDynamicStorage());
      return m_Mt->getStorageDynamic();
   }

   [[nodiscard]] SparseStorage& getStorageTransposed() {
      assert(m_Mt);
      return m_Mt->getStorage();
   }

   [[nodiscard]] const SparseStorage& getStorageTransposed() const {
      assert(m_Mt);
      return m_Mt->getStorage();
   }

   [[nodiscard]] bool hasDynamicStorage() const { return (mStorageDynamic != nullptr); };

   void addNnzPerRow(Vector<int>& nnzVec) const;
   void addNnzPerCol(Vector<int>& nnzVec) const;
   void addNnzPerCol(Vector<int>& nnzVec, int begin_cols, int end_cols) const;

   /** fill vector with absolute minimum/maximum value of each row */
   void getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) const override;

   /** fill vector with absolute minimum/maximum value of each column */
   void getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) const override;

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

   const SparseMatrix& getTranspose() const;
   SparseMatrix& getTranspose();

   void deleteEmptyRowsCols(const Vector<int>& rowNnzVec, const Vector<int>& colNnzVec);

   void deleteEmptyRows(int*& orgIndex);

   void fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset, double* rowsArrayDense,
         int* rowSparsity = nullptr) const;

   void fromGetRowsBlock(int row_start, int n_rows, int array_line_size, int array_line_offset, double* rows_array_dense,
         int* row_sparsity = nullptr) const;

   void
   fromGetColsBlock(const int* colIndices, int nCols, int arrayLineSize, int arrayLineOffset, double* colsArrayDense, int* rowSparsity = nullptr) const;

   void
   fromGetColsBlock(int col_start, int n_cols, int array_line_size, int array_line_offset, double* cols_array_dense, int* row_sparsity = nullptr) const;

   bool hasTransposed() const;

   void freeDynamicStorage();

   int appendRow(const SparseMatrix& matrix_row, int row);
   int appendCol(const SparseMatrix& matrix_col, int col);

   double localRowTimesVec(const SimpleVector<double>& vec, int row) const;
   void axpyWithRowAt(double alpha, SimpleVector<double>& y, int row) const;
   void axpyWithRowAtPosNeg(double alpha, SimpleVector<double>& y_pos, SimpleVector<double>& y_neg, int row) const;

   void removeRow(int row);
   void removeCol(int col);

   void removeRowUsingTransposed(int row, SparseStorageDynamic& mat_trans);

   void removeEntryAtRowCol(int row, int col);
   void removeEntryAtRowColIndex(int row, int col_index);

   void addColToRow(double coeff, int col, int row);

   std::unique_ptr<SparseMatrix> shaveLeft(int n_cols);
   std::unique_ptr<GeneralMatrix> shaveBottom(int n_rows) override;
   void dropNEmptyRowsBottom(int n_rows);
   void dropNEmptyRowsTop(int n_rows);

   ~SparseMatrix() override = default;
};

#endif
