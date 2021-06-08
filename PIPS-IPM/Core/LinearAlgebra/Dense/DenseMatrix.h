/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSEGENMATRIX_H
#define DENSEGENMATRIX_H

#include "DenseStorage.h"
#include <memory>

class SparseMatrix;

/** A class of dense, non-symmetric, possibly non-square, matrices.
 *  @ingroup DenseLinearAlgebra
 */
class DenseMatrix : public GeneralMatrix {
public:
   std::shared_ptr<DenseStorage> mStorage;

   explicit DenseMatrix(int size);
   DenseMatrix(int m, int n);
   DenseMatrix(double A[], int m, int n);

   [[nodiscard]] int is_a(int matType) const override;

   [[nodiscard]] long long n_rows() const override;
   [[nodiscard]] long long n_columns() const override;

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override;

   void atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) override;

   /** Fill a region of this matrix with zeros.
    *
    *  The region starts at (row, col) and extends rowExtent rows
    *  and colExtent columns.
    */
   virtual void atPutZeros(int row, int col, int rowExtent, int colExtent);

   void putZeros();

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;

   void atPutSubmatrix(int destRow, int destCol, const AbstractMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) override;
   void atPutSpRow(int row, const double A[], int lenA, const int jcolA[], int& info) override;

   void putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) override;

   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   virtual void mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;

   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   virtual void transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;

   void matTransDMultMat(const Vector<double>&, SymmetricMatrix**) const override { assert(false && "not implemented"); };
   void matTransDinvMultMat(const Vector<double>&, SymmetricMatrix**) const override { assert(false && "not implemented"); };

   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const override;
   void fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int rowExtent, int& info) const override;

   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void symmetricScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;

   void write_to_stream(std::ostream& out) const override;
   void write_to_streamDense(std::ostream& out) const override;

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int idiag, Vector<double>& v) const override;
   /** Get a row of this matrix. */
   virtual void getRow(int rowIndex, Vector<double>& v_in);

   double* operator[](int index) { return mStorage->M[index]; }

   const double* operator[](int index) const { return mStorage->M[index]; }

   /** Return a pointer to the first element in the matrix */
   double* elements() { return mStorage->M[0]; };
   /** Return mMat, an    */
   double** Mat() { return mStorage->M; };

   DenseStorage& getStorageRef() { return *mStorage; }
   std::shared_ptr<DenseStorage> getStorageHandle() { return mStorage; }

   /* compute beta * res += alpha * this * mat where mat gets multiplied to the submatrix
    * starting at mul_start and the results gets added starting at res_start */
   void multMatAt(int row_start, int row_end, int col_offset_this, double beta, int row_start_res, int col_offset_result, DenseMatrix& res,
         double alpha, const SparseMatrix& mat) const;

   /* adds mat to this starting at row_0 col_0 */
   void addMatAt(const SparseMatrix& mat, int mat_row_start, int mat_row_end, int this_row_0, int this_col_0);
};

#endif
