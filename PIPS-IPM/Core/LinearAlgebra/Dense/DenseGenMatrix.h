/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSEGENMATRIX_H
#define DENSEGENMATRIX_H

#include "DoubleMatrix.h"
#include "DenseStorage.h"
#include <memory>

class SparseGenMatrix;

/** A class of dense, non-symmetric, possibly non-square, matrices.
 *  @ingroup DenseLinearAlgebra
 */
class DenseGenMatrix : public GenMatrix {
public:
   std::shared_ptr<DenseStorage> mStorage;

   DenseGenMatrix(int size);
   DenseGenMatrix(int m, int n);
   DenseGenMatrix(double A[], int m, int n);

   int isKindOf(int matType) const override;

   int getM() const { return mStorage->m; };
   int getN() const { return mStorage->n; };

   void getSize(long long& m, long long& n) const override;
   void getSize(int& m, int& n) const override;

   void atPutDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) override;

   /** Fill a region of this matrix with zeros.
    *
    *  The region starts at (row, col) and extends rowExtent rows
    *  and colExtent columns.
    */
   virtual void atPutZeros(int row, int col, int rowExtent, int colExtent);

   void putZeros();

   void getDiagonal(Vector<double>& vec) override;
   void setToDiagonal(const Vector<double>& vec) override;

   void atPutSubmatrix(int destRow, int destCol, DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) override;
   void atPutSpRow(int row, double A[], int lenA, int jcolA[], int& info) override;

   void putSparseTriple(int irow[], int len, int jcol[], double A[], int& info) override;

   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   virtual void mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;

   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   virtual void transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;

   void matTransDMultMat(Vector<double>&, SymMatrix**) override { assert(false && "not implemented"); };
   void matTransDinvMultMat(Vector<double>&, SymMatrix**) override { assert(false && "not implemented"); };

   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) override;

   void fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int rowExtent, int& info) override;

   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void symmetricScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   double abmaxnorm() const override;
   double abminnormNonZero(double tol = 1e-30) const override;

   void writeToStream(std::ostream& out) const override;
   void writeToStreamDense(std::ostream& out) const override;

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int idiag, Vector<double>& v) override;
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

   /* the following functions added by C.Petra 09/09 */

   /** MatMat product
    *
    * this = alpha* op(A) * op(B) + beta*this
    *
    * op(...) specifies whether to use the matrix or its transpose
    */
   virtual void matMult(double alpha, DenseGenMatrix& A, int transA, DenseGenMatrix& B, int transB, double beta);

   /* compute beta * res += alpha * this * mat where mat gets multiplied to the submatrix
    * starting at mul_start and the results gets added starting at res_start */
   void multMatAt(int row_start, int row_end, int col_offset_this, double beta, int row_start_res, int col_offset_result, DenseGenMatrix& res,
         double alpha, /* const */ SparseGenMatrix& mat) const;

   /* adds mat to this starting at row_0 col_0 */
   void addMatAt(const SparseGenMatrix& mat, int mat_row_start, int mat_row_end, int this_row_0, int this_col_0);
};

#endif
