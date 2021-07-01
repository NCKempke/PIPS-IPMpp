/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSESYMMETRICMATRIX_H
#define DENSESYMMETRICMATRIX_H

#include <memory>
#include "./DenseStorage.h"

class SparseSymmetricMatrix;

class SparseMatrix;

class DenseMatrix;

/** A class representing dense, symmetric matrices
 * @ingroup DenseLinearAlgebra
 */
class DenseSymmetricMatrix : public SymmetricMatrix {
public:
   std::shared_ptr<DenseStorage> mStorage;

   explicit DenseSymmetricMatrix(int size);
   DenseSymmetricMatrix(double Q[], int size);

   [[nodiscard]] int is_a(int matrixType) const override;

   virtual void mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;
   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   virtual void transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;
   void transpose_mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override;
   [[nodiscard]] long long n_rows() const override;
   [[nodiscard]] long long n_columns() const override;

   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;
   void write_to_stream(std::ostream& out) const override;
   void write_to_streamDense(std::ostream& out) const override;

   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const override;

   void fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int rowExtent, int& info) const override;

   void symmetricScale(const Vector<double>& vec) override;
   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   void symAtPutSpRow(int col, const double A[], int lenA, const int irowA[], int& info) override;

   void putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) override;

   virtual void atAddOuterProductOf(int row, int col, double alpha, double* x, int incx, int nx);

   void symAtPutSubmatrix(int destRow, int destCol, const AbstractMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) override;

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int idiag, Vector<double>& v) const override;

   void diagonal_add_constant_from(int from, int length, double value) override;
   void diagonal_set_to_constant_from(int from, int length, double value) override;

   double* operator[](int index) { return mStorage->M[index]; }
   const double* operator[](int index) const { return mStorage->M[index]; }

   /** Return a pointer to the first element in the matrix */
   double* elements() { return mStorage->M[0]; };
   /** Return mMat, an    */
   double** Mat() { return mStorage->M; };

   [[nodiscard]] long long size() const override;

   [[nodiscard]] DenseStorage& getStorage() { return *mStorage; }
   [[nodiscard]] const DenseStorage& getStorage() const { return *mStorage; }

   /* this = alpha * op(A)*op(B)  +   beta * this */
   void matMult(double alpha, GeneralMatrix& A_, int transA, GeneralMatrix& B_, int transB, double beta);
   void symAtPutSubmatrix(int destRow, int destCol, const AbstractMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent, int forceSymUpdate);

   /**
    * Performs a rank-k update. Depending on the value of 'trans', i.e.,
    *   - this=alpha*this + beta*U*U'   if trans=0
    *   - this=alpha*this + beta*U'*U   if trans<>0
    */
   void atRankkUpdate(double alpha, double beta, DenseMatrix& U, int trans);

   [[nodiscard]] int getNumberOfNonZeros() const;

   /** adds matrix to this starting at row_0 col_0 - matrix must lie completely either in lower or upper triangular part */
   void add_matrix_at(const DenseMatrix& matrix, int row_0, int col_0);
   void add_matrix_at_without_diag(const SparseSymmetricMatrix& matrix, int row_0, int col_0);
   void add_matrix_at(const SparseMatrix& matrix, int row_0, int col_0);

};


#endif /* DENSESYMMETRICMATRIX_H */
