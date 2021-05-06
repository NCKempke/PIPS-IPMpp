/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSESYMMATRIX_H
#define DENSESYMMATRIX_H

#include <memory>
#include "DenseStorage.h"
#include "AbstractMatrix.h"
#include "DenseSymMatrixHandle.h"

class SparseSymmetricMatrix;

class SparseMatrix;

class DenseMatrix;

/** A class representing dense, symmetric matrices
 * @ingroup DenseLinearAlgebra
 */
class DenseSymMatrix : public SymmetricMatrix {
public:
   std::shared_ptr<DenseStorage> mStorage;

   explicit DenseSymMatrix(int size);
   DenseSymMatrix(double Q[], int size);

   [[nodiscard]] int is_a(int matrixType) const override;

   virtual void mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;
   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   virtual void transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;
   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   void getSize(long long& m, long long& n) const override;
   void getSize(int& m, int& n) const override;

   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;
   void writeToStream(std::ostream& out) const override;
   void writeToStreamDense(std::ostream& out) const override;

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

   [[nodiscard]] DenseStorage& getStorageRef() { return *mStorage; }
   [[nodiscard]] const DenseStorage& getStorageRef() const { return *mStorage; }
   [[nodiscard]] std::shared_ptr<DenseStorage> getStorageHandle() const { return mStorage; }

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
};


#endif
