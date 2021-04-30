/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSESYMMATRIX_H
#define DENSESYMMATRIX_H

#include <memory>
#include "DenseStorage.h"
#include "DoubleMatrix.h"
#include "DenseSymMatrixHandle.h"

class SparseSymMatrix;

class SparseGenMatrix;

class DenseGenMatrix;

/** A class representing dense, symmetric matrices
 * @ingroup DenseLinearAlgebra
 */
class DenseSymMatrix : public SymMatrix {
public:
   std::shared_ptr<DenseStorage> mStorage;

   DenseSymMatrix(int size);
   DenseSymMatrix(double Q[], int size);

   int isKindOf(int matrixType) const override;

   virtual void mult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;
   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   virtual void transMult(double beta, double y[], int incy, double alpha, const double x[], int incx) const;
   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   void getSize(long long& m, long long& n) const override;
   void getSize(int& m, int& n) const override;

   double abmaxnorm() const override;
   double abminnormNonZero(double tol = 1e-30) const override;
   void writeToStream(std::ostream& out) const override;
   void writeToStreamDense(std::ostream& out) const override;
   void randomizePSD(double* seed) override;

   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) override;

   void fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int rowExtent, int& info) override;

   void symmetricScale(const Vector<double>& vec) override;
   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   void symAtPutSpRow(int col, double A[], int lenA, int irowA[], int& info) override;

   /** Insert the dense array symmetrically (the part that winds up
    *  in the lower triangle of this matrix is significant.)
    */
   virtual void symAtPutDense(int row, int col, double* A, int lda, int rowExtent, int colExtent);
   /** Put a block of zeros into this matrix symmetrically (the part that
    *  winds up in the lower triangle of this matrix is significant.) */
   virtual void symAtPutZeros(int row, int col, int rowExtent, int colExtent);

   void putSparseTriple(int irow[], int len, int jcol[], double A[], int& info) override;

   virtual void atAddOuterProductOf(int row, int col, double alpha, double* x, int incx, int nx);

   void symAtPutSubmatrix(int destRow, int destCol, DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent) override;

   void getDiagonal(Vector<double>& vec) override;
   void setToDiagonal(const Vector<double>& vec) override;

  void atPutDiagonal( int idiag, const Vector<double>& v ) override;
  void atAddDiagonal( int idiag, const Vector<double>& v ) override;
  void fromGetDiagonal( int idiag, Vector<double>& v ) override;

  void diagonal_add_constant_from(int from, int length, double value) override;


   static DenseSymMatrix* randomPSD(int n, double* seed);

   double* operator[](int index) { return mStorage->M[index]; }

   const double* operator[](int index) const { return mStorage->M[index]; }

   /** Return a pointer to the first element in the matrix */
   double* elements() { return mStorage->M[0]; };
   /** Return mMat, an    */
   double** Mat() { return mStorage->M; };

   long long size() const override;

   DenseStorage& getStorageRef() { return *mStorage; }
   const DenseStorage& getStorageRef() const { return *mStorage; }
   std::shared_ptr<DenseStorage> getStorageHandle() { return mStorage; }
   const std::shared_ptr<DenseStorage> getStorageHandle() const { return mStorage; }

   /* this = alpha * op(A)*op(B)  +   beta * this */
   void matMult(double alpha, GenMatrix& A_, int transA, GenMatrix& B_, int transB, double beta);
   virtual void
   symAtPutSubmatrix(int destRow, int destCol, DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent, int forceSymUpdate);

   /**
    * Performs a rank-k update. Depending on the value of 'trans', i.e.,
    *   - this=alpha*this + beta*U*U'   if trans=0
    *   - this=alpha*this + beta*U'*U   if trans<>0
    */
   void atRankkUpdate(double alpha, double beta, DenseGenMatrix& U, int trans);

   int getNumberOfNonZeros() const;
};


#endif
