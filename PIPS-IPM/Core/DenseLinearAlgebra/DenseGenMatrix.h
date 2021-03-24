/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSEGENMATRIX_H
#define DENSEGENMATRIX_H

#include "DoubleMatrix.h"
#include "DenseStorage.h"
#include "DenseGenMatrixHandle.h"

class DoubleLinearSolver;
class SparseGenMatrix;

/** A class of dense, non-symmetric, possibly non-square, matrices.
 *  @ingroup DenseLinearAlgebra
 */
class DenseGenMatrix : public GenMatrix {
public:
  DenseStorageHandle mStorage;

  DenseGenMatrix( int size );
  DenseGenMatrix( int m, int n );
  DenseGenMatrix( double A[], int m, int n );

  int isKindOf( int matType ) const override;

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;

  void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override;

  /** Fill a region of this matrix with zeros.
   *
   *  The region starts at (row, col) and extends rowExtent rows
   *  and colExtent columns.
   */
  virtual void atPutZeros( int row, int col,
			   int rowExtent, int colExtent );

  void getDiagonal( OoqpVector& vec ) override;
  void setToDiagonal( const OoqpVector& vec ) override;

  void atPutSubmatrix( int destRow, int destCol,
			       DoubleMatrix& M,
			       int srcRow, int srcCol,
			       int rowExtent, int colExtent ) override;
  void atPutSpRow( int row, double A[], int lenA, int jcolA[],
			   int& info ) override;

  void putSparseTriple( int irow[], int len, int jcol[], double A[],
				int& info ) override;

  void mult ( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const override;
  virtual void mult ( double beta,  double y[], int incy,
		      double alpha, const double x[], int incx ) const;

  void transMult ( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const override;
  virtual void transMult ( double beta,  double y[], int incy,
			   double alpha, const double x[], int incx ) const;

  void matTransDMultMat( OoqpVector&, SymMatrix** ) override { assert(false && "not implemented"); };
  void matTransDinvMultMat(OoqpVector&, SymMatrix** ) override { assert(false && "not implemented"); };

  void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override;

  void fromGetSpRow( int row, int col,
			     double A[], int lenA, int jcolA[], int& nnz,
			     int rowExtent, int& info ) override;

  void columnScale( const OoqpVector& vec ) override;
  void rowScale( const OoqpVector& vec ) override;
  void symmetricScale( const OoqpVector &vec ) override;
  void scalarMult( double num ) override;

  double abmaxnorm() const override;
  double abminnormNonZero( double tol = 1e-30 ) const override;

  void writeToStream( std::ostream& out ) const override;
  void writeToStreamDense( std::ostream& out ) const override;
  void randomize( double alpha, double beta, double * seed ) override;

  void atPutDiagonal( int idiag, const OoqpVector& v ) override;
  void atAddDiagonal( int idiag, const OoqpVector& v ) override;
  void fromGetDiagonal( int idiag, OoqpVector& v ) override;
  /** Get a row of this matrix. */
  virtual void getRow ( int rowIndex, OoqpVector& v_in);

  double * operator[]( int index ) { return mStorage->M[index]; }

  const double * operator[]( int index ) const
  { return mStorage->M[index]; }

  /** Return a pointer to the first element in the matrix */
  double* elements() { return mStorage->M[0]; };
  /** Return mMat, an    */
  double **Mat() { return mStorage->M; };

  DenseStorage& getStorageRef() { return *mStorage; }
  DenseStorageHandle getStorageHandle() { return mStorage; }

  /* the following functions added by C.Petra 09/09 */

  /** MatMat product
   *
   * this = alpha* op(A) * op(B) + beta*this
   *
   * op(...) specifies whether to use the matrix or its transpose
   */
  virtual void matMult(double alpha,
		       DenseGenMatrix& A, int transA,
		       DenseGenMatrix& B, int transB,
		       double beta);

  /* compute beta * res += alpha * this * mat where mat gets multiplied to the submatrix
   * starting at mul_start and the results gets added starting at res_start */
  void multMatAt( int mul_start, double beta, int res_start, DenseGenMatrix& res, double alpha, /* const */ SparseGenMatrix& mat ) const;

};

#endif
