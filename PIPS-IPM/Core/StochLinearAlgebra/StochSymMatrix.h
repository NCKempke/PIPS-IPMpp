#ifndef STOCHSYMMATRIX_H
#define STOCHSYMMATRIX_H

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"
#include "StringGenMatrix.h"
#include "pipsport.h"

#include <vector>
#include <iostream>
#include <fstream>

#include "mpi.h"

class BorderedSymMatrix;


/*
 * [ Q0  R1^T ... RN^T ]
 * [ R1  Q1            ]
 * [ R2     Q2         ]
 * [ .         .       ]
 * [ .           .     ]
 * [ RN            QN  ]
 */

class StochSymMatrix : public SymMatrix {

private:

  // note: also used for dummy class!
  virtual void deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec, const OoqpVectorBase<int>* linkParent);
  virtual void writeToStreamDenseChild(stringstream& out, int offset) const;

public:
  /** Constructs a matrix with local size 'local_n' having 'local_nnz' local nonzeros
      and set the global size and the id to to 'global_n' and 'id', respectively.
      The parameter 'id' is used for output/debug purposes only.
      The created matrix will have no children.*/
  StochSymMatrix( int id, long long global_n, int local_n, int local_nnz, MPI_Comm mpiComm );
  StochSymMatrix( int id, long long global_n, 
		  int diag_n, int diag_nnz, 
		  int border_n, int border_nnz,
		  MPI_Comm mpiComm_);
  virtual ~StochSymMatrix();

  std::vector<StochSymMatrix*> children;
  SparseSymMatrix* diag;
  SparseGenMatrix* border;
  int id;
  long long n;
  MPI_Comm mpiComm;
  int iAmDistrib;
  
  virtual void AddChild(StochSymMatrix* child);

  virtual StochSymMatrix* clone() const;

  virtual int isKindOf( int type ) const;
  virtual void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent );
  virtual void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent );

  virtual void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info );

  virtual void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info );

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;

  long long size() const override;

  void symAtPutSubmatrix( int destRow, int destCol,
				  DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent ) override;

  void fromGetSpRow( int row, int col,
                             double A[], int lenA, int irowA[], int& nnz,
                             int rowExtent, int& info ) override;

  virtual void atPutZeros( int row, int col,
			   int rowExtent, int colExtent );
  void mult ( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const override;
  void transMult ( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const override;
  
  double abmaxnorm() const override;
  
  void writeToStream(ostream& out) const override;

  void writeToStreamDense(std::ostream& out) const override;

  virtual void randomizePSD(double * seed);
  
  virtual void getDiagonal( OoqpVector& vec );
  void setToDiagonal( const OoqpVector& vec ) override;
  virtual void atPutDiagonal( int idiag, OoqpVector& v );
  virtual void fromGetDiagonal( int idiag, OoqpVector& x );

  virtual void putSparseTriple( int irow[], int len, int jcol[], double A[], 
				int& info );

  void symmetricScale ( const OoqpVector& vec ) override;
  void columnScale ( const OoqpVector& vec ) override;
  void rowScale ( const OoqpVector& vec ) override;
  virtual void scalarMult( double num );

  // note: also used for dummy class!
  virtual void deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec)
  {
     deleteEmptyRowsCols(nnzVec, nullptr);
  }

  // TODO specify border bottom and left..
  virtual BorderedSymMatrix* raiseBorder( int n_vars );

 protected:
  virtual void shaveBorder(int n_vars, StringGenMatrix*& border_vertical);

  StochSymMatrix* parent;
};

/** 
 * Dummy stochastic symmetric matrix
 */

class StochSymDummyMatrix : public StochSymMatrix
{


private:
   void writeToStreamDenseChild(stringstream& out, int offset) const override {};


protected:

public:

  StochSymDummyMatrix(int id_)
    : StochSymMatrix(id_, 0, 0, 0, MPI_COMM_NULL) {};

  virtual ~StochSymDummyMatrix(){};

  StochSymDummyMatrix* clone() const override { return new StochSymDummyMatrix(id); };

  void AddChild(StochSymMatrix* child) override {};

  int isKindOf( int type ) const override;

  void atPutDense( int row, int col, double * A, int lda,
			   int rowExtent, int colExtent ) override {};
  void fromGetDense( int row, int col, double * A, int lda,
			     int rowExtent, int colExtent ) override {};

  void symAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			      int& info ) override {};

  void fsymAtPutSpRow( int row, double A[], int lenA, int jcolA[],
			       int& info ) override {};

  void getSize( long long& m, long long& n ) const override { m = 0; n = 0; }
  void getSize( int& m, int& n ) const override { m = 0; n = 0; }

  long long size() const override { return 0; }

  void symAtPutSubmatrix( int destRow, int destCol,
				  DoubleMatrix& M,
				  int srcRow, int srcCol,
				  int rowExtent, int colExtent ) override {};

  void fromGetSpRow( int row, int col,
                             double A[], int lenA, int irowA[], int& nnz,
                             int rowExtent, int& info ) override {};

  void atPutZeros( int row, int col,
			   int rowExtent, int colExtent ) override {};
  void mult ( double beta,  OoqpVector& y,
		      double alpha, const OoqpVector& x ) const override {};
  void transMult ( double beta,  OoqpVector& y,
			   double alpha, const OoqpVector& x ) const override {};
  
  double abmaxnorm() const override { return 0.0; }
  
  void writeToStream(ostream& out) const override {};
  void writeToStreamDense(std::ostream& out) const override {};

  void randomizePSD(double * seed) override {};
  
  void getDiagonal( OoqpVector& vec ) override {};
  void setToDiagonal( const OoqpVector& vec ) override {};
  void atPutDiagonal( int idiag, OoqpVector& v ) override {};
  void fromGetDiagonal( int idiag, OoqpVector& x ) override {};

  void putSparseTriple( int irow[], int len, int jcol[], double A[],
				int& info ) override {};

  void symmetricScale ( const OoqpVector& vec ) override {};
  void columnScale ( const OoqpVector& vec ) override {};
  void rowScale ( const OoqpVector& vec ) override {};
  void scalarMult( double num ) override {};

  BorderedSymMatrix* raiseBorder( int n_vars ) override { assert(0 && "CANNOT SHAVE BORDER OFF OF A DUMMY MATRIX"); return nullptr; };

 protected:
  void shaveBorder( int n_vars, StringGenMatrix*& border_vertical ) override
     { border_vertical = new StringGenDummyMatrix(); };

};

typedef SmartPointer<StochSymMatrix> StochSymMatrixHandle;

#endif
