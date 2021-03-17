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
  virtual void writeToStreamDenseChild(std::stringstream& out, int offset) const;

public:
  StochSymMatrix( SymMatrix* diag, SparseGenMatrix* border, MPI_Comm mpiComm );

  /** Constructs a matrix with local size 'local_n' having 'local_nnz' local nonzeros
      and set the global size and the id to to 'global_n' and 'id', respectively.
      The parameter 'id' is used for output/debug purposes only.
      The created matrix will have no children.*/
  StochSymMatrix( long long global_n, int local_n, int local_nnz, MPI_Comm mpiComm );

  ~StochSymMatrix() override;

  std::vector<StochSymMatrix*> children;
  SymMatrix* diag{};
  SparseGenMatrix* border{};

  long long n{0};
  MPI_Comm mpiComm{MPI_COMM_NULL};
  int iAmDistrib{0};
  
  void AddChild(StochSymMatrix* child);

  SymMatrix* clone() const override;

  int isKindOf( int type ) const override;

  void fromGetDense( int, int, double*, int, int, int ) override { assert( false && "Not implemented" ); };
  void symAtPutSpRow( int, double[], int, int[], int&) override { assert( false && "Not implemented" ); };

  void getSize( long long& m, long long& n ) const override;
  void getSize( int& m, int& n ) const override;

  long long size() const override;

  void symAtPutSubmatrix( int, int, DoubleMatrix&, int, int, int, int ) override { assert( false && "Not implemented" ); };;
  void fromGetSpRow( int, int, double[], int, int[], int&, int, int& ) override { assert( false && "Not implemented" ); };

  void mult ( double beta,  OoqpVector& y, double alpha, const OoqpVector& x ) const override;
  void transMult ( double beta,  OoqpVector& y, double alpha, const OoqpVector& x ) const override;
  
  double abmaxnorm() const override;
  double abminnormNonZero( double tol = 1e-30 ) const override;

  void writeToStream(std::ostream&) const override { assert( false && "Not implemented" ); };

  void writeToStreamDense(std::ostream& out) const override;

  void randomizePSD( double* ) override { assert( false && "Not implemented" ); };

  void getDiagonal( OoqpVector& vec ) override;
  void setToDiagonal( const OoqpVector& vec ) override;
  void atPutDiagonal( int idiag, OoqpVector& v ) override;
  void fromGetDiagonal( int idiag, OoqpVector& x ) override;

  void putSparseTriple( int[], int, int[], double[], int& ) override { assert( false && "Not implemented" ); };

  void symmetricScale ( const OoqpVector& vec ) override;
  void columnScale ( const OoqpVector& vec ) override;
  void rowScale ( const OoqpVector& vec ) override;

  void scalarMult( double num ) override;

  // note: also used for dummy class!
  virtual void deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec)
  {
     deleteEmptyRowsCols(nnzVec, nullptr);
  }

  // TODO specify border bottom and left..
  virtual BorderedSymMatrix* raiseBorder( int n_vars );
  virtual void splitMatrix( const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms );

  virtual void recomputeSize();
 protected:
  virtual StringGenMatrix* shaveBorder(int n_vars);
  virtual StringGenMatrix* shaveBorder2(int n_vars);

  StochSymMatrix* parent{};
};

/** 
 * Dummy stochastic symmetric matrix
 */

class StochSymDummyMatrix : public StochSymMatrix
{

private:
   void writeToStreamDenseChild( std::stringstream&, int ) const override {};

public:

  StochSymDummyMatrix() : StochSymMatrix(0, 0, 0, MPI_COMM_NULL) {};

  ~StochSymDummyMatrix() override = default;

  SymMatrix* clone() const override { return new StochSymDummyMatrix(); };

  int isKindOf( int type ) const override;

  void fromGetDense( int, int, double*, int, int, int ) override {};
  void symAtPutSpRow( int, double[], int, int[], int&) override {};

  void getSize( long long& m, long long& n ) const override { m = 0; n = 0; }
  void getSize( int& m, int& n ) const override { m = 0; n = 0; }

  long long size() const override { return 0; }

  void symAtPutSubmatrix( int, int, DoubleMatrix&, int, int,int, int ) override {};

  void fromGetSpRow( int, int, double[], int, int[], int&, int, int& ) override {};

  void mult ( double, OoqpVector&, double, const OoqpVector& ) const override {};
  void transMult ( double, OoqpVector&, double, const OoqpVector& ) const override {};
  
  double abmaxnorm() const override { return 0.0; }
  double abminnormNonZero( double ) const override { return std::numeric_limits<double>::infinity(); }

  void writeToStream( std::ostream& ) const override {};
  void writeToStreamDense( std::ostream& ) const override {};

  void randomizePSD( double* ) override {};
  
  void getDiagonal( OoqpVector& ) override {};
  void setToDiagonal( const OoqpVector& ) override {};
  void atPutDiagonal( int, OoqpVector& ) override {};
  void fromGetDiagonal( int, OoqpVector& ) override {};

  void putSparseTriple( int[], int, int[], double[], int& ) override {};

  void symmetricScale ( const OoqpVector& ) override {};
  void columnScale ( const OoqpVector& ) override {};
  void rowScale ( const OoqpVector& ) override {};
  void scalarMult( double ) override {};

  BorderedSymMatrix* raiseBorder( int ) override { assert(0 && "CANNOT SHAVE BORDER OFF OF A DUMMY MATRIX"); return nullptr; };
  void splitMatrix( const std::vector<unsigned int>&, const std::vector<MPI_Comm>& ) override {};
  void recomputeSize() override {};

 protected:
  StringGenMatrix* shaveBorder( int ) override { return new StringGenDummyMatrix(); };
  StringGenMatrix* shaveBorder2( int ) override { return new StringGenDummyMatrix(); };

};

typedef SmartPointer<StochSymMatrix> StochSymMatrixHandle;

#endif
