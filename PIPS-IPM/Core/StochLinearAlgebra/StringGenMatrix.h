/*
 * StringGenMatrix.h
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_STRINGGENMATRIX_H_
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_STRINGGENMATRIX_H_

#include "DoubleMatrix.h"
#include "SparseGenMatrix.h"
#include "DoubleMatrixTypes.h"

#include "mpi.h"

class StringGenMatrix : public GenMatrix
{
   public:
      // TODO : keep public? ...
      // like stoch gen matrix possibly infinite but we will only have one layer of children at all times...
      std::vector<StringGenMatrix*> children;
      SparseGenMatrix* mat; // never null
      SparseGenMatrix* mat_link; // possibly null
      bool is_vertical;

   protected:
      MPI_Comm mpi_comm;
      const bool distributed;
      const int rank;

   public:
//      StringGenMatrix(int id, int A_m, int A_n, int A_nnz,
//           int ALink_m, int ALink_n, int ALink_nnz,
//           MPI_Comm mpiComm_);

      StringGenMatrix(bool is_vertical, SparseGenMatrix* mat, SparseGenMatrix* mat_link, MPI_Comm mpi_comm_);

      virtual ~StringGenMatrix();

      virtual void addChild(StringGenMatrix* child);

      int isKindOf( int matrix ) const override;

      /** y = beta * y + alpha * this * x */
      void mult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override;

      /** y = beta * y + alpha * this^T * x */
      void transMult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override;

      double abmaxnorm() const override;

      void scalarMult( double num ) override;

      // TODO : again not sure what to return here - not implemented for now..
      void getSize( long long& m, long long& n ) const override { assert( "not implemented" && 0 ); };
      void getSize( int& m, int& n ) const override { assert( "not implemented" && 0 ); };

      void writeToStreamDense(std::ostream& out) const override; // TODO : implement

      void getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override; // TODO : implement
      void getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override; // TODO : implement

      void columnScale ( const OoqpVector& vec ) override;
      void rowScale ( const OoqpVector& vec ) override;

      /* methods not needed for Hierarchical approach */
      void atPutDiagonal( int idiag, OoqpVector& x ) override { assert( "not implemented" && 0 ); };
      void fromGetDiagonal( int idiag, OoqpVector& x ) override { assert( "not implemented" && 0 ); };
      void fromGetDense( int row, int col, double * A, int lda, int rowExtent, int colExtent ) override { assert( "not implemented" && 0 ); };
      void fromGetSpRow( int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info ) override { assert( "not implemented" && 0 ); };
      void getDiagonal( OoqpVector& vec ) override { assert( "not implemented" && 0 ); };
      void setToDiagonal( OoqpVector& vec ) override { assert( "not implemented" && 0 ); };
      void matTransDMultMat(OoqpVector& d, SymMatrix** res) override { assert( "not implemented" && 0 ); };
      void matTransDinvMultMat(OoqpVector& d, SymMatrix** res) override { assert( "not implemented" && 0 ); };
      void symmetricScale ( const OoqpVector& vec ) override { assert( "not implemented" && 0 ); };
      void writeToStream(ostream& out) const override { assert( "not implemented" && 0 ); };
      void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent ) override { assert( "not implemented" && 0 ); };
      void atPutDense( int row, int col, double * A, int lda, int rowExtent, int colExtent ) override { assert( "not implemented" && 0 ); };
      void atPutSpRow( int col, double A[], int lenA, int jcolA[], int& info ) override { assert( "not implemented" && 0 ); };
      void putSparseTriple( int irow[], int len, int jcol[], double A[], int& info ) override { assert( "not implemented" && 0 ); };
      void randomize(double alpha, double beta, double * seed) override { assert( "not implemented" && 0 ); };


   protected:
      void multVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;
      void multHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

      void transMultVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;
      void transMultHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

      void columnScaleVertical( const OoqpVector& vec );
      void columnScaleHorizontal( const OoqpVector& vec );

      void rowScaleVertical( const OoqpVector& vec );
      void rowScaleHorizontal( const OoqpVector& vec );
};



/**
 * Dummy version ...
 */
class StringGenDummyMatrix : public StringGenMatrix
{
   public:
      StringGenDummyMatrix(int id);
//       : StringGenMatrix() {};

      virtual ~StringGenDummyMatrix() {};
      void addChild(StringGenMatrix* child) override {};
      int isKindOf( int type ) const override { return type == kStochGenDummyMatrix; };
      void mult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override {};
      void transMult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override {};
      double abmaxnorm() const override {return -std::numeric_limits<double>::infinity();};
      void scalarMult( double num ) override {};
      void writeToStream( std::ostream& out ) const override {};
      void getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override {};
      void getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override {};
      void columnScale( const OoqpVector& vec ) override {};
      void rowScale( const OoqpVector& vec ) override {};
};



#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_STRINGGENMATRIX_H_ */
