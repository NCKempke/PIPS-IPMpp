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
      long long m, n;
      const MPI_Comm mpi_comm;
      const bool distributed;
      const int rank;

   public:
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

      void writeToStreamDense(std::ostream& out) const override; // TODO : implement
      virtual std::string writeToStreamDenseRowChildren(int row) const;

      void getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override;
      void getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override;

      void columnScale ( const OoqpVector& vec ) override;
      void rowScale ( const OoqpVector& vec ) override;

      void addRowSums( OoqpVector& vec ) override;
      void addColSums( OoqpVector& vec ) override;

      void getSize( long long& m_, long long& n_ ) const override { m_ = m; n_ = n; };
      void getSize( int& m_, int& n_ ) const override { m_ = m; n_ = n; };


      /* methods not needed for Hierarchical approach */
      double abminnormNonZero( double tol = 1e-30) const override { assert( false && "TODO: implement" ); return 0.0; };
      void atPutDiagonal( int idiag, OoqpVector& x ) override { assert( "not implemented" && 0 ); };
      void fromGetDiagonal( int idiag, OoqpVector& x ) override { assert( "not implemented" && 0 ); };
      void fromGetDense( int row, int col, double * A, int lda, int rowExtent, int colExtent ) override { assert( "not implemented" && 0 ); };
      void fromGetSpRow( int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info ) override { assert( "not implemented" && 0 ); };
      void getDiagonal( OoqpVector& vec ) override { assert( "not implemented" && 0 ); };
      void setToDiagonal( const OoqpVector& vec ) override { assert( "not implemented" && 0 ); };
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
      StringGenMatrix();

      virtual void multVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;
      virtual void multHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

      virtual void transMultVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;
      virtual void transMultHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const;

      virtual void columnScaleVertical( const OoqpVector& vec );
      virtual void columnScaleHorizontal( const OoqpVector& vec );

      virtual void rowScaleVertical( const OoqpVector& vec );
      virtual void rowScaleHorizontal( const OoqpVector& vec );

      virtual void getRowMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* col_scale, OoqpVector& minmax) const;
      virtual void getRowMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* col_scale, OoqpVector& minmax) const;

      virtual void getColMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* row_scale, OoqpVector& minmax) const;
      virtual void getColMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* row_scale, OoqpVector& minmax) const;

      virtual void addRowSumsVertical( OoqpVector& vec );
      virtual void addRowSumsHorizontal( OoqpVector& vec );

      virtual void addColSumsVertical( OoqpVector& vec );
      virtual void addColSumsHorizontal( OoqpVector& vec );
};

/**
 * Dummy version ...
 */
class StringGenDummyMatrix : public StringGenMatrix
{
   public:
      StringGenDummyMatrix() {};

      virtual ~StringGenDummyMatrix() {};
      void addChild(StringGenMatrix* child) override {};
      int isKindOf( int type ) const override { return type == kStringGenDummyMatrix || type == kStringMatrix || type == kStringGenMatrix; };
      void mult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override {};
      void transMult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override {};
      double abmaxnorm() const override {return -std::numeric_limits<double>::infinity();};
      void scalarMult( double num ) override {};
      void writeToStream( std::ostream& out ) const override {};
      std::string writeToStreamDenseRowChildren(int row) const override { return ""; };

      void getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override {};
      void getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override {};
      void columnScale( const OoqpVector& vec ) override {};
      void rowScale( const OoqpVector& vec ) override {};

      void addRowSums( OoqpVector& vec ) override {};
      void addColSums( OoqpVector& vec ) override {};

   protected:
      void multVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const override {};
      void multHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const override {};

      void transMultVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const override {};
      void transMultHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x) const override {};

      void columnScaleVertical( const OoqpVector& vec ) override {};
      void columnScaleHorizontal( const OoqpVector& vec ) override {};

      void rowScaleVertical( const OoqpVector& vec ) override {};
      void rowScaleHorizontal( const OoqpVector& vec ) override {};

      void getRowMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* col_scale, OoqpVector& minmax) const override {};
      void getRowMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* col_scale, OoqpVector& minmax) const override {};

      void getColMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* row_scale, OoqpVector& minmax) const override {};
      void getColMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* row_scale, OoqpVector& minmax) const override {};

      void addRowSumsVertical( OoqpVector& vec ) override {};
      void addRowSumsHorizontal( OoqpVector& vec ) override {};

      void addColSumsVertical( OoqpVector& vec ) override {};
      void addColSumsHorizontal( OoqpVector& vec ) override {};
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_STRINGGENMATRIX_H_ */
