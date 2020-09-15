/*
 * BorderedGenMatrix.h
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDGENMATRIX_H_
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDGENMATRIX_H_

#include "DoubleMatrix.h"

#include <vector>
#include "mpi.h"

class StringGenMatrix;
class StochGenMatrix;

class BorderedGenMatrix : public GenMatrix
{
   public:
      StochGenMatrix* inner_matrix;
      StringGenMatrix* border_right;
      StringGenMatrix* border_bottom;
	
      // TODO: is StochGenMatrix appropriate? What does this block look like -> it has parts of the diagonals in it for inequality linking constraints and nothing else?
      StochGenMatrix* bottom_right_block;

   protected:

      MPI_Comm mpi_comm;
      const bool distributed;
      const int rank;

   private:
      const long long m; // number rows in border
      const long long n; // number cols in border

   public:
      BorderedGenMatrix(StochGenMatrix* inner_matrix, StringGenMatrix* border_right,
            StringGenMatrix* border_bottom, StochGenMatrix* bottom_right_block, MPI_Comm mpi_comm_);
      virtual ~BorderedGenMatrix();

      int isKindOf( int matrixType ) const override;

      /** y = beta * y + alpha * this * x */
      void mult( double beta,  OoqpVector& y, double alpha, const OoqpVector& x ) const override; // TODO : implement

      /** y = beta * y + alpha * this^T * x */
      void transMult ( double beta,   OoqpVector& y, double alpha,  const OoqpVector& x ) const override; // TODO : implement

      double abmaxnorm() const override;

      // probably needed by scaler ?
      void columnScale ( const OoqpVector& vec ) override; // TODO : implement
      void rowScale ( const OoqpVector& vec ) override; // TODO : implement

      void scalarMult( double num ) override; // TODO : implement

      void getSize( long long& m_, long long& n_ ) const override;
      void getSize( int& m_, int& n_ ) const override;

      // TODO which of the following are necessary?
      void getNnzPerRow(OoqpVectorBase<int>& nnzVec) override; // TODO : implement
      void getNnzPerCol(OoqpVectorBase<int>& nnzVec) override; // TODO : implement
      void addRowSums( OoqpVector& vec ) override; // TODO : implement
      void addColSums( OoqpVector& vec ) override; // TODO : implement
      void getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override; // TODO : implement
      void getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override; // TODO : implement

      /* methods not needed for Hierarchical approach */
      void fromGetDense( int row, int col, double * A, int lda, int rowExtent, int colExtent ) override { assert(0 && "not implemented"); };
      void fromGetSpRow( int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info ) override { assert(0 && "not implemented"); };
      void putSparseTriple( int irow[], int len, int jcol[], double A[], int& info ) override { assert(0 && "not implemented"); };
      void writeToStream( std::ostream& out ) const override { assert(0 && "not implemented"); }; // TODO : implement maybe?
      void writeToStreamDense( std::ostream& out ) const override { assert(0 && "not implemented"); };; // TODO implement maybe?
      void getDiagonal( OoqpVector& vec ) override  { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
      void setToDiagonal( OoqpVector& vec ) override { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
      void atPutDiagonal( int idiag, OoqpVector& x ) override { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
      void fromGetDiagonal( int idiag, OoqpVector& x ) override { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
      void symmetricScale ( const OoqpVector& vec ) override { assert(0 && "not implemented"); };
      void atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent ) override { assert(0 && "not implemented"); };
      void atPutDense( int row, int col, double * A, int lda, int rowExtent, int colExtent ) override { assert(0 && "not implemented"); };
      void atPutSpRow( int col, double A[], int lenA, int jcolA[], int& info ) override { assert(0 && "not implemented"); };
      void randomize(double alpha, double beta, double * seed) override { assert(0 && "not implemented"); };
      void matTransDMultMat(OoqpVector& d, SymMatrix** res) override { assert(0 && "not implemented"); }; // TODO : needed?
      void matTransDinvMultMat(OoqpVector& d, SymMatrix** res) override { assert(0 && "not implemented"); }; // TODO : needed?
   private:

      bool hasVecStructureForBorderedMat( const OoqpVector& vec, bool row_vec ) const;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDGENMATRIX_H_ */
