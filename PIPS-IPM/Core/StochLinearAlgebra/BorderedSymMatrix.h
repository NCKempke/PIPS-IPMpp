/*
 * BorderedSymMatrix.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDSYMMATRIX_H_
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDSYMMATRIX_H_

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"
#include "StochSymMatrix.h"
#include "StringGenMatrix.h"
#include "StochVector.h"
#include <vector>
#include "mpi.h"

class StringGenMatrix;
class StochSymMatrix;

class BorderedSymMatrix : public SymMatrix
{
   public:
      BorderedSymMatrix(int id, StochSymMatrix* inner_matrix, StringGenMatrix* border_vertical,
            SymMatrix* bottom_block, MPI_Comm mpiComm_);

      virtual ~BorderedSymMatrix();

      int isKindOf( int type ) const override;

      /** y = beta * y + alpha * this * x */
      void mult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override;

      /** y = beta * y + alpha * this^T * x */
      void transMult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override;

      void fromGetDiagonal( int idiag, OoqpVector& x ) override;

      double abmaxnorm() const override;
      void scalarMult( double num ) override;
      void getSize( long long& m, long long& n ) const override;
      void getSize( int& m, int& n ) const override;
      long long size() const override;

      // ignored methods
      void writeToStreamDense(ostream& out) const override
         {assert( "Not implemented" && 0 );};
      void symAtPutSpRow( int col, double A[], int lenA, int irowA[], int& info ) override
         {assert( "Not implemented" && 0 );};
      void randomizePSD(double * seed) override
         {assert( "Not implemented" && 0 );};
      void symAtPutSubmatrix( int destRow, int destCol,  DoubleMatrix& M, int srcRow, int srcCol, int rowExtent, int colExtent ) override
         {assert( "Not implemented" && 0 );};
      void atPutDiagonal( int idiag, OoqpVector& x ) override
         {assert( "Not implemented" && 0 );};
      void fromGetDense( int row, int col, double * A, int lda, int rowExtent, int colExtent ) override
         {assert( "Not implemented" && 0 );};
      void fromGetSpRow( int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info ) override
         {assert( "Not implemented" && 0 );};
      void getDiagonal( OoqpVector& vec ) override
         {assert( "Not implemented" && 0 );};
      void setToDiagonal( OoqpVector& vec ) override
         {assert( "Not implemented" && 0 );};
      void symmetricScale ( const OoqpVector& vec ) override
         {assert( "Not implemented" && 0 );};
      void columnScale ( const OoqpVector& vec ) override
         {assert( "Not implemented" && 0 );};
      void rowScale ( const OoqpVector& vec ) override
         {assert( "Not implemented" && 0 );};
      void writeToStream(ostream& out) const override
         {assert( "Not implemented" && 0 );};
      void putSparseTriple( int irow[], int len, int jcol[], double A[], int& info ) override
         {assert( "Not implemented" && 0 );};

      // TODO could be more general..
      StochSymMatrix* inner_matrix;
      StringGenMatrix* border_vertical;

      SymMatrix* bottom_block;

   protected:

      int id;
      MPI_Comm mpiComm;
      const int iAmDistrib;

      long long n;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDSYMMATRIX_H_ */
