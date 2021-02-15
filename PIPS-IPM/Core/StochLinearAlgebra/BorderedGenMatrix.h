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
class SparseGenMatrix;


/*
 * representing a matrix of type
 *
 * [ B  K ]
 * [ C  B']
 *
 * Where K is a StochGenMatrix
 */

// TODO : make more general ? K B B' and C can be any matrices in theory..
class BorderedGenMatrix : public GenMatrix
{
   public:
      StochGenMatrix * const inner_matrix{};
      StringGenMatrix * const border_left{};
      StringGenMatrix * const border_bottom{};
	
      // TODO: is SparseGenMatrix appropriate? What does this block look like -> it has parts of the diagonals in it for inequality linking constraints and nothing else?
      SparseGenMatrix * const bottom_left_block{};

   protected:

      MPI_Comm mpi_comm{MPI_COMM_NULL};
      const bool distributed{false};
      const int rank{-1};

      long long m{0};
      long long n{0};

   public:
      BorderedGenMatrix(StochGenMatrix* inner_matrix, StringGenMatrix* border_left,
            StringGenMatrix* border_bottom, SparseGenMatrix* bottom_left_block, MPI_Comm mpi_comm_);
      virtual ~BorderedGenMatrix();

      int isKindOf( int matrixType ) const override;

      /** y = beta * y + alpha * this * x */
      void mult( double beta,  OoqpVector& y, double alpha, const OoqpVector& x ) const override;
      /** y = beta * y + alpha * this^T * x */
      void transMult( double beta, OoqpVector& y, double alpha,  const OoqpVector& x ) const override;

      double abmaxnorm() const override;
      void columnScale( const OoqpVector& vec ) override;
      void rowScale( const OoqpVector& vec ) override;

      void scalarMult( double num ) override;

      void getSize( long long& m_, long long& n_ ) const override;
      void getSize( int& m_, int& n_ ) const override;

      void getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override;
      void getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override;

      void writeToStreamDense( std::ostream& out ) const;

      void addRowSums( OoqpVector& vec ) const override;
      void addColSums( OoqpVector& vec ) const override;

      /* methods not needed for Hierarchical approach */
      double abminnormNonZero( double ) const override { assert( false && "TODO: implement" ); return 0.0; };
      void writeToStream( std::ostream& ) const override { assert(0 && "not implemented"); }; // TODO : implement maybe?
      void getDiagonal( OoqpVector& ) override  { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
      void setToDiagonal( const OoqpVector& ) override { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
      void atPutDiagonal( int, const OoqpVector& ) override { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
      void fromGetDiagonal( int, OoqpVector& ) override { assert(0 && "not implemented"); }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?

      void matTransDMultMat(OoqpVector&, SymMatrix** ) override { assert(0 && "not implemented"); }; // TODO : needed?
      void matTransDinvMultMat(OoqpVector&, SymMatrix** ) override { assert(0 && "not implemented"); }; // TODO : needed?

      void getNnzPerRow( OoqpVectorBase<int>& ) override { assert( 0 && "not implemented"); };
      void getNnzPerCol( OoqpVectorBase<int>& ) override { assert( 0 && "not implemented"); };
      void fromGetDense( int, int, double*, int, int, int ) override { assert(0 && "not implemented"); };
      void fromGetSpRow( int, int, double*, int, int*, int&, int, int& ) override { assert(0 && "not implemented"); };
      void putSparseTriple( int*, int, int*, double*, int& ) override { assert(0 && "not implemented"); };
      void symmetricScale ( const OoqpVector& ) override { assert(0 && "not implemented"); };
      void atPutSubmatrix( int, int, DoubleMatrix&, int, int, int, int ) override { assert(0 && "not implemented"); };
      void atPutDense( int, int, double*, int, int, int ) override { assert(0 && "not implemented"); };
      void atPutSpRow( int, double*, int, int*, int& ) override { assert(0 && "not implemented"); };
      void randomize(double, double, double* ) override { assert(0 && "not implemented"); };
   private:

      template<typename T>
      bool hasVecStructureForBorderedMat( const OoqpVectorBase<T>& vec, bool row_vec ) const;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDGENMATRIX_H_ */
