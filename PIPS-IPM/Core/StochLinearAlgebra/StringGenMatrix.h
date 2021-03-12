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

class sTreeCallbacks;

class StringGenMatrix : public GenMatrix
{
   public:
      // TODO : keep public? ...
      // like stoch gen matrix possibly infinite but we will only have one layer of children at all times...
      std::vector<StringGenMatrix*> children;
      GenMatrix* mat{}; // never null
      GenMatrix* mat_link{}; // possibly null
      bool is_vertical{false};

   protected:
      long long m{0};
      long long n{0};
      const MPI_Comm mpi_comm{MPI_COMM_NULL};
      const bool distributed{false};
      const int rank{-1};

   public:
      StringGenMatrix() = default;

      StringGenMatrix(bool is_vertical, GenMatrix* mat, GenMatrix* mat_link, MPI_Comm mpi_comm_);

      ~StringGenMatrix() override;

      MPI_Comm getComm() const { return mpi_comm; };

      virtual void addChild(StringGenMatrix* child);

      virtual bool isEmpty() const;
      int isKindOf( int matrix ) const override;

      /** y = beta * y + alpha * this * x */
      void mult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override;
      /** y = beta * y + alpha * this^T * x */
      void transMult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const override;

      double abmaxnorm() const override;
      void scalarMult( double num ) override;

      void writeToStreamDense( std::ostream& ) const override;
      void writeToStreamDenseRow( std::ostream& out, int row ) const override;
      void writeDashedLineToStream( std::ostream& out ) const override;

      void getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec ) override;
      void getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec ) override;

      void columnScale ( const OoqpVector& vec ) override;
      void rowScale ( const OoqpVector& vec ) override;

      void addRowSums( OoqpVector& vec ) const override;
      void addColSums( OoqpVector& vec ) const override;

      void getSize( long long& m_, long long& n_ ) const override { m_ = m; n_ = n; };
      void getSize( int& m_, int& n_ ) const override { m_ = m; n_ = n; };

      /** split the current children according to map_child_subchild: the new StringGenMatrices has one additional layer of StringGenMatrices */
      void combineChildrenInNewChildren( const std::vector<unsigned int>& map_child_subchild, const std::vector<MPI_Comm>& child_comms );
      virtual void splitAlongTree( const sTreeCallbacks& tree );

      GenMatrix* shaveBottom( int n_rows ) override;

      /* methods not needed for Hierarchical approach */
      double abminnormNonZero( double) const override { assert( false && "TODO: implement" ); return 0.0; };
      void atPutDiagonal( int, OoqpVector& ) override { assert( "not implemented" && 0 ); };
      void fromGetDiagonal( int, OoqpVector& ) override { assert( "not implemented" && 0 ); };
      void fromGetDense( int, int, double*, int, int, int ) override { assert( "not implemented" && 0 ); };
      void fromGetSpRow( int, int, double[], int, int[], int&, int, int& ) override { assert( "not implemented" && 0 ); };
      void getDiagonal( OoqpVector& ) override { assert( "not implemented" && 0 ); };
      void setToDiagonal( const OoqpVector& ) override { assert( "not implemented" && 0 ); };
      void matTransDMultMat(OoqpVector&, SymMatrix** ) override { assert( "not implemented" && 0 ); };
      void matTransDinvMultMat(OoqpVector&, SymMatrix** ) override { assert( "not implemented" && 0 ); };
      void symmetricScale ( const OoqpVector& ) override { assert( "not implemented" && 0 ); };
      void writeToStream( std::ostream& ) const override { assert( "not implemented" && 0 ); };
      void atPutSubmatrix( int, int, DoubleMatrix&, int, int, int, int ) override { assert( "not implemented" && 0 ); };
      void atPutDense( int, int, double*, int, int, int ) override { assert( "not implemented" && 0 ); };
      void atPutSpRow( int, double[], int, int[], int& ) override { assert( "not implemented" && 0 ); };
      void putSparseTriple( int[], int, int[], double[], int& ) override { assert( "not implemented" && 0 ); };
      void randomize( double, double, double* ) override { assert( "not implemented" && 0 ); };

   protected:
      virtual void multVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const;
      virtual void multHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x, bool root) const;

      virtual void transMultVertical( double beta, OoqpVector& y, double alpha, const OoqpVector& x, bool root ) const;
      virtual void transMultHorizontal( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const;

      virtual void columnScaleVertical( const OoqpVector& vec );
      virtual void columnScaleHorizontal( const OoqpVector& vec );

      virtual void rowScaleVertical( const OoqpVector& vec );
      virtual void rowScaleHorizontal( const OoqpVector& vec );

      virtual void getRowMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* col_scale, OoqpVector& minmax) const;
      virtual void getRowMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* col_scale, OoqpVector& minmax) const;

      virtual void getColMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* row_scale, OoqpVector& minmax) const;
      virtual void getColMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* row_scale, OoqpVector& minmax) const;

      virtual void addRowSumsVertical( OoqpVector& vec ) const;
      virtual void addRowSumsHorizontal( OoqpVector& vec ) const;

      virtual void addColSumsVertical( OoqpVector& vec ) const;
      virtual void addColSumsHorizontal( OoqpVector& vec ) const;
};

/**
 * Dummy version ...
 */
class StringGenDummyMatrix : public StringGenMatrix
{
   public:
      StringGenDummyMatrix() = default;
      ~StringGenDummyMatrix() override = default;

      void addChild( StringGenMatrix* ) override {};

      bool isEmpty() const override { return true; };
      int isKindOf( int type ) const override { return type == kStringGenDummyMatrix || type == kStringMatrix || type == kStringGenMatrix; };
      void mult( double, OoqpVector&, double, const OoqpVector& ) const override {};
      void transMult( double, OoqpVector&, double, const OoqpVector& ) const override {};
      double abmaxnorm() const override { return -std::numeric_limits<double>::infinity(); };
      void scalarMult( double ) override {};
      void writeToStream( std::ostream& ) const override {};
      void writeToStreamDenseRow( std::ostream&, int ) const override {};
      void writeDashedLineToStream( std::ostream& ) const override {};

      void getRowMinMaxVec( bool, bool, const OoqpVector*, OoqpVector& ) override {};
      void getColMinMaxVec( bool, bool, const OoqpVector*, OoqpVector& ) override {};
      void columnScale( const OoqpVector& ) override {};
      void rowScale( const OoqpVector& ) override {};

      void addRowSums( OoqpVector& ) const override {};
      void addColSums( OoqpVector& ) const override {};

      GenMatrix* shaveBottom( int ) override { return new StringGenDummyMatrix(); };

      void splitAlongTree( const sTreeCallbacks& ) override { assert( false && "should not end up here" ); };

   protected:
      void multVertical( double, OoqpVector&, double, const OoqpVector& ) const override {};
      void multHorizontal( double, OoqpVector&, double, const OoqpVector&, bool ) const override {};

      void transMultVertical( double, OoqpVector&, double, const OoqpVector&, bool ) const override {};
      void transMultHorizontal( double, OoqpVector&, double, const OoqpVector& ) const override {};

      void columnScaleVertical( const OoqpVector& ) override {};
      void columnScaleHorizontal( const OoqpVector& ) override {};

      void rowScaleVertical( const OoqpVector& ) override {};
      void rowScaleHorizontal( const OoqpVector& ) override {};

      void getRowMinMaxVecVertical( bool, bool, const OoqpVector*, OoqpVector& ) const override {};
      void getRowMinMaxVecHorizontal( bool, bool, const OoqpVector*, OoqpVector& ) const override {};

      void getColMinMaxVecVertical( bool, bool, const OoqpVector*, OoqpVector& ) const override {};
      void getColMinMaxVecHorizontal( bool, bool, const OoqpVector*, OoqpVector& ) const override {};

      void addRowSumsVertical( OoqpVector& ) const override {};
      void addRowSumsHorizontal( OoqpVector& ) const override {};

      void addColSumsVertical( OoqpVector& ) const override {};
      void addColSumsHorizontal( OoqpVector& ) const override {};
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_STRINGGENMATRIX_H_ */
