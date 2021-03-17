/*
 * BorderedGenMatrix.C
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */


#include "BorderedGenMatrix.h"
#include "StochVector_fwd.h"
#include "OoqpVector_fwd.h"
#include "StochGenMatrix.h"
#include "SparseGenMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixTypes.h"

#include "pipsdef.h"

#include <vector>
#include <cassert>
#include "StringGenMatrix.h"

BorderedGenMatrix::BorderedGenMatrix(StochGenMatrix* inner_matrix_, StringGenMatrix* border_left_,
            StringGenMatrix* border_bottom_, SparseGenMatrix* bottom_left_block_, MPI_Comm mpi_comm_ ) :
            inner_matrix{inner_matrix_}, border_left{border_left_}, border_bottom{border_bottom_}, bottom_left_block{bottom_left_block_},
            mpi_comm(mpi_comm_), distributed( mpi_comm == MPI_COMM_NULL ), rank( PIPS_MPIgetRank(mpi_comm) )
{
   assert( inner_matrix );
   assert( border_left );
   assert( border_bottom );
   assert( bottom_left_block );

   bottom_left_block->getSize(m, n);

   int m_bottom, n_bottom;
   border_bottom->getSize(m_bottom, n_bottom);

   int m_left, n_left;
   border_left->getSize(m_left, n_left);

   int m_inner, n_inner;
   inner_matrix->getSize(m_inner, n_inner );

   assert( n == n_left );
   assert( m == m_bottom );

   assert( m_inner == m_left );
   assert( n_inner == n_bottom );

   m += m_left;
   n += n_bottom;

}

BorderedGenMatrix::~BorderedGenMatrix()
{
   delete bottom_left_block;
   delete border_bottom;
   delete border_left;
   delete inner_matrix;
}

int BorderedGenMatrix::isKindOf( int type ) const
{
   return type == kBorderedGenMatrix || type == kBorderedMatrix || type == kGenMatrix;
}

void BorderedGenMatrix::mult( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{
   /* x row, y column shaped */
   assert( hasVecStructureForBorderedMat(x_in, true) );
   assert( hasVecStructureForBorderedMat(y_in, false) );

   const StochVector& x = dynamic_cast<const StochVector&>(x_in);
   StochVector& y = dynamic_cast<StochVector&>(y_in);

   border_left->mult(beta, *y.children[0], alpha, *x.vec);
   inner_matrix->mult(1.0, *y.children[0], alpha, *x.children[0]);

   bottom_left_block->mult(beta, *y.vecl, alpha, *x.vec);
   border_bottom->mult(1.0, *y.vecl, alpha, *x.children[0]);
}

/** y = beta * y + alpha * this^T * x */
void BorderedGenMatrix::transMult( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{
   /* x column, y row shaped */
   assert( hasVecStructureForBorderedMat(x_in, false) );
   assert( hasVecStructureForBorderedMat(y_in, true) );

   const StochVector& x = dynamic_cast<const StochVector&>(x_in);
   StochVector& y = dynamic_cast<StochVector&>(y_in);

   border_left->transMult(beta, *y.vec, alpha, *x.children[0]);
   bottom_left_block->transMult(1.0, *y.vec, alpha, *x.vecl);

   inner_matrix->transMult(beta, *y.children[0], alpha, *x.children[0]);
   border_bottom->transMult(1.0, *y.children[0], alpha, *x.vecl);
}

double BorderedGenMatrix::abmaxnorm() const
{
   double norm = -std::numeric_limits<double>::infinity();

   norm = std::max(norm, inner_matrix->abmaxnorm());
   norm = std::max(norm, border_left->abmaxnorm());
   norm = std::max(norm, border_bottom->abmaxnorm());
   norm = std::max(norm, bottom_left_block->abmaxnorm());

   return norm;
}

void BorderedGenMatrix::columnScale( const OoqpVector& vec )
{
   assert( hasVecStructureForBorderedMat(vec, true) );

   const StochVector& svec = dynamic_cast<const StochVector&>(vec);

   border_left->columnScale(*svec.vec);
   bottom_left_block->columnScale(*svec.vec);

   inner_matrix->columnScale(*svec.children[0]);
   border_bottom->columnScale(*svec.children[0]);
}

void BorderedGenMatrix::rowScale ( const OoqpVector& vec )
{
   assert( hasVecStructureForBorderedMat(vec, false) );

   const StochVector& svec = dynamic_cast<const StochVector&>(vec);

   border_left->rowScale(*svec.children[0]);
   inner_matrix->rowScale(*svec.children[0]);

   bottom_left_block->rowScale(*svec.vecl);
   border_bottom->rowScale(*svec.vecl);
}

void BorderedGenMatrix::scalarMult( double num )
{
   inner_matrix->scalarMult(num);
   border_left->scalarMult(num);
   border_bottom->scalarMult(num);
   bottom_left_block->scalarMult(num);
}

void BorderedGenMatrix::getSize( long long& m_, long long& n_ ) const
{
   m_ = m;
   n_ = n;
}

void BorderedGenMatrix::getSize( int& m_, int& n_ ) const
{
   m_ = m;
   n_ = n;
}

void BorderedGenMatrix::getRowMinMaxVec( bool get_min, bool initialize_vec, const OoqpVector* col_scale_in, OoqpVector& minmax_in )
{
   assert( hasVecStructureForBorderedMat(minmax_in, false) );
   const bool has_colscale = (col_scale_in != nullptr);
   if( has_colscale )
      assert( hasVecStructureForBorderedMat(*col_scale_in, true) );

   StochVector& minmax = dynamic_cast<StochVector&>(minmax_in);
   const StochVector* col_scale = has_colscale ? dynamic_cast<const StochVector*>(col_scale_in) : nullptr;

   border_left->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->vec : nullptr, *minmax.children[0]);
   inner_matrix->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->children[0] : nullptr, *minmax.children[0]);

   bottom_left_block->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->vec : nullptr, *minmax.vecl);
   border_bottom->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->children[0] : nullptr, *minmax.vecl);
}

void BorderedGenMatrix::getColMinMaxVec( bool get_min, bool initialize_vec, const OoqpVector* row_scale_in, OoqpVector& minmax_in )
{
   assert( hasVecStructureForBorderedMat(minmax_in, true) );

   const bool has_rowscale = (row_scale_in != nullptr);
   if( has_rowscale )
      assert( hasVecStructureForBorderedMat(*row_scale_in, false) );

   StochVector& minmax = dynamic_cast<StochVector&>(minmax_in);
   const StochVector* row_scale = has_rowscale ? dynamic_cast<const StochVector*>(row_scale_in) : nullptr;

   border_left->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->children[0] : nullptr, *minmax.vec);
   bottom_left_block->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->vecl : nullptr, *minmax.vec);

   inner_matrix->getColMinMaxVec(get_min, initialize_vec, has_rowscale ? row_scale->children[0] : nullptr, *minmax.children[0]);
   border_bottom->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->vecl : nullptr, *minmax.children[0]);
}

void BorderedGenMatrix::addRowSums( OoqpVector& vec_ ) const
{
  assert( hasVecStructureForBorderedMat( vec_, false ) );

  StochVector& vec = dynamic_cast<StochVector&>(vec_);

  border_left->addRowSums( *vec.children[0] );
  inner_matrix->addRowSums( *vec.children[0] );

  bottom_left_block->addRowSums( *vec.vecl );
  border_bottom->addRowSums( *vec.vecl );
}

void BorderedGenMatrix::addColSums( OoqpVector& vec_ ) const
{
   assert( hasVecStructureForBorderedMat( vec_, true ) );

   StochVector& vec = dynamic_cast<StochVector&>(vec_);

   border_left->addColSums( *vec.vec );
   bottom_left_block->addColSums( *vec.vec );

   inner_matrix->addColSums( *vec.children[0] );
   border_bottom->addColSums( *vec.children[0] );
}

template<typename T>
bool BorderedGenMatrix::hasVecStructureForBorderedMat( const OoqpVectorBase<T>& vec, bool row_vec ) const
{
   const StochVectorBase<T>& vecs = dynamic_cast<const StochVectorBase<T>&>(vec);

   if( vecs.children.size() != 1 )
   {
      std::cout << "children" << std::endl;
      return false;
   }

   if( vecs.children[0] == nullptr )
   {
      std::cout << "child[0]" << std::endl;
      return false;
   }

   if( row_vec )
   {
      if( vecs.vecl != nullptr )
      {
         std::cout << "row-vec but root.vecl" << std::endl;
         return false;
      }
      if( vecs.vec == nullptr )
      {
         std::cout << "row-vec but NO root.vec" << std::endl;
         return false;
      }
   }
   else
   {
      if( vecs.vec != nullptr )
      {
         std::cout << "col-vec but root.vec" << std::endl;
         return false;
      }
      if( vecs.vecl == nullptr )
      {
         std::cout << "col-vec but NO root.vecl" << std::endl;
         return false;
      }

   }

   if( row_vec && vecs.length() != n )
   {
      std::cout << "ROW: root.length = " << vecs.length() << " != " << n << " = border.n " << std::endl;
      return false;
   }
   if( !row_vec && vecs.length() != m )
   {
      std::cout << "COL: root.length = " << vecs.length() << " != " << m << " = border.m " << std::endl;
      return false;
   }

   int n_border, m_border;
   border_left->getSize(m_border, n_border);
   if( row_vec && vecs.vec->length() != n_border )
   {
      std::cout << "ROW: root.vec.length = " << vecs.vec->length() << " != " << n_border << " = border.n " << std::endl;
      return false;
   }
   if( !row_vec && vecs.length() != m )
   {
      std::cout << "COL: root.vecl.length = " << vecs.vecl->length() << " != " << m_border << " = border.m " << std::endl;
      return false;
   }

   return true;
}

void BorderedGenMatrix::writeToStreamDense( std::ostream& out ) const
{
   const int my_rank = PIPS_MPIgetRank(mpi_comm);
   const int size = PIPS_MPIgetSize(mpi_comm);
   const bool iAmDistrib = ( size != 0 );

   inner_matrix->writeToStreamDenseBordered( *border_left, out );

   int mL, nL; this->bottom_left_block->getSize(mL, nL);
   if( mL > 0 )
   {
      if( iAmDistrib )
         MPI_Barrier(mpi_comm);

      // for each row r do:
      for( int r = 0; r < mL; r++ )
      {
         MPI_Barrier(mpi_comm);

         // process Zero collects all the information and then prints it.
         if( my_rank == 0 )
         {
            bottom_left_block->writeToStreamDenseRow(out, r);
            out << "|\t";
         }
         border_bottom->writeToStreamDenseRow(out, r);

         if( my_rank == 0 )
            out << "\n";

         MPI_Barrier(mpi_comm);
      }
   }
   if( iAmDistrib )
      MPI_Barrier(mpi_comm);
}

