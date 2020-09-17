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
#include "StochVector.h"
#include "DoubleMatrixTypes.h"

#include "pipsdef.h"

#include <vector>
#include <cassert>
#include "StringGenMatrix.h"

BorderedGenMatrix::BorderedGenMatrix(StochGenMatrix* inner_matrix, StringGenMatrix* border_left,
            StringGenMatrix* border_bottom, StochGenMatrix* bottom_right_block, MPI_Comm mpi_comm_) :
            inner_matrix(inner_matrix), border_left(border_left), border_bottom(border_bottom), bottom_left_block(bottom_right_block),
            mpi_comm(mpi_comm_), distributed( mpi_comm == MPI_COMM_NULL ), rank( PIPS_MPIgetRank(mpi_comm) ),
            m(bottom_right_block->m), n(bottom_right_block->n)
{
   assert( inner_matrix );
   assert( border_left );
   assert( border_bottom );
   assert( bottom_right_block );
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
   assert( hasVecStructureForBorderedMat(y_in, false) );
   assert( hasVecStructureForBorderedMat(x_in, true) );

   const StochVector& x = dynamic_cast<const StochVector&>(x_in);
   StochVector& y = dynamic_cast<StochVector&>(y_in);

   border_left->mult(beta, *y.children[0], alpha, *x.vec);
   inner_matrix->mult(1.0, *y.children[0], alpha, *x.children[0]);
   bottom_left_block->mult(beta, *y.vecl, alpha, *x.vec);
   border_bottom->mult(1.0, *y.vecl, alpha, *x.children[0]);
}

void BorderedGenMatrix::transMult ( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const
{
   assert( 0 && "todo: implement");
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

void BorderedGenMatrix::columnScale ( const OoqpVector& vec )
{
   assert( 0 && "todo: implement");
}

void BorderedGenMatrix::rowScale ( const OoqpVector& vec )
{
   assert( 0 && "todo: implement");
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

void BorderedGenMatrix::getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec )
{
   assert( 0 && "todo: implement");
}

void BorderedGenMatrix::getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec )
{
   assert( 0 && "todo: implement");
}

void BorderedGenMatrix::getNnzPerRow(OoqpVectorBase<int>& nnzVec)
{
   assert( 0 && "todo: implement");
}

void BorderedGenMatrix::getNnzPerCol(OoqpVectorBase<int>& nnzVec)
{
   assert( 0 && "todo: implement");
}

void BorderedGenMatrix::addRowSums( OoqpVector& vec )
{
   assert( 0 && "todo: implement");
}

void BorderedGenMatrix::addColSums( OoqpVector& vec )
{
   assert( 0 && "todo: implement");
}

bool BorderedGenMatrix::hasVecStructureForBorderedMat( const OoqpVector& vec, bool row_vec ) const
{
   const StochVector& vecs = dynamic_cast<const StochVector&>(vec);

   if( vecs.children.size() != 1 )
      return false;

   if( vecs.children[0] == nullptr )
      return false;

   if( row_vec )
   {
      if( vecs.vecl != nullptr )
         return false;
      if( vecs.vec == nullptr )
         return false;
   }
   else
   {
      if( vecs.vec != nullptr )
         return false;
      if( vecs.vecl == nullptr )
         return false;

   }

   if( row_vec && vecs.vec->length() != n )
      return false;
   if( !row_vec && vecs.vecl->length() != m )
      return false;

   return true;
}
