/*
 * BorderedSymMatrix.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#include "BorderedSymMatrix.h"
#include "DoubleMatrixTypes.h"
#include "SimpleVector.h"
#include "pipsdef.h"
#include <algorithm>

BorderedSymMatrix::BorderedSymMatrix(StochSymMatrix* inner_matrix_, StringGenMatrix* border_vertical_, SymMatrix* top_left_block_,
            MPI_Comm mpiComm_) : inner_matrix(inner_matrix_), border_vertical(border_vertical_), top_left_block(top_left_block_), mpiComm( mpiComm_ ),
            iAmDistrib( mpiComm == MPI_COMM_NULL )
{
   assert( inner_matrix );
   assert( border_vertical );
   assert( top_left_block );

   assert( inner_matrix->children.size() == border_vertical->children.size() );
   assert( border_vertical->is_vertical );

   inner_matrix->getSize(n, n);

   int n_border, m_border;
   border_vertical->getSize(m_border, n_border);

   n += n_border;

#ifndef NDEBUG
   int n_bottom;
   top_left_block->getSize(n_bottom, n_bottom);
   int n_inner;
   inner_matrix->getSize(n_inner, n_inner);

   assert( n_inner == m_border );
   assert( n_bottom == n_border );
#endif
}

BorderedSymMatrix::~BorderedSymMatrix()
{
   delete inner_matrix;
   delete border_vertical;
   delete top_left_block;
}

int BorderedSymMatrix::isKindOf( int type ) const
{
   return (type == kBorderedSymMatrix || type == kSymMatrix || type == kBorderedMatrix);
}

/** y = beta * y + alpha * this * x */
void BorderedSymMatrix::mult( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{
   DistributedVector<double>& y = dynamic_cast<DistributedVector<double>&>(y_in);
   const DistributedVector<double>& x = dynamic_cast<const DistributedVector<double>&>(x_in);

   assert(x.children.size() == 1 && y.children.size() == 1);
   assert(x.children[0] && y.children[0]);
   assert( x.first );
   assert( y.first );
   assert( !x.last );
   assert( !y.last );

   top_left_block->mult(beta, *y.first, alpha, *x.first );
   border_vertical->transMult(1.0, *y.first, alpha, *x.children[0] );

   border_vertical->mult( beta, *y.children[0], alpha, *x.first );
   inner_matrix->mult( 1.0, *y.children[0], alpha, *x.children[0] );
}

/** y = beta * y + alpha * this^T * x */
void BorderedSymMatrix::transMult( double beta, OoqpVector& y, double alpha,  const OoqpVector& x ) const
{
   this->mult( beta, y, alpha, x );
}

void BorderedSymMatrix::fromGetDiagonal( int idiag, OoqpVector& x_in )
{
   assert( "The value of the parameter idiag is not supported!" && idiag == 0 );

   DistributedVector<double>& x = dynamic_cast<DistributedVector<double>&>(x_in);
   assert( x.children.size() == 1 );
   assert( x.children[0] );
   assert(!x.last && x.first );

   top_left_block->getDiagonal(*x.first);

   inner_matrix->fromGetDiagonal(idiag, *x.children[0]);
}

double BorderedSymMatrix::abmaxnorm() const
{
   double norm = -std::numeric_limits<double>::infinity();

   norm = std::max(norm, inner_matrix->abmaxnorm());
   norm = std::max(norm, border_vertical->abmaxnorm());
   norm = std::max(norm, top_left_block->abmaxnorm());

   return norm;
}

void BorderedSymMatrix::scalarMult( double num )
{
   inner_matrix->scalarMult( num );
   border_vertical->scalarMult( num );
}

void BorderedSymMatrix::getSize( long long& m_, long long& n_ ) const
{
   m_ = n;
   n_ = n;
}

void BorderedSymMatrix::getSize( int& m_, int& n_ ) const
{
   m_ = n;
   n_ = n;
}

long long BorderedSymMatrix::size() const
{
   return n;
}
