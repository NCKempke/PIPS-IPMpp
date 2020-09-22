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

BorderedSymMatrix::BorderedSymMatrix(int id_, long long n_, StochSymMatrix* inner_matrix_, StringGenMatrix* border_vertical_, SymMatrix* bottom_block_,
            MPI_Comm mpiComm_) : inner_matrix(inner_matrix_), border_vertical(border_vertical_), bottom_block(bottom_block_), id(id_), n(n_), mpiComm( mpiComm_ ),
            iAmDistrib( mpiComm == MPI_COMM_NULL )
{
   assert( inner_matrix );
   assert( border_vertical );
   assert( bottom_block );

   assert(inner_matrix->children.size() == border_vertical->children.size() );
   assert( border_vertical->is_vertical );

#ifndef NDEBUG
   int n, m;
   border_vertical->getSize(m, n);

   assert( n == inner_matrix->n + n);
#endif
}

BorderedSymMatrix::~BorderedSymMatrix()
{
   delete inner_matrix;
   delete border_vertical;
   delete bottom_block;
}

int BorderedSymMatrix::isKindOf( int type ) const
{
   return (type == kBorderedSymMatrix || type == kSymMatrix || type == kBorderedMatrix);
}

/** y = beta * y + alpha * this * x */
void BorderedSymMatrix::mult( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{
   StochVector& y = dynamic_cast<StochVector&>(y_in);
   const StochVector& x = dynamic_cast<const StochVector&>(x_in);

   assert(x.children.size() == 1 && y.children.size() == 1);
   assert(x.children[0] && y.children[0]);
   assert( y.vecl && !x.vecl);
   assert( !y.vec && x.vec );

   border_vertical->mult( beta, *y.children[0], alpha, *x.vec );
   inner_matrix->mult( 1.0, *y.children[0], alpha, *x.children[0] );

   bottom_block->mult( beta, *y.vecl, alpha, *x.vec );
   border_vertical->transMult( 1.0, *y.vecl, alpha, *x.children[0] );
}

/** y = beta * y + alpha * this^T * x */
void BorderedSymMatrix::transMult( double beta, OoqpVector& y, double alpha,  const OoqpVector& x ) const
{
   this->mult( beta, y, alpha, x );
}

void BorderedSymMatrix::fromGetDiagonal( int idiag, OoqpVector& x_in )
{
   assert( "The value of the parameter idiag is not supported!" && idiag == 0 );

   StochVector& x = dynamic_cast<StochVector&>(x_in);
   assert( x.children.size() == 1 );
   assert( x.children[0] );
   assert( !x.vecl && x.vec );

   bottom_block->getDiagonal(*x.vec);

   inner_matrix->fromGetDiagonal(idiag, *x.children[0]);
}

double BorderedSymMatrix::abmaxnorm() const
{
   double norm = -std::numeric_limits<double>::infinity();

   norm = std::max(norm, inner_matrix->abmaxnorm());
   norm = std::max(norm, border_vertical->abmaxnorm());
   norm = std::max(norm, bottom_block->abmaxnorm());

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
