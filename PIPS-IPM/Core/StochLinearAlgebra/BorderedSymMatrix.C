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

BorderedSymMatrix::BorderedSymMatrix(int id_, long long n_, StochSymMatrix* inner_matrix_, StringGenMatrix* border_vertical_, SymMatrix* bottom_right_block_,
            MPI_Comm mpiComm_) : inner_matrix(inner_matrix_), border_vertical(border_vertical_), bottom_right_block(bottom_right_block_), id(id_), n(n_), mpiComm( mpiComm_ ),
            iAmDistrib( mpiComm == MPI_COMM_NULL )
{
   assert( inner_matrix );
   assert( border_vertical );

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
   delete bottom_right_block;
}

int BorderedSymMatrix::isKindOf( int type ) const
{
   return (type == kBorderedSymMatrix || type == kSymMatrix || type == kBorderedMatrix);

}

/** y = beta * y + alpha * this * x */
void BorderedSymMatrix::mult( double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const
{
   assert( 0 && " TODO : implement ..." );
}

/** y = beta * y + alpha * this^T * x */
void BorderedSymMatrix::transMult( double beta, OoqpVector& y, double alpha,  const OoqpVector& x ) const
{
   assert( 0 && " TODO : implement ..." );
}

void BorderedSymMatrix::fromGetDiagonal( int idiag, OoqpVector& x )
{
   assert( 0 && " TODO : implement ..." );
}

double BorderedSymMatrix::abmaxnorm() const
{
   assert( 0 && " TODO : implement ..." );
}

void BorderedSymMatrix::scalarMult( double num )
{
   assert( 0 && " TODO : implement ..." );
}

void BorderedSymMatrix::getSize( long long& m, long long& n ) const
{
   assert( 0 && " TODO : implement ..." );
}

void BorderedSymMatrix::getSize( int& m, int& n ) const
{
   assert( 0 && " TODO : implement ..." );
}

long long BorderedSymMatrix::size() const
{
   assert( 0 && " TODO : implement ..." );
}
