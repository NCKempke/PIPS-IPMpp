/*
 * StringGenMatrix.C
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */

#include "StringGenMatrix.h"

#include "StochVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"

#include "pipsdef.h"
#include <algorithm>

StringGenMatrix::StringGenMatrix(bool is_vertical, SparseGenMatrix* mat, SparseGenMatrix* mat_link, MPI_Comm mpi_comm_)
   : mat(mat), mat_link(mat_link), is_vertical(is_vertical), mpi_comm(mpi_comm_), distributed(mpi_comm == MPI_COMM_NULL), rank(PIPS_MPIgetRank(mpi_comm))
{
   assert(mat);
}

StringGenMatrix::~StringGenMatrix()
{
   for( auto& child : children )
      delete child;

   delete mat;
   delete mat_link;
}

void StringGenMatrix::addChild(StringGenMatrix* child)
{
   children.push_back(child);
}

int StringGenMatrix::isKindOf( int type ) const
{
   return (type == kStringGenMatrix || type == kStringMatrix || type == kGenMatrix);
}

double StringGenMatrix::abmaxnorm() const
{
   double norm = 0.0;

   for( size_t it = 0; it < children.size(); it++ )
      norm = std::max(norm, children[it]->abmaxnorm());

   if( distributed )
      norm = PIPS_MPIgetMax(norm, mpi_comm);

   norm = std::max(norm, mat->abmaxnorm());

   if( mat_link )
      norm = std::max(norm, mat_link->abmaxnorm());

   return norm;
}

void StringGenMatrix::scalarMult( double num )
{
   mat->scalarMult(num);

   if( mat_link )
      mat_link->scalarMult(num);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->scalarMult(num);
}


void StringGenMatrix::mult(double beta, OoqpVector& y, double alpha, const OoqpVector& x ) const
{
   if( is_vertical )
      multVertical(beta, y, alpha, x);
   else
      multHorizontal(beta, y, alpha, x);
}

void StringGenMatrix::transMult(double beta, OoqpVector& y, double alpha, const OoqpVector& x) const
{
   if( is_vertical )
      transMultVertical(beta, y, alpha, x);
   else
      transMultHorizontal(beta, y, alpha, x);
}

void StringGenMatrix::multVertical( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in) const
{
   const SimpleVector& x = dynamic_cast<const SimpleVector&>(x_in);
   StochVector& y = dynamic_cast<StochVector&>(y_in);

   assert(y.children.size() == children.size());
   assert(is_vertical);
   if( y.vecl )
      assert( mat_link );
   else
      assert( mat_link == nullptr );

   mat->mult(beta, *y.vec, alpha, x);

   /* linking matrix present? */
   if( mat_link )
      mat_link->mult(beta, *y.vecl, alpha, x);

   for( size_t i = 0; i < children.size(); ++i )
   {
      if( children[i]->isKindOf(kStringGenDummyMatrix) )
         assert( y.children[i]->isKindOf(kStochDummy) );

      children[i]->multVertical(beta, *y.children[i], alpha, x);
   }
}

void StringGenMatrix::multHorizontal( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{
   const StochVector& x = dynamic_cast<const StochVector&>(x_in);
   SimpleVector& y = dynamic_cast<SimpleVector&>(y_in);

   assert( !is_vertical );
   assert( x.children.size() == children.size() );
   if( children.size() == 0 )
   {
      assert( !distributed );
      assert( rank == 0 );
      assert( mpi_comm == MPI_COMM_NULL );
   }
   assert( ( x.vecl && mat_link ) || ( x.vecl == nullptr && mat_link == nullptr ) );

   if( rank == 0 )
      mat->mult(beta, y, alpha, *x.vec);
   else
      y.setToZero();

   for( size_t i = 0; i < children.size(); ++i )
      children[i]->multHorizontal(1.0, y, alpha, *x.children[i]);

   if( distributed )
      PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);

   if( mat_link )
      mat->mult(1.0, y, alpha, *x.vecl);
}

void StringGenMatrix::transMultVertical ( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{

   const StochVector& x = dynamic_cast<const StochVector&>(x_in);
   SimpleVector& y = dynamic_cast<SimpleVector&>(y_in);

   assert( is_vertical );
   assert(x.children.size() == children.size());
   if( children.size() == 0 )
   {
      assert( !distributed );
      assert( rank == 0 );
      assert( mpi_comm == MPI_COMM_NULL );
   }
   assert( ( x.vecl && mat_link ) || ( x.vecl == nullptr && mat_link == nullptr ) );

   if( rank == 0 )
      mat->transMult(beta, y, alpha, x);
   else
      y.setToZero();

   for( size_t i = 0; i < children.size(); i++ )
      children[i]->transMultVertical(1.0, y, alpha, *x.children[i]);

   if( distributed )
      PIPS_MPIsumArrayInPlace(y.elements(), y.length(), mpi_comm);

   if( mat_link )
      mat_link->transMult(1.0, y, alpha, *x.vecl);
}

void StringGenMatrix::transMultHorizontal ( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{

   const SimpleVector& x = dynamic_cast<const SimpleVector&>(x_in);
   StochVector& y = dynamic_cast<StochVector&>(y_in);

   assert( !is_vertical );
   assert(y.children.size() == children.size());
   if( y.vecl )
      assert( mat_link );
   else
      assert( mat_link == nullptr );

   mat->transMult(beta, *y.vec, alpha, x);

   if( mat_link )
      mat_link->transMult( beta, *y.vecl, alpha, x);

   for( size_t i = 0; i < children.size(); ++i )
   {
      if( children[i]->isKindOf(kStringGenDummyMatrix) )
         assert( y.children[i]->isKindOf(kStochDummy) );

      children[i]->transMultHorizontal(beta, *y.children[i], alpha, x);
   }
}
