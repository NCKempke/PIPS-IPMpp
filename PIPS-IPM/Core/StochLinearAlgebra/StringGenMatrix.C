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

StringGenMatrix::StringGenMatrix() : mat(nullptr), mat_link(nullptr), is_vertical(false), m(0), n(0), id(-1), mpi_comm(MPI_COMM_NULL), distributed(false), rank(-1)
{
}

StringGenMatrix::StringGenMatrix(int id_, bool is_vertical, SparseGenMatrix* mat, SparseGenMatrix* mat_link, MPI_Comm mpi_comm_)
   : mat(mat), mat_link(mat_link), is_vertical(is_vertical), id(id_), mpi_comm(mpi_comm_), distributed( PIPS_MPIgetDistributed(mpi_comm) ), rank( PIPS_MPIgetRank(mpi_comm) )
{
   assert(mat);

   mat->getSize(m, n);

   if( mat_link )
   {
      long long ml, nl;
      mat_link->getSize(ml, nl);

      if( is_vertical )
      {
         assert( n == nl );
         m += ml;
      }
      else
      {
         assert( m == ml );
         n += nl;
      }
   }
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

   long long m_, n_;
   child->getSize(m_, n_);

   assert( child->is_vertical == this->is_vertical || child->isKindOf( kStringGenDummyMatrix ) );

   if( !child->isKindOf( kStringGenDummyMatrix ) )
   {
      if( is_vertical )
      {
         assert( n == n_ );
         m += m_;
      }
      else
      {
         assert( m == m_ );
         n += n_;
      }
   }
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

   assert(is_vertical);
   assert( y.children.size() == children.size() );
   if( y.vecl )
      assert( mat_link );
   else
      assert( mat_link == nullptr );

   mat->mult(beta, *y.vec, alpha, x);

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

void StringGenMatrix::transMultVertical( double beta, OoqpVector& y_in, double alpha, const OoqpVector& x_in ) const
{
   const StochVector& x = dynamic_cast<const StochVector&>(x_in);
   SimpleVector& y = dynamic_cast<SimpleVector&>(y_in);

   assert( is_vertical );
   assert(x.children.size() == children.size());
   if( children.size() == 0 )
   {
      assert( !distributed );
      assert( rank == 0 );
   }

   assert( x.vec && mat );
   assert( ( x.vecl && mat_link ) || ( x.vecl == nullptr && mat_link == nullptr ) );

   if( rank == 0 )
      mat->transMult(beta, y, alpha, *x.vec);
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

void StringGenMatrix::getColMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* row_scale, OoqpVector& minmax_in ) const
{
   assert( !is_vertical );
   StochVector& minmax = dynamic_cast<StochVector&>(minmax_in);

   assert( minmax.vec && mat );
   assert( minmax.children.size() == children.size());
   assert( (minmax.vecl && mat_link) || (minmax.vecl == nullptr && mat_link == nullptr) );

   mat->getColMinMaxVec(get_min, initialize_vec, row_scale, *minmax.vec);

   for( size_t i = 0; i < children.size(); ++i )
   {
      assert( minmax.children[i]);
      if( children[i]->isKindOf(kStringGenDummyMatrix) )
         assert( minmax.children[i]->isKindOf(kStochDummy) );

      children[i]->getColMinMaxVecHorizontal(get_min, initialize_vec, row_scale, *minmax.children[i]);
   }

   if( mat_link )
      mat_link->getColMinMaxVec(get_min, initialize_vec, row_scale, *minmax.vecl);
}

void StringGenMatrix::getColMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* row_scale_in, OoqpVector& minmax ) const
{
   assert( is_vertical );
   const bool has_rowscale = (row_scale_in != nullptr);

   const StochVector* row_scale = dynamic_cast<const StochVector*>(row_scale_in);
   assert( !has_rowscale || row_scale->children.size() == children.size() );
   if( has_rowscale )
      assert( (row_scale->vecl && mat_link) || ( row_scale->vecl == nullptr && mat_link == nullptr) );

   mat->getColMinMaxVec( get_min, initialize_vec, has_rowscale ? row_scale->vec : nullptr, minmax );

   for( size_t i = 0; i < children.size(); i++ )
   {
      if( has_rowscale )
         if( children[i]->isKindOf(kStringGenDummyMatrix) )
            assert( row_scale->children[i]->isKindOf(kStochDummy) );

      children[i]->getColMinMaxVecVertical(get_min, false, has_rowscale ? row_scale->children[i] : nullptr, minmax);
   }

   if( mat_link )
      mat_link->getColMinMaxVec(get_min, false, has_rowscale ? row_scale->vecl : nullptr, minmax);
}

/** StochVector colScaleVec, SimpleVector minmaxVec */
void StringGenMatrix::getRowMinMaxVecHorizontal( bool get_min, bool initialize_vec, const OoqpVector* col_scale_in, OoqpVector& minmax) const
{
   assert( !is_vertical );
   const bool has_colscale = (col_scale_in != nullptr);

   const StochVector* col_scale = dynamic_cast<const StochVector*>(col_scale_in);
   assert( !has_colscale || col_scale->children.size() == children.size() );
   if( has_colscale )
      assert( (col_scale->vecl && mat_link) || ( col_scale->vecl == nullptr && mat_link == nullptr) );

   mat->getRowMinMaxVec(get_min, initialize_vec, has_colscale ? col_scale->vec : nullptr, minmax);

   for( size_t i = 0; i < children.size(); i++ )
   {
      if( has_colscale )
         if( children[i]->isKindOf(kStringGenDummyMatrix) )
            assert( col_scale->children[i]->isKindOf(kStochDummy) );

      children[i]->getRowMinMaxVecHorizontal(get_min, false, has_colscale ? col_scale->children[i] : nullptr, minmax);
   }

   if( mat_link )
      mat_link->getRowMinMaxVec(get_min, false, has_colscale ? col_scale->vecl : nullptr, minmax);
}

/** StochVector minmaxVec, SimpleVector colScaleVec */
void StringGenMatrix::getRowMinMaxVecVertical( bool get_min, bool initialize_vec, const OoqpVector* col_scale, OoqpVector& minmax_in ) const
{
   assert( is_vertical );

   StochVector& minmax = dynamic_cast<StochVector&>(minmax_in);

   assert( minmax.vec && mat );
   assert( minmax.children.size() == children.size());
   assert( (minmax.vecl && mat_link) || (minmax.vecl == nullptr && mat_link == nullptr) );

   mat->getRowMinMaxVec(get_min, initialize_vec, col_scale, *minmax.vec);

   for( size_t i = 0; i < children.size(); ++i )
   {
      assert( minmax.children[i]);
      if( children[i]->isKindOf(kStringGenDummyMatrix) )
         assert( minmax.children[i]->isKindOf(kStochDummy) );

      children[i]->getRowMinMaxVecVertical(get_min, initialize_vec, col_scale, *minmax.children[i]);
   }

   if( mat_link )
      mat_link->getRowMinMaxVec(get_min, initialize_vec, col_scale, *minmax.vecl);
}

void StringGenMatrix::getRowMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* colScaleVec, OoqpVector& minmaxVec )
{
   if( is_vertical )
      getRowMinMaxVecVertical(getMin, initializeVec, colScaleVec, minmaxVec);
   else
      getRowMinMaxVecHorizontal(getMin, initializeVec, colScaleVec, minmaxVec);
}


void StringGenMatrix::getColMinMaxVec( bool getMin, bool initializeVec, const OoqpVector* rowScaleVec, OoqpVector& minmaxVec )
{
   if( is_vertical )
      getColMinMaxVecVertical(getMin, initializeVec, rowScaleVec, minmaxVec);
   else
      getColMinMaxVecHorizontal(getMin, initializeVec, rowScaleVec, minmaxVec);
}

void StringGenMatrix::columnScaleVertical( const OoqpVector& vec )
{
   assert( is_vertical );

   mat->columnScale( vec );

   for( size_t i = 0; i < children.size(); i++ )
      children[i]->columnScaleVertical( vec );

   if( mat_link )
      mat_link->columnScale( vec );
}

void StringGenMatrix::columnScaleHorizontal( const OoqpVector& vec_in )
{
   assert( !is_vertical );

   const StochVector& vec = dynamic_cast<const StochVector&>(vec_in);

   assert(vec.vec);
   assert(vec.children.size() == children.size());
   assert( (vec.vecl && mat_link) || (vec.vecl == nullptr && mat_link == nullptr) );

   mat->columnScale(*vec.vec);

   for( size_t i = 0; i < children.size(); ++i )
   {
      assert(vec.children[i]);
      if( children[i]->isKindOf(kStochGenDummyMatrix) )
         assert( vec.children[i]->isKindOf(kStochDummy) );

      children[i]->columnScaleHorizontal(*vec.children[i]);
   }

   if( mat_link )
      mat_link->columnScale(*vec.vecl);
}

void StringGenMatrix::rowScaleVertical( const OoqpVector& vec_in )
{
   assert( is_vertical );

   const StochVector& vec = dynamic_cast<const StochVector&>(vec_in);

   assert(vec.vec);
   assert(vec.children.size() == children.size());
   assert( (vec.vecl && mat_link) || (vec.vecl == nullptr && mat_link == nullptr) );

   mat->rowScale(*vec.vec);

   for( size_t i = 0; i < children.size(); i++ )
   {
      assert( vec.children[i] );
      if( children[i]->isKindOf(kStringGenDummyMatrix) )
         assert( vec.children[i]->isKindOf(kStochDummy) );

      children[i]->rowScaleVertical(*vec.children[i]);
   }

   if( mat_link )
      mat_link->rowScale(*vec.vecl);
}

void StringGenMatrix::rowScaleHorizontal( const OoqpVector& vec )
{
   assert( !is_vertical );
   mat->rowScale( vec );

   if( mat_link )
      mat_link->rowScale( vec );

   for( size_t i = 0; i < children.size(); ++i )
      children[i]->rowScaleHorizontal(vec);
}

void StringGenMatrix::columnScale ( const OoqpVector& vec )
{
   if( is_vertical )
      columnScaleVertical(vec);
   else
      columnScaleHorizontal(vec);
}

void StringGenMatrix::rowScale ( const OoqpVector& vec )
{
   if( is_vertical )
      rowScaleVertical(vec);
   else
      rowScaleHorizontal(vec);
}

void StringGenMatrix::writeToStreamDense(std::ostream& out) const
{
   assert( 0 && "TODO: implement...");
}

std::string StringGenMatrix::writeToStreamDenseRowChildren(int row) const
{
   std::string str_all;
   for( size_t it = 0; it < children.size(); it++ )
   {
      if( children[it]->isKindOf(kStringGenDummyMatrix) )
         continue;
      else
         assert( children[it]->mat );

      std::string str = children[it]->mat->writeToStreamDenseRow(row);
      str_all.append(str);
   }
   return str_all;
}

