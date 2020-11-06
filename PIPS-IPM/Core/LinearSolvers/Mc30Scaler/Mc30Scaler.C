/*
 * Mc30Scaler.C
 *
 *  Created on: Nov 6, 2020
 *      Author: bzfkempk
 */
#include "Mc30Scaler.h"
#include "SimpleVector.h"

Mc30Scaler::Mc30Scaler() : scaling_output_control(1), scaling_error(0)
{
}

Mc30Scaler::~Mc30Scaler()
{
}

void Mc30Scaler::getFortranIndex(const int* rowM, const int* colM, int length, int max_index)
{
   assert( length > 0 );
   if( static_cast<size_t>(length) > rowM_ft_indexed.size() )
   {
      rowM_ft_indexed.resize(length);
      colM_ft_indexed.resize(length);
   }

   for( int i = 0; i < length; ++i )
   {
      assert( rowM[i] + 1 <= max_index );
      assert( colM[i] + 1 <= max_index );
      assert( 0 <= rowM[i] );
      assert( 0 <= colM[i] );

      rowM_ft_indexed[i] = rowM[i] + 1;
      colM_ft_indexed[i] = colM[i] + 1;
   }
}

void Mc30Scaler::scaleMatrixTripletFormat( int n, int nnz, double* M, const int* rowM, const int* colM, bool fortran_indexed )
{
   assert( n > 0 );
   assert( nnz > 0);

   if( !fortran_indexed )
      getFortranIndex(rowM, colM, nnz, n);

   /* set working arrays for scaling */
   if( scaling_factors.size() < static_cast<size_t>(n) )
   {
      scaling_factors.resize(n);
      scaling_workspace.resize(4 * n);
   }

   /* compute the log scaling factors */
   FNAME(mc30ad)( &n, &nnz, M,
         fortran_indexed ? rowM : rowM_ft_indexed.data(),
         fortran_indexed ? colM : colM_ft_indexed.data(),
         scaling_factors.data(), scaling_workspace.data(), &scaling_output_control, &scaling_error);
   assert( scaling_error == 0 );

   /* compute the scaling factors */
   for( size_t i = 0; i < scaling_factors.size(); ++i )
      scaling_factors[i] = std::exp(scaling_factors[i]);

   /* do the actual scaling */
   if( fortran_indexed )
   {
      for( int i = 0; i < nnz; ++i )
      {
         const double scale_ai = scaling_factors[rowM[i] - 1] * scaling_factors[colM[i] - 1];
         M[i] *= scale_ai;
      }
   }
   else
   {
      for( int i = 0; i < nnz; ++i )
      {
         const double scale_ai = scaling_factors[rowM[i]] * scaling_factors[colM[i]];
         M[i] *= scale_ai;
      }
   }
}

void Mc30Scaler::scaleVector( OoqpVector& vec_in ) const
{
   SimpleVector& vec = dynamic_cast<SimpleVector&>(vec_in);
   assert( static_cast<size_t>(vec.length()) == scaling_factors.size() );

   for( size_t i = 0; i < scaling_factors.size(); ++i )
      vec[i] *= scaling_factors[i];
}

void Mc30Scaler::unscaleVector( OoqpVector& vec_in ) const
{
   SimpleVector& vec = dynamic_cast<SimpleVector&>(vec_in);
   assert( static_cast<size_t>(vec.length()) == scaling_factors.size() );

   for( size_t i = 0; i < scaling_factors.size(); ++i )
      vec[i] /= scaling_factors[i];
}

