/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP 
 *
 * 10/2/2009 MatMult function added by C. Petra
 */
 
#include <cassert>
#include <cmath>

#include "DenseGenMatrix.h"
#include "DenseSymMatrix.h"
#include "OoqpBlas.h"
#include "SimpleVector.h"

#include "DoubleMatrixTypes.h"
#include "SparseGenMatrix.h"

int DenseGenMatrix::isKindOf( int type ) const
{
  return type == kDenseGenMatrix || type == kGenMatrix;
}


DenseGenMatrix::DenseGenMatrix( int size )
{
  mStorage = DenseStorageHandle( new DenseStorage( size, size ) );
}


DenseGenMatrix::DenseGenMatrix( double A[], int m, int n )
{
  mStorage = DenseStorageHandle( new DenseStorage( A, m, n ) );
}


DenseGenMatrix::DenseGenMatrix( int m, int n  )
{
  mStorage = DenseStorageHandle( new DenseStorage( m, n ) );
}


// Delegate these methods to the storage
void DenseGenMatrix::atPutDense( int row, int col, double * A, int lda,
				     int rowExtent, int colExtent )
{
  mStorage->atPutDense( row, col, A, lda, rowExtent, colExtent );
}


void DenseGenMatrix::atPutZeros( int row, int col,
				     int rowExtent, int colExtent )
{
  mStorage->atPutZeros( row, col, rowExtent, colExtent );
}

void DenseGenMatrix::putZeros()
{
   mStorage->putZeros();
}

void DenseGenMatrix::putSparseTriple( int irow[], int len,
					  int jcol[], double A[], 
					  int& info )
{
  mStorage->putSparseTriple( irow, len, jcol, A, info );
}


void DenseGenMatrix::fromGetSpRow( int row, int col,
				       double A[], int lenA,
				       int jcolA[], int& nnz,
				       int colExtent, int& info )
{
  mStorage->fromGetSpRow( row, col, A, lenA, jcolA, nnz, colExtent, info );
}


void DenseGenMatrix::atPutSpRow( int  row, double * A,
				     int lenA , int * jcolA,
				     int& info )
{
  mStorage->atPutSpRow( row, A, lenA, jcolA, info );
}


void DenseGenMatrix::getDiagonal( OoqpVector& vec )
{
  mStorage->getDiagonal( vec );
}


void DenseGenMatrix::setToDiagonal( const OoqpVector& vec )
{
  mStorage->setToDiagonal( vec );
}

void DenseGenMatrix::getSize( long long& m, long long& n ) const
{
  m = mStorage->m;
  n = mStorage->n;
}

void DenseGenMatrix::getSize( int& m, int& n ) const
{
  m = mStorage->m;
  n = mStorage->n;
}

void DenseGenMatrix::atPutSubmatrix( int destRow, int destCol,
					 DoubleMatrix & Mat,
					 int srcRow, int srcCol,
					 int rowExtent, int colExtent )
{
  int m = mStorage->m, n = mStorage->n;
  double ** M = mStorage->M;

  assert( destRow >= 0 && destRow + rowExtent <= m );
  assert( destCol >= 0 && destCol + colExtent <= n );

  // If assertions are turned off, clip to the actual size of this matrix
  destRow = ( destRow >= 0 ) ? destRow : 0;
  destCol = ( destCol >= 0 ) ? destCol : 0;
  rowExtent = ( destRow + rowExtent <= m ) ?  rowExtent : m - destRow;
  colExtent = ( destCol + colExtent <= n ) ?  colExtent : n - destCol;

  Mat.fromGetDense( srcRow, srcCol, &M[destRow][destCol], n,
		     rowExtent, colExtent );
}


void DenseGenMatrix::mult ( double beta,  double y[], int incy,
				double alpha, const double x[], int incx ) const
{
  char fortranTrans = 'T';
  int n = mStorage->n, m = mStorage->m;
  
  dgemv_( &fortranTrans, &n, &m, &alpha, &mStorage->M[0][0], &n,
	  x, &incx, &beta, y, &incy );
}


void DenseGenMatrix::mult ( double beta,  OoqpVector& y_in,
				double alpha, const OoqpVector& x_in ) const
{
  char fortranTrans = 'T';
  int n = mStorage->n, m = mStorage->m;
  double ** M = mStorage->M;
  int incx = 1, incy = 1;

  const SimpleVector & x = dynamic_cast<const SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);

  if( n != 0 && m != 0 ) {
    dgemv_( &fortranTrans, &n, &m, &alpha, &M[0][0], &n,
	    &x[0], &incx, &beta, &y[0], &incy );
  } else {
    if( m != 0 ) y.scale( beta );
  }
}


void DenseGenMatrix::transMult ( double beta,  double y[], int incy,
				     double alpha, const double x[], int incx ) const
{
  char fortranTrans = 'N';
  int n = mStorage->n, m = mStorage->m; 
  double ** M = mStorage->M;
 
  dgemv_( &fortranTrans, &n, &m, &alpha, &M[0][0], &n,
	  x, &incx, &beta, y, &incy );
}


void DenseGenMatrix::transMult ( double beta,  OoqpVector& y_in,
				     double alpha, const OoqpVector& x_in ) const
{
  char fortranTrans = 'N';
  int n = mStorage->n, m = mStorage->m;
  double ** M = mStorage->M;
  const SimpleVector & x = dynamic_cast<const SimpleVector &>(x_in);
  SimpleVector & y = dynamic_cast<SimpleVector &>(y_in);
  int incx = 1, incy = 1;
  
  if( m != 0 && n != 0 ) {
    dgemv_( &fortranTrans, &n, &m, &alpha, &M[0][0], &n,
	    &x[0], &incx, &beta, &y[0], &incy );
  } else {
    if( n != 0 ) {
	y.scale( beta );
    }
  }
}

double DenseGenMatrix::abmaxnorm() const
{
   assert( mStorage.notNil() );
   return mStorage->abmaxnorm();
}
double DenseGenMatrix::abminnormNonZero( double tol ) const
{
   assert( mStorage.notNil() );
   return mStorage->abminnormNonZero(tol);
}

void DenseGenMatrix::writeToStream( std::ostream& out ) const
{
   for( int i = 0; i < mStorage->m; i++ )
   {
      for( int j = 0; j < mStorage->n; j++ )
         out << mStorage->M[i][j] << "\t";

      out << std::endl;
   }
}

void DenseGenMatrix::writeToStreamDense( std::ostream& out ) const
{
   writeToStream(out);
}


void DenseGenMatrix::randomize( double alpha, double beta, double * seed )
{
  int m = mStorage->m, n = mStorage->n;
  double ** M = mStorage->M;
  double drand(double *);
  int i, j;

  double scale = beta - alpha;
  double shift = alpha/scale;

  for( i = 0; i < m; i++ ) {
    for( j = 0; j < n; j++ ) {
      M[i][j] = scale * (drand(seed) + shift);
    }
  }
}


void DenseGenMatrix::fromGetDense( int row, int col, double * A,
				       int lda,
				       int rowExtent, int colExtent )
{
  int m = mStorage->m, n = mStorage->n;

  assert( row >= 0 && row + rowExtent <= m );
  assert( col >= 0 && col + colExtent <= n );

  // If assertions are turned off, clip to the actual size of this matrix
  row = ( row >= 0 ) ? row : 0;
  col = ( col >= 0 ) ? col : 0;
  int lrow = ( row + rowExtent <= m ) ?  rowExtent : m - row;
  int lcol = ( col + colExtent <= n ) ?  colExtent : n - col;
  
  mStorage->fromGetDense( row, col, A, lda, lrow, lcol );
}

void DenseGenMatrix::atPutDiagonal( int idiag, OoqpVector& v )
{
  mStorage->atPutDiagonal( idiag, v );
}

void DenseGenMatrix::fromGetDiagonal( int idiag, OoqpVector& v )
{
  mStorage->fromGetDiagonal( idiag, v );
}

void DenseGenMatrix::getRow ( int rowIndex, OoqpVector& v_in)
{
  assert (rowIndex >= 0 && rowIndex <= mStorage->m);
  SimpleVector & v = dynamic_cast<SimpleVector &>(v_in);

  mStorage->fromGetDense(rowIndex, 0, &v[0], 1, 1, mStorage->n);
}

void DenseGenMatrix::columnScale( const OoqpVector& vec )
{
  mStorage->columnScale( vec );
}

void DenseGenMatrix::symmetricScale( const OoqpVector& vec )
{
  mStorage->symmetricScale( vec );
}

void DenseGenMatrix::rowScale( const OoqpVector& vec )
{
  mStorage->columnScale( vec );
}

void DenseGenMatrix::scalarMult( double num )
{
  mStorage->scalarMult( num );
}

void DenseGenMatrix::matMult(double alpha, 
			     DenseGenMatrix& A, int transA, 
			     DenseGenMatrix& B, int transB,
			     double beta)
{
  assert(0);
  char fortranTransA = (transA==0?'N':'T');
  char fortranTransB = (transB==0?'N':'T');

  DenseGenMatrix& C = *this;

  int m,n,k,tmp;


  if(transA) {
    // A^T op(B)
    A.getSize(k,m);
  } else {
    // A op(B)
    A.getSize(m,k);
  }

  if(transB) {
    //op(A) B^T
    B.getSize(n,tmp);
  } else { 
    // op(A) B
    B.getSize(tmp, n);
  }

  assert(m == C.mStorage->m);
  assert(n == C.mStorage->n);
  assert(k == tmp);

  double ** CC = C.mStorage->M;
  double ** AA = A.mStorage->M;
  double ** BB = B.mStorage->M;

  dgemm_(&fortranTransA, &fortranTransB, 
	 &m,&n,&k,
	 &alpha, 
	 &AA[0][0], &m,
	 &BB[0][0], &k,
	 &beta,
	 &CC[0][0], &m);

  //  if( n != 0 && m != 0 ) {
  //  dgemv_( &fortranTrans, &n, &m, &alpha, &C[0][0], &n,
  //	    &x[0], &incx, &beta, &y[0], &incy );

}

// TODO : probably move to some utility class..
/* compute beta * res += alpha * this * mat where mat gets multiplied to the submatrix
 * starting at mul_start and the results gets added starting at res_start */
void DenseGenMatrix::multMatAt( int row_start, int row_end, int col_offset_this, double beta, int row_start_res, int col_offset_result, DenseGenMatrix& res, double alpha, /* const */ SparseGenMatrix& mat ) const
{
   assert( 0 <= col_offset_this );
   assert( 0 <= row_start );
   assert( row_start <= row_end );
   assert( row_end <= mStorage->m );

   int mat_m, mat_n; mat.getSize( mat_m, mat_n );

   assert( col_offset_this <= mStorage->n && col_offset_this + mat_m <= mStorage->n );

   const int n_rows = row_end - row_start;
   assert( col_offset_result >= 0 );
   assert( col_offset_result <= res.mStorage->n && col_offset_result + mat_n <= res.mStorage->n );
   assert( n_rows + row_start_res <= res.mStorage->m );

   SparseStorage& mat_tp = mat.getTranspose().getStorageRef();
   for( int j = 0; j < n_rows; ++j )
   {
      for( int i = 0; i < mat_n; ++i )
      {
         if( beta != 1.0 )
            res[row_start + j][col_offset_result + i] *= beta;

         const int col_start = mat_tp.krowM[i];
         const int col_end = mat_tp.krowM[i + 1];

         for( int k = col_start; k < col_end; ++k )
         {
            const int row = mat_tp.jcolM[k];
            const double val = mat_tp.M[k];

            assert( col_offset_result + i < res.mStorage->n );
            assert( row < mat_m );
            assert( row + col_offset_this < mStorage->n );
            res[row_start_res + j][col_offset_result + i] += mStorage->M[row_start + j][row + col_offset_this] * val * alpha;
         }
      }
   }
}

void DenseGenMatrix::addMatAt( const SparseGenMatrix& mat, int mat_row_start, int mat_row_end, int this_row_0, int this_col_0 )
{
   int mmat, nmat; mat.getSize( mmat, nmat );
   if( mmat <= 0 || nmat <= 0 )
      return;

   assert( 0 <= mat_row_start && mat_row_start <= mat_row_end && mat_row_end - 1 < mmat );

   const int n_rows = mat_row_end - mat_row_start;
   assert( 0 <= this_row_0 && this_row_0 + n_rows - 1 < mStorage->m );
   assert( 0 <= this_col_0 && this_col_0 + nmat - 1 < mStorage->n );

   for( int row = 0; row < n_rows; ++row )
   {
      const int row_in_mat = mat_row_start + row;
      const int row_start = mat.krowM()[row_in_mat];
      const int row_end = mat.krowM()[row_in_mat + 1];

      for( int j = row_start; j < row_end; ++j )
      {
         const int col = mat.jcolM()[j];
         const double val = mat.M()[j];

         assert( col < nmat );
         assert( this_col_0 + col < mStorage->n );
         assert( this_row_0 + row < mStorage->m );
         this->operator [](this_row_0 + row)[this_col_0 + col] += val;
      }
   }
}
