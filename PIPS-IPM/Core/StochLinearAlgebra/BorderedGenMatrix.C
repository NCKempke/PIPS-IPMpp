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

BorderedGenMatrix::BorderedGenMatrix(int id_, StochGenMatrix* inner_matrix, StringGenMatrix* border_left,
            StringGenMatrix* border_bottom, SparseGenMatrix* bottom_left_block, MPI_Comm mpi_comm_) :
            inner_matrix(inner_matrix), border_left(border_left), border_bottom(border_bottom), bottom_left_block(bottom_left_block),
            id(id_), mpi_comm(mpi_comm_), distributed( mpi_comm == MPI_COMM_NULL ), rank( PIPS_MPIgetRank(mpi_comm) )
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
   MPI_Status status;

   const int my_rank = PIPS_MPIgetRank(mpi_comm);
   const int size = PIPS_MPIgetSize(mpi_comm);
   const bool iAmDistrib = ( size != 0 );
   this->inner_matrix->writeToStreamDenseBordered( *border_left, out );

   int mL, nL; this->bottom_left_block->getSize(mL, nL);
   if( mL > 0 )
   {
      if( iAmDistrib )
         MPI_Barrier(mpi_comm);

      // for each row r do:
      for( int r = 0; r < mL; r++ )
      {
         if( iAmDistrib )
         {
            MPI_Barrier(mpi_comm);

            // process Zero collects all the information and then prints it.
            if( my_rank == 0 )
            {
               out << bottom_left_block->writeToStreamDenseRow(r);
               out << "|\t";

               out << this->border_bottom->mat->writeToStreamDenseRow(r);

               out << border_bottom->writeToStreamDenseRowChildren(r);

               for( int p = 1; p < size; p++ )
               {
                  int l;
                  MPI_Probe(p, r + 1, mpi_comm, &status);
                  MPI_Get_count(&status, MPI_CHAR, &l);
                  char *buf = new char[l];
                  MPI_Recv(buf, l, MPI_CHAR, p, r + 1, mpi_comm, &status);
                  std::string rowPartFromP(buf, l);
                  out << rowPartFromP;
                  delete[] buf;
               }
               out << std::endl;

            }
            else // rank != 0
            {
               std::string str = border_bottom->writeToStreamDenseRowChildren(r);
               MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, r + 1, mpi_comm);
            }
         }
         else // not distributed
         {
            std::stringstream sout;
            bottom_left_block->writeToStreamDenseRow(sout, r);
            sout << "|\t";

            this->border_bottom->mat->writeToStreamDenseRow(sout, r);

            for( size_t it = 0; it < border_bottom->children.size(); it++ )
               border_bottom->children[it]->mat->writeToStreamDenseRow(sout, r);

            out << sout.rdbuf() << std::endl;
         }
      }
   }
   if( iAmDistrib )
      MPI_Barrier(mpi_comm);
}

