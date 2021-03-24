#include "StochGenMatrix.h"

#include "DoubleMatrixTypes.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "pipsdef.h"
#include "pipsport.h"

#include <limits>
#include <algorithm>
#include <numeric>
#include <memory>

StochGenMatrix::StochGenMatrix( GenMatrix* Amat, GenMatrix* Bmat, GenMatrix* Blmat, MPI_Comm mpiComm_, bool inner_leaf, bool inner_root )
   : Amat{Amat}, Bmat{Bmat}, Blmat{Blmat}, mpiComm{mpiComm_}, iAmDistrib{ PIPS_MPIgetDistributed(mpiComm) },
     inner_leaf{inner_leaf}, inner_root{inner_root}
{
   assert( Amat );
   assert( Bmat );
   assert( Blmat );

#ifndef NDEBUG
   int mA, nA;
   Amat->getSize(mA, nA);
   int mB, nB;
   Bmat->getSize(mB, nB);
   int mBl, nBl;
   Blmat->getSize(mBl, nBl);

   if( nA == 0 && mA == 0 )
      assert( nBl == nB );
   else
   {
      assert( mA == mB );
      assert( nB == nBl );
   }
#endif

   recomputeSize();
}

StochGenMatrix::StochGenMatrix(long long global_m, long long global_n, int A_m,
      int A_n, int A_nnz, int B_m, int B_n, int B_nnz, MPI_Comm mpiComm_ ) :
      m(global_m), n(global_n), mpiComm(mpiComm_), iAmDistrib(PIPS_MPIgetDistributed(mpiComm))
{
   Amat = new SparseGenMatrix(A_m, A_n, A_nnz);
   Bmat = new SparseGenMatrix(B_m, B_n, B_nnz);
   Blmat = new SparseGenMatrix(0, B_n, 0);
}

StochGenMatrix::StochGenMatrix(long long global_m, long long global_n, int A_m,
      int A_n, int A_nnz, int B_m, int B_n, int B_nnz, int Bl_m, int Bl_n,
      int Bl_nnz, MPI_Comm mpiComm_ ) :
      m(global_m), n(global_n), mpiComm(mpiComm_), iAmDistrib(PIPS_MPIgetDistributed(mpiComm))
{
   Amat = new SparseGenMatrix(A_m, A_n, A_nnz);
   Bmat = new SparseGenMatrix(B_m, B_n, B_nnz);
   Blmat = new SparseGenMatrix(Bl_m, Bl_n, Bl_nnz);
}

StochGenMatrix::StochGenMatrix(long long global_m, long long global_n,
      MPI_Comm mpiComm_) :
      m(global_m), n(global_n), mpiComm(mpiComm_), iAmDistrib(PIPS_MPIgetDistributed(mpiComm))
{
}

StochGenMatrix::~StochGenMatrix()
{
  for(size_t it = 0; it < children.size(); it++)
    delete children[it];

  delete Amat;
  delete Bmat;
  delete Blmat;
}


bool StochGenMatrix::amatEmpty() const
{
   int mA,nA;
   Amat->getSize(mA, nA);
   return mA <= 0 || nA <= 0;
}

bool StochGenMatrix::hasSparseMatrices() const
{
   return Amat->isKindOf(kSparseGenMatrix) && Bmat->isKindOf(kSparseGenMatrix)
      && Blmat->isKindOf(kSparseGenMatrix);
}
GenMatrix* StochGenMatrix::cloneFull(bool switchToDynamicStorage) const
{
   StochGenMatrix* clone = new StochGenMatrix(m, n, mpiComm);
   assert( hasSparseMatrices() );

   // clone submatrices
   clone->Amat = Amat->cloneFull(switchToDynamicStorage);
   clone->Bmat = Bmat->cloneFull(switchToDynamicStorage);
   clone->Blmat = Blmat->cloneFull(switchToDynamicStorage);

   for( size_t it = 0; it < children.size(); it++ )
      clone->children.push_back( dynamic_cast<StochGenMatrix*>( children[it]->cloneFull(switchToDynamicStorage) ) );

   return clone;	
}

/* creates an empty copy of the matrix with n = 0 for all submatrices and m (cols) as before */
GenMatrix* StochGenMatrix::cloneEmptyRows(bool switchToDynamicStorage) const
{
   StochGenMatrix* clone = new StochGenMatrix(m, n, mpiComm);
   assert( hasSparseMatrices() );

   // clone submatrices
   clone->Amat = Amat->cloneEmptyRows(switchToDynamicStorage);
   clone->Bmat = Bmat->cloneEmptyRows(switchToDynamicStorage);
   clone->Blmat = Blmat->cloneEmptyRows(switchToDynamicStorage);

   for( size_t it = 0; it < children.size(); it++ )
      clone->children.push_back( dynamic_cast<StochGenMatrix*>( children[it]->cloneEmptyRows(switchToDynamicStorage) ) );

   return clone;
}

void StochGenMatrix::AddChild(StochGenMatrix* child)
{
   children.push_back(child);
}

int StochGenMatrix::isKindOf( int type ) const
{
  return type == kStochGenMatrix || type == kGenMatrix;
}

int StochGenDummyMatrix::isKindOf( int type ) const
{
  return type == kStochGenDummyMatrix || type == kGenMatrix || type == kStochGenMatrix;
}

void StochGenMatrix::getSize( long long& m_out, long long& n_out ) const
{
  m_out = m; n_out = n;
}

void StochGenMatrix::getSize( int& m_out, int& n_out ) const
{
  m_out = m; n_out = n;
}

void StochGenMatrix::columnScale2( const OoqpVector& vec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);
   assert(scalevec.children.size() == 0 && children.size() == 0);

   Bmat->columnScale( *scalevec.vec );

   if( !amatEmpty() )
      Amat->columnScale( *scalevec.getLinkingVecNotHierarchicalTop() );

   Blmat->columnScale( *scalevec.vec );
}

void StochGenMatrix::columnScale( const OoqpVector& vec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);

   assert( amatEmpty() );

   Bmat->columnScale( *scalevec.getLinkingVecNotHierarchicalTop() );
   Blmat->columnScale( *scalevec.getLinkingVecNotHierarchicalTop() );

   assert(children.size() == scalevec.children.size());
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->columnScale2( *(scalevec.children[it]) );
}

void StochGenMatrix::rowScale2( const OoqpVector& vec, const OoqpVector* linkingvec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);

   assert(scalevec.children.size() == 0 && children.size() == 0);

   if( !amatEmpty() )
      Amat->rowScale(*scalevec.vec);

   Bmat->rowScale(*scalevec.vec);

   if( linkingvec )
      Blmat->rowScale(*linkingvec);
}

void StochGenMatrix::rowScale( const OoqpVector& vec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);

   Bmat->rowScale(*scalevec.vec);
   if( scalevec.vecl )
      Blmat->rowScale(*scalevec.vecl);

   assert(children.size() == scalevec.children.size());

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->rowScale2(*(scalevec.children[it]), scalevec.vecl);
}

void StochGenMatrix::scalarMult( double num)
{
  Amat->scalarMult(num);
  Bmat->scalarMult(num);
  Blmat->scalarMult(num);

  for (size_t it = 0; it < children.size(); it++) 
    children[it]->scalarMult(num);
}

void StochGenMatrix::getDiagonal( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);

  Bmat->getDiagonal(*vec.vec);

  assert( children.size() == vec.children.size() );

  for( size_t it = 0; it < children.size(); it++ )
    children[it]->getDiagonal(*vec.children[it]);
}
 
void StochGenMatrix::setToDiagonal( const OoqpVector& vec_ )
{
  const StochVector& vec = dynamic_cast<const StochVector&>(vec_);

  Bmat->setToDiagonal(*vec.vec);

  assert(children.size() == vec.children.size());

  for(size_t it = 0; it < children.size(); it++)
    children[it]->setToDiagonal(*vec.children[it]);
}

/* y = beta * y + alpha * this * x */
void StochGenMatrix::mult( double beta, OoqpVector& y_,
			   double alpha, const OoqpVector& x_ ) const
{
   if( 0.0 == alpha )
   {
      y_.scale( beta );
      return;
   }

   const StochVector & x = dynamic_cast<const StochVector&>(x_);
   StochVector& y = dynamic_cast<StochVector&>(y_);

   assert( amatEmpty() );
   Bmat->mult(beta, *y.vec, alpha, *x.getLinkingVecNotHierarchicalTop() );

   if( y.vecl )
   {
      if( iAmSpecial(iAmDistrib, mpiComm ) )
         Blmat->mult(beta, *y.vecl, alpha, *x.getLinkingVecNotHierarchicalTop() );
      else
         y.vecl->setToZero();
   }

   assert( y.children.size() == children.size() );
   assert( x.children.size() == children.size() );

   for(size_t it = 0; it < children.size(); it++)
      children[it]->mult2(beta, *y.children[it], alpha, *x.children[it], y.vecl);

   if( iAmDistrib && y.vecl )
      PIPS_MPIsumArrayInPlace( dynamic_cast<SimpleVector*>(y.vecl)->elements(), y.vecl->length(), mpiComm);
}


/* mult method for children; needed only for linking constraints */
void StochGenMatrix::mult2( double beta,  StochVector& y,
			   double alpha, StochVector& x, OoqpVector* yparentl_ )
{
   assert( alpha != 0.0 );
   assert( children.size() == 0 );
   assert( y.children.size() == children.size() );
   assert( x.children.size() == children.size() );
   assert( x.vec );
   assert( y.vec );

   Bmat->mult(beta, *y.vec, alpha, *x.vec);

   if( yparentl_ )
   {
      Blmat->mult(1.0, *yparentl_, alpha, *x.vec);
      if( !iAmSpecial(iAmDistrib, mpiComm) )
         yparentl_->setToZero();
   }

   if( !amatEmpty() )
   {
      const OoqpVector* link_vec = x.getLinkingVecNotHierarchicalTop();
      assert( link_vec != x.vec );
      Amat->mult(1.0, *y.vec, alpha, *link_vec);
   }
}


void StochGenMatrix::transMult ( double beta, OoqpVector& y_,
				 double alpha, const OoqpVector& x_ ) const
{
   if( 0.0 == alpha )
   {
      y_.scale(beta);
      return;
   }

   const StochVector& x = dynamic_cast<const StochVector&>(x_);
   StochVector& y = dynamic_cast<StochVector&>(y_);

   const bool at_root = y.vec == y.getLinkingVecNotHierarchicalTop();
   assert( y.vec );
   assert( x.vec );

   if( iAmSpecial(iAmDistrib, mpiComm) )
   {
      Bmat->transMult(at_root ? beta : 1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.vec);

      if( x.vecl )
         Blmat->transMult(1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.vecl);
   }
   else if( at_root )
      y.vec->setToZero();

   assert(y.children.size() == children.size());
   assert(x.children.size() == children.size());

   for(size_t it = 0; it < children.size(); it++)
      children[it]->transMult2(beta, *y.children[it], alpha, *x.children[it], x.vecl);

   if( iAmDistrib && y.vec == y.getLinkingVecNotHierarchicalTop() )
      PIPS_MPIsumArrayInPlace( dynamic_cast<SimpleVector*>(y.vec)->elements(), y.vec->length(), mpiComm);
}

void StochGenMatrix::transMult2 ( double beta, StochVector& y,
				  double alpha, StochVector& x, const OoqpVector* xvecl) const
{
   assert( alpha != 0.0 );
   assert( x.vec );
   assert( y.vec );
   assert( y.children.size() - children.size() == 0 );
   assert( x.children.size() - children.size() == 0 );
   assert( children.size() == 0 );

   if( !amatEmpty() )
      Amat->transMult(1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.vec);

   Bmat->transMult(beta, *y.vec, alpha, *x.vec);
   if( xvecl )
      Blmat->transMult(1.0, *y.vec, alpha, *xvecl );
}

double StochGenMatrix::abmaxnorm() const
{
   double nrm = 0.0;
  
   for(size_t it = 0; it < children.size(); it++)
      nrm = std::max(nrm, children[it]->abmaxnorm());

   if(iAmDistrib)
      PIPS_MPIgetMaxInPlace( nrm, mpiComm );

   nrm = std::max(nrm, std::max(Amat->abmaxnorm(), Bmat->abmaxnorm()));
   nrm = std::max(nrm, Blmat->abmaxnorm());

   return nrm;
}

double StochGenMatrix::abminnormNonZero( double tol ) const
{
  double nrm = std::numeric_limits<double>::infinity();

  for(size_t it = 0; it < children.size(); it++)
    nrm = std::min(nrm, children[it]->abminnormNonZero( tol ) );

  if( iAmDistrib )
     PIPS_MPIgetMinInPlace( nrm, mpiComm );

  nrm = std::min(nrm, std::min(Amat->abminnormNonZero( tol ), Bmat->abminnormNonZero( tol )));
  nrm = std::min(nrm, Blmat->abminnormNonZero( tol ));

  return nrm;
}

void StochGenMatrix::getLinkVarsNnz(std::vector<int>& vec) const
{
   assert( hasSparseMatrices() );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->getLinkVarsNnzChild(vec);

   if( iAmDistrib )
   {
      int* buffer = new int[vec.size()];
      MPI_Allreduce(&vec[0], buffer, vec.size(), MPI_INT, MPI_SUM, mpiComm);

      std::memcpy(&vec[0], buffer, vec.size() * sizeof(int));

      delete[] buffer;
   }
}

void StochGenMatrix::getLinkVarsNnzChild(std::vector<int>& vec) const
{
   assert(children.size() == 0);
   assert( hasSparseMatrices() );

   dynamic_cast<const SparseGenMatrix*>(Amat)->getLinkVarsNnz(vec);
}

void StochGenMatrix::writeToStreamDenseBorderedChild( const StringGenMatrix& border_left, std::ostream& out, int offset ) const
{
   assert( border_left.children.size() == this->children.size() );

   if( Bmat->isKindOf(kStochGenMatrix) )
   {
      assert( border_left.mat->isKindOf(kStringGenMatrix) );
      dynamic_cast<StochGenMatrix*>(Bmat)->writeToStreamDenseBordered( dynamic_cast<StringGenMatrix&>(*border_left.mat), out, offset);
   }
   /// Border.mat | Amat | offset | Bmat ///
   else
   {
      assert( border_left.mat->isKindOf(kSparseGenMatrix) );
      assert( hasSparseMatrices() );
      assert( PIPS_MPIgetRank(mpiComm) == 0 );

      int nB, mB;
      Bmat->getSize(mB, nB);
      int nA, mA;
      Amat->getSize(mA, nA);
      int nBd, mBd;
      border_left.mat->getSize(mBd, nBd);

      assert( mB == mA && mA == mBd );

      for( int row = 0; row < mB; ++row )
      {
         border_left.mat->writeToStreamDenseRow( out, row );
         out << "|\t";
         Amat->writeToStreamDenseRow(out, row);

         for( int i = 0; i < offset; ++i )
            out << "\t";
         Bmat->writeToStreamDenseRow(out, row);
         out << "\n";
      }
   }
}


void StochGenMatrix::writeToStreamDenseBordered( const StringGenMatrix& border_left, std::ostream& out, int offset ) const
{
   assert( border_left.children.size() == this->children.size() );
   assert( hasSparseMatrices() );
   assert( children.size() != 0 );

   const int original_offset = offset;
   const int my_rank = PIPS_MPIgetRank(mpiComm);

   if( iAmDistrib )
      MPI_Barrier(mpiComm);

   int mBmat, nBmat; this->Bmat->getSize(mBmat, nBmat);
   int mBd, nBd; border_left.mat->getSize(mBd, nBd);
   assert( mBmat == mBd );

   assert( Bmat->isKindOf(kSparseGenMatrix) );
   assert( border_left.mat->isKindOf(kSparseGenMatrix) );

   /// Border.mat | Bmat ///
   if( my_rank == 0 )
   {
      for( int i = 0; i < mBmat; ++i )
      {
         dynamic_cast<const SparseGenMatrix&>(*border_left.mat).writeToStreamDenseRow(out, i);
         out << "|\t";
         dynamic_cast<const SparseGenMatrix*>(this->Bmat)->writeToStreamDenseRow(out, i);
         out << "\n";
      }
   }

   /// write children ///
   /// BorderChild | offset | Child ///
   std::ostringstream child_stream{};
   if( !amatEmpty() )
      assert( children.size() ==  0 );

   for( size_t it = 0; it < children.size(); it++ )
   {
      MPI_Barrier(mpiComm);
      const StringGenMatrix& border_child = *border_left.children[it];

      children[it]->writeToStreamDenseBorderedChild( border_child, child_stream, offset );

      int mChild, nChild;
      children[it]->getSize(mChild, nChild);
      offset += PIPS_MPIgetSum( nChild, mpiComm );

      MPI_Barrier(mpiComm);
   }

   const std::string children_string = child_stream.str();
   const std::string all_children = PIPS_MPIallgatherString( children_string, mpiComm );
   if( my_rank == 0 )
      out << all_children;

   /// border.bl_mat | Blmat | offset | children ///
   int mlink, nlink;
   this->Blmat->getSize(mlink, nlink);
   if( mlink > 0 )
   {
      assert( border_left.mat_link );
      int mBdl, nBdl; border_left.mat_link->getSize(mBdl, nBdl);
      assert( mBdl == mlink );

      // for each row r do:
      for( int r = 0; r < mlink; r++ )
      {
         MPI_Barrier(mpiComm);
         std::ostringstream link_row_stream;

         if( my_rank == 0 )
         {
            dynamic_cast<const SparseGenMatrix&>(*border_left.mat_link).writeToStreamDenseRow(link_row_stream, r);
            link_row_stream << "|\t";
            dynamic_cast<const SparseGenMatrix*>(this->Blmat)->writeToStreamDenseRow(link_row_stream, r);
            for( int i = 0; i < original_offset; ++i )
               link_row_stream << "\t";
         }

         for( size_t it = 0; it < children.size(); it++ )
            children[it]->Blmat->writeToStreamDenseRow( link_row_stream, r );

         const std::string children_link_row = link_row_stream.str();
         const std::string link_row = PIPS_MPIallgatherString( children_link_row, mpiComm );

         if( my_rank == 0 )
         {
            out << link_row;
            out << "\n";
         }
         MPI_Barrier(mpiComm);
      }

      std::ostringstream dashed_line_stream;
      if( my_rank == 0 )
      {
         dynamic_cast<const SparseGenMatrix&>(*border_left.mat_link).writeDashedLineToStream(dashed_line_stream);
         dashed_line_stream << "|\t";
      }
      writeDashedLineToStream( dashed_line_stream, original_offset );

      if( my_rank == 0 )
      {
         const std::string dashed_line = dashed_line_stream.str();
         out << dashed_line << "\n";
      }
   }
   if( iAmDistrib )
      MPI_Barrier(mpiComm);
}

void StochGenMatrix::writeDashedLineToStream( std::ostream& out, int offset ) const
{
   const int my_rank = PIPS_MPIgetRank(mpiComm);

   MPI_Barrier(mpiComm);
   std::ostringstream link_row_stream;

   if( my_rank == 0 )
   {
      dynamic_cast<const SparseGenMatrix*>(this->Blmat)->writeDashedLineToStream(link_row_stream);
      for( int i = 0; i < offset; ++i )
         link_row_stream << "-\t";
   }

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->Blmat->writeDashedLineToStream( link_row_stream );

   const std::string children_link_row = link_row_stream.str();
   const std::string link_row = PIPS_MPIallgatherString( children_link_row, mpiComm );

   if( my_rank == 0 )
   {
      out << link_row;
      out << "\n";
   }
   MPI_Barrier(mpiComm);

}

void StochGenMatrix::writeToStreamDense(std::ostream& out, int offset) const
{
   assert( hasSparseMatrices() );
   assert( children.size() != 0 );

   const int original_offset = offset;
   const int my_rank = PIPS_MPIgetRank(mpiComm);

   if( iAmDistrib )
      MPI_Barrier(mpiComm);

   int mBmat, nBmat; Bmat->getSize(mBmat, nBmat);

   assert( Bmat->isKindOf(kSparseGenMatrix) );

   /// Bmat ///
   if( my_rank == 0 )
   {
      for( int i = 0; i < mBmat; ++i )
      {
         for( int j = 0; j < offset; ++j )
            out << "\t";
         dynamic_cast<const SparseGenMatrix*>(this->Bmat)->writeToStreamDenseRow(out, i);
         out << "\n";
      }
   }

   /// write children ///
   /// offset | Child ///
   std::ostringstream child_stream{};
   if( !amatEmpty() )
      assert( children.size() ==  0 );

   for( size_t it = 0; it < children.size(); it++ )
   {
      MPI_Barrier(mpiComm);
      children[it]->writeToStreamDenseChild( child_stream, offset );

      int mChild, nChild;
      children[it]->getSize(mChild, nChild);
      offset += PIPS_MPIgetSum( nChild, mpiComm );

      MPI_Barrier(mpiComm);
   }

   const std::string children_string = child_stream.str();
   const std::string all_children = PIPS_MPIallgatherString( children_string, mpiComm );
   if( my_rank == 0 )
      out << all_children;

   /// Blmat | offset | children ///
   int mlink, nlink;
   this->Blmat->getSize(mlink, nlink);
   if( mlink > 0 )
   {
      // for each row r do:
      for( int r = 0; r < mlink; r++ )
      {
         MPI_Barrier(mpiComm);
         std::ostringstream link_row_stream;

         if( my_rank == 0 )
         {
            dynamic_cast<const SparseGenMatrix*>(this->Blmat)->writeToStreamDenseRow(link_row_stream, r);
            for( int i = 0; i < original_offset; ++i )
               link_row_stream << "\t";
         }

         for( size_t it = 0; it < children.size(); it++ )
            children[it]->Blmat->writeToStreamDenseRow( link_row_stream, r );

         const std::string children_link_row = link_row_stream.str();
         const std::string link_row = PIPS_MPIallgatherString( children_link_row, mpiComm );

         if( my_rank == 0 )
         {
            out << link_row;
            out << "\n";
         }
         MPI_Barrier(mpiComm);
      }
   }
   if( iAmDistrib )
      MPI_Barrier(mpiComm);
}

/** writes child matrix blocks, offset indicates the offset between A and B block. */
void StochGenMatrix::writeToStreamDenseChild(std::ostream& out, int offset) const
{
   if( Bmat->isKindOf(kStochGenMatrix) )
      dynamic_cast<StochGenMatrix*>(Bmat)->writeToStreamDense( out, offset);
   /// Border.mat | Amat | offset | Bmat ///
   else
   {
      assert( hasSparseMatrices() );
      assert( PIPS_MPIgetRank(mpiComm) == 0 );

      int nB, mB;
      Bmat->getSize(mB, nB);
      int nA, mA;
      Amat->getSize(mA, nA);
      assert( mB == mA );

      for( int row = 0; row < mB; ++row )
      {
         Amat->writeToStreamDenseRow(out, row);
         for( int i = 0; i < offset; ++i )
            out << "\t";
         Bmat->writeToStreamDenseRow(out, row);
         out << "\n";
      }
   }
}

/** returns a string containing the linking-row rowidx of the children. */
void StochGenMatrix::writeToStreamDenseRowLink(std::ostream& out, int rowidx) const
{
   assert( !Blmat->isKindOf(kStringGenMatrix) );
   dynamic_cast<const SparseGenMatrix&>(*Blmat).writeToStreamDenseRow( out, rowidx );
}

void StochGenMatrix::writeMPSformatRows(std::ostream& out, int rowType, OoqpVector* irhs) const
{
   assert( hasSparseMatrices() );

   const int myRank = PIPS_MPIgetRank(mpiComm);
   std::string rt;
   if( rowType == 0 )
      rt = "E";
   else if( rowType == 1 )
      rt = "L";
   else if( rowType == 2 )
      rt = "G";
   else
      assert(0);

   StochVector* irhsStoch;
   if( irhs )
      irhsStoch = dynamic_cast<StochVector*>(irhs);
   else
      irhsStoch = nullptr;

   int m, n;
   if( myRank == 0)
   {
      // A_0 block:
      this->Bmat->getSize(m, n);
      for(int i=0; i<m; i++)
      {
         if( !irhs || (irhsStoch && dynamic_cast<SimpleVector*>(irhsStoch->vec)->elements()[i] != 0.0) )
            out << " " << rt << " row_" << rt << "_" << "R" << "_" << i << "\n";
      }
      // linking rows:
      if( Blmat )
      {
         this->Blmat->getSize(m, n);
         for(int i=0; i<m; i++)
         {
            if( !irhs || (irhsStoch && dynamic_cast<SimpleVector*>(irhsStoch->vecl)->elements()[i] != 0.0) )
               out << " " << rt << " row_" << rt << "_" << "L" << "_" << i << "\n";
         }
      }
   }
   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->Amat->getSize(m, n);
      for(int i = 0; i < m; i++)
      {
         if( !irhs || (irhsStoch && dynamic_cast<SimpleVector*>(irhsStoch->children[it]->vec)->elements()[i] != 0.0) )
            out << " "<< rt << " row_" << rt << "_" << it << "_" << i << "\n";
      }
   }
}

void StochGenMatrix::initTransposed(bool dynamic)
{
   assert( hasSparseMatrices() );
   dynamic_cast<SparseGenMatrix*>(Bmat)->initTransposed(dynamic);
   dynamic_cast<SparseGenMatrix*>(Blmat)->initTransposed(dynamic);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->initTransposedChild(dynamic);
}

void StochGenMatrix::deleteTransposed()
{
   assert( hasSparseMatrices() );

   dynamic_cast<SparseGenMatrix*>(Amat)->deleteTransposed();
   dynamic_cast<SparseGenMatrix*>(Bmat)->deleteTransposed();
   dynamic_cast<SparseGenMatrix*>(Blmat)->deleteTransposed();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->deleteTransposed();
}

void StochGenMatrix::initTransposedChild(bool dynamic)
{
   assert( hasSparseMatrices() );
   dynamic_cast<SparseGenMatrix*>(Amat)->initTransposed(dynamic);
   dynamic_cast<SparseGenMatrix*>(Bmat)->initTransposed(dynamic);

   if( Blmat != nullptr )
      dynamic_cast<SparseGenMatrix*>(Blmat)->initTransposed(dynamic);
}

int StochGenMatrix::numberOfNonZeros() const
{
   assert( hasSparseMatrices() );
   unsigned int nnz = 0;

   for(size_t it = 0; it < children.size(); it++)
      nnz += children[it]->numberOfNonZeros();

   if( iAmDistrib )
      PIPS_MPIgetSumInPlace(nnz, mpiComm);

   nnz += Amat->numberOfNonZeros() + Bmat->numberOfNonZeros() + Blmat->numberOfNonZeros();

   return nnz;
}

void StochGenMatrix::getNnzPerRow(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent)
{
   assert( hasSparseMatrices() );
   StochVectorBase<int>& nnzVecStoch = dynamic_cast<StochVectorBase<int>&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   SimpleVectorBase<int>* nnzvecl = nullptr;

   dynamic_cast<SparseGenMatrix*>(Bmat)->addNnzPerRow(*(nnzVecStoch.vec));

   if( linkParent != nullptr )
      dynamic_cast<SparseGenMatrix*>(Amat)->addNnzPerRow(*(nnzVecStoch.vec));

   /* with linking constraints? */
   if( nnzVecStoch.vecl || linkParent )
   {
      assert(nnzVecStoch.vecl == nullptr || linkParent == nullptr);

      if( linkParent )
         nnzvecl = dynamic_cast<SimpleVectorBase<int>*>(linkParent);
      else
         nnzvecl = dynamic_cast<SimpleVectorBase<int>*>(nnzVecStoch.vecl);

      if( linkParent != nullptr || iAmSpecial(iAmDistrib, mpiComm) )
         dynamic_cast<SparseGenMatrix*>(Blmat)->addNnzPerRow(*nnzvecl);
   }


   for( size_t it = 0; it < children.size(); it++ )
     children[it]->getNnzPerRow(*(nnzVecStoch.children[it]), nnzvecl);

   // distributed, with linking constraints, and at root?
   if( iAmDistrib && nnzVecStoch.vecl != nullptr && linkParent == nullptr )
   {
      PIPS_MPIsumArrayInPlace(nnzvecl->elements(), nnzvecl->length(), mpiComm);
   }
}

void StochGenMatrix::getNnzPerCol(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent)
{
   assert( hasSparseMatrices() );
   StochVectorBase<int>& nnzVecStoch = dynamic_cast<StochVectorBase<int>&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   SimpleVectorBase<int>* vec = dynamic_cast<SimpleVectorBase<int>*>(nnzVecStoch.vec);

   if( iAmSpecial(iAmDistrib, mpiComm) || linkParent != nullptr )
   {
      dynamic_cast<SparseGenMatrix*>(Bmat)->addNnzPerCol(*(vec));

      int blm, bln;
      Blmat->getSize(blm, bln);

      /* with linking constraints? */
      if( blm > 0 )
         dynamic_cast<SparseGenMatrix*>(Blmat)->addNnzPerCol(*vec);
   }

   // not at root?
   if( linkParent != nullptr )
      dynamic_cast<SparseGenMatrix*>(Amat)->addNnzPerCol(*linkParent);
   else
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->getNnzPerCol(*(nnzVecStoch.children[it]), vec);
   }

   // distributed and at root?
   if( iAmDistrib && linkParent == nullptr )
   {
      PIPS_MPIsumArrayInPlace(vec->elements(), vec->length(), mpiComm);
   }
}

void StochGenMatrix::getRowMinMaxVec(bool getMin, bool initializeVec, const OoqpVector* col_scale_,
      OoqpVector& minmax_ )
{
   assert( amatEmpty() );

   StochVector& minmax = dynamic_cast<StochVector&>(minmax_);

   const bool scale = col_scale_;
   const bool has_linking = minmax.vecl;

   const StochVector* const col_scale = scale ? dynamic_cast<const StochVector*>(col_scale_) : nullptr;
   const OoqpVector* const col_scale_vec = scale ? col_scale->getLinkingVecNotHierarchicalTop() : nullptr;

   Bmat->getRowMinMaxVec(getMin, initializeVec, col_scale_vec, *minmax.vec);

   if( has_linking )
   {
      if( initializeVec )
      {
         if( getMin )
            minmax.vecl->setToConstant(std::numeric_limits<double>::max());
         else
            minmax.vecl->setToZero();
      }

      if( iAmSpecial(iAmDistrib, mpiComm) )
         Blmat->getRowMinMaxVec(getMin, false, col_scale_vec, *minmax.vecl);
   }

   assert( minmax.children.size() == children.size() );

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->getRowMinMaxVecChild(getMin, initializeVec, scale ? col_scale->children[it] : nullptr, *(minmax.children[it]), minmax.vecl);

   if( iAmDistrib )
   {
      if( getMin )
        PIPS_MPIminArrayInPlace(dynamic_cast<SimpleVector&>(*minmax.vecl).elements(), minmax.vecl->length(), mpiComm);
      else
        PIPS_MPImaxArrayInPlace(dynamic_cast<SimpleVector&>(*minmax.vecl).elements(), minmax.vecl->length(), mpiComm);
   }
}

void StochGenMatrix::getRowMinMaxVecChild(bool getMin, bool initializeVec, const OoqpVector* col_scale_,
      OoqpVector& minmax_, OoqpVector* minmax_linking_cons)
{
   assert( children.empty() );
   StochVector& minmax = dynamic_cast<StochVector&>(minmax_);

   const bool scale = col_scale_;
   const bool has_linking = minmax_linking_cons;

   const StochVector* const col_scale = dynamic_cast<const StochVector*>(col_scale_);

   const OoqpVector* const col_scale_vec = scale ? col_scale->vec : nullptr;
   const OoqpVector* const col_scale_linkingvar_vec = scale ? col_scale->getLinkingVecNotHierarchicalTop() : nullptr;

   Bmat->getRowMinMaxVec(getMin, initializeVec, col_scale_vec, *(minmax.vec));

   if( !amatEmpty() )
   {
      assert( !Bmat->isKindOf(kStochGenMatrix) );
      Amat->getRowMinMaxVec(getMin, false, col_scale_linkingvar_vec, *(minmax.vec));
   }

   if( has_linking )
      Blmat->getRowMinMaxVec(getMin, false, col_scale_vec, *minmax_linking_cons);
}

void StochGenMatrix::getColMinMaxVec(bool getMin, bool initializeVec, const OoqpVector* rowScaleVec_, OoqpVector& minmaxVec_ )
{
   assert( amatEmpty() );
   StochVector& minmaxVec = dynamic_cast<StochVector&>(minmaxVec_);
   const StochVector* rowScaleVec = dynamic_cast<const StochVector*>(rowScaleVec_);

   int blm, bln;
   Blmat->getSize(blm, bln);
   const bool scale = rowScaleVec;
   const bool has_linking = blm > 0;

   const OoqpVector* row_scale_vec = scale ? rowScaleVec->vec : nullptr;
   const OoqpVector* row_scale_link = scale ? rowScaleVec->vecl : nullptr;

   if( minmaxVec.vec == minmaxVec.getLinkingVecNotHierarchicalTop() )
      Bmat->getColMinMaxVec(getMin, initializeVec, row_scale_vec, *minmaxVec.getLinkingVecNotHierarchicalTop() );
   else
      Bmat->getColMinMaxVec(getMin, false, row_scale_vec, *minmaxVec.getLinkingVecNotHierarchicalTop() );

   if( has_linking )
      Blmat->getColMinMaxVec(getMin, false, row_scale_link, *minmaxVec.getLinkingVecNotHierarchicalTop() );

   assert(minmaxVec.children.size() == children.size());

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->getColMinMaxVecChild(getMin, initializeVec, scale ? rowScaleVec->children[it] : nullptr, row_scale_link, *(minmaxVec.children[it]) );

   if( iAmDistrib )
   {
      SimpleVector& mvec = dynamic_cast<SimpleVector&>(*minmaxVec.vec);
      if( getMin )
         PIPS_MPIminArrayInPlace( mvec.elements(), mvec.length(), mpiComm );
      else
         PIPS_MPImaxArrayInPlace( mvec.elements(), mvec.length(), mpiComm );
   }
}

void StochGenMatrix::getColMinMaxVecChild( bool getMin, bool initializeVec, const OoqpVector* rowScale_, const OoqpVector* rowScaleParent,
      OoqpVector& minmaxVec_ )
{
   assert( children.empty() );
   StochVector& minmaxVec = dynamic_cast<StochVector&>(minmaxVec_);

   int blm, bln;
   Blmat->getSize(blm, bln);
   const bool scale = rowScale_;
   const bool has_linking = blm > 0;

   const StochVector* rowScale = dynamic_cast<const StochVector*>(rowScale_);
   const OoqpVector* row_scale_vec = scale ? rowScale->vec : nullptr;

   Bmat->getColMinMaxVec(getMin, initializeVec, row_scale_vec, *minmaxVec.vec);

   if( !amatEmpty() )
   {
      assert( !Bmat->isKindOf(kStochGenMatrix) );
      Amat->getColMinMaxVec(getMin, false, row_scale_vec, *minmaxVec.getLinkingVecNotHierarchicalTop() );
   }

   if( has_linking )
      Blmat->getColMinMaxVec(getMin, false, rowScaleParent, *minmaxVec.vec );
}



void StochGenMatrix::addRowSums( OoqpVector& sumVec, OoqpVector* linkParent ) const
{
   assert( false && "TODO : hierarchical version");
   assert( hasSparseMatrices() );

   StochVector& sumVecStoch = dynamic_cast<StochVector&>(sumVec);
   SimpleVector* mvecl = nullptr;

   // assert tree compatibility
   assert(sumVecStoch.children.size() == children.size());

   Bmat->addRowSums(*sumVecStoch.vec);

   // not at root?
   if( linkParent != nullptr )
      Amat->addRowSums(*sumVecStoch.vec);

   /* with linking constraints? */
   if( sumVecStoch.vecl || linkParent )
   {
      assert(sumVecStoch.vecl == nullptr || linkParent == nullptr);

      // at root?
      if( linkParent == nullptr )
         mvecl = dynamic_cast<SimpleVector*>(sumVecStoch.vecl);
      else
         mvecl = dynamic_cast<SimpleVector*>(linkParent);

      if( linkParent != nullptr || iAmSpecial(iAmDistrib, mpiComm) )
         Blmat->addRowSums(*mvecl);
   }

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->addRowSums(*(sumVecStoch.children[it]), mvecl);

   // distributed, with linking constraints, and at root?
   if( iAmDistrib && sumVecStoch.vecl != nullptr && linkParent == nullptr )
   {
      assert(mvecl != nullptr);

      // sum up linking constraints vectors
      const int locn = mvecl->length();
      double* buffer = new double[locn];

      MPI_Allreduce(mvecl->elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      mvecl->copyFromArray(buffer);

      delete[] buffer;
   }
}

void StochGenMatrix::addColSums( OoqpVector& sumVec, OoqpVector* linkParent ) const
{
   assert( false && "TODO : hierarchical version");
   assert( hasSparseMatrices() );

   StochVector& sumVecStoch = dynamic_cast<StochVector&>(sumVec);

   // assert tree compatibility
   assert(sumVecStoch.children.size() == children.size());

   SimpleVector* const mvec = dynamic_cast<SimpleVector*>(sumVecStoch.vec);

   if( iAmSpecial(iAmDistrib, mpiComm) || linkParent != nullptr )
      Bmat->addColSums(*mvec);

   int blm, bln;
   Blmat->getSize(blm, bln);

   /* with linking constraints? */
   if( blm > 0 && (iAmSpecial(iAmDistrib, mpiComm) || linkParent != nullptr) )
      Blmat->addColSums(*mvec);

   // not at root?
   if( linkParent != nullptr )
      Amat->addColSums(*linkParent);
   else
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->addColSums(*(sumVecStoch.children[it]), mvec);
   }

   // distributed and at root?
   if( iAmDistrib && linkParent == nullptr )
   {
      const int locn = mvec->length();
      double* const entries = mvec->elements();
      double* buffer = new double[locn];

      MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      mvec->copyFromArray(buffer);

      delete[] buffer;
   }
}

void StochGenMatrix::initStaticStorageFromDynamic(const OoqpVectorBase<int>& rowNnzVec, const OoqpVectorBase<int>& colNnzVec,
  const OoqpVectorBase<int>* rowLinkVec, const OoqpVectorBase<int>* colParentVec)
{
   assert( hasSparseMatrices() );

   const StochVectorBase<int>& rowNnzVecStoch = dynamic_cast<const StochVectorBase<int>&>(rowNnzVec);
   const StochVectorBase<int>& colNnzVecStoch = dynamic_cast<const StochVectorBase<int>&>(colNnzVec);

   assert(rowNnzVecStoch.children.size() == colNnzVecStoch.children.size());

   const SimpleVectorBase<int>* const rowvec = dynamic_cast<const SimpleVectorBase<int>*>(rowNnzVecStoch.vec);
   const SimpleVectorBase<int>* const colvec = dynamic_cast<const SimpleVectorBase<int>*>(colNnzVecStoch.vec);

   const SimpleVectorBase<int>* const rowlink = dynamic_cast<const SimpleVectorBase<int>*>(rowNnzVecStoch.vecl);
   assert(rowvec); assert(colvec);

   dynamic_cast<SparseGenMatrix*>(Amat)->initStaticStorageFromDynamic(*rowvec, colParentVec); // initialized with colVec == nullptr for parent
   dynamic_cast<SparseGenMatrix*>(Bmat)->initStaticStorageFromDynamic(*rowvec, colvec);

   // at root?
   if( colParentVec == nullptr )
   {
      assert(rowLinkVec == nullptr);

      if( rowlink != nullptr )
         dynamic_cast<SparseGenMatrix*>(Blmat)->initStaticStorageFromDynamic(*rowlink, colvec);

      for( size_t it = 0; it < children.size(); it++ )
         children[it]->initStaticStorageFromDynamic(*(rowNnzVecStoch.children[it]), *(colNnzVecStoch.children[it]), rowlink, colvec);
   }
   else
   {
      assert( children.size() == 0 );
      if( rowLinkVec != nullptr)
         dynamic_cast<SparseGenMatrix*>(Blmat)->initStaticStorageFromDynamic(*rowLinkVec, colvec);
   }

}

void StochGenMatrix::freeDynamicStorage()
{
   assert( hasSparseMatrices() );

   dynamic_cast<SparseGenMatrix*>(Amat)->freeDynamicStorage();
   dynamic_cast<SparseGenMatrix*>(Bmat)->freeDynamicStorage();
   dynamic_cast<SparseGenMatrix*>(Blmat)->freeDynamicStorage();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->freeDynamicStorage();
}

void StochGenMatrix::recomputeSize( StochGenMatrix* parent )
{
   m = 0;
   n = 0;

   if( Bmat->isKindOf(kStochGenMatrix) )
   {
      assert( children.empty() );
      assert( inner_leaf );

      dynamic_cast<StochGenMatrix*>(Bmat)->recomputeSize();
   }

   if( !inner_root )
      Bmat->getSize(m, n);

   assert( m >= 0 );
   assert( n >= 0 );

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->recomputeSize( this );

      int m_child, n_child;
      children[it]->getSize(m_child, n_child);

      m += m_child;
      n += n_child;
   }

   if( !parent )
   {
      int bl_mat_m = 0; int bl_mat_n = 0;
      Blmat->getSize(bl_mat_m, bl_mat_n);
      m += bl_mat_m;
   }
}

void StochGenMatrix::updateKLinkConsCount(std::vector<int>& linkCount) const
{
   assert( hasSparseMatrices() );

   if( Blmat == nullptr )
      return;

   int m, n;

   Blmat->getSize(m, n);
   assert(m > 0);
   assert(linkCount.size() == size_t(m));

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         assert(children[it]->Blmat);
         dynamic_cast<SparseGenMatrix*>(children[it]->Blmat)->updateNonEmptyRowsCount(linkCount);
      }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &linkCount[0], m, MPI_INT, MPI_SUM, mpiComm);
}

void StochGenMatrix::updateKLinkVarsCount(std::vector<int>& link_block_count) const
{
   assert(hasSparseMatrices());
   int m, n;
   Bmat->getSize(m, n);

   if( n == 0 )
      return;

   assert(link_block_count.size() == size_t(n));

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         dynamic_cast<SparseGenMatrix*>(children[it]->Amat)->getTranspose().updateNonEmptyRowsCount(link_block_count);
         dynamic_cast<SparseGenMatrix*>(children[it]->Amat)->deleteTransposed();
      }

   if( iAmDistrib )
      PIPS_MPIsumArrayInPlace( link_block_count, mpiComm );
}

void StochGenMatrix::get2LinkStartBlocksAndCountsNew(std::vector<int>& block_start, std::vector<int>& block_count) const
{
   assert(hasSparseMatrices());
   block_start.clear();
   block_count.clear();

   if( Blmat == nullptr )
      return;

   int m, n; Blmat->getSize(m, n);

   if( m == 0 )
      return;
   assert( m > 0 );

   const int n_blocks = children.size();
   block_count.resize(m);
   std::fill(block_count.begin(), block_count.end(), 0);

   /* init with max + 1 and max - 1 for allreduce later */
   block_start.resize(m);
   std::fill(block_start.begin(), block_start.end(), n_blocks);
   std::vector<int> block_end(m, -1);

   std::vector<bool> is_2_link(m, false);

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         assert(children[it]->Blmat);
         dynamic_cast<const SparseGenMatrix*>(children[it]->Blmat)->updateNonEmptyRowsCountNew(it, block_count, block_start, block_end);
      }

   if( iAmDistrib )
   {
      // TODO : one can filter the non-candidates locally first on all processes
      PIPS_MPIminArrayInPlace(block_start, mpiComm);
      PIPS_MPImaxArrayInPlace(block_end, mpiComm);
      PIPS_MPIsumArrayInPlace(block_count, mpiComm);
   }

   for( int it = 0; it < m; ++it )
   {
      const int start = block_start[it];
      const int end = block_end[it];

      assert( start == n_blocks || (0 <= start && start < n_blocks) );
      if( start == n_blocks )
         assert( end == -1 );
      else
         assert( start <= end && end < n_blocks );

      if( end == start + 1 )
         assert( block_count[it] == 2 );
      else
      {
         /* not a consecutive 2 link */
         block_start[it] = -1;
      }
   }
}

std::vector<int> StochGenMatrix::get2LinkStartBlocks() const
{
   assert(hasSparseMatrices());
   if( Blmat == nullptr )
      return std::vector<int>();

   int m, n;
   Blmat->getSize(m, n);

   if( m == 0 )
      return std::vector<int>();
   assert(m > 0);

   std::vector<int> linkBlockCount(m, 0);
   std::vector<int> linkBlockStart(m, -1);
   std::vector<int> linkBlockEnd(m, -1);
   std::vector<bool> is2link(m, false);

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         assert(children[it]->Blmat);
         dynamic_cast<const SparseGenMatrix*>(children[it]->Blmat)->updateNonEmptyRowsCount(it, linkBlockCount, linkBlockStart, linkBlockEnd);
      }

   if( iAmDistrib )
      PIPS_MPIsumArrayInPlace(linkBlockCount, mpiComm);

   /* filter out process local two links already */
   for( int i = 0; i < m; i++ )
   {
      assert(linkBlockEnd[i] == -1 || linkBlockStart[i] <= linkBlockEnd[i]);

      if( linkBlockCount[i] == 2 && (linkBlockEnd[i] - linkBlockStart[i]) == 1 )
      {
         assert(linkBlockStart[i] >= 0 && linkBlockEnd[i] >= 0);
         is2link[i] = true;
      }
   }


   if( iAmDistrib )
   {
      // find 2-links between different processes

      const int size = PIPS_MPIgetSize(mpiComm);
      assert(size > 0);

      // 1. allgather number of local 2-link candidates
      std::vector<int> localCandsIdx;
      std::vector<int> localCandsBlock;
      std::vector<int> candsPerProc(size, -1);

      /* a local candidate is a linking row that appears in exactly two blocks, starts on this process but where the second block is is not stored on this process */
      for( int i = 0; i < m; i++ )
         if( linkBlockCount[i] == 2 && linkBlockStart[i] >= 0 && linkBlockEnd[i] == -1 )
         {
            assert(!is2link[i]);

            localCandsIdx.push_back(i);
            assert(unsigned(linkBlockStart[i]) < children.size());

            localCandsBlock.push_back(linkBlockStart[i]);
         }

      const int localcount = localCandsIdx.size();

      PIPS_MPIallgather(&localcount, 1, &candsPerProc[0], 1, mpiComm);

#ifndef NDEBUG
      for( size_t i = 0; i < candsPerProc.size(); i++ )
         assert(candsPerProc[i] >= 0);
#endif

      // 2. allgatherv 2-link candidates
      std::vector<int> displacements(size + 1, 0);
      for( int i = 1; i <= size; i++ )
         displacements[i] = candsPerProc[i - 1] + displacements[i - 1];

      const int nAllCands = displacements[size];

      std::vector<int> allCandsRow(nAllCands, -1);
      std::vector<int> allCandsBlock(nAllCands, -1);

      MPI_Allgatherv(&localCandsIdx[0], localcount, MPI_INT, &allCandsRow[0], &candsPerProc[0],
            &displacements[0], MPI_INT, mpiComm);

      MPI_Allgatherv(&localCandsBlock[0], localcount, MPI_INT, &allCandsBlock[0], &candsPerProc[0],
            &displacements[0], MPI_INT, mpiComm);

#ifndef NDEBUG
      for( size_t i = 0; i < allCandsRow.size(); i++ )
         assert(allCandsRow[i] >= 0 && allCandsRow[i] < m && allCandsBlock[i] >= 0 && allCandsBlock[i] < static_cast<int>(children.size()));
#endif


      // 3. check which candidates are indeed 2-links
      std::vector<int> blocksHash(m, -1);
      for( int i = 0; i < size - 1; i++ )
      {
         // hash
         for( int j = displacements[i]; j < displacements[i + 1]; j++ )
            blocksHash[allCandsRow[j]] = allCandsBlock[j];

         // compare with next
         for( int j = displacements[i + 1]; j < displacements[i + 2]; j++ )
         {
            assert(allCandsBlock[j] > 0);
            const int candRow = allCandsRow[j];
            const int candBlock = allCandsBlock[j];
            if( blocksHash[candRow] >= 0 )
            {
               assert(blocksHash[candRow] != candBlock);

               const int startBlock = std::min(blocksHash[candRow], candBlock);
               const int endBlock = std::max(blocksHash[candRow], candBlock);

               assert(startBlock >= 0 && endBlock >= 0);

               if( endBlock - startBlock != 1 )
                  continue;

               assert( !is2link[candRow] );
               is2link[candRow] = true;

               // start block owned by this MPI process?
               if( !children[startBlock]->isKindOf(kStochGenDummyMatrix) )
               {
                  assert(children[endBlock]->isKindOf(kStochGenDummyMatrix));
                  linkBlockStart[candRow] = startBlock;
               }
               else
               {
                  assert(children[startBlock]->isKindOf(kStochGenDummyMatrix));
                  linkBlockStart[candRow] = -1;
               }
            }
         }

         // un-hash
         for( int j = displacements[i]; j < displacements[i + 1]; j++ )
            blocksHash[allCandsRow[j]] = -1;
      }
   }

   // correct block identifier
   for( int i = 0; i < m; i++ )
      if( !is2link[i] )
         linkBlockStart[i] = -1;

   if( iAmDistrib )
      PIPS_MPImaxArrayInPlace(linkBlockStart, mpiComm);

   return linkBlockStart;
}


void StochGenMatrix::permuteLinkingVars(const std::vector<unsigned int>& permvec)
{
   assert(hasSparseMatrices());
   if( Blmat )
      dynamic_cast<SparseGenMatrix*>(Blmat)->permuteCols(permvec);

   dynamic_cast<SparseGenMatrix*>(Bmat)->permuteCols(permvec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->permuteLinkingVarsChild(permvec);
}

void StochGenMatrix::permuteLinkingVarsChild(const std::vector<unsigned int>& permvec)
{
   assert(hasSparseMatrices());

   dynamic_cast<SparseGenMatrix*>(Amat)->permuteCols(permvec);

   assert(children.size() == 0);
}

void StochGenMatrix::permuteLinkingCons(const std::vector<unsigned int>& permvec)
{
   assert(hasSparseMatrices());
   if( Blmat )
      dynamic_cast<SparseGenMatrix*>(Blmat)->permuteRows(permvec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->permuteLinkingCons(permvec);
}


void StochGenMatrix::updateTransposed()
{
   assert(hasSparseMatrices());
   dynamic_cast<SparseGenMatrix*>(Amat)->updateTransposed();
   dynamic_cast<SparseGenMatrix*>(Bmat)->updateTransposed();
   dynamic_cast<SparseGenMatrix*>(Blmat)->updateTransposed();

  for( size_t it = 0; it < children.size(); it++ )
     children[it]->updateTransposed();
}

/* check whether root node date is same in all processes
 *
 * todo: check this
 * todo: make better use of std::vector and iterators
 * root node data is Amat (empty), Bmat and Blmat of root node. Children not checked.
 */
bool StochGenMatrix::isRootNodeInSync() const
{
   assert(hasSparseMatrices());
   bool in_sync = true;

   assert( Amat );
   assert( Bmat );
   assert( Blmat );

   /* no need to check if not distributed or not in root node */
   if( !iAmDistrib || children.size() == 0)
      return in_sync;

   const int my_rank = PIPS_MPIgetRank(mpiComm);

   const SparseGenMatrix& amat_sp = dynamic_cast<const SparseGenMatrix&>(*Amat);
   const SparseGenMatrix& bmat_sp = dynamic_cast<const SparseGenMatrix&>(*Bmat);
   const SparseGenMatrix& blmat_sp = dynamic_cast<const SparseGenMatrix&>(*Blmat);
   /* if matrix has static storage */
   if( amat_sp.getStorageHandle().notNil() || bmat_sp.getStorageHandle().notNil() ||
	    blmat_sp.getStorageHandle().notNil() )
   {
      assert( amat_sp.getStorageHandle().notNil() );
      assert( bmat_sp.getStorageHandle().notNil() );
      assert( blmat_sp.getStorageHandle().notNil() );

      /* since we are in root node Amat should look as follows */
      assert( amat_sp.getStorageRef().len == 0 );
      assert( amat_sp.getStorageRef().n == -1 );

      /* static storage */
      const int lenght_entries_bmat = bmat_sp.getStorageRef().len;
      const int length_columns_bmat = bmat_sp.getStorageRef().len;
      const int lenght_rowoffest_bmat = bmat_sp.getStorageRef().m + 1;

      const int lenght_entries_blmat = blmat_sp.getStorageRef().len;
      const int length_columns_blmat = blmat_sp.getStorageRef().len;
      const int lenght_rowoffest_blmat = blmat_sp.getStorageRef().m + 1;

      const long long count_row_cols = length_columns_bmat + lenght_rowoffest_bmat + length_columns_blmat + lenght_rowoffest_blmat;
      const long long count_entries = lenght_entries_bmat + lenght_entries_blmat;

      assert( count_row_cols < std::numeric_limits<int>::max());
      assert( count_entries < std::numeric_limits<int>::max());

      std::vector<double> sendbuf_entries(count_entries, 0.0);
      std::vector<double> recvbuf_entries(count_entries, 0.0);

      std::vector<int> sendbuf_row_col(count_row_cols, 0);
      std::vector<int> recvbuf_row_col(count_row_cols, 0);

      /* fill Bmat into send buffers */
      const double * M = bmat_sp.getStorageRef().M;
      const int * krowM = bmat_sp.getStorageRef().krowM;
      const int * jColM = bmat_sp.getStorageRef().jcolM;

      std::copy(M, M + lenght_entries_bmat, sendbuf_entries.begin());

      std::copy(krowM, krowM + lenght_rowoffest_bmat, sendbuf_row_col.begin());
      std::copy(jColM, jColM + lenght_entries_bmat, sendbuf_row_col.begin() + lenght_rowoffest_bmat);

      /* fill Blmat into send buffers */
      const double * Ml = blmat_sp.getStorageRef().M;
      const int * krowMl = blmat_sp.getStorageRef().krowM;
      const int * jColMl = blmat_sp.getStorageRef().jcolM;

      std::copy(Ml, Ml + lenght_entries_blmat, sendbuf_entries.begin() + lenght_entries_bmat);
      std::copy(krowMl, krowMl + lenght_rowoffest_blmat, sendbuf_row_col.begin() + lenght_rowoffest_bmat + lenght_entries_bmat);
      std::copy(jColMl, jColMl + lenght_entries_blmat, sendbuf_row_col.begin() + lenght_rowoffest_bmat + lenght_entries_bmat + lenght_rowoffest_blmat);

      /* Reduce Bmat and Blmat buffers */
      MPI_Allreduce(&sendbuf_entries[0], &recvbuf_entries[0], static_cast<int>(count_entries), MPI_DOUBLE, MPI_MAX, mpiComm);

      MPI_Allreduce(&sendbuf_row_col[0], &recvbuf_row_col[0], static_cast<int>(count_row_cols), MPI_INT, MPI_MAX, mpiComm);

      /* check recvbuf_entries */
      for( int i = 0; i < count_entries; ++i )
      {
         if( !PIPSisEQ(sendbuf_entries[i], recvbuf_entries[i]) )
         {
            /* someone else had a higher value here */
            if(my_rank == 0)
               std::cout << "matrix entries out of sync\n";
            in_sync = false;
            break;
         }
      }

      for( int i = 0; i < count_row_cols; ++i )
      {
         if( !PIPSisEQ(sendbuf_row_col[i], recvbuf_row_col[i]) )
         {
            /* someone else had a higher value here */
            if(my_rank == 0)
               std::cout << "matrix indices (col or row) out of sync\n";
            in_sync = false;
         }
      }
   }

   /* if stoch mat has dynamic storage also check that */
   if(bmat_sp.hasDynamicStorage() || blmat_sp.hasDynamicStorage())
   {
      assert(bmat_sp.hasDynamicStorage());
      assert(blmat_sp.hasDynamicStorage());

      const SparseStorageDynamic& Bmat_dyn = bmat_sp.getStorageDynamicRef();
      const SparseStorageDynamic& Blmat_dyn = blmat_sp.getStorageDynamicRef();

      /* dynamic storage */
      int bmat_dyn_len = 0;
      for(int i = 0; i < Bmat_dyn.getM(); ++i)
         bmat_dyn_len += (Bmat_dyn.getRowPtr(i).end - Bmat_dyn.getRowPtr(i).start);

      int blmat_dyn_len = 0;
      for(int i = 0; i < Blmat_dyn.getM(); ++i)
         blmat_dyn_len += (Blmat_dyn.getRowPtr(i).end - Blmat_dyn.getRowPtr(i).start);

      const int lenght_entries_bmat_dynamic = bmat_dyn_len;
      const int length_columns_bmat_dynamic = bmat_dyn_len;
      const int lenght_rowoffest_bmat_dynamic = Bmat_dyn.getM() + 1;

      const int lenght_entries_blmat_dynamic = blmat_dyn_len;
      const int length_columns_blmat_dynamic = blmat_dyn_len;
      const int lenght_rowoffest_blmat_dynamic = Blmat_dyn.getM() + 1;

      const long long count_row_cols_dyn = length_columns_bmat_dynamic + 2 * lenght_rowoffest_bmat_dynamic + length_columns_blmat_dynamic + 2 * lenght_rowoffest_blmat_dynamic;
      const long long count_entries_dyn = lenght_entries_bmat_dynamic + lenght_entries_blmat_dynamic;

      assert( count_row_cols_dyn < std::numeric_limits<int>::max());
      assert( count_entries_dyn < std::numeric_limits<int>::max());

      std::vector<double> sendbuf_entries_dynamic(count_entries_dyn, 0.0);
      std::vector<double> recvbuf_entries_dynamic(count_entries_dyn, 0.0);

      std::vector<int> sendbuf_row_col_dynamic(count_row_cols_dyn, 0);
      std::vector<int> recvbuf_row_coldynamic(count_row_cols_dyn, 0);;

      /* fill Bmat into send buffers */
      const double * M = Bmat_dyn.getMat();
      const int * jColM = Bmat_dyn.getJcolM();

      int count_entries = 0;
      int count_row_col = 0;

      /* entries Bmat into double array */
      for(int i = 0; i < Bmat_dyn.getM(); ++i)
      {
         for(int j = Bmat_dyn.getRowPtr(i).start; j < Bmat_dyn.getRowPtr(i).end; ++j)
         {
            sendbuf_entries_dynamic[count_entries] = M[j];
            count_entries++;
         }
      }
      assert(count_entries == lenght_entries_bmat_dynamic);

      /* row pointers Bmat into int array */
      for(int i = 0; i < lenght_rowoffest_bmat_dynamic; ++i)
      {
         sendbuf_row_col_dynamic[count_row_col] = Bmat_dyn.getRowPtr()->start;
         sendbuf_row_col_dynamic[count_row_col + 1] = Bmat_dyn.getRowPtr()->end;
         count_row_col += 2;
      }
      assert(count_row_col == 2 * lenght_rowoffest_bmat_dynamic);

      /* col indices of Bmat into int array */
      for(int i = 0; i < Bmat_dyn.getM(); ++i)
      {
         for(int j = Bmat_dyn.getRowPtr(i).start; j < Bmat_dyn.getRowPtr(i).end; ++j)
         {
            sendbuf_row_col_dynamic[count_row_col] = jColM[j];
            count_row_col++;
         }
      }
      assert(count_row_col == 2 * lenght_rowoffest_bmat_dynamic + length_columns_bmat_dynamic);

      /* fill Blmat into send buffers */
      const double * Ml = Blmat_dyn.getMat();
      const int * jColMl = Blmat_dyn.getJcolM();

      /* entries Blmat into double array */
      for(int i = 0; i < Blmat_dyn.getM(); ++i)
      {
         for(int j = Blmat_dyn.getRowPtr(i).start; j < Blmat_dyn.getRowPtr(i).end; ++j)
         {
            sendbuf_entries_dynamic[count_entries] = Ml[j];
            count_entries++;
         }
      }
      assert(count_entries == lenght_entries_bmat_dynamic + lenght_entries_blmat_dynamic);

      /* row pointers Blmat into int array */
      for(int i = 0; i < lenght_rowoffest_blmat_dynamic; ++i)
      {
         assert(2 * lenght_rowoffest_bmat_dynamic + lenght_entries_bmat_dynamic + 2 * i + 1 < count_row_cols_dyn);
         sendbuf_row_col_dynamic[count_row_col] = Blmat_dyn.getRowPtr()->start;
         sendbuf_row_col_dynamic[count_row_col + 1] = Blmat_dyn.getRowPtr()->end;
         count_row_col += 2;
      }
      assert(count_row_col == 2 * lenght_rowoffest_bmat_dynamic + length_columns_bmat_dynamic + 2 * lenght_rowoffest_blmat_dynamic);

      /* col indices of Bmat into int array */
      for(int i = 0; i < Blmat_dyn.getM(); ++i)
      {
         for(int j = Blmat_dyn.getRowPtr(i).start; j < Blmat_dyn.getRowPtr(i).end; ++j)
         {
            sendbuf_row_col_dynamic[count_row_col] = jColMl[j];
            count_row_col++;
         }
      }
      assert(count_row_col == 2 * lenght_rowoffest_bmat_dynamic + length_columns_bmat_dynamic + 2 * lenght_rowoffest_blmat_dynamic + length_columns_blmat_dynamic);

      /* Reduce Bmat and Blmat buffers */
      MPI_Allreduce(&sendbuf_entries_dynamic[0], &recvbuf_entries_dynamic[0], static_cast<int>(count_entries_dyn), MPI_DOUBLE, MPI_MAX, mpiComm);

      MPI_Allreduce(&sendbuf_row_col_dynamic[0], &recvbuf_row_coldynamic[0], static_cast<int>(count_row_cols_dyn), MPI_INT, MPI_MAX, mpiComm);

      /* check recvbuf_entries */
      for( int i = 0; i < count_entries_dyn; ++i )
      {
         if( !PIPSisEQ(sendbuf_entries_dynamic[i], recvbuf_entries_dynamic[i]) )
         {
            /* someone else had a higher value here */
            if(my_rank == 0)
               std::cout << "matrix entries in dynamic storage out of sync\n";
            in_sync = false;
         }
      }
      for( int i = 0; i < count_row_cols_dyn; ++i )
      {
         if( !PIPSisEQ(sendbuf_row_col_dynamic[i], recvbuf_row_coldynamic[i]) )
         {
            /* someone else had a higher value here */
            if(my_rank == 0)
               std::cout << "matrix indices in dynamic storage out of sync\n";
            in_sync = false;
         }
      }
   }


   return in_sync;
}

/* Find correct matrices to append row to
 *  Can only be called in root node
 *
 *  Child -1 is parent
 *
 * @return rowindex (in specified block row) of newly appended row
 */
int StochGenMatrix::appendRow( const StochGenMatrix& matrix_row, int child, int row, bool linking ) 
{
   assert(hasSparseMatrices());
  // todo: check that matrix is in correct format
  assert( matrix_row.children.size() == children.size() );
  assert( children.size() != 0 );
  assert( -1 <= child && child <= (int) children.size() );

  int index_row;

  // append row to all matrices necessary
  // todo maybe this can be done nicer - maybe we can just recursively call some method also on the dummies 
  if(linking)
  {
    index_row = dynamic_cast<SparseGenMatrix*>(Blmat)->appendRow( dynamic_cast<const SparseGenMatrix&>(*matrix_row.Blmat), row );

    for(unsigned int i = 0; i < children.size(); ++i)
    { 
      if( !children[i]->isKindOf(kStochGenDummyMatrix) )
      {
        assert( !matrix_row.children[i]->isKindOf(kStochGenDummyMatrix) );
        dynamic_cast<SparseGenMatrix*>(children[i]->Blmat)->appendRow( dynamic_cast<const SparseGenMatrix&>(*matrix_row.children[i]->Blmat), row);
      }
    }
  }
  else
  {
    if(child != -1)
    {
      index_row = dynamic_cast<SparseGenMatrix&>(*children[child]->Amat).appendRow( dynamic_cast<const SparseGenMatrix&>(*matrix_row.children[child]->Amat), row );
#ifndef NDEBUG
      const int index_row1 = dynamic_cast<SparseGenMatrix&>(*children[child]->Bmat).appendRow(
            dynamic_cast<const SparseGenMatrix&>(*matrix_row.children[child]->Bmat), row );
#else
      dynamic_cast<SparseGenMatrix&>(*children[child]->Bmat).appendRow( dynamic_cast<SparseGenMatrix&>(*matrix_row.children[child]->Bmat), row );
#endif
      assert(index_row1 == index_row);
    }
    else
      index_row = dynamic_cast<SparseGenMatrix&>(*Bmat).appendRow( dynamic_cast<const SparseGenMatrix&>(*matrix_row.Bmat), row );
  }

  return index_row;
};

/* y += alpha RowAt(child, row, linking) */
void StochGenMatrix::axpyWithRowAt( double alpha, StochVector* y, SimpleVector* y_linking, int child, int row, bool linking) const
{
   assert(hasSparseMatrices());
   assert( y );
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(y->children.size() == children.size());

   /* go through all available children and calculate y += alpha * rowAt(row) */
   if( linking )
   {
      assert( Blmat );
      if( y_linking )
         dynamic_cast<const SparseGenMatrix*>(Blmat)->axpyWithRowAt(alpha, *y_linking, row);
      else
      {
         assert( y->vec );
         dynamic_cast<const SparseGenMatrix*>(Blmat)->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->vec), row);
      }

      for( unsigned int i = 0; i < children.size(); ++i )
      {
         if( !children[i]->isKindOf(kStochGenDummyMatrix) )
         {
            assert(children[i]->Blmat);
            assert(y->children[i]->vec);
            dynamic_cast<const SparseGenMatrix*>(children[i]->Blmat)->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->children[i]->vec), row);
         }
      }
   }
   else
   {
      if( child == -1 )
      {
         assert(Bmat);
         if( y_linking )
            dynamic_cast<const SparseGenMatrix*>(Bmat)->axpyWithRowAt(alpha, *y_linking, row);
         else
         {
            assert(y->vec);
            dynamic_cast<const SparseGenMatrix*>(Bmat)->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->vec), row);
         }
      }
      else
      {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);

         assert(y->children[child]->vec);
         dynamic_cast<const SparseGenMatrix*>(children[child]->Bmat)->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->children[child]->vec), row);

         if( y_linking )
            dynamic_cast<const SparseGenMatrix*>(children[child]->Amat)->axpyWithRowAt(alpha, *y_linking, row);
         else
         {
            assert(y->vec);
            dynamic_cast<const SparseGenMatrix*>(children[child]->Amat)->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->vec), row);
         }
      }
   }
}

void StochGenMatrix::axpyWithRowAtPosNeg( double alpha, StochVector* y_pos, SimpleVector* y_link_pos,
      StochVector* y_neg, SimpleVector* y_link_neg, int child, int row, bool linking ) const
{
   assert(hasSparseMatrices());
   assert( y_pos && y_neg );
   assert( (y_link_neg && y_link_pos) || (!y_link_neg && !y_link_pos) );
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(y_neg->children.size() == children.size());
   assert(y_pos->children.size() == children.size());

   /* go through all available children and calculate y += alpha * rowAt(row) */
   if( linking )
   {
      assert( Blmat );
      if( y_link_pos )
         dynamic_cast<const SparseGenMatrix*>(Blmat)->axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
      else
      {
         assert( y_pos->vec );
         assert( y_neg->vec );
         dynamic_cast<const SparseGenMatrix*>(Blmat)->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->vec), dynamic_cast<SimpleVector&>(*y_neg->vec), row);
      }

      for( unsigned int i = 0; i < children.size(); ++i )
      {
         if( !children[i]->isKindOf(kStochGenDummyMatrix) )
         {
            assert(children[i]->Blmat);
            assert(y_pos->children[i]->vec);
            assert(y_neg->children[i]->vec);
            dynamic_cast<const SparseGenMatrix*>(children[i]->Blmat)->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->children[i]->vec), dynamic_cast<SimpleVector&>(*y_neg->children[i]->vec), row);
         }
      }
   }
   else
   {
      if( child == -1 )
      {
         assert(Bmat);
         if( y_link_pos )
            dynamic_cast<const SparseGenMatrix*>(Bmat)->axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
         else
         {
            assert(y_pos->vec);
            assert(y_neg->vec);
            dynamic_cast<const SparseGenMatrix*>(Bmat)->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->vec), dynamic_cast<SimpleVector&>(*y_neg->vec), row);
         }
      }
      else
      {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);

         assert(y_pos->children[child]->vec);
         assert(y_neg->children[child]->vec);
         dynamic_cast<const SparseGenMatrix*>(children[child]->Bmat)->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->children[child]->vec), dynamic_cast<SimpleVector&>(*y_neg->children[child]->vec), row);

         if( y_link_pos )
            dynamic_cast<const SparseGenMatrix*>(children[child]->Amat)->axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
         else
         {
            assert(y_pos->vec);
            assert(y_neg->vec);
            dynamic_cast<const SparseGenMatrix*>(children[child]->Amat)->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->vec), dynamic_cast<SimpleVector&>(*y_neg->vec), row);
         }
      }
   }
}

double StochGenMatrix::localRowTimesVec(const StochVector &vec, int child, int row, bool linking) const
{
   assert(hasSparseMatrices());
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(vec.children.size() == children.size());

   double res = 0.0;

   /* go through all available children and multiply the vec times row in submatrix */
   if( linking )
   {
      assert(Blmat);
      assert(vec.vec);
      res += dynamic_cast<const SparseGenMatrix*>(Blmat)->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);

      for( unsigned int i = 0; i < children.size(); ++i )
      {
         if( !children[i]->isKindOf(kStochGenDummyMatrix) )
         {
            assert(children[i]->Blmat);
            assert(vec.children[i]->vec);
            res += dynamic_cast<const SparseGenMatrix*>(children[i]->Blmat)->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.children[i]->vec), row);
         }
      }
   }
   else
   {
      if( child == -1 )
      {
         assert(Bmat);
         assert(vec.vec);
         res += dynamic_cast<const SparseGenMatrix*>(Bmat)->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);
      }
      else
      {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);
         assert(vec.vec);
         assert(vec.children[child]->vec);
         res += dynamic_cast<const SparseGenMatrix*>(children[child]->Amat)->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);
         res += dynamic_cast<const SparseGenMatrix*>(children[child]->Bmat)->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.children[child]->vec), row);
      }
   }

   return res;
}

// TODO specify border and left from sData...
BorderedGenMatrix* StochGenMatrix::raiseBorder( int m_conss, int n_vars )
{
#ifndef NDEBUG
   int m_link, n_link;
   Blmat->getSize(m_link, n_link);
   assert(m_conss <= m_link && n_vars <= n_link);
#endif

   SparseGenMatrix* const A_left = dynamic_cast<SparseGenMatrix*>(Bmat)->shaveLeft(n_vars);

   SparseGenMatrix* const Bl_left_top = dynamic_cast<SparseGenMatrix*>(Blmat)->shaveLeft(n_vars);
   SparseGenMatrix* const bottom_left_block = dynamic_cast<SparseGenMatrix*>(Bl_left_top->shaveBottom(m_conss));

   SparseGenMatrix* const Bl_right_bottom = dynamic_cast<SparseGenMatrix*>(Blmat->shaveBottom(m_conss));

   StringGenMatrix* const border_bottom = new StringGenMatrix(false, Bl_right_bottom, nullptr, mpiComm);
   StringGenMatrix* const border_left = new StringGenMatrix(true, A_left, Bl_left_top, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->shaveBorder(m_conss, n_vars, border_left, border_bottom);

   m -= m_conss;
   n -= n_vars;

   BorderedGenMatrix* const bordered_matrix = new BorderedGenMatrix(this, border_left, border_bottom, bottom_left_block, mpiComm);
   StochGenMatrix* me = this;
   IotrAddRef( &me );

   assert(m >= 0 && n >= 0);

   return bordered_matrix;
}

void StochGenMatrix::shaveBorder( int m_conss, int n_vars, StringGenMatrix* border_left, StringGenMatrix* border_bottom )
{
   if( Bmat->isKindOf(kStochGenMatrix) )
   {
      assert( amatEmpty() );
      border_left->addChild( new StringGenMatrix( true, dynamic_cast<StringGenMatrix*>(dynamic_cast<StochGenMatrix*>(Bmat)->shaveLeftBorder( n_vars ) ), nullptr, mpiComm ) );
      border_bottom->addChild( new StringGenMatrix( false, dynamic_cast<StringGenMatrix*>(Blmat->shaveBottom( m_conss )), nullptr, mpiComm) );
   }
   else
   {
      assert( hasSparseMatrices() );
      assert( children.size() == 0 );
      assert( PIPS_MPIgetSize( mpiComm ) == 1 );

      SparseGenMatrix* const border_a_mat = dynamic_cast<SparseGenMatrix*>(Amat)->shaveLeft(n_vars);
      SparseGenMatrix* const border_bl_mat = dynamic_cast<SparseGenMatrix*>(dynamic_cast<SparseGenMatrix*>(Blmat)->shaveBottom(m_conss));

      StringGenMatrix* border_left_child = new StringGenMatrix(true, border_a_mat, nullptr, mpiComm);
      StringGenMatrix* border_bottom_child = new StringGenMatrix(false, border_bl_mat, nullptr, mpiComm);

      border_left->addChild(border_left_child);
      border_bottom->addChild(border_bottom_child);
   }
}

StringGenMatrix* StochGenMatrix::shaveLeftBorder( int n_vars )
{
   assert( children.size() > 0 );
   assert( hasSparseMatrices() );
   assert( amatEmpty() );

   SparseGenMatrix* const border_b_mat = dynamic_cast<SparseGenMatrix*>(Bmat)->shaveLeft(n_vars);
   SparseGenMatrix* const border_bl_mat = dynamic_cast<SparseGenMatrix*>(Blmat)->shaveLeft(n_vars);

   StringGenMatrix* border = new StringGenMatrix( true, border_b_mat, border_bl_mat, mpiComm );

   for( auto& child : children )
      border->addChild( dynamic_cast<StringGenMatrix*>( child->shaveLeftBorderChild(n_vars) ) );

   return border;
}

StringGenMatrix* StochGenMatrix::shaveLeftBorderChild( int n_vars )
{
   assert( children.empty() );

   if( Bmat->isKindOf(kStochGenMatrix) )
   {
      assert(amatEmpty());
      return new StringGenMatrix( true, dynamic_cast<StochGenMatrix*>(Bmat)->shaveLeftBorder(n_vars), nullptr, mpiComm );
   }
   else
      return new StringGenMatrix( true, dynamic_cast<SparseGenMatrix*>(Amat)->shaveLeft(n_vars), nullptr, mpiComm );
}

StringGenMatrix* StochGenMatrix::shaveLinkingConstraints( unsigned int n_conss )
{
   assert( hasSparseMatrices() );

   SparseGenMatrix* border_bl_mat = dynamic_cast<SparseGenMatrix*>(Blmat->shaveBottom(n_conss));
   StringGenMatrix* border = new StringGenMatrix(false, border_bl_mat, nullptr, mpiComm);

   if( children.size() == 0 )
      assert( PIPS_MPIgetSize(mpiComm) == 1 );
   for( auto& child : children )
   {
      StringGenMatrix* border_child = child->shaveLinkingConstraints(n_conss);
      border->addChild(border_child);
   }

   return border;
}

void StochGenMatrix::splitMatrix( const std::vector<int>& twolinks_start_in_block, const std::vector<unsigned int>& map_blocks_children, unsigned int n_links_in_root,
      const std::vector<MPI_Comm>& child_comms )
{
   const unsigned int n_curr_children = children.size();

   assert( hasSparseMatrices() );
   assert( n_curr_children == map_blocks_children.size() );
   assert( n_curr_children == twolinks_start_in_block.size() );
   assert( twolinks_start_in_block.back() == 0 );

   int nBl, m_links_left;
   Blmat->getSize(m_links_left, nBl);
   assert( std::accumulate( twolinks_start_in_block.begin(), twolinks_start_in_block.end(), 0 ) <= m_links_left );

   const unsigned int n_new_children = getNDistinctValues(map_blocks_children);
   std::vector<StochGenMatrix*> new_children(n_new_children);

   StringGenMatrix* Blmat_new = shaveLinkingConstraints( n_links_in_root );

   Blmat_new->combineChildrenInNewChildren( map_blocks_children, child_comms );

   SparseGenMatrix* Blmat_leftover = dynamic_cast<SparseGenMatrix*>(Blmat);

#ifndef NDEBUG
   int n_child_links_sum{0};
   const unsigned int n_links_orig = m_links_left;
#endif
   m_links_left -= n_links_in_root;
   /* for each future new child collect its children and add them to the new child
    * then shave off the linking constraints that stay at the new child's level
    */
   unsigned int m_links_so_far{0};
   unsigned int begin_curr_child_blocks{0};
   unsigned int end_curr_child_blocks{0};
   for( unsigned int i = 0; i < n_new_children; ++i )
   {
      while( end_curr_child_blocks != (n_curr_children - 1) &&
            map_blocks_children[end_curr_child_blocks] == map_blocks_children[end_curr_child_blocks + 1] )
         ++end_curr_child_blocks;

      const int n_links_for_child = std::accumulate( twolinks_start_in_block.begin() + begin_curr_child_blocks,
            twolinks_start_in_block.begin() + end_curr_child_blocks, 0 );
      const unsigned int n_blocks_for_child = end_curr_child_blocks - begin_curr_child_blocks + 1;

#ifndef NDEBUG
      n_child_links_sum += n_links_for_child;
#endif
      /* combine children in new StochGenMatrix Bmat */
      /* create root node with only Blmat */
      SparseGenMatrix* Blmat_child = Blmat_leftover;
      Blmat_leftover = dynamic_cast<SparseGenMatrix*>(Blmat_child->shaveBottom(m_links_left - n_links_for_child));

      if( child_comms[i] == MPI_COMM_NULL )
      {
         delete Blmat_child;
         Blmat_child = nullptr;
      }

      StochGenMatrix* Bmat = (child_comms[i] == MPI_COMM_NULL) ? nullptr :
            new StochGenMatrix( new SparseGenMatrix(0, 0, 0), new SparseGenMatrix(0, nBl, 0), Blmat_child, child_comms[i], false, true );

      /* shave off empty two link part from respective children and add them to the new root/remove them from the old root */
      for( unsigned int j = 0; j < n_blocks_for_child; ++j )
      {
         assert( m_links_left >= n_links_for_child );

         StochGenMatrix* child = children.front();
         children.erase(children.begin());

         if( child_comms[i] == MPI_COMM_NULL )
            assert( child->mpiComm == MPI_COMM_NULL );

         if( child->mpiComm != MPI_COMM_NULL )
         {
            dynamic_cast<SparseGenMatrix*>(child->Blmat)->dropNEmptyRowsTop( m_links_so_far );
            dynamic_cast<SparseGenMatrix*>(child->Blmat)->dropNEmptyRowsBottom( m_links_left - n_links_for_child );
         }

#ifndef NDEBUG
         if( child->mpiComm != MPI_COMM_NULL )
         {
            int blm,bln;
            child->Blmat->getSize(blm, bln);
            assert( blm == n_links_for_child );
         }
#endif
         if( Bmat )
            Bmat->AddChild(child);
//         else
//            delete child;
      }
      if( Bmat )
         Bmat->recomputeSize();

      /* create child holding the new Bmat and it's Blmat part */
      if( child_comms[i] == MPI_COMM_NULL )
      {
         assert( Blmat_new->children[i]->isKindOf(kStringGenDummyMatrix));
         delete Blmat_new->children[i];
         Blmat_new->children[i] = nullptr;
      }
      else
         assert( Blmat_new->children[i]->isKindOf(kStringGenMatrix) );

      new_children[i] = (child_comms[i] != MPI_COMM_NULL) ? new StochGenMatrix( new SparseGenMatrix(0, 0, 0), Bmat, Blmat_new->children[i], child_comms[i], true, false) :
            new StochGenDummyMatrix();

      ++end_curr_child_blocks;
      begin_curr_child_blocks = end_curr_child_blocks;
      m_links_left -= n_links_for_child;
      m_links_so_far += n_links_for_child;
   }

   assert( n_child_links_sum + n_links_in_root == n_links_orig );
   assert( children.size() == 0 );
   assert( m_links_left == 0 );

   /* exchange children and recompute sizes */
   children.insert(children.end(), new_children.begin(), new_children.end());

   Blmat_new->children.clear();
   Blmat = Blmat_new->mat;
   Blmat_new->mat = nullptr;
   delete Blmat_new;
   delete Blmat_leftover;

   recomputeSize();
}
