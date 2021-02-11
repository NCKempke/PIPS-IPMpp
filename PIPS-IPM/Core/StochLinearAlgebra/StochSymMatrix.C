#include "StochSymMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixTypes.h"
#include "BorderedSymMatrix.h"
#include "StringGenMatrix.h"
#include "pipsport.h"

#include <cassert>

StochSymMatrix::StochSymMatrix( SymMatrix* diag_, SparseGenMatrix* border_, MPI_Comm mpiComm_ ) :
   diag{ diag_ }, border{ border_ }, mpiComm{ mpiComm_ }, iAmDistrib( PIPS_MPIgetDistributed(mpiComm) )
{
   assert( diag );
   recomputeSize();
}

/**
 * This the constructor is usually called for the root node. In this case
 * parent is set to nullptr; the cross Hessian does not exist, so
 * border is set up to be an empty matrix. 
 *
 * If it is called to create a child, the calling code should call 
 *   this->AddChild(c)
 * 'AddChild' method correctly sets the parent and (re)creates an EMPTY
 * border with correct sizes.
 */
StochSymMatrix::StochSymMatrix(long long global_n, int local_n, int local_nnz,
			       MPI_Comm mpiComm_)
  : n(global_n), mpiComm(mpiComm_), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) )
{
  diag = new SparseSymMatrix(local_n, local_nnz);
  // the cross Hessian is nullptr for the root node; it may be also nullptr for 
  // children in the case when the Hessian does not have cross terms and
  // the children are created with this constructor. The border will be 
  // set up to correct sizes later for this case.
}

void StochSymMatrix::AddChild(StochSymMatrix* child)
{
   child->parent = this;
   assert( this->border == nullptr );

   if( !child->border )
      child->border = new SparseGenMatrix(child->diag->size(), this->diag->size(), 0);

   children.push_back(child);
}

StochSymMatrix::~StochSymMatrix()
{
   for(size_t it = 0; it < children.size(); it++)
      delete children[it];

   delete diag;
   delete border;
}

SymMatrix* StochSymMatrix::clone() const
{
   SymMatrix* diag_clone = diag->clone();
   SparseGenMatrix* border_clone = border ? dynamic_cast<SparseGenMatrix*>(border->cloneFull()) : nullptr;

   StochSymMatrix* clone = new StochSymMatrix( diag_clone, border_clone, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
   {
      StochSymMatrix* child = dynamic_cast<StochSymMatrix*>(children[it]->clone());
      clone->AddChild(child);
      clone->n += child->n;
   }

   assert( size() == clone->size() );
   return clone;
}

void StochSymMatrix::recomputeSize()
{
   assert( diag );
   n = 0;
   for( auto& child : children )
   {
      child->recomputeSize();

      n += child->size();
   }

   if( diag->isKindOf(kStochSymMatrix) )
      dynamic_cast<StochSymMatrix*>(diag)->recomputeSize();

   n += diag->size();
}


int StochSymMatrix::isKindOf( int type ) const
{
  return type == kStochSymMatrix || type == kSymMatrix;
}

void StochSymMatrix::getSize( long long& m_, long long& n_ ) const
{
  m_=n; n_=n;
}

void StochSymMatrix::getSize( int& m_, int& n_ ) const
{
  m_=n; n_=n;
}


long long StochSymMatrix::size() const
{
  return n;
}

/** y = beta * y + alpha * this * x 
 * 
 *           [ Q0*x0+ sum(Ri^T*xi) ]
 *           [        .            ]
 * this * x =[        .            ]
 *           [        .            ]
 *           [   Ri*x0 + Qi*xi     ]
 *
 * Here Qi are diagonal blocks, Ri are left bordering blocks
 */
void StochSymMatrix::mult( double beta,  OoqpVector& y_,
			    double alpha, const OoqpVector& x_ ) const
{
   const StochVector & x = dynamic_cast<const StochVector&>(x_);
   StochVector & y = dynamic_cast<StochVector&>(y_);

   if( 0.0 == alpha)
   {
      y.vec->scale( beta );
      return;
   }

   if( !parent )
   {
      if( PIPS_MPIiAmSpecial(iAmDistrib, mpiComm) )
         diag->mult( beta, *y.vec, alpha, *x.vec );
      else
         y.vec->setToZero();
   }
   else
      diag->mult( beta, *y.vec, alpha, *x.vec );

   // y0 = y0 + alpha * border^T * xi
   if( border )
   {
      assert( parent );
      assert( x.parent );
      assert( y.parent );

      border->transMult( 1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.vec );
      border->mult(1.0, *y.vec, alpha, *x.getLinkingVecNotHierarchicalTop() );
   }

   assert(y.children.size() == children.size() );
   assert(x.children.size() == children.size() );

   // recursively multiply the children
   for (size_t it = 0; it < children.size(); it++)
      children[it]->mult(beta, *(y.children[it]), alpha, *(x.children[it]));

   if( iAmDistrib && !parent )
      PIPS_MPIsumArrayInPlace( dynamic_cast<SimpleVector*>(y.vec)->elements(), y.vec->length(), mpiComm );
}

/** y = beta * y + alpha * this^T * x */
void StochSymMatrix::transMult ( double beta,  OoqpVector& y_,
				 double alpha, const OoqpVector& x_) const
{
  // We are symmetric, this^T = this, therefore call 'mult' method
  this->mult(beta, y_, alpha, x_);
}
  
/** the magnitude of the element in this matrix with largest absolute value.
   */
double StochSymMatrix::abmaxnorm() const
{
  double maxNorm = 0.0;

  for(size_t it = 0; it < children.size(); it++)
    maxNorm = std::max( maxNorm, children[it]->abmaxnorm() );

  if( iAmDistrib )
     PIPS_MPIgetMaxInPlace( maxNorm, mpiComm );

  maxNorm = std::max( maxNorm, diag->abmaxnorm() );
  if( border )
     maxNorm = std::max( maxNorm, border->abmaxnorm() );
  return maxNorm;
}

double StochSymMatrix::abminnormNonZero( double tol ) const
{
  double min = std::numeric_limits<double>::infinity();

  for(size_t it = 0; it < children.size(); it++)
    min = std::min( min, children[it]->abminnormNonZero(tol) );

  if( iAmDistrib )
     PIPS_MPIgetMinInPlace( min, mpiComm );

  min = std::min( min, diag->abminnormNonZero(tol) );
  if( border )
     min = std::min( min, diag->abminnormNonZero(tol) );
  return min;
}

void StochSymMatrix::writeToStreamDense(std::ostream& out) const
{
   const int rank = PIPS_MPIgetRank(mpiComm);
   const int world_size = PIPS_MPIgetSize(mpiComm);

   int m, n;
   int offset = 0;
   std::stringstream sout;
   MPI_Status status;
   int l;

   /* this is at the root node - thus there is no border */
   assert(this->border == nullptr );

   if( iAmDistrib )
      MPI_Barrier(mpiComm);

   if( iAmDistrib && rank > 0 )  // receive offset from previous process
      MPI_Recv(&offset, 1, MPI_INT, (rank - 1), 0, mpiComm, MPI_STATUS_IGNORE);
   else  //  !iAmDistrib || (iAmDistrib && rank == 0)
      this->diag->writeToStreamDense(sout);

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->writeToStreamDenseChild(sout, offset);
      children[it]->diag->getSize(m, n);
      offset += n;
   }

   if( iAmDistrib && rank > 0 )
   {
      std::string str = sout.str();
      // send string to rank ZERO to print it there:
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, rank, mpiComm);
      // send offset to next process:
      if( rank < world_size - 1 )
         MPI_Ssend(&offset, 1, MPI_INT, rank + 1, 0, mpiComm);
   }
   else if( !iAmDistrib )
      out << sout.str();
   else if( iAmDistrib && rank == 0 )
   {
      out << sout.str();
      MPI_Ssend(&offset, 1, MPI_INT, rank + 1, 0, mpiComm);

      for( int p = 1; p < world_size; p++ )
      {
         MPI_Probe(p, p, mpiComm, &status);
         MPI_Get_count(&status, MPI_CHAR, &l);
         char *buf = new char[l];
         MPI_Recv(buf, l, MPI_CHAR, p, p, mpiComm, &status);
         std::string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }
   }

   if( iAmDistrib )
      MPI_Barrier(mpiComm);
   std::cout << " done\n";
}

void StochSymMatrix::writeToStreamDenseChild(std::stringstream& out, int offset) const
{
   if( diag->isKindOf(kSparseSymMatrix) )
   {
      for( int r = 0; r < diag->size(); r++)
      {
         if( border )
            border->writeToStreamDenseRow(out, r);

         for(int i = 0; i < offset; i++)
            out <<'\t';

         dynamic_cast<const SparseSymMatrix&>(*diag).writeToStreamDenseRow(out, r);
         out << "\n";
      }
   }
   else
      dynamic_cast<const StochSymMatrix&>(*diag).writeToStreamDense(out);
}



void StochSymMatrix::getDiagonal( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  diag->getDiagonal( *vec.vec);

  for(size_t it=0; it<children.size(); it++)
    children[it]->getDiagonal(*vec.children[it]);
}

void StochSymMatrix::setToDiagonal( const OoqpVector& vec_ )
{
  const StochVector& vec = dynamic_cast<const StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  diag->setToDiagonal(*vec.vec);

  for(size_t it=0; it<children.size(); it++)
    children[it]->setToDiagonal(*vec.children[it]);
}

void StochSymMatrix::atPutDiagonal( int idiag, OoqpVector& v_ )
{
  StochVector& v = dynamic_cast<StochVector&>(v_);

  //check the tree compatibility
  int nChildren = children.size();
  assert(v.children.size() - nChildren==0);

  //check the node size compatibility
  assert(this->diag->size() == v.vec->length());

  diag->atPutDiagonal ( idiag, *v.vec);

  for (int it=0; it<nChildren; it++) 
    children[it]->atPutDiagonal( idiag, *v.children[it]);
}

void StochSymMatrix::fromGetDiagonal( int, OoqpVector& x_ )
{
   assert( "The value of the parameter is not supported!" );

   StochVector& x = dynamic_cast<StochVector&>(x_);
   assert(x.children.size() == children.size());

   diag->getDiagonal(*x.vec);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->getDiagonal(*x.children[it]);
}

void StochSymMatrix::symmetricScale( const OoqpVector& vec_ )
{
  const StochVector& vec = dynamic_cast<const StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  diag->symmetricScale(*vec.vec);

  for (size_t it=0; it<children.size(); it++) 
    children[it]->symmetricScale(*vec.children[it]);
}

void StochSymMatrix::columnScale( const OoqpVector& vec_ )
{
   const StochVector& vec = dynamic_cast<const StochVector&>(vec_);
   assert(children.size() == vec.children.size());

   diag->columnScale(*vec.vec);

   if( border )
      border->columnScale( *vec.getLinkingVecNotHierarchicalTop() );

   for (size_t it = 0; it < children.size(); it++)
      children[it]->columnScale(*vec.children[it]);
}

void StochSymMatrix::rowScale ( const OoqpVector& vec_ )
{
   const StochVector& vec = dynamic_cast<const StochVector&>(vec_);
   assert(children.size() == vec.children.size());

   diag->rowScale( *vec.vec );
   if( border )
      border->rowScale( *vec.vec );

   for (size_t it = 0; it < children.size(); it++)
      children[it]->rowScale(*vec.children[it]);
}


void StochSymMatrix::scalarMult( double num )
{
   diag->scalarMult(num);
   if( border )
      border->scalarMult(num);
   for (size_t it = 0; it < children.size(); it++)
      children[it]->scalarMult(num);
}


void StochSymMatrix::deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec, const OoqpVectorBase<int>* linkParent)
{
   const StochVectorBase<int>& nnzVecStoch = dynamic_cast<const StochVectorBase<int>&>(nnzVec);
   assert( children.size() == nnzVecStoch.children.size() );

   const SimpleVectorBase<int>* vec = dynamic_cast<const SimpleVectorBase<int>*>( nnzVecStoch.vec );
   assert( vec );

   const int n_old = n;
   for( size_t it = 0; it < children.size(); it++ )
      children[it]->deleteEmptyRowsCols(*nnzVecStoch.children[it], vec);

   if( linkParent != nullptr )
   {
      assert( border );
      assert( children.size() == 0 );
     // adapt border
      assert(dynamic_cast<const SimpleVectorBase<int>*>(linkParent));
      border->deleteEmptyRowsCols(*vec, *linkParent);
   }
   else
      assert( !border );

   assert( diag->isKindOf(kSparseSymMatrix) );
   dynamic_cast<SparseSymMatrix&>(*diag).deleteEmptyRowsCols(*vec);

   const int nnzs = vec->getNnzs();
   const int length = vec->length();

   assert( nnzs <= length );
   n -= (length - nnzs);

   assert( n <= n_old );
   assert( n >= 0 );

   const int n_changes = n_old - n;

   if( linkParent != nullptr )
   {
      assert( parent != nullptr );
      parent->n -= n_changes;
   }
}


int StochSymDummyMatrix::isKindOf( int type ) const
{ 
  return type == kStochSymDummyMatrix || type == kStochSymMatrix;
}

void StochSymMatrix::splitMatrix( const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms )
{
   const unsigned int n_curr_children = children.size();
   assert( n_curr_children == map_blocks_children.size() );

   const unsigned int n_new_children = getNDistinctValues(map_blocks_children);
   std::vector<StochSymMatrix*> new_children(n_new_children);

   unsigned int begin_curr_child_blocks{0};
   unsigned int end_curr_child_blocks{0};
   for( unsigned int i = 0; i < n_new_children; ++i )
   {
      while( end_curr_child_blocks != (n_curr_children - 1) &&
            map_blocks_children[end_curr_child_blocks] == map_blocks_children[end_curr_child_blocks + 1] )
         ++end_curr_child_blocks;

      const unsigned int n_blocks_for_child = end_curr_child_blocks - begin_curr_child_blocks + 1;

      StochSymMatrix* diag = (child_comms[i] == MPI_COMM_NULL) ? nullptr :
            new StochSymMatrix( 0, 0, 0, child_comms[i] );

      /* shave off empty two link part from respective children and add them to the new root/remove them from the old root */
      for( unsigned int j = 0; j < n_blocks_for_child; ++j )
      {
         StochSymMatrix* child = children.front();
         children.erase(children.begin());

         if( child_comms[i] == MPI_COMM_NULL )
            assert( child->mpiComm == MPI_COMM_NULL );

         if( diag )
            diag->AddChild(child);
         else
            delete child;
      }
      if( diag )
         diag->recomputeSize();

      /* create child holding the new Bmat and it's Blmat part */
      new_children[i] = (child_comms[i] != MPI_COMM_NULL) ? new StochSymMatrix( diag, nullptr, child_comms[i] ) : new StochSymDummyMatrix();
      if( child_comms[i] != MPI_COMM_NULL )
         dynamic_cast<StochSymMatrix*>(new_children[i]->diag)->parent = new_children[i];

      ++end_curr_child_blocks;
      begin_curr_child_blocks = end_curr_child_blocks;
   }
   assert( children.size() == 0 );

   /* exchange children and recompute sizes */
   children.insert(children.end(), new_children.begin(), new_children.end());

   for( auto& child : children )
   {
      child->recomputeSize();
      child->parent = this;
   }
   this->recomputeSize();
}

BorderedSymMatrix* StochSymMatrix::raiseBorder(int n_vars)
{
   assert( parent == nullptr );
   assert( border == nullptr );

   StringGenMatrix* border_vertical = shaveBorder(n_vars);

   BorderedSymMatrix* const border_layer = new BorderedSymMatrix(this, border_vertical, new SparseSymMatrix(n_vars, 0, false), mpiComm);

   StochSymMatrix* me = this;
   IotrAddRef(&me);

   assert(n >= 0);

   return border_layer;
}

StringGenMatrix* StochSymMatrix::shaveBorder( int n_vars )
{
   assert( diag->isKindOf(kSparseSymMatrix) );

   SparseGenMatrix* const border_top_left = parent ? new SparseGenMatrix(0,n_vars,0) : dynamic_cast<SparseSymMatrix*>(diag)->shaveSymLeftBottom(n_vars);
   StringGenMatrix* const border_vertical = new StringGenMatrix(true, border_top_left, nullptr, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
      border_vertical->addChild( children[it]->shaveBorder2(n_vars) );

   n -= n_vars;
   return border_vertical;
}

StringGenMatrix* StochSymMatrix::shaveBorder2(int n_vars)
{
   n -= n_vars;

   if( border )
   {
      SparseGenMatrix* const border_block = border->shaveLeft(n_vars);
      return new StringGenMatrix(true, border_block, nullptr, mpiComm);
   }
   else
   {
      assert( children.size() == 0 );
      assert( diag->isKindOf(kStochSymMatrix) );

      StochSymMatrix& diags = dynamic_cast<StochSymMatrix&>(*diag);
      return new StringGenMatrix(true, diags.shaveBorder(n_vars), nullptr, mpiComm);
   }
}
