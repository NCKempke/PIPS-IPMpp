#include "StochSymMatrix.h"
#include "StochVector.h"
#include "DoubleMatrixTypes.h"
#include "BorderedSymMatrix.h"
#include "StringGenMatrix.h"
#include "pipsport.h"

#include <cassert>

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
  :n(global_n), mpiComm(mpiComm_), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) ), parent(nullptr)
{
  diag = new SparseSymMatrix(local_n, local_nnz);
  // the cross Hessian is nullptr for the root node; it may be also nullptr for 
  // children in the case when the Hessian does not have cross terms and
  // the children are created with this constructor. The border will be 
  // set up to correct sizes later for this case.
  border = nullptr;
}

StochSymMatrix::StochSymMatrix( long long global_n,
				int nrows, int diag_nnz, 
				int nbordercols, int border_nnz,
				MPI_Comm mpiComm_)
  :n(global_n), mpiComm(mpiComm_), iAmDistrib( PIPS_MPIgetDistributed(mpiComm) ), parent(nullptr)
{
   diag = new SparseSymMatrix(nrows, diag_nnz);
   border = new SparseGenMatrix(nrows, nbordercols, border_nnz);
}

void StochSymMatrix::AddChild(StochSymMatrix* child)
{
  child->parent=this;

  if( !border )
    child->border = new SparseGenMatrix(child->diag->size(), this->diag->size(), 0);

  children.push_back(child);
}

StochSymMatrix::~StochSymMatrix()
{
   for(size_t it = 0; it < children.size(); it++)
      delete children[it];

   if( diag )
      delete diag;
   if( border )
      delete border;
}

// TODO : does not clone the border..
StochSymMatrix* StochSymMatrix::clone() const
{
   const int local_n = diag->getStorageRef().n;
   const int local_nnz = diag->getStorageRef().len;

   StochSymMatrix* clone = new StochSymMatrix(n, local_n, local_nnz, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
   {
      StochSymMatrix* child = children[it]->clone();
      clone->AddChild(child);
   }

   return clone;
}

int StochSymMatrix::isKindOf( int type ) const
{
  return type == kStochSymMatrix || type == kSymMatrix;
}

void StochSymMatrix::atPutDense( int /* row */, int /* col */,
				 double * /* A */, int /* lda */,
				 int /* rowExtent */,
				 int /* colExtent */ )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::fromGetDense( int row, int col, double * A, int lda,
		   int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::symAtPutSpRow( int row, 
				    double A[], int lenA, int jcolA[],
				    int& info )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::fsymAtPutSpRow( int row, 
				     double A[], int lenA, int jcolA[],
				     int& info )
{
  assert( "Not implemented" && 0 );
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

void StochSymMatrix::symAtPutSubmatrix( int destRow, int destCol,
					DoubleMatrix& M,
					int srcRow, int srcCol,
					int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::fromGetSpRow( int row, int col,
				   double A[], int lenA, int irowA[], int& nnz,
				   int rowExtent, int& info )
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::atPutZeros( int row, int col, int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
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

   size_t nChildren = children.size();
   assert(y.children.size() == nChildren );
   assert(x.children.size() == nChildren );

   assert(this->diag->size() == y.vec->length());
   assert(this->diag->size() == x.vec->length());

   SimpleVector & yvec = dynamic_cast<SimpleVector&>(*y.vec);
   const SimpleVector & xvec = dynamic_cast<const SimpleVector&>(*x.vec);

   if(0.0 == alpha)
   {
      y.vec->scale( beta );
      return;
   }

   const bool iAmRoot = (parent == nullptr);
   const int my_rank = PIPS_MPIgetRank( mpiComm );
   const bool iAmSpecial = (my_rank == 0);

   if( iAmRoot )
   {
      // y0 = beta * y0 + alpha * Q0 * x0
      if(iAmSpecial)
         diag->mult( beta, yvec, alpha, xvec );
      else
         yvec.setToZero();
   }
   else
      // yi = beta * yi + alpha * Qi * xi
      diag->mult( beta, yvec, alpha, xvec );

   // y0 = y0 + alpha * border^T * xi
   if( border )
   {
      assert( !iAmRoot );
      assert( x.parent );
      assert( y.parent );

      border->transMult( 1.0, *x.parent->vec, alpha, xvec );
      border->mult(1.0, yvec, alpha, *x.parent->vec);
   }

   // recursively multiply the children
   for (size_t it = 0; it < nChildren; it++)
      children[it]->mult(beta, *(y.children[it]), alpha, *(x.children[it]));

   if( iAmDistrib && nChildren > 0 )
      PIPS_MPIsumArrayInPlace( yvec.elements(), yvec.length(), mpiComm );

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

void StochSymMatrix::writeToStream(ostream& out) const
{
  assert( "Not implemented" && 0 );
}

void StochSymMatrix::randomizePSD(double * seed)
{
  assert( "Not implemented" && 0 );
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
   std::cout << " done " << std::endl;
}

void StochSymMatrix::writeToStreamDenseChild(stringstream& out, int offset) const
{
   assert( border != nullptr );
   int m_diag, m_border, n;
   this->diag->getSize(m_diag, n);
   this->border->getSize(m_border, n);

   assert( m_diag == m_border );

   for(int r = 0; r < m_diag; r++)
   {
      this->border->writeToStreamDenseRow(out, r);

      for(int i = 0; i < offset; i++)
         out <<'\t';

      this->diag->writeToStreamDenseRow(out, r);
      out << std::endl;
   }
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

void StochSymMatrix::fromGetDiagonal( int idiag, OoqpVector& x_ )
{
   assert("The value of the parameter is not supported!" && idiag==0);

   StochVector& x = dynamic_cast<StochVector&>(x_);
   assert(x.children.size() == children.size());

   diag->getDiagonal(*x.vec);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->getDiagonal(*x.children[it]);
}

void StochSymMatrix::putSparseTriple( int irow[], int len, int jcol[], 
				      double A[], int& info )
{
  assert("Not implemented!" && 0);
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

  for (size_t it = 0; it < children.size(); it++)
    children[it]->columnScale(*vec.children[it]);
}

void StochSymMatrix::rowScale ( const OoqpVector& vec_ )
{
  const StochVector& vec = dynamic_cast<const StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  diag->rowScale(*vec.vec);

  for (size_t it=0; it<children.size(); it++) 
    children[it]->rowScale(*vec.children[it]);
}


void StochSymMatrix::scalarMult( double num )
{
  diag->scalarMult(num);
  for (size_t it=0; it<children.size(); it++) 
    children[it]->scalarMult(num);
}


void StochSymMatrix::deleteEmptyRowsCols(const OoqpVectorBase<int>& nnzVec, const OoqpVectorBase<int>* linkParent)
{
   const StochVectorBase<int>& nnzVecStoch = dynamic_cast<const StochVectorBase<int>&>(nnzVec);
   assert(children.size() == nnzVecStoch.children.size());

   const SimpleVectorBase<int>* const vec = dynamic_cast<const SimpleVectorBase<int>*>(nnzVecStoch.vec);
   assert(vec);

   const int n_old = n;

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->deleteEmptyRowsCols(*nnzVecStoch.children[it], vec);

   if( linkParent != nullptr )
   {
      assert( children.size() == 0 );
     // adapt border
      assert(dynamic_cast<const SimpleVectorBase<int>*>(linkParent));
      border->deleteEmptyRowsCols(*vec, *linkParent);
   }
   else
      assert( border == nullptr );

   diag->deleteEmptyRowsCols(*vec);

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
  return type==kStochSymDummyMatrix;
}

/*StochSymMatrix::StochSymMatrix(const vector<StochSymMatrix*> &blocks) 
{
  mpiComm = blocks[0]->mpiComm;
  n = blocks[0]->n;
  id = blocks[0]->id;

  vector<SparseSymMatrix*> v(blocks.size());
  for(size_t i = 0; i < blocks.size(); i++) 
    v[i] = blocks[i]->mat;

  mat = new SparseSymMatrix(v);

  border = nullptr;

}
*/

BorderedSymMatrix* StochSymMatrix::raiseBorder(int n_vars)
{
   assert( parent == nullptr );
   assert( border == nullptr );

   SparseGenMatrix* const border_top_left = diag->shaveSymLeftBottom(n_vars);
   StringGenMatrix* const border_vertical = new StringGenMatrix(true, border_top_left, nullptr, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
   {
      StringGenMatrix* border_vertical_child;
      children[it]->shaveBorder(n_vars, border_vertical_child);
      border_vertical->addChild(border_vertical_child);
   }

   n -= n_vars;

   BorderedSymMatrix* const border_layer = new BorderedSymMatrix(this, border_vertical, new SparseSymMatrix(n_vars, 0, false), mpiComm);

   StochSymMatrix* me = this;
   IotrAddRef(&me);

   assert(n >= 0);

   return border_layer;
}

void StochSymMatrix::shaveBorder(int n_vars, StringGenMatrix*& border_vertical)
{
   assert( border );
   SparseGenMatrix* const border_block = border->shaveLeft(n_vars);
   border_vertical = new StringGenMatrix(true, border_block, nullptr, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
   {
      assert("deactivated for now" && 0);

      StringGenMatrix* border_child;
      children[it]->shaveBorder(n_vars, border_child);

      border_vertical->addChild(border_child);
   }
}
