#include "StochGenMatrix.h"
#include "DoubleMatrixTypes.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "pipsdef.h"
#include "pipsport.h"
#include <limits>
#include <algorithm>

StochGenMatrix::StochGenMatrix(long long global_m, long long global_n, int A_m,
      int A_n, int A_nnz, int B_m, int B_n, int B_nnz, MPI_Comm mpiComm_) :
      m(global_m), n(global_n), mpiComm(mpiComm_), iAmDistrib(PIPS_MPIgetDistributed(mpiComm))
{
   Amat = new SparseGenMatrix(A_m, A_n, A_nnz);
   Bmat = new SparseGenMatrix(B_m, B_n, B_nnz);
   Blmat = new SparseGenMatrix(0, 0, 0);
}

StochGenMatrix::StochGenMatrix(long long global_m, long long global_n, int A_m,
      int A_n, int A_nnz, int B_m, int B_n, int B_nnz, int Bl_m, int Bl_n,
      int Bl_nnz, MPI_Comm mpiComm_) :
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
   Amat = nullptr;
   Bmat = nullptr;
   Blmat = nullptr;
}

StochGenMatrix::~StochGenMatrix()
{
  for(size_t it = 0; it < children.size(); it++)
    delete children[it];

  if (Amat)
    delete Amat;

  if (Bmat)
    delete Bmat;

  if (Blmat)
    delete Blmat;
}


StochGenMatrix* StochGenMatrix::cloneFull(bool switchToDynamicStorage) const
{
   StochGenMatrix* clone = new StochGenMatrix(m, n, mpiComm);

   // clone submatrices
   clone->Amat = Amat->cloneFull(switchToDynamicStorage);
   clone->Bmat = Bmat->cloneFull(switchToDynamicStorage);
   clone->Blmat = Blmat->cloneFull(switchToDynamicStorage);

   for( size_t it = 0; it < children.size(); it++ )
      clone->children.push_back(children[it]->cloneFull(switchToDynamicStorage));

   return clone;	
}

/* creates an empty copy of the matrix with n = 0 for all submatrices and m (cols) as before */
StochGenMatrix* StochGenMatrix::cloneEmptyRows(bool switchToDynamicStorage) const
{
   StochGenMatrix* clone = new StochGenMatrix(m, n, mpiComm);

   // clone submatrices
   clone->Amat = Amat->cloneEmptyRows(switchToDynamicStorage);
   clone->Bmat = Bmat->cloneEmptyRows(switchToDynamicStorage);
	 clone->Blmat = Blmat->cloneEmptyRows(switchToDynamicStorage);

   for( size_t it = 0; it < children.size(); it++ )
      clone->children.push_back(children[it]->cloneEmptyRows(switchToDynamicStorage));

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
  return type == kStochGenDummyMatrix;
}
void StochGenMatrix::getSize( long long& m_out, long long& n_out ) const
{
  m_out = m; n_out=n;
}
void StochGenMatrix::getSize( int& m_out, int& n_out ) const
{
  m_out = m; n_out=n;
}


void StochGenMatrix::atPutDense( int row, int col, double * A, int lda,
				 int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::fromGetDense( int row, int col, double * A, int lda,
				   int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::columnScale2( const OoqpVector& vec, const OoqpVector& parentvec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);
   const SimpleVector& scalevecparent = dynamic_cast<const SimpleVector&>(parentvec);

   assert(scalevec.children.size() == 0 && children.size() == 0);

   Amat->columnScale(scalevecparent);
   Bmat->columnScale(*scalevec.vec);
   Blmat->columnScale(*scalevec.vec);
}

void StochGenMatrix::columnScale( const OoqpVector& vec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);

   assert(children.size() == scalevec.children.size());

   Bmat->columnScale(*scalevec.vec);
   Blmat->columnScale(*scalevec.vec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->columnScale2(*(scalevec.children[it]), *scalevec.vec);
}

void StochGenMatrix::rowScale2( const OoqpVector& vec, const OoqpVector* linkingvec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);

   assert(scalevec.children.size() == 0 && children.size() == 0);

   Amat->rowScale(*scalevec.vec);
   Bmat->rowScale(*scalevec.vec);

   if( linkingvec )
   {
      const SimpleVector* vecl = dynamic_cast<const SimpleVector*>(linkingvec);
      Blmat->rowScale(*vecl);
   }
}

void StochGenMatrix::rowScale( const OoqpVector& vec )
{
   const StochVector& scalevec = dynamic_cast<const StochVector&>(vec);

   assert(children.size() == scalevec.children.size());

   Bmat->rowScale(*scalevec.vec);

   const SimpleVector* vecl = dynamic_cast<const SimpleVector*>(scalevec.vecl);

   if( vecl )
      Blmat->rowScale(*vecl);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->rowScale2(*(scalevec.children[it]), vecl);
}

void StochGenMatrix::symmetricScale( const OoqpVector &vec)
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::scalarMult( double num)
{
  Amat->scalarMult(num);
  Bmat->scalarMult(num);
  Blmat->scalarMult(num);

  for (size_t it = 0; it < children.size(); it++) 
    children[it]->scalarMult(num);
}

void StochGenMatrix::fromGetSpRow( int row, int col,
				   double A[], int lenA, int jcolA[], int& nnz,
				   int colExtent, int& info )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::atPutSubmatrix( int destRow, int destCol, DoubleMatrix& M,
				     int srcRow, int srcCol,
				     int rowExtent, int colExtent )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::atPutSpRow( int col, 
				 double A[], int lenA, int jcolA[],
				 int& info )
{
  assert( "Not implemented" && 0 );
}


void StochGenMatrix::putSparseTriple( int irow[], int len, int jcol[], double A[], 
				      int& info )
{
  assert( "Not implemented" && 0 );
}

void StochGenMatrix::getDiagonal( OoqpVector& vec_ )
{
  StochVector& vec = dynamic_cast<StochVector&>(vec_);

  //check compatibility
  assert( children.size() == vec.children.size() );

  // only local B gives the diagonal
  Bmat->getDiagonal(*vec.vec);

  //do it recursively
  for(size_t it=0; it<children.size(); it++)
    children[it]->getDiagonal(*vec.children[it]);
}
 
void StochGenMatrix::setToDiagonal( const OoqpVector& vec_ )
{
  const StochVector& vec = dynamic_cast<const StochVector&>(vec_);
  assert(children.size() == vec.children.size());

  Bmat->setToDiagonal( *vec.vec);

  for(size_t it=0; it<children.size(); it++)
    children[it]->setToDiagonal(*vec.children[it]);
}

/* y = beta * y + alpha * this * x */
void StochGenMatrix::mult( double beta,  OoqpVector& y_,
			   double alpha, const OoqpVector& x_ ) const
{
  const StochVector & x = dynamic_cast<const StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  const SimpleVector& xvec = dynamic_cast<const SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

  if (0.0 == alpha) {
    y.vec->scale( beta );

    for(size_t it = 0; it < children.size(); it++)
      children[it]->mult(beta, *y.children[it], alpha, *x.children[it]);

    return;
  } else {
    //if( alpha != 1.0 || beta != 1.0 ) {
    //  y.vec->scale( beta/alpha );
    //}

    Bmat->mult(beta, yvec, alpha, xvec);

    long long mA, nA; 
    Amat->getSize(mA,nA);

    // not at root?
    if(nA>0) {
      Amat->mult(1.0, yvec, alpha, *x.parent->vec);

    }

    //if( 1.0 != alpha ) {
    //  y.vec->scale(alpha);
    //}
  }

  int blm, bln;
  Blmat->getSize(blm, bln);

  /* linking constraints present? */
  if( blm > 0 )
  {
	SimpleVector& yvecl = dynamic_cast<SimpleVector&>(*y.vecl);

	if( iAmSpecial(iAmDistrib, mpiComm) )
	  Blmat->mult(beta, yvecl, alpha, xvec);
	else
	  yvecl.setToZero();

    for(size_t it = 0; it < children.size(); it++)
	   children[it]->mult2(beta, *y.children[it], alpha, *x.children[it], yvecl);

	if(iAmDistrib) {
	  // sum up linking constraints vectors
	  int locn = yvecl.length();
	  double* buffer = new double[locn];

	  MPI_Allreduce(yvecl.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
	  yvecl.copyFromArray(buffer);

	  delete[] buffer;
    }
  }
  else
  {
    for(size_t it=0; it<children.size(); it++)
      children[it]->mult(beta, *y.children[it], alpha, *x.children[it]);
  }
}


/* mult method for children; needed only for linking constraints */
void StochGenMatrix::mult2( double beta,  OoqpVector& y_,
			   double alpha, OoqpVector& x_, OoqpVector& yparentl_ )
{
  StochVector & x = dynamic_cast<StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  if( 0.0 == alpha ) {
    y.vec->scale( beta );
    return;
  }

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

  Bmat->mult(beta, yvec, alpha, xvec);
  Blmat->mult(1.0, yparentl_, alpha, xvec);
  Amat->mult(1.0, *y.vec, alpha, *x.parent->vec);

  // not implemented
  assert(children.size() == 0);
}


void StochGenMatrix::transMult ( double beta,   OoqpVector& y_,
				 double alpha,  const OoqpVector& x_ ) const
{
  const StochVector & x = dynamic_cast<const StochVector&>(x_);
  StochVector & y = dynamic_cast<StochVector&>(y_);

  // assert tree compatibility
  assert(y.children.size() == children.size());
  assert(x.children.size() == children.size());

  const SimpleVector& xvec = dynamic_cast<const SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

  int blm, bln;
  Blmat->getSize(blm, bln);

  // with linking constraints?
  if( blm > 0 )
  {
    assert(x.vecl);
    const SimpleVector& xvecl = dynamic_cast<const SimpleVector&>(*x.vecl);

    if( iAmSpecial(iAmDistrib, mpiComm) )
    {
      //y_i = beta* y_i  +  alpha* B_i^T* x_i
      Bmat->transMult(beta, yvec, alpha, xvec);

	   //y_i = y_i  +  alpha* Bl_0^T* xl_i
	   Blmat->transMult(1.0, yvec, alpha, xvecl);
    }
    else
    {
      yvec.setToZero();
    }

    //!opt alloc buffer here and send it through the tree to be used by
    //!children when MPI_Allreduce
    //let the children compute their contribution
    for(size_t it = 0; it < children.size(); it++) {
      children[it]->transMult2(beta, *y.children[it], alpha, *x.children[it], yvec, xvecl);
    }
  }
  else // no linking constraints
  {
    if( iAmSpecial(iAmDistrib, mpiComm) )
      Bmat->transMult(beta, yvec, alpha, xvec);
    else
      yvec.setToZero();

    for(size_t it=0; it<children.size(); it++) {
      children[it]->transMult2(beta, *y.children[it], alpha, *x.children[it], yvec);
    }
  }

  if(iAmDistrib) {
    int locn=yvec.length();
    double* buffer = new double[locn];

    MPI_Allreduce(yvec.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
    yvec.copyFromArray(buffer);

    delete[] buffer;
  }
}


void StochGenMatrix::transMult2 ( double beta,   StochVector& y,
				  double alpha,  StochVector& x,
				  OoqpVector& yvecParent, const OoqpVector& xvecl) const
{
  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  const SimpleVector& xvec = dynamic_cast<const SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

#ifdef STOCH_TESTING
  int nA, mA;
  Amat->getSize(mA, nA);
  // this should NOT be the root
  assert(nA>0);
#endif

  //do A_i^T x_i and add it to yvecParent which already contains B_0^T x_0
  Amat->transMult(1.0, yvecParent, alpha, *x.vec);

  if( iAmSpecial(iAmDistrib, mpiComm) )
  {
    Bmat->transMult(beta, yvec, alpha, xvec);

	 //y_i = y_i  +  alpha* Bl_i^T* xl_i
	 Blmat->transMult(1.0, yvec, alpha, xvecl);
  }
  else
    yvec.setToZero();
  
  assert(children.size() == 0);
#if 0
  for(size_t it=0; it<children.size(); it++)
    children[it]->transMult2(beta, *y.children[it], 
			     alpha, *x.children[it],
			     yvec);
#endif


  if(iAmDistrib) {
    int locn=yvec.length();
    double* buffer = new double[locn];

    MPI_Allreduce(yvec.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
    yvec.copyFromArray(buffer);

    delete[] buffer;
  }
}

void StochGenMatrix::transMult2 ( double beta,   StochVector& y,
				  double alpha,  StochVector& x,
				  OoqpVector& yvecParent)
{
  //assert tree compatibility
  assert(y.children.size() - children.size() == 0);
  assert(x.children.size() - children.size() == 0);

  SimpleVector& xvec = dynamic_cast<SimpleVector&>(*x.vec);
  SimpleVector& yvec = dynamic_cast<SimpleVector&>(*y.vec);

#ifdef STOCH_TESTING
  int nA, mA;
  Amat->getSize(mA, nA);
  // this should NOT be the root
  assert(nA>0);
#endif

  //do A_i^T x_i and add it to yvecParent which already contains B_0^T x_0
  Amat->transMult(1.0, yvecParent, alpha, *x.vec);

  if( iAmSpecial(iAmDistrib, mpiComm) )
  {
    //do A_i^T x_i and add it to yvecParent which already contains
    //B_0^T x_0
    Bmat->transMult(beta, yvec, alpha, xvec);
  }
  else
    yvec.setToZero();

  for(size_t it=0; it<children.size(); it++)
    children[it]->transMult2(beta, *y.children[it],
			     alpha, *x.children[it],
			     yvec);

  if(iAmDistrib) {
    int locn=yvec.length();
    double* buffer = new double[locn];

    MPI_Allreduce(yvec.elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);
    yvec.copyFromArray(buffer);

    delete[] buffer;
  }
}


double StochGenMatrix::abmaxnorm() const
{
  double nrm = 0.0;
  
  for(size_t it = 0; it < children.size(); it++)
    nrm = std::max(nrm, children[it]->abmaxnorm());

  if(iAmDistrib)
  {
     double nrmG = 0;
     MPI_Allreduce(&nrm, &nrmG, 1, MPI_DOUBLE, MPI_MAX, mpiComm);
     nrm = nrmG;
  }

  nrm = std::max(nrm, max(Amat->abmaxnorm(), Bmat->abmaxnorm()));
  nrm = std::max(nrm, Blmat->abmaxnorm());

  nrm = std::max(nrm, max(Amat->abmaxnorm(), Bmat->abmaxnorm()));
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

   Amat->getLinkVarsNnz(vec);
}

void StochGenMatrix::writeToStream(ostream& out) const
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::writeToStreamDenseBordered( const StringGenMatrix& border_left, std::ostream& out ) const
{
   assert( border_left.children.size() == this->children.size() );

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   const int world_size = PIPS_MPIgetSize(mpiComm);

   int offset = 0;
   stringstream sout;
   MPI_Status status;
   int l;

   if( iAmDistrib )
      MPI_Barrier(mpiComm);

   if( iAmDistrib && my_rank > 0 )  // receive offset from previous process
      MPI_Recv(&offset, 1, MPI_INT, (my_rank - 1), 0, mpiComm, MPI_STATUS_IGNORE);
   else  //  !iAmDistrib || (iAmDistrib && rank == 0)
   {
      assert( border_left.mat );
      assert( this->Bmat );

      int mBmat, nBmat; this->Bmat->getSize(mBmat, nBmat);
      int mBd, nBd; border_left.mat->getSize(mBd, nBd);

      assert( mBmat == mBd );
      for( int i = 0; i < mBmat; ++i )
      {
         border_left.mat->writeToStreamDenseRow(sout, i);
         sout << "|\t";
         this->Bmat->writeToStreamDenseRow(sout, i);
         sout << endl;
      }
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      const StringGenMatrix& border_child = *border_left.children[it];

      children[it]->writeToStreamDenseChildBordered(sout, offset, *border_child.mat);
      int mBmat, nBmat; children[it]->Bmat->getSize(mBmat, nBmat);
      offset += nBmat;
   }

   if( iAmDistrib && my_rank > 0 )
   {
      std::string str = sout.str();
      // send string to rank ZERO to print it there:
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, my_rank, mpiComm);
      // send offset to next process:
      if( my_rank < world_size - 1 )
         MPI_Ssend(&offset, 1, MPI_INT, my_rank + 1, 0, mpiComm);
   }
   else if( !iAmDistrib )
      out << sout.str();
   else if( iAmDistrib && my_rank == 0 )
   {
      out << sout.str();
      MPI_Ssend(&offset, 1, MPI_INT, my_rank + 1, 0, mpiComm);

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

   // linking contraints ?
   int mlink, nlink;
   this->Blmat->getSize(mlink, nlink);
   if( mlink > 0 )
   {
      if( iAmDistrib )
         MPI_Barrier(mpiComm);

      assert( border_left.mat_link );
      int mBdl, nBdl; border_left.mat_link->getSize(mBdl, nBdl);
      assert( mBdl == mlink );

      // for each row r do:
      for( int r = 0; r < mlink; r++ )
      {
         if( iAmDistrib )
         {
            MPI_Barrier(mpiComm);

            // process Zero collects all the information and then prints it.
            if( my_rank == 0 )
            {
               out << border_left.mat_link->writeToStreamDenseRow(r);
               out << "|\t";
               out << this->Blmat->writeToStreamDenseRow(r);

               out << writeToStreamDenseRowLink(r);

               for( int p = 1; p < world_size; p++ )
               {
                  MPI_Probe(p, r + 1, mpiComm, &status);
                  MPI_Get_count(&status, MPI_CHAR, &l);
                  char *buf = new char[l];
                  MPI_Recv(buf, l, MPI_CHAR, p, r + 1, mpiComm, &status);
                  std::string rowPartFromP(buf, l);
                  out << rowPartFromP;
                  delete[] buf;
               }
               out << std::endl;

            }
            else // rank != 0
            {
               std::string str = writeToStreamDenseRowLink(r);
               MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, r + 1, mpiComm);
            }
         }
         else // not distributed
         {
            std::stringstream sout;
            border_left.mat_link->writeToStreamDenseRow(sout, r);
            sout << "|\t";
            this->Blmat->writeToStreamDenseRow(sout, r);

            for( size_t it = 0; it < children.size(); it++ )
               children[it]->Blmat->writeToStreamDenseRow(sout, r);

            out << sout.rdbuf() << std::endl;
         }
      }
   }

   assert( this->Bmat );
   int mBd, nBd; border_left.getSize(mBd, nBd);
   int mB0mat, nB0mat; this->Bmat->getSize(mB0mat, nB0mat);

   if( iAmDistrib )
   {
      // send offset to next process:
      if( my_rank == world_size - 1 )
         MPI_Ssend(&offset, 1, MPI_INT, 0, 0, mpiComm);
      else if( my_rank == 0 )
      {
         MPI_Recv(&offset, 1, MPI_INT, world_size - 1, 0, mpiComm, MPI_STATUS_IGNORE);

         for( int i = 0; i < offset + nBd + nB0mat + 1; ++i )
            out << "_\t";
         out << std::endl;
      }
   }
   else
   {
      for( int i = 0; i < offset + nBd + nB0mat + 1; ++i )
         out << "_\t";
      out << std::endl;
   }

   if( iAmDistrib )
      MPI_Barrier(mpiComm);
}

void StochGenMatrix::writeToStreamDense(ostream& out) const
{
   const int my_rank = PIPS_MPIgetRank(mpiComm);
   const int world_size = PIPS_MPIgetSize(mpiComm);

   int m, n;
   int offset = 0;
   stringstream sout;
   MPI_Status status;
   int l;

   if( iAmDistrib )
      MPI_Barrier(mpiComm);

   if( iAmDistrib && my_rank > 0 )  // receive offset from previous process
      MPI_Recv(&offset, 1, MPI_INT, (my_rank - 1), 0, mpiComm, MPI_STATUS_IGNORE);

   else  //  !iAmDistrib || (iAmDistrib && rank == 0)
      this->Bmat->writeToStreamDense(out);

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->writeToStreamDenseChild(sout, offset);
      children[it]->Bmat->getSize(m, n);
      offset += n;
   }

   if( iAmDistrib && my_rank > 0 )
   {
      std::string str = sout.str();
      // send string to rank ZERO to print it there:
      MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, my_rank, mpiComm);
      // send offset to next process:
      if( my_rank < world_size - 1 )
         MPI_Ssend(&offset, 1, MPI_INT, my_rank + 1, 0, mpiComm);
   }

   else if( !iAmDistrib )
      out << sout.str();

   else if( iAmDistrib && my_rank == 0 )
   {
      out << sout.str();
      MPI_Ssend(&offset, 1, MPI_INT, my_rank + 1, 0, mpiComm);

      for( int p = 1; p < world_size; p++ )
      {
         MPI_Probe(p, p, mpiComm, &status);
         MPI_Get_count(&status, MPI_CHAR, &l);
         char *buf = new char[l];
         MPI_Recv(buf, l, MPI_CHAR, p, p, mpiComm, &status);
         string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }

   }

   // linking contraints ?
   int mlink, nlink;
   this->Blmat->getSize(mlink, nlink);
   if( mlink > 0 )
   {
      if( iAmDistrib )
         MPI_Barrier(mpiComm);

      // for each row r do:
      for( int r = 0; r < mlink; r++ )
      {
         if( iAmDistrib )
         {
            MPI_Barrier(mpiComm);

            // process Zero collects all the information and then prints it.
            if( my_rank == 0 )
            {
               out << this->Blmat->writeToStreamDenseRow(r);

               out << writeToStreamDenseRowLink(r);

               for( int p = 1; p < world_size; p++ )
               {
                  MPI_Probe(p, r + 1, mpiComm, &status);
                  MPI_Get_count(&status, MPI_CHAR, &l);
                  char *buf = new char[l];
                  MPI_Recv(buf, l, MPI_CHAR, p, r + 1, mpiComm, &status);
                  string rowPartFromP(buf, l);
                  out << rowPartFromP;
                  delete[] buf;
               }
               out << endl;

            }
            else // rank != 0
            {
               std::string str = writeToStreamDenseRowLink(r);
               MPI_Ssend(str.c_str(), str.length(), MPI_CHAR, 0, r + 1, mpiComm);
            }
         }
         else // not distributed
         {
            stringstream sout;
            this->Blmat->writeToStreamDenseRow(sout, r);
            for( size_t it = 0; it < children.size(); it++ )
               children[it]->Blmat->writeToStreamDenseRow(sout, r);

            out << sout.rdbuf() << endl;
         }
      }
   }
   if( iAmDistrib )
      MPI_Barrier(mpiComm);
}

/** writes child matrix blocks, offset indicates the offset between A and B block. */
void StochGenMatrix::writeToStreamDenseChild(stringstream& out, int offset) const
{
   int mA, mB, n;
   this->Amat->getSize(mA, n);
   this->Bmat->getSize(mB, n);
   assert(mA == mB );
   for(int r=0; r < mA; r++)
   {
      this->Amat->writeToStreamDenseRow(out, r);
      for(int i=0; i<offset; i++)
         out <<'\t';
      this->Bmat->writeToStreamDenseRow(out, r);
      out << endl;
   }
}

void StochGenMatrix::writeToStreamDenseChildBordered(stringstream& out, int offset, const SparseGenMatrix& border) const
{
   int mA, mB, mBd, n;
   this->Amat->getSize(mA, n);
   this->Bmat->getSize(mB, n);
   border.getSize(mBd, n);

   assert( mBd == mA );
   assert( mA == mB );

   for(int r = 0; r < mA; r++)
   {
      border.writeToStreamDenseRow(out, r);
      out << "|\t";

      this->Amat->writeToStreamDenseRow(out, r);
      for(int i = 0; i < offset; i++)
         out <<'\t';

      this->Bmat->writeToStreamDenseRow(out, r);
      out << std::endl;
   }
}

/** returns a string containing the linking-row rowidx of the children. */
std::string StochGenMatrix::writeToStreamDenseRowLink(int rowidx) const
{
   std::string str_all;
   for( size_t it = 0; it < children.size(); it++ )
   {
      std::string str = children[it]->Blmat->writeToStreamDenseRow(rowidx);
      str_all.append(str);
   }
   return str_all;
}

void StochGenMatrix::writeMPSformatRows(ostream& out, int rowType, OoqpVector* irhs) const
{
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   string rt;
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
            out<< " "<<rt<<" row_"<<rt<<"_"<<"R" <<"_"<<i <<endl;
      }
      // linking rows:
      if( Blmat )
      {
         this->Blmat->getSize(m, n);
         for(int i=0; i<m; i++)
         {
            if( !irhs || (irhsStoch && dynamic_cast<SimpleVector*>(irhsStoch->vecl)->elements()[i] != 0.0) )
               out<<" "<< rt<<" row_"<<rt<<"_"<<"L" <<"_"<<i <<endl;
         }
      }
   }
   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->Amat->getSize(m, n);
      for(int i=0; i<m; i++)
      {
         if( !irhs || (irhsStoch && dynamic_cast<SimpleVector*>(irhsStoch->children[it]->vec)->elements()[i] != 0.0) )
            out<<" "<< rt<<" row_"<<rt<<"_"<<it <<"_"<<i <<endl;
      }
   }
}

/* Make the elements in this matrix symmetric. The elements of interest
 *  must be in the lower triangle, and the upper triangle must be empty.
 *  @param info zero if the operation succeeded. Otherwise, insufficient
 *  space was allocated to symmetrize the matrix.
 */
void StochGenMatrix::symmetrize( int& info )
{
  assert( "Has not been yet implemented" && 0 );
}


void StochGenMatrix::randomize( double alpha, double beta, double * seed )
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::initTransposed(bool dynamic)
{
   Bmat->initTransposed(dynamic);
   Blmat->initTransposed(dynamic);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->initTransposedChild(dynamic);
}

void StochGenMatrix::deleteTransposed()
{
   Amat->deleteTransposed();
   Bmat->deleteTransposed();
   Blmat->deleteTransposed();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->deleteTransposed();
}

void StochGenMatrix::initTransposedChild(bool dynamic)
{
   Amat->initTransposed(dynamic);
   Bmat->initTransposed(dynamic);

   if( Blmat != nullptr )
      Blmat->initTransposed(dynamic);
}


void StochGenMatrix::atPutDiagonal( int idiag, OoqpVector& v )
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::fromGetDiagonal( int idiag, OoqpVector& v )
{
  assert( "Has not been yet implemented" && 0 );
}

int StochGenMatrix::numberOfNonZeros()
{
  int nnz = 0;

  for(size_t it=0; it<children.size(); it++)
    nnz += children[it]->numberOfNonZeros();

  if(iAmDistrib) {
    int nnzG = 0;
    MPI_Allreduce(&nnz, &nnzG, 1, MPI_INT, MPI_SUM, mpiComm);
    nnz=nnzG;
  }

  nnz += Amat->numberOfNonZeros() + Bmat->numberOfNonZeros() + Blmat->numberOfNonZeros();

  return nnz;
}

void StochGenMatrix::matTransDMultMat(OoqpVector& d, SymMatrix** res)
{
  assert( "Has not been yet implemented" && 0 );
}
void StochGenMatrix::matTransDinvMultMat(OoqpVector& d, SymMatrix** res)
{
  assert( "Has not been yet implemented" && 0 );
}

void StochGenMatrix::getNnzPerRow(OoqpVectorBase<int>& nnzVec, OoqpVectorBase<int>* linkParent)
{
   StochVectorBase<int>& nnzVecStoch = dynamic_cast<StochVectorBase<int>&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   SimpleVectorBase<int>* nnzvecl = nullptr;

   Bmat->addNnzPerRow(*(nnzVecStoch.vec));

   if( linkParent != nullptr )
      Amat->addNnzPerRow(*(nnzVecStoch.vec));

   /* with linking constraints? */
   if( nnzVecStoch.vecl || linkParent )
   {
      assert(nnzVecStoch.vecl == nullptr || linkParent == nullptr);

      if( linkParent )
         nnzvecl = dynamic_cast<SimpleVectorBase<int>*>(linkParent);
      else
         nnzvecl = dynamic_cast<SimpleVectorBase<int>*>(nnzVecStoch.vecl);

      if( linkParent != nullptr || iAmSpecial(iAmDistrib, mpiComm) )
         Blmat->addNnzPerRow(*nnzvecl);
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
   StochVectorBase<int>& nnzVecStoch = dynamic_cast<StochVectorBase<int>&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   SimpleVectorBase<int>* vec = dynamic_cast<SimpleVectorBase<int>*>(nnzVecStoch.vec);

   if( iAmSpecial(iAmDistrib, mpiComm) || linkParent != nullptr )
   {
      Bmat->addNnzPerCol(*(vec));

      int blm, bln;
      Blmat->getSize(blm, bln);

      /* with linking constraints? */
      if( blm > 0 )
         Blmat->addNnzPerCol(*vec);
   }

   // not at root?
   if( linkParent != nullptr )
      Amat->addNnzPerCol(*linkParent);
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

void StochGenMatrix::getRowMinMaxVec(bool getMin, bool initializeVec,
      const OoqpVector* colScaleVec, const OoqpVector* colScaleParent, OoqpVector& minmaxVec, OoqpVector* linkParent)
{
   StochVector& minmaxVecStoch = dynamic_cast<StochVector&>(minmaxVec);
   const StochVector* const colScaleVecStoch = dynamic_cast<const StochVector*>(colScaleVec);

   SimpleVector* mvecl = nullptr;
   const SimpleVector* const covecparent = dynamic_cast<const SimpleVector*>(colScaleParent);
   const SimpleVector* const covec = colScaleVecStoch != nullptr ?
         dynamic_cast<SimpleVector*>(colScaleVecStoch->vec) : nullptr;

   // assert tree compatibility
   assert(minmaxVecStoch.children.size() == children.size());

   Bmat->getRowMinMaxVec(getMin, initializeVec, covec, *(minmaxVecStoch.vec));

   // not at root?
   if( linkParent != nullptr )
      Amat->getRowMinMaxVec(getMin, false, covecparent, *(minmaxVecStoch.vec));

   /* with linking constraints? */
   if( minmaxVecStoch.vecl || linkParent )
   {
      assert(minmaxVecStoch.vecl == nullptr || linkParent == nullptr);

      // at root?
      if( linkParent == nullptr )
      {
         mvecl = dynamic_cast<SimpleVector*>(minmaxVecStoch.vecl);

         if( initializeVec )
         {
            if( getMin )
               mvecl->setToConstant(std::numeric_limits<double>::max());
            else
               mvecl->setToZero();
         }
      }
      else
      {
         mvecl = dynamic_cast<SimpleVector*>(linkParent);
      }

      if( linkParent != nullptr || iAmSpecial(iAmDistrib, mpiComm) )
         Blmat->getRowMinMaxVec(getMin, false, covec, *mvecl);
   }

   if( colScaleVec != nullptr )
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->getRowMinMaxVec(getMin, initializeVec, colScaleVecStoch->children[it], covec,
               *(minmaxVecStoch.children[it]), mvecl);
   }
   else
   {
      for( size_t it = 0; it < children.size(); it++ )
         children[it]->getRowMinMaxVec(getMin, initializeVec, nullptr, nullptr,
               *(minmaxVecStoch.children[it]), mvecl);
   }

   // distributed, with linking constraints, and at root?
   if( iAmDistrib && minmaxVecStoch.vecl != nullptr && linkParent == nullptr )
   {
      assert(mvecl != nullptr);

      // sum up linking constraints vectors
      if( getMin )
        PIPS_MPIminArrayInPlace(mvecl->elements(), mvecl->length(), mpiComm);
      else
        PIPS_MPImaxArrayInPlace(mvecl->elements(), mvecl->length(), mpiComm);
   }
}

void StochGenMatrix::getColMinMaxVec(bool getMin, bool initializeVec,
        const OoqpVector* rowScaleVec, const OoqpVector* rowScaleLink, OoqpVector& minmaxVec, OoqpVector* minmaxParent)
{
   StochVector& minmaxVecStoch = dynamic_cast<StochVector&>(minmaxVec);
   const StochVector* rowScaleVecStoch = dynamic_cast<const StochVector*>(rowScaleVec);

   // assert tree compatibility
   assert(minmaxVecStoch.children.size() == children.size());

   SimpleVector* const mvec = dynamic_cast<SimpleVector*>(minmaxVecStoch.vec);
   const SimpleVector* const covec = rowScaleVecStoch != nullptr ?
         dynamic_cast<SimpleVector*>(rowScaleVecStoch->vec) : nullptr;

   Bmat->getColMinMaxVec(getMin, initializeVec, covec, *(mvec));

   int blm, bln;
   Blmat->getSize(blm, bln);

   const SimpleVector* covecl = nullptr;

   /* with linking constraints? */
   if( blm > 0 )
   {
      // with rowScale vector?
      if( rowScaleVecStoch != nullptr )
      {
         // at root?
         if( minmaxParent == nullptr )
            covecl = dynamic_cast<SimpleVector*>(rowScaleVecStoch->vecl);
         else
            covecl = dynamic_cast<const SimpleVector*>(rowScaleLink);

         assert(covecl != nullptr);
      }

      if( iAmSpecial(iAmDistrib, mpiComm) || minmaxParent != nullptr )
         Blmat->getColMinMaxVec(getMin, false, covecl, *mvec);
   }

   // not at root?
   if( minmaxParent != nullptr )
      Amat->getColMinMaxVec(getMin, false, covec, *(minmaxParent));
   else
   {
      if( rowScaleVecStoch )
      {
         for( size_t it = 0; it < children.size(); it++ )
            children[it]->getColMinMaxVec(getMin, initializeVec, rowScaleVecStoch->children[it], covecl, *(minmaxVecStoch.children[it]), mvec);
      }
      else
      {
         for( size_t it = 0; it < children.size(); it++ )
            children[it]->getColMinMaxVec(getMin, initializeVec, nullptr, nullptr, *(minmaxVecStoch.children[it]), mvec);
      }
   }

   // distributed and at root?
   if( iAmDistrib && minmaxParent == nullptr )
   {
      const int locn = mvec->length();
      double* const entries = mvec->elements();
      double* buffer = new double[locn];

      if( getMin )
         MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_MIN, mpiComm);
      else
         MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_MAX, mpiComm);

      mvec->copyFromArray(buffer);

      delete[] buffer;
   }
}


void StochGenMatrix::addRowSums( OoqpVector& sumVec, OoqpVector* linkParent ) const
{
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
   const StochVectorBase<int>& rowNnzVecStoch = dynamic_cast<const StochVectorBase<int>&>(rowNnzVec);
   const StochVectorBase<int>& colNnzVecStoch = dynamic_cast<const StochVectorBase<int>&>(colNnzVec);

   assert(rowNnzVecStoch.children.size() == colNnzVecStoch.children.size());

   const SimpleVectorBase<int>* const rowvec = dynamic_cast<const SimpleVectorBase<int>*>(rowNnzVecStoch.vec);
   const SimpleVectorBase<int>* const colvec = dynamic_cast<const SimpleVectorBase<int>*>(colNnzVecStoch.vec);

   const SimpleVectorBase<int>* const rowlink = dynamic_cast<const SimpleVectorBase<int>*>(rowNnzVecStoch.vecl);
   assert(rowvec); assert(colvec);

   Amat->initStaticStorageFromDynamic(*rowvec, colParentVec); // initialized with colVec == nullptr for parent
   Bmat->initStaticStorageFromDynamic(*rowvec, colvec);

   // at root?
   if( colParentVec == nullptr )
   {
      assert(rowLinkVec == nullptr);

      if( rowlink != nullptr )
         Blmat->initStaticStorageFromDynamic(*rowlink, colvec);

      for( size_t it = 0; it < children.size(); it++ )
         children[it]->initStaticStorageFromDynamic(*(rowNnzVecStoch.children[it]), *(colNnzVecStoch.children[it]), rowlink, colvec);
   }
   else
   {
      assert( children.size() == 0 );
      if( rowLinkVec != nullptr)
         Blmat->initStaticStorageFromDynamic(*rowLinkVec, colvec);
   }

}

void StochGenMatrix::freeDynamicStorage()
{
   Amat->freeDynamicStorage();
   Bmat->freeDynamicStorage();
   Blmat->freeDynamicStorage();

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->freeDynamicStorage();
}

void StochGenMatrix::recomputeSize( StochGenMatrix* parent )
{
   m = 0;
   n = 0;

   if( parent )
      Amat->getSize(m,n);

   assert( m >= 0 );
   assert( n >= 0 );

   int b_mat_m = 0;
   int b_mat_n = 0;
   Bmat->getSize(b_mat_m, b_mat_n);
   assert( b_mat_m >= 0 );
   assert( b_mat_n >= 0 );

   if( !parent )
      m += b_mat_m;
   n += b_mat_n;

   if( parent != nullptr )
   {
      parent->m += b_mat_m;
      parent->n += b_mat_n;
   }

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->recomputeSize( this );

   int bl_mat_m = 0; int bl_mat_n = 0;
   Blmat->getSize(bl_mat_m, bl_mat_n);

   m += bl_mat_m;
}

void StochGenMatrix::updateKLinkConsCount(std::vector<int>& linkCount) const
{
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
         children[it]->Blmat->updateNonEmptyRowsCount(linkCount);
      }

   if( iAmDistrib )
      MPI_Allreduce(MPI_IN_PLACE, &linkCount[0], m, MPI_INT, MPI_SUM, mpiComm);
}

void StochGenMatrix::updateKLinkVarsCount(std::vector<int>& link_block_count) const
{
   int m, n;
   Bmat->getSize(m, n);

   if( n == 0 )
      return;

   assert(link_block_count.size() == size_t(n));

   for( size_t it = 0; it < children.size(); it++ )
      if( !(children[it]->isKindOf(kStochGenDummyMatrix)) )
      {
         children[it]->Amat->getTranspose().updateNonEmptyRowsCount(link_block_count);
         children[it]->Amat->deleteTransposed();
      }

   if( iAmDistrib )
      PIPS_MPIsumArrayInPlace( link_block_count, mpiComm );
}

void StochGenMatrix::get2LinkStartBlocksAndCountsNew(std::vector<int>& block_start, std::vector<int>& block_count) const
{
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
         children[it]->Blmat->updateNonEmptyRowsCountNew(it, block_count, block_start, block_end);
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
         children[it]->Blmat->updateNonEmptyRowsCount(it, linkBlockCount, linkBlockStart, linkBlockEnd);
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
   if( Blmat )
      Blmat->permuteCols(permvec);

   Bmat->permuteCols(permvec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->permuteLinkingVarsChild(permvec);
}

void StochGenMatrix::permuteLinkingVarsChild(const std::vector<unsigned int>& permvec)
{
   Amat->permuteCols(permvec);

   assert(children.size() == 0);
}

void StochGenMatrix::permuteLinkingCons(const std::vector<unsigned int>& permvec)
{
   if( Blmat )
      Blmat->permuteRows(permvec);

   for( size_t it = 0; it < children.size(); it++ )
      children[it]->permuteLinkingCons(permvec);
}


void StochGenMatrix::updateTransposed()
{
  Amat->updateTransposed();
  Bmat->updateTransposed();
  Blmat->updateTransposed();

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
   bool in_sync = true;

   assert( Amat );
   assert( Bmat );
   assert( Blmat );

   /* no need to check if not distributed or not in root node */
   if( !iAmDistrib || children.size() == 0)
      return in_sync;


   int my_rank, world_size;
   MPI_Comm_rank(mpiComm, &my_rank);
   MPI_Comm_size(mpiComm, &world_size);

   /* if matrix has static storage */
   if( Amat->getStorageHandle().notNil() || Bmat->getStorageHandle().notNil() ||
	    Blmat->getStorageHandle().notNil() )
   {
		  assert( Amat->getStorageHandle().notNil() );
      assert( Bmat->getStorageHandle().notNil() );
      assert( Blmat->getStorageHandle().notNil() );

      /* since we are in root node Amat should look as follows */
      assert( Amat->getStorageRef().len == 0 );
      assert( Amat->getStorageRef().n == -1 );

      /* static storage */
      const int lenght_entries_bmat = Bmat->getStorageRef().len;
      const int length_columns_bmat = Bmat->getStorageRef().len;
      const int lenght_rowoffest_bmat = Bmat->getStorageRef().m + 1;

      const int lenght_entries_blmat = Blmat->getStorageRef().len;
      const int length_columns_blmat = Blmat->getStorageRef().len;
      const int lenght_rowoffest_blmat = Blmat->getStorageRef().m + 1;

      const long long count_row_cols = length_columns_bmat + lenght_rowoffest_bmat + length_columns_blmat + lenght_rowoffest_blmat;
      const long long count_entries = lenght_entries_bmat + lenght_entries_blmat;

      assert( count_row_cols < std::numeric_limits<int>::max());
      assert( count_entries < std::numeric_limits<int>::max());

      std::vector<double> sendbuf_entries(count_entries, 0.0);
      std::vector<double> recvbuf_entries(count_entries, 0.0);

      std::vector<int> sendbuf_row_col(count_row_cols, 0);
      std::vector<int> recvbuf_row_col(count_row_cols, 0);

      /* fill Bmat into send buffers */
      const double * M = Bmat->getStorageRef().M;
      const int * krowM = Bmat->getStorageRef().krowM;
      const int * jColM = Bmat->getStorageRef().jcolM;

      std::copy(M, M + lenght_entries_bmat, sendbuf_entries.begin());

      std::copy(krowM, krowM + lenght_rowoffest_bmat, sendbuf_row_col.begin());
      std::copy(jColM, jColM + lenght_entries_bmat, sendbuf_row_col.begin() + lenght_rowoffest_bmat);

      /* fill Blmat into send buffers */
      const double * Ml = Blmat->getStorageRef().M;
      const int * krowMl = Blmat->getStorageRef().krowM;
      const int * jColMl = Blmat->getStorageRef().jcolM;

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
               std::cout << "matrix entries out of sync" << std::endl;
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
               std::cout << "matrix indices (col or row) out of sync" << std::endl;
            in_sync = false;
         }
      }
   }

   /* if stoch mat has dynamic storage also check that */
   if(Bmat->hasDynamicStorage() || Blmat->hasDynamicStorage())
   {
      assert(Bmat->hasDynamicStorage());
      assert(Blmat->hasDynamicStorage());

      SparseStorageDynamic& Bmat_dyn = Bmat->getStorageDynamicRef();
      SparseStorageDynamic& Blmat_dyn = Blmat->getStorageDynamicRef();

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
               std::cout << "matrix entries in dynamic storage out of sync" << std::endl;
            in_sync = false;
         }
      }
      for( int i = 0; i < count_row_cols_dyn; ++i )
      {
         if( !PIPSisEQ(sendbuf_row_col_dynamic[i], recvbuf_row_coldynamic[i]) )
         {
            /* someone else had a higher value here */
            if(my_rank == 0)
               std::cout << "matrix indices in dynamic storage out of sync" << std::endl;
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
  // todo: check that matrix is in correct format
  assert( matrix_row.children.size() == children.size() );
  assert( children.size() != 0 );
  assert( -1 <= child && child <= (int) children.size() );

  int index_row;

  // append row to all matrices necessary
  // todo maybe this can be done nicer - maybe we can just recursively call some method also on the dummies 
  if(linking)
  {
    index_row = Blmat->appendRow( *matrix_row.Blmat, row );

    for(unsigned int i = 0; i < children.size(); ++i)
    { 
      if( !children[i]->isKindOf(kStochGenDummyMatrix) )
      {
        assert( !matrix_row.children[i]->isKindOf(kStochGenDummyMatrix) );
        children[i]->Blmat->appendRow( *matrix_row.children[i]->Blmat, row);
      }
    }
  }
  else
  {
    if(child != -1)
    {
      index_row = children[child]->Amat->appendRow( *matrix_row.children[child]->Amat, row );
#ifndef NDEBUG
      const int index_row1 = children[child]->Bmat->appendRow( *matrix_row.children[child]->Bmat, row );
#else
      children[child]->Bmat->appendRow( *matrix_row.children[child]->Bmat, row );
#endif
      assert(index_row1 == index_row);
    }
    else
    {
      index_row = Bmat->appendRow( *matrix_row.Bmat, row );
    }
  }

  return index_row;
};

/* y += alpha RowAt(child, row, linking) */
void StochGenMatrix::axpyWithRowAt( double alpha, StochVector* y, SimpleVector* y_linking, int child, int row, bool linking) const
{
   assert( y );
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(y->children.size() == children.size());

   /* go through all available children and calculate y += alpha * rowAt(row) */
   if( linking )
   {
      assert( Blmat );
      if( y_linking )
         Blmat->axpyWithRowAt(alpha, *y_linking, row);
      else
      {
         assert( y->vec );
         Blmat->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->vec), row);
      }

      for( unsigned int i = 0; i < children.size(); ++i )
      {
         if( !children[i]->isKindOf(kStochGenDummyMatrix) )
         {
            assert(children[i]->Blmat);
            assert(y->children[i]->vec);
            children[i]->Blmat->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->children[i]->vec), row);
         }
      }
   }
   else
   {
      if( child == -1 )
      {
         assert(Bmat);
         if( y_linking )
            Bmat->axpyWithRowAt(alpha, *y_linking, row);
         else
         {
            assert(y->vec);
            Bmat->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->vec), row);
         }
      }
      else
      {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);

         assert(y->children[child]->vec);
         children[child]->Bmat->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->children[child]->vec), row);

         if( y_linking )
            children[child]->Amat->axpyWithRowAt(alpha, *y_linking, row);
         else
         {
            assert(y->vec);
            children[child]->Amat->axpyWithRowAt(alpha, dynamic_cast<SimpleVector&>(*y->vec), row);
         }
      }
   }
}

void StochGenMatrix::axpyWithRowAtPosNeg( double alpha, StochVector* y_pos, SimpleVector* y_link_pos,
      StochVector* y_neg, SimpleVector* y_link_neg, int child, int row, bool linking ) const
{
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
         Blmat->axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
      else
      {
         assert( y_pos->vec );
         assert( y_neg->vec );
         Blmat->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->vec), dynamic_cast<SimpleVector&>(*y_neg->vec), row);
      }

      for( unsigned int i = 0; i < children.size(); ++i )
      {
         if( !children[i]->isKindOf(kStochGenDummyMatrix) )
         {
            assert(children[i]->Blmat);
            assert(y_pos->children[i]->vec);
            assert(y_neg->children[i]->vec);
            children[i]->Blmat->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->children[i]->vec), dynamic_cast<SimpleVector&>(*y_neg->children[i]->vec), row);
         }
      }
   }
   else
   {
      if( child == -1 )
      {
         assert(Bmat);
         if( y_link_pos )
            Bmat->axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
         else
         {
            assert(y_pos->vec);
            assert(y_neg->vec);
            Bmat->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->vec), dynamic_cast<SimpleVector&>(*y_neg->vec), row);
         }
      }
      else
      {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);

         assert(y_pos->children[child]->vec);
         assert(y_neg->children[child]->vec);
         children[child]->Bmat->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->children[child]->vec), dynamic_cast<SimpleVector&>(*y_neg->children[child]->vec), row);

         if( y_link_pos )
            children[child]->Amat->axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
         else
         {
            assert(y_pos->vec);
            assert(y_neg->vec);
            children[child]->Amat->axpyWithRowAtPosNeg(alpha, dynamic_cast<SimpleVector&>(*y_pos->vec), dynamic_cast<SimpleVector&>(*y_neg->vec), row);
         }
      }
   }
}

double StochGenMatrix::localRowTimesVec(const StochVector &vec, int child, int row, bool linking) const
{
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(vec.children.size() == children.size());

   double res = 0.0;

   /* go through all available children and multiply the vec times row in submatrix */
   if( linking )
   {
      assert(Blmat);
      assert(vec.vec);
      res += Blmat->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);

      for( unsigned int i = 0; i < children.size(); ++i )
      {
         if( !children[i]->isKindOf(kStochGenDummyMatrix) )
         {
            assert(children[i]->Blmat);
            assert(vec.children[i]->vec);
            res += children[i]->Blmat->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.children[i]->vec), row);
         }
      }
   }
   else
   {
      if( child == -1 )
      {
         assert(Bmat);
         assert(vec.vec);
         res += Bmat->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);
      }
      else
      {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);
         assert(vec.vec);
         assert(vec.children[child]->vec);
         res += children[child]->Amat->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.vec), row);
         res += children[child]->Bmat->localRowTimesVec(dynamic_cast<const SimpleVector&>(*vec.children[child]->vec), row);
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

   SparseGenMatrix* const A_left = Bmat->shaveLeft(n_vars);

   SparseGenMatrix* const Bl_left_top = Blmat->shaveLeft(n_vars);
   SparseGenMatrix* const bottom_left_block = Bl_left_top->shaveBottom(m_conss);

   SparseGenMatrix* const Bl_right_bottom = Blmat->shaveBottom(m_conss);

   StringGenMatrix* const border_bottom = new StringGenMatrix(false, Bl_right_bottom, nullptr, mpiComm);
   StringGenMatrix* const border_left = new StringGenMatrix(true, A_left, Bl_left_top, mpiComm);

   for( size_t it = 0; it < children.size(); it++ )
   {
      StringGenMatrix* border_left_child = nullptr;
      StringGenMatrix* border_bottom_child = nullptr;

      children[it]->shaveBorder(m_conss, n_vars, border_left_child, border_bottom_child);

      border_left->addChild(border_left_child);
      border_bottom->addChild(border_bottom_child);
   }

   m -= m_conss;
   n -= n_vars;

   BorderedGenMatrix* const bordered_matrix = new BorderedGenMatrix(this, border_left, border_bottom, bottom_left_block, mpiComm);
   StochGenMatrix* me = this;
   IotrAddRef( &me );

   assert(m >= 0 && n >= 0);

   return bordered_matrix;
}

void StochGenMatrix::shaveBorder( int m_conss, int n_vars, StringGenMatrix*& border_left, StringGenMatrix*& border_bottom )
{
   SparseGenMatrix* const border_a_mat = Amat->shaveLeft(n_vars);
   SparseGenMatrix* const border_bl_mat = Blmat->shaveBottom(m_conss);

   border_bottom = new StringGenMatrix(false, border_bl_mat, nullptr, mpiComm);
   border_left = new StringGenMatrix(true, border_a_mat, nullptr, mpiComm);

   if( children.size() == 0 )
      assert( PIPS_MPIgetSize( mpiComm ) == 1 );

   for( size_t it = 0; it < children.size(); it++ )
   {
      assert(" should not end up here! : todo implement?");

      StringGenMatrix* border_left_child;
      StringGenMatrix* border_bottom_child;

      children[it]->shaveBorder(m_conss, n_vars, border_left_child, border_bottom_child);

      border_left->addChild(border_left_child);
      border_bottom->addChild(border_bottom_child);
   }
}
