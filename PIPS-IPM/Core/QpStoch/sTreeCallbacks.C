/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "sTreeCallbacks.h"
#include "sData.h"

#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "SimpleVector.h"
#include "DoubleMatrixTypes.h"

#include <cmath>
#include <algorithm>    // std::swap
#include <numeric>

#include "pipsport.h"

#ifndef UCTRANS // see note in smlParDriver.C
#define UCTRANS
#endif

sTreeCallbacks::~sTreeCallbacks()
{
}

sTreeCallbacks::sTreeCallbacks() 
   : N_INACTIVE(-1), MY_INACTIVE(-1), MZ_INACTIVE(-1), MYL_INACTIVE(-1), MZL_INACTIVE(-1),
    nx_active(0),  my_active(0),  mz_active(0),  myl_active(0),  mzl_active(0),
    nx_inactive(-1),  my_inactive(-1),  mz_inactive(-1),  myl_inactive(-1),  mzl_inactive(-1),
    isDataPresolved(false), hasPresolvedData(false), data(nullptr)

{
   if( -1 == rankMe ) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
   if( -1 == numProcs ) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

sTreeCallbacks::sTreeCallbacks(StochInputTree* inputTree)
  : sTree(),
    N_INACTIVE(-1), MY_INACTIVE(-1), MZ_INACTIVE(-1), MYL_INACTIVE(-1), MZL_INACTIVE(-1),
    nx_active(0),  my_active(0),  mz_active(0),  myl_active(0),  mzl_active(0),
    nx_inactive(-1),  my_inactive(-1),  mz_inactive(-1),  myl_inactive(-1),  mzl_inactive(-1),
    isDataPresolved(false), hasPresolvedData(false)
{
   if( -1 == rankMe )
      MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
   if( -1 == numProcs )
      MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  data = inputTree->nodeInput;

  for(size_t it = 0; it < inputTree->children.size(); it++)
     children.push_back(new sTreeCallbacks(inputTree->children[it]));
}

sTreeCallbacks::sTreeCallbacks(InputNode* data_)
  : sTree(), 
    N_INACTIVE(-1), MY_INACTIVE(-1), MZ_INACTIVE(-1), MYL_INACTIVE(-1), MZL_INACTIVE(-1),
    nx_active(data_->n),  my_active(data_->my),  mz_active(data_->mz),  myl_active(data_->myl),  mzl_active(data_->mzl),
    nx_inactive(-1),  my_inactive(-1),  mz_inactive(-1),  myl_inactive(-1),  mzl_inactive(-1),
    isDataPresolved(false), hasPresolvedData(false), data(data_)
{
   assert( 0 && "Not used currently" );
   if( -1 == rankMe ) MPI_Comm_rank(MPI_COMM_WORLD, &rankMe);
   if( -1 == numProcs ) MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
}

void sTreeCallbacks::loadLocalSizes()
{
  //already in data structure
}

void sTreeCallbacks::switchToPresolvedData()
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );
   assert(!isDataPresolved);
   assert(hasPresolvedData);

   std::swap(N_INACTIVE, N);
   std::swap(MY_INACTIVE, MY);
   std::swap(MZ_INACTIVE, MZ);
   std::swap(MYL_INACTIVE, MYL);
   std::swap(MZL_INACTIVE, MZL);

   std::swap(nx_active, nx_inactive);
   std::swap(my_active, my_inactive);
   std::swap(mz_active, mz_inactive);
   std::swap(myl_active, myl_inactive);
   std::swap(mzl_active, mzl_inactive);

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->np = this->nx_active;
      dynamic_cast<sTreeCallbacks*>(children[it])->switchToPresolvedData();
   }

   isDataPresolved = true;
}


void sTreeCallbacks::switchToOriginalData()
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );
   assert(isDataPresolved);

   std::swap(N_INACTIVE, N);
   std::swap(MY_INACTIVE, MY);
   std::swap(MZ_INACTIVE, MZ);
   std::swap(MYL_INACTIVE, MYL);
   std::swap(MZL_INACTIVE, MZL);

   std::swap(nx_active, nx_inactive);
   std::swap(my_active, my_inactive);
   std::swap(mz_active, mz_inactive);
   std::swap(myl_active, myl_inactive);
   std::swap(mzl_active, mzl_inactive);

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->np = this->nx_active;
      dynamic_cast<sTreeCallbacks*>(children[it])->switchToOriginalData();
   }

   isDataPresolved = false;
}


bool sTreeCallbacks::isPresolved()
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );
   return isDataPresolved;
}

bool sTreeCallbacks::hasPresolved()
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );
   return hasPresolvedData;
}

void sTreeCallbacks::initPresolvedData( const sData& presolved_data )
{
   const StochSymMatrix& Q = dynamic_cast<const StochSymMatrix&>(*presolved_data.Q);
   const StochGenMatrix& A = dynamic_cast<const StochGenMatrix&>(*presolved_data.A);
   const StochGenMatrix& C = dynamic_cast<const StochGenMatrix&>(*presolved_data.C);

   const StochVector& g = dynamic_cast<const StochVector&>(*presolved_data.g);
   const StochVector& b = dynamic_cast<const StochVector&>(*presolved_data.bA);

   assert( presolved_data.icupp.notNil() || presolved_data.icupp.notNil() );
   const StochVector& ic = presolved_data.iclow.notNil() ? dynamic_cast<const StochVector&>(*presolved_data.iclow) :
         dynamic_cast<const StochVector&>(*presolved_data.icupp);

   initPresolvedData(Q, A, C, g, b, ic, -1, -1);
}

void sTreeCallbacks::initPresolvedData(const StochSymMatrix& Q, const StochGenMatrix& A, const StochGenMatrix& C,
      const StochVector& nxVec, const StochVector& myVec, const StochVector& mzVec, int mylParent, int mzlParent)
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );
   assert(!hasPresolvedData);

   assert(nxVec.children.size() == children.size());
   assert(myVec.children.size() == children.size());
   assert(mzVec.children.size() == children.size());
   assert(A.children.size() == children.size());
   assert(C.children.size() == children.size());

   const SimpleVector& nxVecSimple = dynamic_cast<const SimpleVector&>(*nxVec.vec);
   const SimpleVector& myVecSimple = dynamic_cast<const SimpleVector&>(*myVec.vec);
   const SimpleVector& mzVecSimple = dynamic_cast<const SimpleVector&>(*mzVec.vec);
   const SimpleVector* const myVecSimpleLink = dynamic_cast<const SimpleVector*>(myVec.vecl);
   const SimpleVector* const mzVecSimpleLink = dynamic_cast<const SimpleVector*>(mzVec.vecl);

   N_INACTIVE = nxVecSimple.n;
   MY_INACTIVE = myVecSimple.n;
   MZ_INACTIVE = mzVecSimple.n;

   nx_inactive = N_INACTIVE;
   my_inactive = MY_INACTIVE;
   mz_inactive = MZ_INACTIVE;

   if( myVecSimpleLink != nullptr )
   {
      assert(np == -1);
      MYL_INACTIVE = myVecSimpleLink->n;
      myl_inactive = MYL_INACTIVE;
   }
   else
   {
      myl_inactive = mylParent;
   }

   if( mzVecSimpleLink != nullptr )
   {
      assert(np == -1);
      MZL_INACTIVE = mzVecSimpleLink->n;
      mzl_inactive = MZL_INACTIVE;
   }
   else
   {
      mzl_inactive = mzlParent;
   }

   // empty child?
   if( N_INACTIVE == 0 && MY_INACTIVE == 0 && MZ_INACTIVE == 0 && np != -1 )
   {
      myl_inactive = 0;
      mzl_inactive = 0;
   }

   // are we at the root?
   if( np == -1 )
   {
      assert(mylParent == -1);
      assert(mzlParent == -1);
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      assert(children[it]->np == this->nx_active);
      sTreeCallbacks* sTreeCallbacksChild = dynamic_cast<sTreeCallbacks*>(children[it]);

      sTreeCallbacksChild->initPresolvedData(*Q.children[it], *A.children[it], *C.children[it], *nxVec.children[it],
            *myVec.children[it], *mzVec.children[it], myl_inactive, mzl_inactive);

      N_INACTIVE += sTreeCallbacksChild->N_INACTIVE;
      MY_INACTIVE += sTreeCallbacksChild->MY_INACTIVE;
      MZ_INACTIVE += sTreeCallbacksChild->MZ_INACTIVE;
      MYL_INACTIVE += sTreeCallbacksChild->MYL_INACTIVE;
      MZL_INACTIVE += sTreeCallbacksChild->MZL_INACTIVE;
   }

   hasPresolvedData = true;
}


void sTreeCallbacks::writeSizes( std::ostream& sout ) const
{
   const int myRank = PIPS_MPIgetRank(commWrkrs);

   MPI_Barrier(commWrkrs);
   if( myRank == 0 )
   {
      sout << "N          : "  <<  N           << "\n";
      sout << "MY         : "  <<  MY          << "\n";
      sout << "MZ         : "  <<  MZ          << "\n";
      sout << "MYL        : "  <<  MYL         << "\n";
      sout << "MZL        : "  <<  MZL         << "\n";
      sout << "nx_active  : "  <<  nx_active    << "\n";
      sout << "my_active  : "  <<  my_active    << "\n";
      sout << "mz_active  : "  <<  mz_active    << "\n";
      sout << "myl_active : "  <<  myl_active   << "\n";
      sout << "mzl_active : "  <<  mzl_active   << "\n";
   }


   for( size_t it = 0; it < children.size(); it++ )
   {
      if( children[it]->commWrkrs != MPI_COMM_NULL )
      {
         sout << "child " << it << ": \n\n";
         dynamic_cast<sTreeCallbacks*>(children[it])->writeSizes(sout);
      }
      MPI_Barrier(commWrkrs);
   }

   if( myRank == 0 )
      std::cout << std::endl;
}


// this is usually called after assigning processes
void sTreeCallbacks::computeGlobalSizes()
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );
   assert( !isDataPresolved );

   const int my_rank = PIPS_MPIgetRank(MPI_COMM_WORLD);

   if( data && sTree::isInVector( my_rank, myProcs ) )
   {
      // callback must be used for sizes
      assert( data->nCall );
      assert( data->myCall );
      assert(data->mzCall );

      /* load sizes from callbacks */
      data->nCall(data->user_data, data->id, &data->n);
      data->myCall(data->user_data, data->id, &data->my);
      data->mzCall(data->user_data, data->id, &data->mz);

      if( data->mylCall )
         data->mylCall(data->user_data, data->id, &data->myl);
      else
         data->myl = 0;

      if( data->mzlCall )
         data->mzlCall(data->user_data, data->id, &data->mzl);
      else
         data->myl = 0;

      N = data->n;
      MY = data->my;
      MZ = data->mz;
      MYL = data->myl;
      MZL = data->mzl;

      nx_active = data->n;
      my_active = data->my;
      mz_active = data->mz;
      myl_active = data->myl;
      mzl_active = data->mzl;
   }
   else
   {
      nx_active = my_active = mz_active = myl_active = mzl_active = 0;
      N = MY = MZ = MYL = MZL = 0;
   }

   for( size_t it = 0; it < children.size(); it++ )
   {
      children[it]->np = this->data->n;
      children[it]->computeGlobalSizes();
      N += children[it]->N;
      MY += children[it]->MY;
      MZ += children[it]->MZ;
   }
}

StochSymMatrix* sTreeCallbacks::createQ() const
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );

   //is this node a dead-end for this process?
   if( commWrkrs == MPI_COMM_NULL )
      return new StochSymDummyMatrix();
  
   if( data->nnzQ < 0 )
      data->fnnzQ(data->user_data, data->id, &data->nnzQ);

   StochSymMatrix *Q = new StochSymMatrix(N, data->n, data->nnzQ,
         commWrkrs);

   data->fQ(data->user_data, data->id, Q->diag->krowM(), Q->diag->jcolM(),
         Q->diag->M());

   for( size_t it = 0; it < children.size(); it++ )
   {
      StochSymMatrix *child = children[it]->createQ();
      Q->AddChild(child);
   }
   return Q;
}

StochGenMatrix* sTreeCallbacks::createMatrix( DATA_INT m_ABmat, DATA_INT n_Mat,
      DATA_INT nnzAmat, DATA_NNZ fnnzAmat, DATA_MAT Amat, DATA_INT nnzBmat,
      DATA_NNZ fnnzBmat, DATA_MAT Bmat, DATA_INT m_Blmat, DATA_INT nnzBlmat,
      DATA_NNZ fnnzBlmat, DATA_MAT Blmat ) const
{
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );

   if( commWrkrs == MPI_COMM_NULL )
      return new StochGenDummyMatrix();

   const bool root = (np == -1);
   const bool has_linking = (data->*fnnzBlmat != nullptr);

   if( has_linking && data->*nnzBlmat < 0 )
      (data->*fnnzBlmat)(data->user_data, data->id, &(data->*nnzBlmat));

   if ( data->*nnzAmat < 0 )
      (data->*fnnzAmat)(data->user_data, data->id, &(data->*nnzAmat));

   StochGenMatrix* A = nullptr;

   if( root )
   {
      data->*nnzBmat = 0;

      if( data->*fnnzBlmat > 0 )
      {
         // populate B with A's data B_0 is the A_0 from the theoretical form; also fill Bl
         // (i.e. the first block of linking constraints)
         A = new StochGenMatrix(MY + MYL, N,
               data->*m_ABmat, np, data->*nnzBmat,
               data->*m_ABmat, data->*n_Mat, data->*nnzAmat,
               data->*m_Blmat, data->*n_Mat, data->*nnzBlmat,
               commWrkrs);
      }
      else
      {
         // populate B with A's data B_0 is the A_0 from the theoretical form
         A = new StochGenMatrix(MY + MYL, N,
               data->*m_ABmat, np, data->*nnzBmat,
               data->*m_ABmat, data->*n_Mat,  data->*nnzAmat,
               commWrkrs);
      }

      //populate submatrix B
      (data->*Amat)(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());

      printf("root  -- my=%d  myl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n",
            data->*m_ABmat, data->*m_Blmat, data->*n_Mat, np, data->*nnzAmat, data->*nnzBmat, data->*nnzBlmat);
   }
   else
   {
      if( data->*nnzBmat < 0 )
         (data->*fnnzBmat)(data->user_data, data->id, &(data->*nnzBmat));

      if( data->fnnzBl > 0 )
      {
         A = new StochGenMatrix(MY + MYL, N,
               data->*m_ABmat, np, data->*nnzAmat,
               data->*m_ABmat, data->*n_Mat, data->*nnzBmat,
               data->*m_Blmat, data->*n_Mat, data->*nnzBlmat,
               commWrkrs);
      }
      else
      {
         A = new StochGenMatrix(MY + MYL, N,
               data->*m_ABmat, np, data->*nnzAmat,
               data->*m_ABmat, data->*n_Mat,  data->*nnzBmat,
               commWrkrs);
      }

      //populate the submatrices A, B
      (data->*Amat)(data->user_data, data->id, A->Amat->krowM(), A->Amat->jcolM(), A->Amat->M());
      (data->*Bmat)(data->user_data, data->id, A->Bmat->krowM(), A->Bmat->jcolM(), A->Bmat->M());

      printf("  -- my=%d  myl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n",
            data->*m_ABmat, data->*m_Blmat, data->*n_Mat, np, data->*nnzAmat, data->*nnzBmat, data->*nnzBlmat);
   }

   // populate Bl if existent
   if( data->*Blmat )
      (data->*Blmat)(data->user_data, data->id, A->Blmat->krowM(), A->Blmat->jcolM(), A->Blmat->M());

   for(size_t it = 0; it < children.size(); it++)
   {
      StochGenMatrix* child = dynamic_cast<sTreeCallbacks*>(children[it])->createMatrix( m_ABmat, n_Mat, nnzAmat, fnnzAmat, Amat,
            nnzBmat, fnnzBmat, Bmat, m_Blmat, nnzBlmat, fnnzBlmat, Blmat );
      A->AddChild(child);
   }
   return A;
}

StochGenMatrix* sTreeCallbacks::createA() const
{
   DATA_INT m_ABmat = &InputNode::my;
   DATA_INT n_Mat = &InputNode::n;
   DATA_INT nnzAmat = &InputNode::nnzA;
   DATA_NNZ fnnzAmat = &InputNode::fnnzA;
   DATA_MAT Amat = &InputNode::fA;

   DATA_INT nnzBmat = &InputNode::nnzB;
   DATA_NNZ fnnzBmat = &InputNode::fnnzB;
   DATA_MAT Bmat = &InputNode::fB;

   DATA_INT m_Blmat = &InputNode::myl;
   DATA_INT nnzBlmat = &InputNode::nnzBl;
   DATA_NNZ fnnzBlmat = &InputNode::fnnzBl;
   DATA_MAT Blmat = &InputNode::fBl;

   return createMatrix( m_ABmat, n_Mat, nnzAmat, fnnzAmat, Amat, nnzBmat,
         fnnzBmat, Bmat, m_Blmat, nnzBlmat, fnnzBlmat, Blmat );
}

StochGenMatrix* sTreeCallbacks::createC() const
{
   DATA_INT m_CDmat = &InputNode::mz;
   DATA_INT n_Mat = &InputNode::n;
   DATA_INT nnzCmat = &InputNode::nnzC;
   DATA_NNZ fnnzCmat = &InputNode::fnnzC;
   DATA_MAT Cmat = &InputNode::fC;

   DATA_INT nnzDmat = &InputNode::nnzD;
   DATA_NNZ fnnzDmat = &InputNode::fnnzD;
   DATA_MAT Dmat = &InputNode::fD;

   DATA_INT m_Dlmat = &InputNode::mzl;
   DATA_INT nnzDlmat = &InputNode::nnzDl;
   DATA_NNZ fnnzDlmat = &InputNode::fnnzDl;
   DATA_MAT Dlmat = &InputNode::fDl;

   return createMatrix( m_CDmat, n_Mat, nnzCmat, fnnzCmat, Cmat, nnzDmat,
         fnnzDmat, Dmat, m_Dlmat, nnzDlmat, fnnzDlmat, Dlmat );
}

int sTreeCallbacks::nx() const
{
   return nx_active;
}

int sTreeCallbacks::my() const
{
   return my_active;
}

int sTreeCallbacks::myl() const
{
   return myl_active;
}

int sTreeCallbacks::mz() const
{
   return mz_active;
}

int sTreeCallbacks::mzl() const
{
   return mzl_active;
}

int sTreeCallbacks::id() const
{
   return data->id;
}

StochVector* sTreeCallbacks::createVector( DATA_INT n_vec, DATA_VEC vec, DATA_INT n_linking_vec, DATA_VEC linking_vec ) const
{
   assert( n_vec );
   assert( vec );

   assert( !(is_hierarchical_root || is_hierarchical_inner || is_hierarchical_leaf)
         || (false && "cannot be used with hierarchical data") );

   if( commWrkrs == MPI_COMM_NULL )
      return new StochDummyVector();

   const int nlinking = (np == - 1 && linking_vec != nullptr) ? data->*n_linking_vec : -1;

   StochVector* svec = new StochVector( data->*n_vec, nlinking, commWrkrs );

   assert( svec->vec );
   double* elems = dynamic_cast<SimpleVector*>(svec->vec)->elements();
   double* elems_link = (nlinking != -1) ? dynamic_cast<SimpleVector*>(svec->vecl)->elements() : nullptr;

   (data->*vec)(data->user_data, data->id, elems, data->*n_vec);

   // at root and with linking constraints?
   if( nlinking != -1 )
   {
      assert( n_linking_vec );
      assert( linking_vec );
      (data->*linking_vec)(data->user_data, data->id, elems_link, data->*n_linking_vec);
   }

   for( const sTree* child_tree : children )
   {
      StochVector* child_vec = dynamic_cast<const sTreeCallbacks*>(child_tree)->createVector( n_vec, vec, n_linking_vec, linking_vec );
      svec->AddChild( child_vec );
    }

   return svec;
}

StochVector* sTreeCallbacks::createc() const
{
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fc;

   return createVector( n_func, func, nullptr, nullptr );
}

StochVector* sTreeCallbacks::createxlow() const
{
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fxlow;

   return createVector( n_func, func, nullptr, nullptr );
}

StochVector* sTreeCallbacks::createixlow() const
{
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fixlow;

   return createVector( n_func, func, nullptr, nullptr );
}

StochVector* sTreeCallbacks::createxupp() const
{
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fxupp;

   return createVector( n_func, func, nullptr, nullptr );
}

StochVector* sTreeCallbacks::createixupp() const
{
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fixupp;

   return createVector( n_func, func, nullptr, nullptr );
}

StochVector* sTreeCallbacks::createb() const
{
   DATA_INT n_func = &InputNode::my;
   DATA_VEC func = &InputNode::fb;

   DATA_INT n_func_link = &InputNode::myl;
   DATA_VEC func_link = &InputNode::fbl;

   return createVector( n_func, func, n_func_link, func_link );
}

StochVector* sTreeCallbacks::createclow() const
{
   DATA_INT n_func = &InputNode::mz;
   DATA_VEC func = &InputNode::fclow;

   DATA_INT n_func_link = &InputNode::mzl;
   DATA_VEC func_link = &InputNode::fdllow;

   return createVector( n_func, func, n_func_link, func_link );
}

StochVector* sTreeCallbacks::createiclow() const
{
   DATA_VEC func = &InputNode::ficlow;
   DATA_INT n_func = &InputNode::mz;

   DATA_VEC func_link = &InputNode::fidllow;
   DATA_INT n_func_link = &InputNode::mzl;

   return createVector( n_func, func, n_func_link, func_link );
}

StochVector* sTreeCallbacks::createcupp() const
{
   DATA_INT n_func = &InputNode::mz;
   DATA_VEC func = &InputNode::fcupp;

   DATA_INT n_func_link = &InputNode::mzl;
   DATA_VEC func_link = &InputNode::fdlupp;

   return createVector( n_func, func, n_func_link, func_link );
}

StochVector* sTreeCallbacks::createicupp() const
{
   DATA_INT n_func = &InputNode::mz;
   DATA_VEC func = &InputNode::ficupp;

   DATA_INT n_func_link = &InputNode::mzl;
   DATA_VEC func_link = &InputNode::fidlupp;

   return createVector( n_func, func, n_func_link, func_link );
   assert(!is_hierarchical_root || ( false && "cannot be used with hierarchical data" ) );
}

sTree* sTreeCallbacks::shaveDenseBorder( int nx_to_shave, int myl_to_shave, int mzl_to_shave )
{
   sTreeCallbacks* top_layer = new sTreeCallbacks();

   /* sTree members */
   top_layer->commWrkrs = commWrkrs;
   top_layer->myProcs = myProcs;

   /* this must happen at MPI_COMM_WORLD level */
   assert( rankMe == PIPS_MPIgetRank() );
   assert( numProcs == PIPS_MPIgetSize() );

   top_layer->N = N;
   this->N -= nx_to_shave;

   top_layer->MY = MY;
   top_layer->MYL = MYL;
   this->MYL -= myl_to_shave;

   top_layer->MZ = MZ;
   top_layer->MZL = MZL;
   this->MZL -= mzl_to_shave;
   top_layer->np = -1;

   assert( IPMIterExecTIME == -1 );
   top_layer->children.push_back(this);
   top_layer->numProcs = numProcs;

   // TODO: not sure about the ressources monitors..: resMon, iterMon..
   top_layer->is_hierarchical_root = true;

   /* sTreeCallbacks members */
   top_layer->nx_active = nx_to_shave;
   this->nx_active -= nx_to_shave;

   top_layer->my_active = -1;
   top_layer->mz_active = -1;

   top_layer->myl_active = myl_to_shave;
   this->myl_active -= myl_to_shave;

   top_layer->mzl_active = mzl_to_shave;
   this->mzl_active -= mzl_to_shave;

   return top_layer;
}

void sTreeCallbacks::createSubcommunicatorsAndChildren( std::vector<unsigned int>& map_my_procs_to_sub_comm, std::vector<unsigned int>& map_child_to_sub_comm )
{
   // TODO : maybe adjust the split at some point - the first n_leftovers procs have one additional block assigned to them .. this will lead to some inbalance in the tree here..
   const int size_comm = PIPS_MPIgetSize(commWrkrs);
   assert( size_comm == static_cast<int>(myProcs.size()) );

   if( size_comm < 3 )
      assert( false && "need at least 4 MPI_COMMS to split a root linsys");

   /* we will split the current communicator into sqrt(size) (rounded) new ones */
   const int n_split_comms = std::floor(std::sqrt(size_comm));
   const int leftover_comms = size_comm % n_split_comms;

   if( PIPS_MPIgetRank(commWrkrs) == 0 )
      std::cout << "Splitting communicator of size " << size_comm << " into " << n_split_comms << " new ones" << std::endl;

   int my_comm = -1;
   int current_comm_until = 0;

   map_my_procs_to_sub_comm.clear();
   map_my_procs_to_sub_comm.resize(size_comm);

   const int my_rank = PIPS_MPIgetRank(commWrkrs);
   for( int i = 0; i < n_split_comms; ++i )
   {
      const int current_comm_start = current_comm_until;

      if( current_comm_until < leftover_comms * (n_split_comms + 1) )
         current_comm_until += (n_split_comms + 1);
      else
         current_comm_until += n_split_comms;

      assert( current_comm_until <= size_comm );

      for( int j = current_comm_start; j < current_comm_until; ++j )
      {
         map_my_procs_to_sub_comm[j] = i;

         if( my_rank == j )
            my_comm = i;
      }
   }
   assert( 0 <= my_comm && my_comm < n_split_comms );

   MPI_Comm comm_my_subsystem;
   MPI_Comm_split(commWrkrs, my_comm, my_rank, &comm_my_subsystem);


   std::vector<sTreeCallbacks*> sub_roots(n_split_comms);

   /* create new sub-roots */
   for( int i = 0; i < n_split_comms; ++i )
   {
      sub_roots[i] = new sTreeCallbacks();

      if( i == my_comm )
         sub_roots[i]->commWrkrs = comm_my_subsystem;
      else
         sub_roots[i]->commWrkrs = MPI_COMM_NULL;
   }

   map_child_to_sub_comm.resize(children.size());

   /* assign children to sub_comm systems and assign processes to them  */
   for( size_t i = 0; i < children.size(); ++i )
   {
      assert( children[i]->myProcs.size() == 1 );
      const int child_proc = children[i]->myProcs[0];
      const unsigned int assigned_sub_comm = map_my_procs_to_sub_comm[ child_proc ];
      sub_roots[assigned_sub_comm]->children.push_back( children[i] );
      map_child_to_sub_comm[i] = assigned_sub_comm;

      // TODO : inefficient..
      if( !isInVector(child_proc, sub_roots[assigned_sub_comm]->myProcs) )
         sub_roots[assigned_sub_comm]->myProcs.push_back( child_proc );
      assert( isInVector( child_proc, sub_roots[assigned_sub_comm]->myProcs ) );
   }

   /* add sub_roots as this new children */
   children.clear();
   children.insert( this->children.begin(), sub_roots.begin(), sub_roots.end() );
   assert( children[my_comm]->myProcs.size() == static_cast<const size_t>(PIPS_MPIgetSize(comm_my_subsystem)) );
}


void sTreeCallbacks::splitTreeSquareRoot( const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC )
{
   assert( commWrkrs != MPI_COMM_NULL );
   assert( !is_hierarchical_root );

   assert( twoLinksStartBlockA.size() == twoLinksStartBlockC.size() );
   assert( twoLinksStartBlockA.size() == children.size() );

   assert( rankMe == PIPS_MPIgetRank(commWrkrs) );

   // assert commWrkrs is MPI_COMM_WORLD
   assert( numProcs == PIPS_MPIgetSize(commWrkrs) );

   const size_t n_leafs = children.size();

   createSubcommunicatorsAndChildren( map_proc_subcomm, map_block_subcomm );
   /* children should now be created and an additional layer in the tree should exist */
   assert( map_proc_subcomm.size() == this->myProcs.size() );
   assert( map_block_subcomm.size() == n_leafs );

   this->is_hierarchical_inner = true;

   /* count how many 2 links will stay at this */
   int two_links_root_eq = 0;
   int two_links_root_ineq = 0;
   int two_links_children_eq_sum = 0;
   int two_links_children_ineq_sum = 0;
   std::vector<int> two_links_children_eq(children.size(), 0);
   std::vector<int> two_links_children_ineq(children.size(), 0);

   /* count two links for all new sub-communicators */
   size_t block = 0;
   for( size_t i = 0; i < children.size(); ++i )
   {
      sTreeCallbacks& child = dynamic_cast<sTreeCallbacks&>( *children[i] );

      const std::vector<int>& child_procs = child.myProcs;
      assert( is_sorted(child_procs.begin(), child_procs.end()) );

      for( size_t j = 0; j < child.children.size(); ++j )
      {
         sTreeCallbacks& childchild = dynamic_cast<sTreeCallbacks&>(*child.children[i]);

         /* do not allow splitting of inner nodes yet */
         assert( childchild.children.size() == 0 );

         if( childchild.children.size() == 0 )
         {
            /* assert - either my child or not assigned to me */
            assert( childchild.commWrkrs == MPI_COMM_SELF || childchild.commWrkrs == MPI_COMM_NULL );

            if( childchild.commWrkrs == MPI_COMM_NULL )
            {
               assert( childchild.MYL == 0 );
               assert( childchild.MZL == 0 );
            }
            childchild.is_hierarchical_leaf = true;
         }

         assert( block < twoLinksStartBlockA.size() );
         /* the two links from our last child will stay linking in the root node */
         if( j == child.children.size() - 1 )
         {
            two_links_root_eq += twoLinksStartBlockA[block];
            two_links_root_ineq += twoLinksStartBlockC[block];
         }
         else
         {
            two_links_children_eq_sum += twoLinksStartBlockA[block];
            two_links_children_ineq_sum += twoLinksStartBlockC[block];
            two_links_children_eq[i] += twoLinksStartBlockA[block];
            two_links_children_ineq[i] += twoLinksStartBlockC[block];
         }

         ++block;
      }
   }

   assert( std::accumulate(two_links_children_eq.begin(), two_links_children_eq.end(),
         decltype(two_links_children_eq)::value_type(0)) == two_links_children_eq_sum );
   assert( std::accumulate(two_links_children_ineq.begin(), two_links_children_ineq.end(),
         decltype(two_links_children_ineq)::value_type(0)) == two_links_children_ineq_sum );

   if( rankMe == 0 )
   {
      std::cout << "Splitting " << two_links_children_eq_sum + two_links_root_eq <<
            " equality two-links into " << two_links_root_eq << " root and " << two_links_children_eq_sum << " child links" << std::endl;
      std::cout << "Splitting " << two_links_children_ineq_sum + two_links_root_ineq <<
            " inequality two-links into " << two_links_root_ineq << " root and " << two_links_children_ineq_sum << " child links" << std::endl;
   }

   assert( block == n_leafs );

   /* now set the linking sizes of all children and this */
   this->myl_active -= two_links_children_eq_sum;
   this->mzl_active -= two_links_children_ineq_sum;

   /* recompute sizes for the children */
   for( size_t i = 0; i < children.size(); ++i )
   {
      sTreeCallbacks& child = dynamic_cast<sTreeCallbacks&>(*children[i]);

      child.MYL = two_links_children_eq[i];
      child.myl_active = two_links_children_eq[i];

      child.MZL = two_links_children_ineq[i];
      child.mzl_active = two_links_children_ineq[i];

      child.nx_active = 0;
      child.my_active = 0;
      child.mz_active = 0;
      child.np = this->nx_active;

      /* if we are part of the child we compute the sizes - else they get set to 0 */
      if( isInVector(rankMe, child.myProcs) )
      {
         int nx = 0;
         int my = 0;
         int mz = 0;

         for( size_t j = 0; j < child.children.size(); ++j )
         {

            sTreeCallbacks& childchild = dynamic_cast<sTreeCallbacks&>(*child.children[j]);
            assert( childchild.myProcs.size() == 1 );
            assert( isInVector( childchild.myProcs[0], child.myProcs ) );

            if( childchild.commWrkrs == MPI_COMM_SELF )
            {
               assert( isInVector(rankMe, childchild.myProcs ) );
               nx += childchild.N;
               my += childchild.MY;
               mz += childchild.MZ;
               childchild.MYL = child.MYL;
               childchild.myl_active = child.myl_active;
               childchild.MZL = child.MZL;
               childchild.mzl_active = child.mzl_active;
            }
            else
            {
               assert( childchild.MYL == 0 );
               assert( childchild.myl_active == 0 );
               assert( childchild.MZL == 0 );
               assert( childchild.mzl_active == 0 );
            }
         }
         child.N = nx;
         child.MY = my;
         child.MZ = mz;

      }
      else
      {
         child.N = child.MY = child.MZ = child.MYL = child.MZL = 0;
         child.nx_active = child.my_active = child.mz_active = child.myl_active = child.mzl_active = 0;
      }

   }
#ifndef NEDBUG
   int nx = 0;
   int my = 0;
   int mz = 0;

   for( size_t i = 0; i < children.size(); ++i )
   {
      sTreeCallbacks& child = dynamic_cast<sTreeCallbacks&>(*children[i]);
      nx += child.N;
      my += child.MY;
      mz += child.MZ;
   }

   assert( this->N == nx + this->nx_active );
   assert( this->MY == my + this->my_active );
   assert( this->MZ == mz + this->mz_active );
#endif
}


sTree* sTreeCallbacks::switchToHierarchicalTree( int nx_to_shave, int myl_to_shave, int mzl_to_shave,
      const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC )
{
   assert( !is_hierarchical_root );
   assert( np == -1 );

   assert( nx_to_shave >= 0 );
   assert( myl_to_shave >=0 );
   assert( mzl_to_shave >= 0 );

   assert( nx_to_shave <= nx_active );
   assert( myl_to_shave <= myl_active );
   assert( mzl_to_shave <= mzl_active );

   /* distributed preconditioner must be deactivated */
   assert( rankZeroW == 0 );
   assert( rankPrcnd == -1 );
   assert( commP2ZeroW == MPI_COMM_NULL );

//   this->splitTreeSquareRoot( twoLinksStartBlockA, twoLinksStartBlockC );

   sTreeCallbacks* top_layer = dynamic_cast<sTreeCallbacks*>( shaveDenseBorder( nx_to_shave, myl_to_shave, mzl_to_shave ) );

   return top_layer;
}

void sTreeCallbacks::splitMatrixAccordingToTree( StochSymMatrix& mat ) const
{

   assert( false && "TODO : implement " );

}

void sTreeCallbacks::splitMatrixAccordingToTree( StochGenMatrix& mat ) const
{
   assert( false && "TODO : implement " );

}

void sTreeCallbacks::splitVectorAccordingToTree( StochVector& vec ) const
{
   assert( map_block_subcomm.size() == vec.children.size() );

   std::vector<StochVector*> sub_roots(this->children.size());

   for( size_t i = 0; i < this->children.size(); ++i )
   {
//      sub_roots = new StochVector( , , this->children[i]->commWrkrs);
   }

   assert( false && "TODO : implement " );

}

void sTreeCallbacks::splitDataAccordingToTree( sData& data ) const
{
   splitVectorAccordingToTree( dynamic_cast<StochVector&>(*data.g) );
   splitVectorAccordingToTree( dynamic_cast<StochVector&>(*data.bux) );
   splitVectorAccordingToTree( dynamic_cast<StochVector&>(*data.ixupp) );
   splitVectorAccordingToTree( dynamic_cast<StochVector&>(*data.blx) );
   splitVectorAccordingToTree( dynamic_cast<StochVector&>(*data.ixlow) );

   //   /* we ordered global linking vars first and global linking rows to the end */

   splitMatrixAccordingToTree(dynamic_cast<StochSymMatrix&>(*data.Q));
   assert( false && "TODO : implement " );
//
//   StochVector* bA_hier = dynamic_cast<StochVector&>(*bA).raiseBorder(n_global_eq_linking_conss, true, false);
//
//   StochVector* bu_hier = dynamic_cast<StochVector&>(*bu).raiseBorder(n_global_ineq_linking_conss, true, false);
//   StochVector* icupp_hier = dynamic_cast<StochVector&>(*icupp).raiseBorder(n_global_ineq_linking_conss, true, false);
//   StochVector* bl_hier = dynamic_cast<StochVector&>(*bl).raiseBorder(n_global_ineq_linking_conss, true, false);
//   StochVector* iclow_hier = dynamic_cast<StochVector&>(*iclow).raiseBorder(n_global_ineq_linking_conss, true, false);
//

}
