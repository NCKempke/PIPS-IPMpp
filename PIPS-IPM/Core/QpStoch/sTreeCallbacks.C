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

sTreeCallbacks::sTreeCallbacks()
   : print_tree_sizes_on_reading{ pips_options::getBoolParameter( "PRINT_TREESIZES_ON_READ" ) }
{
   if( -1 == rankMe )
      rankMe = PIPS_MPIgetRank();
   if( -1 == numProcs )
      numProcs = PIPS_MPIgetSize();
}

sTreeCallbacks::sTreeCallbacks(StochInputTree* inputTree)
  : sTree(), print_tree_sizes_on_reading{ pips_options::getBoolParameter( "PRINT_TREESIZES_ON_READ" ) },
    data{ inputTree->nodeInput }
{
   if( -1 == rankMe )
      rankMe = PIPS_MPIgetRank();
   if( -1 == numProcs )
      numProcs = PIPS_MPIgetSize();

  for(size_t it = 0; it < inputTree->children.size(); it++)
     children.push_back(new sTreeCallbacks(inputTree->children[it]));
}

sTreeCallbacks::sTreeCallbacks(InputNode* data_)
  : sTree(), nx_active(data_->n), my_active(data_->my), mz_active(data_->mz),
    myl_active(data_->myl), mzl_active(data_->mzl),
    print_tree_sizes_on_reading{ pips_options::getBoolParameter( "PRINT_TREESIZES_ON_READ" ) },
    data(data_)
{
   assert( false && "Not used currently" );
   if( -1 == rankMe )
      rankMe = PIPS_MPIgetRank();
   if( -1 == numProcs )
      numProcs = PIPS_MPIgetSize();
}

void sTreeCallbacks::addChild( sTreeCallbacks* child )
{
   N += child->N;

   if( child->MY > 0 )
      MY += child->MY;

   if( child->MZ > 0 )
      MZ += child->MZ;

   child->np = nx_active;

   MYL = child->MYL;
   MZL = child->MZL;
   myl_active = child->myl_active;
   mzl_active = child->mzl_active;

   if( child->commWrkrs != MPI_COMM_NULL )
   {
      assert( child->myl_active <= child->MYL );
      assert( child->mzl_active <= child->MZL );
   }

   this->children.push_back(child);
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
   assertTreeStructureCorrect();
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
   assertTreeStructureCorrect();
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

   N_INACTIVE = nxVecSimple.length();
   MY_INACTIVE = myVecSimple.length();
   MZ_INACTIVE = mzVecSimple.length();

   nx_inactive = N_INACTIVE;
   my_inactive = MY_INACTIVE;
   mz_inactive = MZ_INACTIVE;

   if( myVecSimpleLink != nullptr )
   {
      assert(np == -1);
      MYL_INACTIVE = myVecSimpleLink->length();
      myl_inactive = MYL_INACTIVE;
   }
   else
   {
      MYL_INACTIVE =  mylParent;
      myl_inactive = mylParent;
   }

   if( mzVecSimpleLink != nullptr )
   {
      assert(np == -1);
      MZL_INACTIVE = mzVecSimpleLink->length();
      mzl_inactive = MZL_INACTIVE;
   }
   else
   {
      MZL_INACTIVE = mzlParent;
      mzl_inactive = mzlParent;
   }

   // empty child?
   if( N_INACTIVE == 0 && MY_INACTIVE == 0 && MZ_INACTIVE == 0 && np != -1 )
   {
      MYL_INACTIVE = 0;
      MZL_INACTIVE = 0;
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
   }

   hasPresolvedData = true;

}


void sTreeCallbacks::writeSizes( std::ostream& sout ) const
{
   const int myRank = PIPS_MPIgetRank(commWrkrs);

   MPI_Barrier(commWrkrs);
   if( myRank == 0 )
   {
      sout << "N          : " << N << "\n";
      sout << "MY         : " << MY << "\n";
      sout << "MZ         : " << MZ << "\n";
      sout << "MYL        : " << MYL << "\n";
      sout << "MZL        : " << MZL << "\n";
      sout << "nx_active  : " << nx_active << "\n";
      sout << "my_active  : " << my_active << "\n";
      sout << "mz_active  : " << mz_active << "\n";
      sout << "myl_active : " << myl_active << "\n";
      sout << "mzl_active : " << mzl_active << "\n";
   }

   if( sub_root )
   {
      sout << "subroot: \n";
      dynamic_cast<sTreeCallbacks*>(sub_root)->writeSizes(sout);
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

   if( data && sTree::isInVector( rankMe, myProcs ) )
   {
      // callbacks can be used for sizes
      assert( data->nCall || data->n != -1 );
      assert( data->myCall || data->my != -1 );
      assert( data->mzCall || data->mz != -1 );

      /* load sizes from callbacks */
      if( data->n == -1 )
         data->nCall(data->user_data, data->id, &data->n);
      if( data->my == -1 )
         data->myCall(data->user_data, data->id, &data->my);
      if( data->mz == -1 )
         data->mzCall(data->user_data, data->id, &data->mz);

      if( data->mylCall )
         data->mylCall(data->user_data, data->id, &data->myl);
      else if( data->myl == -1 )
         data->myl = 0;

      if( data->mzlCall )
         data->mzlCall(data->user_data, data->id, &data->mzl);
      else if( data->mzl == -1 )
         data->mzl = 0;

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
      children[it]->np = nx_active;
      children[it]->computeGlobalSizes();
      N += children[it]->N;
      MY += children[it]->MY;
      MZ += children[it]->MZ;
   }
   assertTreeStructureCorrect();
}

void sTreeCallbacks::assertTreeStructureChildren() const
{
   for( const sTree* child : children )
   {
      dynamic_cast<const sTreeCallbacks*>(child)->assertTreeStructureCorrect();

      if( child->commWrkrs != MPI_COMM_NULL )
      {
         assert( sTree::isInVector(rankMe, child->myProcs) );
         if( !is_hierarchical_root )
            assert( child->np == nx_active );
      }

      assert( std::includes(myProcs.begin(), myProcs.end(), child->myProcs.begin(), child->myProcs.end()) );
   }
}

void sTreeCallbacks::assertSubRoot() const
{
   if( sub_root )
   {
      assert( children.size() == 0 );
      dynamic_cast<const sTreeCallbacks*>(sub_root)->assertTreeStructureCorrect();
   }
}

void sTreeCallbacks::assertTreeStructureIsNotMyNode() const
{
   assert( !sTree::isInVector(rankMe, myProcs) );

   assert( nx_active == 0 );
   assert( my_active == 0 );
   assert( mz_active == 0 );
   assert( myl_active == 0 );
   assert( mzl_active == 0 );
   assert( N == 0 );
   assert( MY == 0 );
   assert( MZ == 0 );
   assert( MYL == 0 );
   assert( MZL == 0 );
}

void sTreeCallbacks::assertTreeStructureIsMyNodeChildren() const
{
   assert( !sub_root );

   int NX_children{0};
   int MY_children{0};
   int MZ_children{0};

   int MYL_children{0};
   int MZL_children{0};

   for( const sTree* child_ : children )
   {
      const sTreeCallbacks* child = dynamic_cast<const sTreeCallbacks*>(child_);

      if( sTree::isInVector(rankMe, child->myProcs) )
      {
         NX_children += child->N;
         MY_children += child->MY;
         MZ_children += child->MZ;

         assert( child->MYL <= MYL );
         assert( child->MZL <= MZL );
         if( is_hierarchical_inner_root )
         {
            assert( child->MYL >= myl_active );
            assert( child->MZL >= mzl_active );

            MYL_children += child->MYL - myl_active;
            MZL_children += child->MZL - mzl_active;
         }
         else if( !is_hierarchical_root )
         {
            if( child->is_hierarchical_inner_leaf )
            {
               assert( MZL >= child->MZL );
               assert( MYL >= child->MYL );
            }
            else
            {
               assert( MZL == child->MZL );
               assert( MYL == child->MYL );
            }


            assert( myl_active <= child->MYL );
            assert( mzl_active <= child->MZL );
            assert( child->myl_active == myl_active );
            assert( child->mzl_active == mzl_active );
         }
         else
         {
            assert( child->MYL <= MYL );
            assert( child->MZL <= MZL );

            MYL_children += child->MYL;
            MZL_children += child->MZL;
         }
      }
   }

   assert( N == NX_children + nx_active );
   assert( MY == MY_children + my_active );
   assert( MZ == MZ_children + mz_active );
   assert( MYL == MYL_children + myl_active );
   assert( MZL == MZL_children + mzl_active );
}

void sTreeCallbacks::assertTreeStructureIsMyNodeSubRoot() const
{
   assert( sub_root );
   assert( sTree::isInVector(rankMe, myProcs) );
   assert( sub_root->np = -1 );
   assert( is_hierarchical_inner_leaf );

   assert( MYL == myl_active + sub_root->MYL );
   assert( MZL == mzl_active + sub_root->MZL );
   assert( nx_active == sub_root->N );
   assert( my_active == sub_root->MY );
   assert( mz_active == sub_root->MZ );
}

void sTreeCallbacks::assertTreeStructureIsMyNode() const
{
   assert( sTree::isInVector(rankMe, myProcs) );

   if( children.size() != 0 )
      assertTreeStructureIsMyNodeChildren();
   else if( sub_root )
      assertTreeStructureIsMyNodeSubRoot();
}

void sTreeCallbacks::assertTreeStructureCorrect() const
{
   assert( myProcs.size() <= children.size() || children.size() == 0 );
   assert( is_sorted(myProcs.begin(), myProcs.end()) );

   assertTreeStructureChildren();
   assertSubRoot();

   if( commWrkrs == MPI_COMM_NULL )
      assertTreeStructureIsNotMyNode();
   else
      assertTreeStructureIsMyNode();
}

StochSymMatrix* sTreeCallbacks::createQ() const
{
   assert(!is_hierarchical_root && !is_hierarchical_inner_root && !is_hierarchical_inner_leaf );

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

StochGenMatrix* sTreeCallbacks::createMatrix( TREE_SIZE MY, TREE_SIZE MYL, DATA_INT m_ABmat, DATA_INT n_Mat,
      DATA_INT nnzAmat, DATA_NNZ fnnzAmat, DATA_MAT Amat, DATA_INT nnzBmat,
      DATA_NNZ fnnzBmat, DATA_MAT Bmat, DATA_INT m_Blmat, DATA_INT nnzBlmat,
      DATA_NNZ fnnzBlmat, DATA_MAT Blmat ) const
{
   assert(!is_hierarchical_root && !is_hierarchical_inner_root && !is_hierarchical_inner_leaf );

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

      if( data->*fnnzBlmat != nullptr )
      {
         // populate B with A's data B_0 is the A_0 from the theoretical form; also fill Bl
         // (i.e. the first block of linking constraints)
         A = new StochGenMatrix(this->*MY + this->*MYL, N,
               data->*m_ABmat, np, data->*nnzBmat,
               data->*m_ABmat, data->*n_Mat, data->*nnzAmat,
               data->*m_Blmat, data->*n_Mat, data->*nnzBlmat,
               commWrkrs);
      }
      else
      {
         // populate B with A's data B_0 is the A_0 from the theoretical form
         A = new StochGenMatrix(this->*MY + this->*MYL, N,
               data->*m_ABmat, np, data->*nnzBmat,
               data->*m_ABmat, data->*n_Mat,  data->*nnzAmat,
               commWrkrs);
      }

      //populate submatrix B
      (data->*Amat)(data->user_data, data->id, dynamic_cast<SparseGenMatrix*>(A->Bmat)->krowM(),
            dynamic_cast<SparseGenMatrix*>(A->Bmat)->jcolM(), dynamic_cast<SparseGenMatrix*>(A->Bmat)->M());

      if( print_tree_sizes_on_reading )
         printf("root  -- my=%d  myl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n",
               data->*m_ABmat, data->*m_Blmat, data->*n_Mat, np, data->*nnzAmat, data->*nnzBmat, data->*nnzBlmat);
   }
   else
   {
      if( data->*nnzBmat < 0 )
         (data->*fnnzBmat)(data->user_data, data->id, &(data->*nnzBmat));

      if( data->fnnzBl != nullptr )
      {
         A = new StochGenMatrix(this->*MY, N,
               data->*m_ABmat, np, data->*nnzAmat,
               data->*m_ABmat, data->*n_Mat, data->*nnzBmat,
               data->*m_Blmat, data->*n_Mat, data->*nnzBlmat,
               commWrkrs);
      }
      else
      {
         A = new StochGenMatrix(this->*MY, N,
               data->*m_ABmat, np, data->*nnzAmat,
               data->*m_ABmat, data->*n_Mat,  data->*nnzBmat,
               commWrkrs);
      }

      //populate the submatrices A, B
      (data->*Amat)(data->user_data, data->id, dynamic_cast<SparseGenMatrix*>(A->Amat)->krowM(),
            dynamic_cast<SparseGenMatrix*>(A->Amat)->jcolM(), dynamic_cast<SparseGenMatrix*>(A->Amat)->M());
      (data->*Bmat)(data->user_data, data->id, dynamic_cast<SparseGenMatrix*>(A->Bmat)->krowM(),
            dynamic_cast<SparseGenMatrix*>(A->Bmat)->jcolM(), dynamic_cast<SparseGenMatrix*>(A->Bmat)->M());

      if( print_tree_sizes_on_reading )
         printf("  -- my=%d  myl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n",
               data->*m_ABmat, data->*m_Blmat, data->*n_Mat, np, data->*nnzAmat, data->*nnzBmat, data->*nnzBlmat);
   }

   // populate Bl if existent
   if( data->*Blmat )
      (data->*Blmat)(data->user_data, data->id, dynamic_cast<SparseGenMatrix*>(A->Blmat)->krowM(),
            dynamic_cast<SparseGenMatrix*>(A->Blmat)->jcolM(), dynamic_cast<SparseGenMatrix*>(A->Blmat)->M());

   for(size_t it = 0; it < children.size(); it++)
   {
      StochGenMatrix* child = dynamic_cast<sTreeCallbacks*>(children[it])->createMatrix( MY, MYL, m_ABmat, n_Mat, nnzAmat, fnnzAmat, Amat,
            nnzBmat, fnnzBmat, Bmat, m_Blmat, nnzBlmat, fnnzBlmat, Blmat );
      A->AddChild(child);
   }
   return A;
}

StochGenMatrix* sTreeCallbacks::createA() const
{
   TREE_SIZE MY = &sTree::MY;
   TREE_SIZE MYL = &sTree::MYL;

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

   return createMatrix( MY, MYL, m_ABmat, n_Mat, nnzAmat, fnnzAmat, Amat, nnzBmat,
         fnnzBmat, Bmat, m_Blmat, nnzBlmat, fnnzBlmat, Blmat );
}

StochGenMatrix* sTreeCallbacks::createC() const
{
   TREE_SIZE MZ= &sTree::MZ;
   TREE_SIZE MZL = &sTree::MZL;

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

   return createMatrix( MZ, MZL, m_CDmat, n_Mat, nnzCmat, fnnzCmat, Cmat, nnzDmat,
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

   assert( !(is_hierarchical_root || is_hierarchical_inner_root || is_hierarchical_inner_leaf )
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
   if( PIPS_MPIgetRank() == 0 && !pips_options::getBoolParameter("SILENT") )
      std::cout << "Trimming " << nx_to_shave << " vars, " << myl_to_shave << " dense equalities, and " <<
         mzl_to_shave << " inequalities for the border\n";

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

   top_layer->MZ = MZ;
   top_layer->MZL = MZL;
   top_layer->np = -1;

   assert( IPMIterExecTIME == -1 );
   top_layer->children.push_back(this);
   top_layer->numProcs = numProcs;

   // TODO: not sure about the ressources monitors..: resMon, iterMon..
   top_layer->is_hierarchical_root = true;

   /* sTreeCallbacks members */
   top_layer->nx_active = nx_to_shave;
   this->nx_active -= nx_to_shave;

   top_layer->my_active = 0;
   top_layer->mz_active = 0;

   top_layer->myl_active = myl_to_shave;
   this->adjustActiveMylBy(-myl_to_shave);

   top_layer->mzl_active = mzl_to_shave;
   this->adjustActiveMzlBy(-mzl_to_shave);

   for( auto& child : children )
      dynamic_cast<sTreeCallbacks*>(child)->np = nx_active;

   assert( myl_active >= 0 );
   assert( mzl_active >= 0 );
   top_layer->assertTreeStructureCorrect();
   return top_layer;
}

void sTreeCallbacks::mapChildrenToNSubTrees( std::vector<unsigned int>& map_child_to_sub_tree, unsigned int n_children, unsigned int n_subtrees )
{
   map_child_to_sub_tree.clear();
   map_child_to_sub_tree.reserve(n_children);

   if( n_subtrees == 0 )
      return;

   const unsigned int everyone_gets = std::floor(n_children / n_subtrees);
   const unsigned int n_leftovers = n_children % n_subtrees;

   std::vector<unsigned int> children_per_tree( n_subtrees, everyone_gets );

   if( n_leftovers > 0 )
   {
      const unsigned int free_in_leftover_row = n_subtrees - n_leftovers;

      const unsigned int free_after_each_leftover = std::floor(free_in_leftover_row / n_leftovers);
      const unsigned int additional_frees = free_in_leftover_row % n_leftovers;

      unsigned int additional_frees_assigned = 0;
      for( unsigned int i = 0; i < n_subtrees; i += free_after_each_leftover )
      {
         ++children_per_tree[i];
         ++i;

         if( additional_frees != additional_frees_assigned )
         {
            ++i;
            ++additional_frees_assigned;
         }
      }
      assert( additional_frees_assigned == additional_frees );
   }

   assert( std::accumulate( children_per_tree.begin(), children_per_tree.end(),
         decltype(children_per_tree)::value_type(0) ) == n_children );

   for( unsigned int i = 0; i < children_per_tree.size(); ++i )
   {
      for( unsigned int j = 0; j < children_per_tree[i]; ++j )
         map_child_to_sub_tree.push_back( i );
   }
}

unsigned int sTreeCallbacks::getMapChildrenToSqrtNSubTrees( std::vector<unsigned int>& map_child_to_sub_tree, unsigned int n_children )
{
   const unsigned int n_new_roots = std::round( std::sqrt( n_children ) );

   mapChildrenToNSubTrees(map_child_to_sub_tree, n_children, n_new_roots);
   return n_new_roots;
}

void sTreeCallbacks::createSubcommunicatorsAndChildren( std::vector<unsigned int>& map_child_to_sub_tree )
{
   // TODO : maybe adjust the split at some point - the first n_leftovers procs have one additional block assigned to them .. this will lead to some imbalance in the tree here..
   PIPS_MPIabortIf( children.size() < 3, "Need at lest 4 child tree nodes to split a root node");

   const unsigned int n_new_roots = getMapChildrenToSqrtNSubTrees( map_child_to_sub_tree, children.size() );
   assert( map_child_to_sub_tree.size() == children.size() );

   /* create new sub-roots */
   std::vector<sTreeCallbacks*> new_leafs(n_new_roots);
   for( auto& leaf : new_leafs )
   {
      leaf = new sTreeCallbacks();
      leaf->setHierarchicalInnerLeaf();

      leaf->sub_root = new sTreeCallbacks();
      leaf->sub_root->setHierarchicalInnerRoot();
   }

   /* determine processes for each subtree and assign children */
   for( size_t child = 0; child < children.size(); ++child )
   {
      const auto& child_procs = children[child]->myProcs;
      assert( child_procs.size() >= 1 );

      const unsigned int assigned_sub_root_for_child = map_child_to_sub_tree[ child ];
      assert( assigned_sub_root_for_child < new_leafs.size() );
      sTreeCallbacks* assigned_leaf = new_leafs[ assigned_sub_root_for_child ];

      for( int process : child_procs )
      {
         // assuming sorted..
         if( assigned_leaf->myProcs.empty() || assigned_leaf->myProcs.back() != process )
         {
            if( !assigned_leaf->myProcs.empty() )
            {
               assert( process > assigned_leaf->myProcs.back() );
               assert( !isInVector( process, assigned_leaf->myProcs ) );
            }

            assigned_leaf->myProcs.push_back( process );
         }
      }

      dynamic_cast<sTreeCallbacks*>(assigned_leaf->sub_root)->addChild( dynamic_cast<sTreeCallbacks*>(children[child]) );
   }

   /* create all sub-communicators */
   for( auto new_leaf : new_leafs )
   {
      new_leaf->commWrkrs = PIPS_MPIcreateGroupFromRanks( new_leaf->myProcs );

      if( !isInVector( rankMe, new_leaf->myProcs) )
         assert( new_leaf->commWrkrs == MPI_COMM_NULL );

      new_leaf->sub_root->commWrkrs = new_leaf->commWrkrs;
      new_leaf->sub_root->myProcs = new_leaf->myProcs;
   }

#ifndef NDEBUG
   for( size_t child = 0; child < children.size(); ++child )
   {
      if( isInVector(rankMe, children[child]->myProcs ) )
      {
         auto child_new_root = new_leafs[map_child_to_sub_tree[child]];
         assert( isInVector(rankMe, child_new_root->myProcs) );
         assert( child_new_root->commWrkrs != MPI_COMM_NULL );
      }
   }
#endif

   /* add sub_roots as this new children */
   children.clear();
   children.insert( this->children.begin(), new_leafs.begin(), new_leafs.end() );
   is_hierarchical_inner_root = true;
}

void sTreeCallbacks::countTwoLinksForChildTrees(const std::vector<int>& two_links_start_in_child_A, const std::vector<int>& two_links_start_in_child_C,
      std::vector<unsigned int>& two_links_children_eq, std::vector<unsigned int>& two_links_children_ineq,
      unsigned int& two_links_root_eq, unsigned int& two_links_root_ineq ) const
{
   const unsigned int n_children = children.size();
#ifndef NDEBUG
   const unsigned int n_leafs = map_node_sub_root.size();
#endif
   assert( is_hierarchical_inner_root );
   assert( n_leafs == two_links_start_in_child_A.size() );
   assert( n_leafs == two_links_start_in_child_C.size() );

   two_links_children_eq.clear();
   two_links_children_eq.resize(n_children);
   two_links_children_ineq.clear();
   two_links_children_ineq.resize(n_children);

   /* count two links for all new sub-communicators */
   unsigned int leaf = 0;
   for( unsigned int i = 0; i < n_children; ++i )
   {
      const sTreeCallbacks& sub_root = dynamic_cast<const sTreeCallbacks&>( *children[i]->sub_root );

      for( size_t j = 0; j < sub_root.children.size(); ++j )
      {
         assert( leaf < n_leafs );

         /* the two links of each last child will stay in the root node */
         if( j == sub_root.children.size() - 1 )
         {
            two_links_root_eq += two_links_start_in_child_A[leaf];
            two_links_root_ineq += two_links_start_in_child_C[leaf];
         }
         else
         {
            two_links_children_eq[i] += two_links_start_in_child_A[leaf];
            two_links_children_ineq[i] += two_links_start_in_child_C[leaf];
         }

         ++leaf;
      }
   }
   assert( leaf == n_leafs );
}

void sTreeCallbacks::adjustActiveMylBy( int adjustment )
{
   assert( !is_hierarchical_inner_leaf );
   myl_active += adjustment;
   MYL += adjustment;

   assert( myl_active >= 0 );
   assert( MYL >= 0 );

   for( sTree* child_ : children )
   {
      sTreeCallbacks* child = dynamic_cast<sTreeCallbacks*>(child_);

      assert( child->children.size() == 0 );
      child->myl_active += adjustment;
      child->MYL += adjustment;
      assert( child->myl_active >= 0 );
      assert( child->MYL >= 0 );
   }
}

void sTreeCallbacks::adjustActiveMzlBy( int adjustment )
{
   assert( !is_hierarchical_inner_leaf );
   mzl_active += adjustment;
   MZL += adjustment;

   for( auto& child_ : children )
   {
      sTreeCallbacks* child = dynamic_cast<sTreeCallbacks*>(child_);

      assert( child->children.size() == 0 );
      child->mzl_active += adjustment;
      child->MZL += adjustment;
      assert( child->mzl_active >= 0 );
      assert( child->MZL >= 0 );
   }
}

void sTreeCallbacks::adjustSizesAfterSplit( const std::vector<unsigned int>& two_links_children_eq,
      const std::vector<unsigned int>& two_links_children_ineq, unsigned int two_links_children_eq_sum, unsigned int two_links_children_ineq_sum)
{
   assert( is_hierarchical_inner_root );
   assert( static_cast<unsigned int>(myl_active) >= two_links_children_eq_sum );
   assert( static_cast<unsigned int>(mzl_active) >= two_links_children_ineq_sum );

   myl_active -= two_links_children_eq_sum;
   mzl_active -= two_links_children_ineq_sum;

   /* recompute sizes for the children */
   for( size_t i = 0; i < children.size(); ++i )
   {
      sTreeCallbacks& inner_leaf = dynamic_cast<sTreeCallbacks&>(*children[i]);
      sTreeCallbacks& sub_root = dynamic_cast<sTreeCallbacks&>(*inner_leaf.sub_root);
      assert( inner_leaf.is_hierarchical_inner_leaf );

      if( isInVector(rankMe, inner_leaf.myProcs) )
      {
         inner_leaf.MYL = myl_active + two_links_children_eq[i];
         inner_leaf.myl_active = myl_active;
         if( MYL >= 0 )
            sub_root.adjustActiveMylBy( -myl_active - (two_links_children_eq_sum - two_links_children_eq[i]) );

         inner_leaf.MZL = mzl_active + two_links_children_ineq[i];
         inner_leaf.mzl_active = mzl_active;
         if( MZL >= 0 )
            sub_root.adjustActiveMzlBy( -mzl_active - (two_links_children_ineq_sum - two_links_children_ineq[i]) );

         inner_leaf.np = nx_active;

         inner_leaf.N = sub_root.N;
         inner_leaf.nx_active = sub_root.N;

         inner_leaf.MY = sub_root.MY;
         inner_leaf.my_active = sub_root.MY;

         inner_leaf.MZ = sub_root.MZ;
         inner_leaf.mz_active = sub_root.MZ;
      }
      else
      {
         inner_leaf.N = inner_leaf.MY = inner_leaf.MZ = inner_leaf.MYL = inner_leaf.MZL = 0;
         inner_leaf.nx_active = inner_leaf.my_active = inner_leaf.mz_active = inner_leaf.myl_active = inner_leaf.mzl_active = 0;
      }

   }
}


void sTreeCallbacks::splitTreeSquareRoot( const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC )
{
   assert( commWrkrs != MPI_COMM_NULL );
   assert( !is_hierarchical_root );

#ifndef NDEBUG
   const size_t n_old_leafs = children.size();
#endif

   createSubcommunicatorsAndChildren( map_node_sub_root );
   assert( map_node_sub_root.size() == n_old_leafs );

   std::vector<unsigned int> two_links_children_eq, two_links_children_ineq;
   unsigned int two_links_root_eq{0}; unsigned int two_links_root_ineq{0};

   countTwoLinksForChildTrees(twoLinksStartBlockA, twoLinksStartBlockC, two_links_children_eq,
         two_links_children_ineq, two_links_root_eq, two_links_root_ineq );

   const unsigned int two_links_children_eq_sum = std::accumulate( two_links_children_eq.begin(), two_links_children_eq.end(), unsigned(0) );
   const unsigned int two_links_children_ineq_sum = std::accumulate( two_links_children_ineq.begin(), two_links_children_ineq.end(), unsigned(0) );

   if( rankMe == 0 && !pips_options::getBoolParameter("SILENT") )
   {
      std::cout << "Splitting node into " << this->children.size() << " subroots\n";
      std::cout << "Splitting " << two_links_children_eq_sum + two_links_root_eq << " equality two-links into " << two_links_root_eq
            << " root and " << two_links_children_eq_sum << " child links\n";
      std::cout << "Splitting " << two_links_children_ineq_sum + two_links_root_ineq << " inequality two-links into " << two_links_root_ineq
            << " root and " << two_links_children_ineq_sum << " child links\n";
   }

   adjustSizesAfterSplit(two_links_children_eq, two_links_children_ineq, two_links_children_eq_sum, two_links_children_ineq_sum);

   assertTreeStructureCorrect();
}


sTree* sTreeCallbacks::switchToHierarchicalTree( int nx_to_shave, int myl_to_shave, int mzl_to_shave,
      const std::vector<int>& twoLinksStartBlockA, const std::vector<int>& twoLinksStartBlockC )
{
   assert( !is_hierarchical_root );
   assert( np == -1 );
   assertTreeStructureCorrect();

   assert( nx_to_shave >= 0 );
   assert( myl_to_shave >=0 );
   assert( mzl_to_shave >= 0 );

   assert( nx_to_shave <= nx_active );
   assert( myl_to_shave <= myl_active );
   assert( mzl_to_shave <= mzl_active );

   /* distributed preconditioner must be deactivated */
   assert( !distributedPreconditionerActive() );

   this->splitTreeSquareRoot( twoLinksStartBlockA, twoLinksStartBlockC );

   sTreeCallbacks* top_layer = dynamic_cast<sTreeCallbacks*>( shaveDenseBorder( nx_to_shave, myl_to_shave, mzl_to_shave ) );

   return top_layer;
}

sTree* sTreeCallbacks::collapseDenseBorder()
{
   /* this must happen at MPI_COMM_WORLD level */
   assert( is_hierarchical_root );
   assert( rankMe == PIPS_MPIgetRank() );
   assert( numProcs == PIPS_MPIgetSize() );
   assert( children.size() == 1 );

   sTreeCallbacks* new_top = dynamic_cast<sTreeCallbacks*>(children[0]);
   children.clear();

   commWrkrs = MPI_COMM_NULL;
   myProcs.clear();

   new_top->N = N;
   N = -1;
   MY = -1;
   new_top->MYL = MYL;
   MYL = -1;
   MZ = -1;
   new_top->MZL = MZL;
   MZL = -1;

   numProcs = -1;

   new_top->nx_active += nx_active;
   nx_active = -1;
   assert( my_active == - 1 );
   assert( mz_active == - 1 );

   new_top->myl_active += myl_active;
   myl_active = -1;
   new_top->mzl_active += mzl_active;
   mzl_active = -1;

   return new_top;
}

sTree* sTreeCallbacks::collapseHierarchicalTree()
{
   return collapseDenseBorder();
}

std::vector<MPI_Comm> sTreeCallbacks::getChildComms() const
{
   std::vector<MPI_Comm> comms(children.size(), MPI_COMM_NULL);

   for( unsigned int i = 0; i < children.size(); ++i )
      comms[i] = children[i]->commWrkrs;

   return comms;
}

void sTreeCallbacks::splitMatrixAccordingToTree( StochSymMatrix& /*mat*/ ) const
{

   assert( false && "TODO : implement " );

}

void sTreeCallbacks::splitMatrixAccordingToTree( StochGenMatrix& /*mat*/ ) const
{
   assert( false && "TODO : implement " );

}

void sTreeCallbacks::splitVectorAccordingToTree( StochVector& vec ) const
{
   assert( map_node_sub_root.size() == vec.children.size() );

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
