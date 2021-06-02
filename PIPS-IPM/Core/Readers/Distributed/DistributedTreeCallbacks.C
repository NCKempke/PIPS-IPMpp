/* PIPS-IPM                                                           *
 * Author:  Cosmin G. Petra                                           *
 * (C) 2012 Argonne National Laboratory. See Copyright Notification.  */

#include "DistributedTreeCallbacks.h"
#include "DistributedQP.hpp"
#include "DistributedSymmetricMatrix.h"
#include "DistributedMatrix.h"
#include "DistributedVector.h"
#include "SimpleVector.h"
#include <cmath>
#include <algorithm>    // std::swap
#include <numeric>


DistributedTreeCallbacks::DistributedTreeCallbacks(const DistributedTreeCallbacks& other) : DistributedTree(other), N_INACTIVE{other.N_INACTIVE}, MY_INACTIVE{other.MY_INACTIVE},
      MZ_INACTIVE{other.MZ_INACTIVE}, MYL_INACTIVE{other.MYL_INACTIVE}, MZL_INACTIVE{other.MZL_INACTIVE}, nx_active{other.nx_active},
      my_active{other.my_active}, mz_active{other.mz_active}, myl_active{other.myl_active}, mzl_active{other.mzl_active},
      nx_inactive{other.nx_inactive}, my_inactive{other.my_inactive}, mz_inactive{other.mz_inactive}, myl_inactive{other.myl_inactive},
      mzl_inactive{other.mzl_inactive}, is_data_presolved{other.is_data_presolved}, has_presolved_data{other.has_presolved_data},
      print_tree_sizes_on_reading{other.print_tree_sizes_on_reading},
      map_node_sub_root(other.map_node_sub_root.begin(), other.map_node_sub_root.end()), data{other.data} {
}

DistributedTree* DistributedTreeCallbacks::clone() const {
   return new DistributedTreeCallbacks(*this);
}

DistributedTreeCallbacks::DistributedTreeCallbacks() : print_tree_sizes_on_reading{pipsipmpp_options::get_bool_parameter("PRINT_TREESIZES_ON_READ")} {
   if (-1 == rankMe)
      rankMe = PIPS_MPIgetRank();
   if (-1 == numProcs)
      numProcs = PIPS_MPIgetSize();
}

DistributedTreeCallbacks::DistributedTreeCallbacks(DistributedInputTree* inputTree)
      : DistributedTree(), print_tree_sizes_on_reading{pipsipmpp_options::get_bool_parameter("PRINT_TREESIZES_ON_READ")}, data{inputTree->nodeInput} {
   if (-1 == rankMe)
      rankMe = PIPS_MPIgetRank();
   if (-1 == numProcs)
      numProcs = PIPS_MPIgetSize();

   for (auto & it : inputTree->children)
      children.push_back(new DistributedTreeCallbacks(it));
}

void DistributedTreeCallbacks::addChild(DistributedTreeCallbacks* child) {
   N += child->N;

   if (child->MY > 0)
      MY += child->MY;

   if (child->MZ > 0)
      MZ += child->MZ;

   child->np = nx_active;

   if (child->commWrkrs != MPI_COMM_NULL) {
      MYL = child->MYL;
      MZL = child->MZL;
      myl_active = child->myl_active;
      mzl_active = child->mzl_active;
   }

   if (child->commWrkrs != MPI_COMM_NULL) {
      assert(containsSorted(child->myProcs, myProcs));
      assert(child->myl_active <= child->MYL);
      assert(child->mzl_active <= child->MZL);
   }

   this->children.push_back(child);
}


void DistributedTreeCallbacks::switchToPresolvedData() {
   assert(!is_hierarchical_root || (false && "cannot be used with hierarchical data"));
   assert(!is_data_presolved);
   assert(has_presolved_data);

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

   for (auto & it : children) {
      it->np = this->nx_active;
      dynamic_cast<DistributedTreeCallbacks*>(it)->switchToPresolvedData();
   }

   is_data_presolved = true;
   assertTreeStructureCorrect();
}


void DistributedTreeCallbacks::switchToOriginalData() {
   assert(!is_hierarchical_root || (false && "cannot be used with hierarchical data"));
   assert(is_data_presolved);

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

   for (auto & it : children) {
      it->np = this->nx_active;
      dynamic_cast<DistributedTreeCallbacks*>(it)->switchToOriginalData();
   }

   is_data_presolved = false;
   assertTreeStructureCorrect();
}


bool DistributedTreeCallbacks::isPresolved() {
   assert(!is_hierarchical_root || (false && "cannot be used with hierarchical data"));
   return is_data_presolved;
}

bool DistributedTreeCallbacks::hasPresolved() {
   assert(!is_hierarchical_root || (false && "cannot be used with hierarchical data"));
   return has_presolved_data;
}

void DistributedTreeCallbacks::initPresolvedData(const DistributedQP& presolved_data) {
   const auto& Q = dynamic_cast<const DistributedSymmetricMatrix&>(*presolved_data.Q);
   const auto& A = dynamic_cast<const DistributedMatrix&>(*presolved_data.A);
   const auto& C = dynamic_cast<const DistributedMatrix&>(*presolved_data.C);

   const auto& g = dynamic_cast<const DistributedVector<double>&>(*presolved_data.g);
   const auto& b = dynamic_cast<const DistributedVector<double>&>(*presolved_data.bA);

   assert(presolved_data.icupp || presolved_data.icupp);
   const DistributedVector<double>& ic = presolved_data.iclow ? dynamic_cast<const DistributedVector<double>&>(*presolved_data.iclow)
                                                                       : dynamic_cast<const DistributedVector<double>&>(*presolved_data.icupp);

   initPresolvedData(Q, A, C, g, b, ic, -1, -1);
}

void
DistributedTreeCallbacks::initPresolvedData(const DistributedSymmetricMatrix& Q, const DistributedMatrix& A, const DistributedMatrix& C, const DistributedVector<double>& nxVec,
      const DistributedVector<double>& myVec, const DistributedVector<double>& mzVec, int mylParent, int mzlParent) {
   assert(!is_hierarchical_root || (false && "cannot be used with hierarchical data"));
   assert(!has_presolved_data);

   assert(nxVec.children.size() == children.size());
   assert(myVec.children.size() == children.size());
   assert(mzVec.children.size() == children.size());
   assert(A.children.size() == children.size());
   assert(C.children.size() == children.size());

   const auto& nxVecSimple = dynamic_cast<const SimpleVector<double>&>(*nxVec.first);
   const auto& myVecSimple = dynamic_cast<const SimpleVector<double>&>(*myVec.first);
   const auto& mzVecSimple = dynamic_cast<const SimpleVector<double>&>(*mzVec.first);
   const auto* const myVecSimpleLink = dynamic_cast<const SimpleVector<double>*>(myVec.last.get());
   const auto* const mzVecSimpleLink = dynamic_cast<const SimpleVector<double>*>(mzVec.last.get());

   N_INACTIVE = nxVecSimple.length();
   MY_INACTIVE = myVecSimple.length();
   MZ_INACTIVE = mzVecSimple.length();

   nx_inactive = N_INACTIVE;
   my_inactive = MY_INACTIVE;
   mz_inactive = MZ_INACTIVE;

   if (myVecSimpleLink) {
      assert(np == -1);
      MYL_INACTIVE = myVecSimpleLink->length();
      myl_inactive = MYL_INACTIVE;
   }
   else {
      MYL_INACTIVE = mylParent;
      myl_inactive = mylParent;
   }

   if (mzVecSimpleLink) {
      assert(np == -1);
      MZL_INACTIVE = mzVecSimpleLink->length();
      mzl_inactive = MZL_INACTIVE;
   }
   else {
      MZL_INACTIVE = mzlParent;
      mzl_inactive = mzlParent;
   }

   // empty child?
   if (N_INACTIVE == 0 && MY_INACTIVE == 0 && MZ_INACTIVE == 0 && np != -1) {
      MYL_INACTIVE = 0;
      MZL_INACTIVE = 0;
      myl_inactive = 0;
      mzl_inactive = 0;
   }

   // are we at the root?
   if (np == -1) {
      assert(mylParent == -1);
      assert(mzlParent == -1);
   }

   for (size_t it = 0; it < children.size(); it++) {
      assert(children[it]->np == this->nx_active);
      auto* DistributedTreeCallbacksChild = dynamic_cast<DistributedTreeCallbacks*>(children[it]);

      DistributedTreeCallbacksChild->initPresolvedData(*Q.children[it], *A.children[it], *C.children[it], *nxVec.children[it], *myVec.children[it],
            *mzVec.children[it], myl_inactive, mzl_inactive);

      N_INACTIVE += DistributedTreeCallbacksChild->N_INACTIVE;
      MY_INACTIVE += DistributedTreeCallbacksChild->MY_INACTIVE;
      MZ_INACTIVE += DistributedTreeCallbacksChild->MZ_INACTIVE;
   }

   has_presolved_data = true;

}


void DistributedTreeCallbacks::writeSizes(std::ostream& sout) const {
   const int myRank = PIPS_MPIgetRank(commWrkrs);

   MPI_Barrier(commWrkrs);
   if (myRank == 0) {
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

   if (sub_root) {
      sout << "subroot: \n";
      dynamic_cast<DistributedTreeCallbacks*>(sub_root)->writeSizes(sout);
   }

   for (size_t it = 0; it < children.size(); it++) {
      if (children[it]->commWrkrs != MPI_COMM_NULL) {
         sout << "child " << it << ": \n\n";
         dynamic_cast<DistributedTreeCallbacks*>(children[it])->writeSizes(sout);
      }
      MPI_Barrier(commWrkrs);
   }

   if (myRank == 0)
      std::cout << "\n";
}


// this is usually called after assigning processes
void DistributedTreeCallbacks::computeGlobalSizes() {
   assert(!is_hierarchical_root || (false && "cannot be used with hierarchical data"));
   assert(!is_data_presolved);

   if (data && isInVector(rankMe, myProcs)) {
      // callbacks can be used for sizes
      assert(data->nCall || data->n != -1);
      assert(data->myCall || data->my != -1);
      assert(data->mzCall || data->mz != -1);

      /* load sizes from callbacks */
      if (data->n == -1)
         data->nCall(data->user_data, data->id, &data->n);
      if (data->my == -1)
         data->myCall(data->user_data, data->id, &data->my);
      if (data->mz == -1)
         data->mzCall(data->user_data, data->id, &data->mz);

      if (data->mylCall)
         data->mylCall(data->user_data, data->id, &data->myl);
      else if (data->myl == -1)
         data->myl = 0;

      if (data->mzlCall)
         data->mzlCall(data->user_data, data->id, &data->mzl);
      else if (data->mzl == -1)
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
   else {
      nx_active = my_active = mz_active = myl_active = mzl_active = 0;
      N = MY = MZ = MYL = MZL = 0;
   }

   for (auto & it : children) {
      it->np = nx_active;
      it->computeGlobalSizes();
      N += it->N;
      MY += it->MY;
      MZ += it->MZ;
   }
   assertTreeStructureCorrect();
}

void DistributedTreeCallbacks::assertTreeStructureChildren() const {
   for (const DistributedTree* child : children) {
      dynamic_cast<const DistributedTreeCallbacks*>(child)->assertTreeStructureCorrect();

      if (child->commWrkrs != MPI_COMM_NULL) {
         assert(isInVector(rankMe, child->myProcs));
         if (!is_hierarchical_root)
            assert(child->np == nx_active);
      }

      assert(std::includes(myProcs.begin(), myProcs.end(), child->myProcs.begin(), child->myProcs.end()));
   }
}

void DistributedTreeCallbacks::assertSubRoot() const {
   if (sub_root) {
      assert(children.empty());
      dynamic_cast<const DistributedTreeCallbacks*>(sub_root)->assertTreeStructureCorrect();
   }
}

void DistributedTreeCallbacks::assertTreeStructureIsNotMyNode() const {
   assert(!isInVector(rankMe, myProcs));

   assert(nx_active == 0);
   assert(my_active == 0);
   assert(mz_active == 0);
   assert(myl_active == 0);
   assert(mzl_active == 0);
   assert(N == 0);
   assert(MY == 0);
   assert(MZ == 0);

// might be non-zero in case of hierarchical approach
//   assert( MYL == 0 );
//   assert( MZL == 0 );
}

void DistributedTreeCallbacks::assertTreeStructureIsMyNodeChildren() const {
   assert(!sub_root);

   int NX_children{0};
   int MY_children{0};
   int MZ_children{0};

   int MYL_children{0};
   int MZL_children{0};

   for (const DistributedTree* child_ : children) {
      const auto* child = dynamic_cast<const DistributedTreeCallbacks*>(child_);

      if (isInVector(rankMe, child->myProcs)) {
         NX_children += child->N;
         MY_children += child->MY;
         MZ_children += child->MZ;

         assert(child->MYL <= MYL);
         assert(child->MZL <= MZL);

         if (child->sub_root) {
            assert(child->MYL >= myl_active);
            assert(child->MZL >= mzl_active);

            MYL_children += child->MYL - myl_active;
            MZL_children += child->MZL - mzl_active;
         }
         else if (!is_hierarchical_root) {
            if (child->is_hierarchical_inner_leaf) {
               assert(MZL >= child->MZL);
               assert(MYL >= child->MYL);
            }
            else {
               assert(MZL == child->MZL);
               assert(MYL == child->MYL);
            }


            assert(myl_active <= child->MYL);
            assert(mzl_active <= child->MZL);
            assert(child->myl_active == myl_active);
            assert(child->mzl_active == mzl_active);
         }
         else {
            assert(child->MYL <= MYL);
            assert(child->MZL <= MZL);

            MYL_children += child->MYL;
            MZL_children += child->MZL;
         }
      }
   }

   assert(N == NX_children + nx_active);
   assert(MY == MY_children + std::max(0, my_active));
   assert(MZ == MZ_children + std::max(0, mz_active));
   assert(MYL == MYL_children + myl_active);
   assert(MZL == MZL_children + mzl_active);
}

void DistributedTreeCallbacks::assertTreeStructureIsMyNodeSubRoot() const {
   assert(sub_root);
   assert(isInVector(rankMe, myProcs));
   assert(sub_root->np == -1);
   assert(is_hierarchical_inner_leaf);

   assert(MYL == myl_active + sub_root->MYL);
   assert(MZL == mzl_active + sub_root->MZL);
   assert(nx_active == sub_root->N);

   assert(my_active == sub_root->MY);
   assert(mz_active == sub_root->MZ);
}

void DistributedTreeCallbacks::assertTreeStructureIsMyNode() const {
   assert(isInVector(rankMe, myProcs));

   if (!children.empty())
      assertTreeStructureIsMyNodeChildren();
   else if (sub_root)
      assertTreeStructureIsMyNodeSubRoot();
}

void DistributedTreeCallbacks::assertTreeStructureCorrect() const {
   assert(std::is_sorted(myProcs.begin(), myProcs.end()));

   assertTreeStructureChildren();
   assertSubRoot();

   if (commWrkrs == MPI_COMM_NULL)
      assertTreeStructureIsNotMyNode();
   else
      assertTreeStructureIsMyNode();
}

DistributedSymmetricMatrix* DistributedTreeCallbacks::createQ() const {
   assert(!is_hierarchical_root && !is_hierarchical_inner_root && !is_hierarchical_inner_leaf);

   //is this node a dead-end for this process?
   if (commWrkrs == MPI_COMM_NULL)
      return new StochSymDummyMatrix();

   if (data->nnzQ < 0)
      data->fnnzQ(data->user_data, data->id, &data->nnzQ);

   auto* Q = new DistributedSymmetricMatrix(N, data->n, data->nnzQ, commWrkrs);

   data->fQ(data->user_data, data->id, dynamic_cast<SparseSymmetricMatrix&>(*Q->diag).krowM(), dynamic_cast<SparseSymmetricMatrix&>(*Q->diag).jcolM(),
         dynamic_cast<SparseSymmetricMatrix&>(*Q->diag).M());

   for (auto it : children) {
      std::shared_ptr<DistributedSymmetricMatrix> child{it->createQ()};
      Q->AddChild(child);
   }
   return Q;
}

DistributedMatrix*
DistributedTreeCallbacks::createMatrix(TREE_SIZE MY, TREE_SIZE MYL, DATA_INT m_ABmat, DATA_INT n_Mat, DATA_INT nnzAmat, DATA_NNZ fnnzAmat, DATA_MAT Amat,
      DATA_INT nnzBmat, DATA_NNZ fnnzBmat, DATA_MAT Bmat, DATA_INT m_Blmat, DATA_INT nnzBlmat, DATA_NNZ fnnzBlmat, DATA_MAT Blmat,
      const std::string& prefix_for_print) const {
   assert(!is_hierarchical_root && !is_hierarchical_inner_root && !is_hierarchical_inner_leaf);

   if (commWrkrs == MPI_COMM_NULL)
      return new StochGenDummyMatrix();

   const bool root = (np == -1);
   const bool has_linking = (data->*fnnzBlmat != nullptr);

   if (has_linking && data->*nnzBlmat < 0)
      (data->*fnnzBlmat)(data->user_data, data->id, &(data->*nnzBlmat));

   if (data->*nnzAmat < 0)
      (data->*fnnzAmat)(data->user_data, data->id, &(data->*nnzAmat));

   DistributedMatrix* A = nullptr;

   if (root) {
      data->*nnzBmat = 0;

      if (data->*fnnzBlmat != nullptr) {
         // populate B with A's data B_0 is the A_0 from the theoretical form; also fill Bl
         // (i.e. the first block of linking constraints)
         A = new DistributedMatrix(this->*MY + this->*MYL, N, data->*m_ABmat, np, data->*nnzBmat, data->*m_ABmat, data->*n_Mat, data->*nnzAmat,
               data->*m_Blmat, data->*n_Mat, data->*nnzBlmat, commWrkrs);
      }
      else {
         // populate B with A's data B_0 is the A_0 from the theoretical form
         A = new DistributedMatrix(this->*MY + this->*MYL, N, data->*m_ABmat, np, data->*nnzBmat, data->*m_ABmat, data->*n_Mat, data->*nnzAmat,
               commWrkrs);
      }

      //populate submatrix B
      (data->*Amat)(data->user_data, data->id, dynamic_cast<SparseMatrix&>(*A->Bmat).krowM(), dynamic_cast<SparseMatrix&>(*A->Bmat).jcolM(),
            dynamic_cast<SparseMatrix&>(*A->Bmat).M());

      if (print_tree_sizes_on_reading)
         printf("root  -- m%s=%d  m%sl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n", prefix_for_print.c_str(), data->*m_ABmat,
               prefix_for_print.c_str(), data->*m_Blmat, data->*n_Mat, np, data->*nnzAmat, data->*nnzBmat, data->*nnzBlmat);
   }
   else {
      assert(commWrkrs == MPI_COMM_SELF);

      if (data->*nnzBmat < 0)
         (data->*fnnzBmat)(data->user_data, data->id, &(data->*nnzBmat));

      if (data->fnnzBl != nullptr) {
         A = new DistributedMatrix(this->*MY, N, data->*m_ABmat, np, data->*nnzAmat, data->*m_ABmat, data->*n_Mat, data->*nnzBmat, data->*m_Blmat,
               data->*n_Mat, data->*nnzBlmat, commWrkrs);
      }
      else {
         A = new DistributedMatrix(this->*MY, N, data->*m_ABmat, np, data->*nnzAmat, data->*m_ABmat, data->*n_Mat, data->*nnzBmat, commWrkrs);
      }

      //populate the submatrices A, B
      (data->*Amat)(data->user_data, data->id, dynamic_cast<SparseMatrix&>(*A->Amat).krowM(), dynamic_cast<SparseMatrix&>(*A->Amat).jcolM(),
            dynamic_cast<SparseMatrix&>(*A->Amat).M());
      (data->*Bmat)(data->user_data, data->id, dynamic_cast<SparseMatrix&>(*A->Bmat).krowM(), dynamic_cast<SparseMatrix&>(*A->Bmat).jcolM(),
            dynamic_cast<SparseMatrix&>(*A->Bmat).M());

      if (print_tree_sizes_on_reading)
         printf("  -- m%s=%d  m%sl=%d nx=%d   1st stg nx=%d nnzA=%d nnzB=%d, nnzBl=%d\n", prefix_for_print.c_str(), data->*m_ABmat,
               prefix_for_print.c_str(), data->*m_Blmat, data->*n_Mat, np, data->*nnzAmat, data->*nnzBmat, data->*nnzBlmat);
   }

   // populate Bl if existent
   if (data->*Blmat)
      (data->*Blmat)(data->user_data, data->id, dynamic_cast<SparseMatrix&>(*A->Blmat).krowM(), dynamic_cast<SparseMatrix&>(*A->Blmat).jcolM(),
            dynamic_cast<SparseMatrix&>(*A->Blmat).M());

   for (auto it : children) {
      std::shared_ptr<DistributedMatrix> child{dynamic_cast<DistributedTreeCallbacks*>(it)->createMatrix(MY, MYL, m_ABmat, n_Mat, nnzAmat, fnnzAmat, Amat, nnzBmat,
            fnnzBmat, Bmat, m_Blmat, nnzBlmat, fnnzBlmat, Blmat, prefix_for_print)};
      A->AddChild(child);
   }
   return A;
}

DistributedMatrix* DistributedTreeCallbacks::createA() const {
   TREE_SIZE MY = &DistributedTree::MY;
   TREE_SIZE MYL = &DistributedTree::MYL;

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

   const std::string prefix = "y";
   return createMatrix(MY, MYL, m_ABmat, n_Mat, nnzAmat, fnnzAmat, Amat, nnzBmat, fnnzBmat, Bmat, m_Blmat, nnzBlmat, fnnzBlmat, Blmat, prefix);
}

DistributedMatrix* DistributedTreeCallbacks::createC() const {
   TREE_SIZE MZ = &DistributedTree::MZ;
   TREE_SIZE MZL = &DistributedTree::MZL;

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

   const std::string prefix = "z";
   return createMatrix(MZ, MZL, m_CDmat, n_Mat, nnzCmat, fnnzCmat, Cmat, nnzDmat, fnnzDmat, Dmat, m_Dlmat, nnzDlmat, fnnzDlmat, Dlmat, prefix);
}

int DistributedTreeCallbacks::nx() const {
   return nx_active;
}

int DistributedTreeCallbacks::my() const {
   return my_active;
}

int DistributedTreeCallbacks::myl() const {
   return myl_active;
}

int DistributedTreeCallbacks::mz() const {
   return mz_active;
}

int DistributedTreeCallbacks::mzl() const {
   return mzl_active;
}

int DistributedTreeCallbacks::id() const {
   return data->id;
}

DistributedVector<double>* DistributedTreeCallbacks::createVector(DATA_INT n_vec, DATA_VEC vec, DATA_INT n_linking_vec, DATA_VEC linking_vec) const {
   assert(n_vec);
   assert(vec);

   assert(!(is_hierarchical_root || is_hierarchical_inner_root || is_hierarchical_inner_leaf) || (false && "cannot be used with hierarchical data"));

   if (commWrkrs == MPI_COMM_NULL)
      return new DistributedDummyVector<double>();

   const int nlinking = (np == -1 && linking_vec != nullptr) ? data->*n_linking_vec : -1;

   auto* svec = new DistributedVector<double>(data->*n_vec, nlinking, commWrkrs);

   assert(svec->first);
   double* elems = dynamic_cast<SimpleVector<double>&>(*svec->first).elements();
   double* elems_link = (nlinking != -1) ? dynamic_cast<SimpleVector<double>&>(*svec->last).elements() : nullptr;

   (data->*vec)(data->user_data, data->id, elems, data->*n_vec);

   // at root and with linking constraints?
   if (nlinking != -1) {
      assert(n_linking_vec);
      assert(linking_vec);
      (data->*linking_vec)(data->user_data, data->id, elems_link, data->*n_linking_vec);
   }

   for (const DistributedTree* child_tree : children) {
      std::shared_ptr<DistributedVector<double>> child_vec{dynamic_cast<const DistributedTreeCallbacks*>(child_tree)->createVector(n_vec, vec, n_linking_vec, linking_vec)};
      svec->AddChild(child_vec);
   }

   return svec;
}

DistributedVector<double>* DistributedTreeCallbacks::createc() const {
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fc;

   return createVector(n_func, func, nullptr, nullptr);
}

DistributedVector<double>* DistributedTreeCallbacks::createxlow() const {
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fxlow;

   return createVector(n_func, func, nullptr, nullptr);
}

DistributedVector<double>* DistributedTreeCallbacks::createixlow() const {
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fixlow;

   return createVector(n_func, func, nullptr, nullptr);
}

DistributedVector<double>* DistributedTreeCallbacks::createxupp() const {
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fxupp;

   return createVector(n_func, func, nullptr, nullptr);
}

DistributedVector<double>* DistributedTreeCallbacks::createixupp() const {
   DATA_INT n_func = &InputNode::n;
   DATA_VEC func = &InputNode::fixupp;

   return createVector(n_func, func, nullptr, nullptr);
}

DistributedVector<double>* DistributedTreeCallbacks::createb() const {
   DATA_INT n_func = &InputNode::my;
   DATA_VEC func = &InputNode::fb;

   DATA_INT n_func_link = &InputNode::myl;
   DATA_VEC func_link = &InputNode::fbl;

   return createVector(n_func, func, n_func_link, func_link);
}

DistributedVector<double>* DistributedTreeCallbacks::createclow() const {
   DATA_INT n_func = &InputNode::mz;
   DATA_VEC func = &InputNode::fclow;

   DATA_INT n_func_link = &InputNode::mzl;
   DATA_VEC func_link = &InputNode::fdllow;

   return createVector(n_func, func, n_func_link, func_link);
}

DistributedVector<double>* DistributedTreeCallbacks::createiclow() const {
   DATA_VEC func = &InputNode::ficlow;
   DATA_INT n_func = &InputNode::mz;

   DATA_VEC func_link = &InputNode::fidllow;
   DATA_INT n_func_link = &InputNode::mzl;

   return createVector(n_func, func, n_func_link, func_link);
}

DistributedVector<double>* DistributedTreeCallbacks::createcupp() const {
   DATA_INT n_func = &InputNode::mz;
   DATA_VEC func = &InputNode::fcupp;

   DATA_INT n_func_link = &InputNode::mzl;
   DATA_VEC func_link = &InputNode::fdlupp;

   return createVector(n_func, func, n_func_link, func_link);
}

DistributedVector<double>* DistributedTreeCallbacks::createicupp() const {
   DATA_INT n_func = &InputNode::mz;
   DATA_VEC func = &InputNode::ficupp;

   DATA_INT n_func_link = &InputNode::mzl;
   DATA_VEC func_link = &InputNode::fidlupp;

   return createVector(n_func, func, n_func_link, func_link);
   assert(!is_hierarchical_root || (false && "cannot be used with hierarchical data"));
}

DistributedTree* DistributedTreeCallbacks::shaveDenseBorder(int nx_to_shave, int myl_to_shave, int mzl_to_shave) {
   assertTreeStructureCorrect();
   if (PIPS_MPIgetRank() == 0 && !pipsipmpp_options::get_bool_parameter("SILENT"))
      std::cout << "Trimming " << nx_to_shave << " vars, " << myl_to_shave << " dense equalities, and " << mzl_to_shave
                << " inequalities for the border\n";

   auto* top_layer = new DistributedTreeCallbacks();

   /* sTree members */
   top_layer->commWrkrs = commWrkrs;
   top_layer->myProcs = myProcs;

   /* this must happen at MPI_COMM_WORLD level */
   assert(rankMe == PIPS_MPIgetRank());
   assert(numProcs == PIPS_MPIgetSize());

   top_layer->N = N;
   this->N -= nx_to_shave;

   top_layer->MY = MY;
   top_layer->MYL = MYL;

   top_layer->MZ = MZ;
   top_layer->MZL = MZL;
   top_layer->np = -1;

   assert(IPMIterExecTIME == -1);
   top_layer->children.push_back(this);

   // TODO: not sure about the ressources monitors..: resMon, iterMon..
   top_layer->is_hierarchical_root = true;

   /* DistributedTreeCallbacks members */
   top_layer->nx_active = nx_to_shave;
   this->nx_active -= nx_to_shave;

   /* we move all of A0 to the border layer - else we have structural rank deficiency in the schur complement later */
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL_MOVE_A0_TO_DENSE_LAYER")) {
      top_layer->my_active = this->my_active;
      this->MY -= this->my_active;
      this->my_active = -1;
      this->was_a0_moved_to_border = true;
   } else {
      top_layer->my_active = -1;
   }
   top_layer->mz_active = -1;

   top_layer->myl_active = myl_to_shave;
   this->adjust_active_size_by(&DistributedTreeCallbacks::myl_active, &DistributedTree::MYL, -myl_to_shave);

   top_layer->mzl_active = mzl_to_shave;
   this->adjust_active_size_by(&DistributedTreeCallbacks::mzl_active, &DistributedTree::MZL, -mzl_to_shave);

   for (auto& child : children)
      dynamic_cast<DistributedTreeCallbacks*>(child)->np = nx_active;

   assert(myl_active >= 0);
   assert(mzl_active >= 0);
   top_layer->assertTreeStructureCorrect();
   return top_layer;
}

unsigned int
DistributedTreeCallbacks::getMapChildrenToNthRootSubTrees(int& take_nth_root, std::vector<unsigned int>& map_child_to_sub_tree, unsigned int n_children,
      unsigned int n_procs, const std::vector<unsigned int>& child_procs) {
   assert(take_nth_root > 1);
   assert(n_procs <= n_children);
   assert(child_procs.size() == n_children);

   // old way - results in a lot of wait at Allreduces..
//   const unsigned int n_new_roots = std::round( std::sqrt( n_children ) );
   unsigned int n_new_roots = 0;
   while (n_new_roots <= 1 && take_nth_root > 1) {
      n_new_roots = std::floor(std::pow(n_children, 1.0 / take_nth_root));

      if (n_procs <= n_new_roots && ((n_new_roots % n_procs) != 0)) {
         while (n_new_roots % n_procs != 0)
            --n_new_roots;
      }
      else if (n_procs > n_new_roots && n_procs % n_new_roots != 0) {
         while (n_procs % n_new_roots != 0)
            --n_new_roots;
      }

      if (n_new_roots <= 1) {
         if (PIPS_MPIgetRank(commWrkrs) == 0)
            std::cout << "Too many layers for hierarchical split specified - the number of blocks and MPI processes does not allow for "
                      << take_nth_root << " layers - decreasing amount of layers\n";
         take_nth_root -= 1;
      }
   }

   if (take_nth_root <= 1)
      return n_new_roots;

   assert(n_new_roots > 1);

   /* now map new roots to procs or procs to new_roots */
   if (n_procs >= n_new_roots) {
      std::vector<unsigned int> map_proc_to_sub_tree;
      mapChildrenToNSubTrees(map_proc_to_sub_tree, n_procs, n_new_roots);

      map_child_to_sub_tree.clear();
      map_child_to_sub_tree.resize(n_children);

      for (size_t i = 0; i < child_procs.size(); ++i)
         map_child_to_sub_tree[i] = map_proc_to_sub_tree[child_procs[i]];
   }
   else {
      std::vector<unsigned int> map_sub_tree_to_proc;
      mapChildrenToNSubTrees(map_sub_tree_to_proc, n_new_roots, n_procs);

      mapChildrenToNSubTrees(map_child_to_sub_tree, n_children, n_new_roots);

      for (size_t i = 0; i < child_procs.size(); ++i) {
         if (child_procs[i] != map_sub_tree_to_proc[map_child_to_sub_tree[i]]) {
            assert(child_procs[i] == map_sub_tree_to_proc[map_child_to_sub_tree[i] + 1]);
            map_child_to_sub_tree[i] += 1;
         }
      }
   }

   assert(std::is_sorted(map_child_to_sub_tree.begin(), map_child_to_sub_tree.end()));
   return n_new_roots;
}

void DistributedTreeCallbacks::createSubcommunicatorsAndChildren(int& take_nth_root, std::vector<unsigned int>& map_child_to_sub_tree) {
   assert(children.size() > 1);
   assert(myProcs.size() == getNDistinctValues(myProcs));

   std::vector<unsigned int> my_child_procs;
   my_child_procs.reserve(children.size());

   assert(children.front()->myProcs.size() == 1);
   my_child_procs.push_back(0);
   unsigned int proc = 0;
   for (size_t i = 1; i < children.size(); ++i) {
      assert(children[i]->myProcs.size() == 1);
      if (children[i - 1]->myProcs != children[i]->myProcs)
         ++proc;

      my_child_procs.push_back(proc);
   }

   const unsigned int n_new_roots = getMapChildrenToNthRootSubTrees(take_nth_root, map_child_to_sub_tree, children.size(), myProcs.size(),
         my_child_procs);

   if (n_new_roots <= 1) {
      assert(take_nth_root == 1);
      return;
   }
   assert(map_child_to_sub_tree.size() == children.size());

   /* create new sub-roots */
   std::vector<DistributedTreeCallbacks*> new_leafs(n_new_roots);
   for (auto& leaf : new_leafs) {
      leaf = new DistributedTreeCallbacks();
      leaf->setHierarchicalInnerLeaf();

      leaf->sub_root = new DistributedTreeCallbacks();
      leaf->sub_root->setHierarchicalInnerLeaf();
   }

   /* determine processes for each subtree and assign children */
   for (size_t child = 0; child < children.size(); ++child) {
      const auto& child_procs = children[child]->myProcs;

      assert(!child_procs.empty());

      const unsigned int assigned_sub_root_for_child = map_child_to_sub_tree[child];
      assert(assigned_sub_root_for_child < new_leafs.size());
      DistributedTreeCallbacks* assigned_leaf = new_leafs[assigned_sub_root_for_child];

      for (int process : child_procs) {
         // assuming sorted..
         if (assigned_leaf->myProcs.empty() || assigned_leaf->myProcs.back() != process) {
            if (!assigned_leaf->myProcs.empty()) {
               assert(process > assigned_leaf->myProcs.back());
               assert(!isInVector(process, assigned_leaf->myProcs));
            }

            assigned_leaf->myProcs.push_back(process);
            assigned_leaf->sub_root->myProcs.push_back(process);
         }
      }

      dynamic_cast<DistributedTreeCallbacks*>(assigned_leaf->sub_root)->addChild(dynamic_cast<DistributedTreeCallbacks*>(children[child]));
   }

   /* create all sub-communicators */
   for (auto new_leaf : new_leafs) {
      new_leaf->commWrkrs = PIPS_MPIcreateGroupFromRanks(new_leaf->myProcs, commWrkrs);

      if (!isInVector(rankMe, new_leaf->myProcs))
         assert(new_leaf->commWrkrs == MPI_COMM_NULL);

      new_leaf->sub_root->commWrkrs = new_leaf->commWrkrs;
      new_leaf->sub_root->myProcs = new_leaf->myProcs;
   }

#ifndef NDEBUG
   for (size_t child = 0; child < children.size(); ++child) {
      if (isInVector(rankMe, children[child]->myProcs)) {
         auto child_new_root = new_leafs[map_child_to_sub_tree[child]];
         assert(child_new_root->commWrkrs != MPI_COMM_NULL);
      }
   }
#endif

   /* add sub_roots as this new children */
   children.clear();
   children.insert(children.begin(), new_leafs.begin(), new_leafs.end());

   if (!is_hierarchical_inner_leaf)
      is_hierarchical_inner_root = true;
}

void
DistributedTreeCallbacks::countTwoLinksForChildTrees(const std::vector<int>& two_links_start_in_child_A, const std::vector<int>& two_links_start_in_child_C,
      std::vector<unsigned int>& two_links_children_eq, std::vector<unsigned int>& two_links_children_ineq, unsigned int& two_links_root_eq,
      unsigned int& two_links_root_ineq) const {
   const unsigned int n_children = children.size();
#ifndef NDEBUG
   const unsigned int n_leafs = map_node_sub_root.size();
#endif
   assert(n_leafs == two_links_start_in_child_A.size());
   assert(n_leafs == two_links_start_in_child_C.size());

   two_links_children_eq.clear();
   two_links_children_eq.resize(n_children);
   two_links_children_ineq.clear();
   two_links_children_ineq.resize(n_children);

   /* count two links for all new sub-communicators */
   unsigned int leaf = 0;
   for (unsigned int i = 0; i < n_children; ++i) {
      const auto& sub_root = dynamic_cast<const DistributedTreeCallbacks&>( *children[i]->sub_root );

      for (size_t j = 0; j < sub_root.children.size(); ++j) {
         assert(leaf < n_leafs);

         /* the two links of each last child will stay in the root node */
         if (j == sub_root.children.size() - 1) {
            two_links_root_eq += two_links_start_in_child_A[leaf];
            two_links_root_ineq += two_links_start_in_child_C[leaf];
         }
         else {
            two_links_children_eq[i] += two_links_start_in_child_A[leaf];
            two_links_children_ineq[i] += two_links_start_in_child_C[leaf];
         }

         ++leaf;
      }
   }
   assert(leaf == n_leafs);
}

void DistributedTreeCallbacks::adjust_active_size_by(int DistributedTreeCallbacks::* active_size_to_adjust, long long DistributedTree::* glob_size_to_adjust, int adjustment)
{
   assert(isInVector(rankMe, myProcs));

   this->*active_size_to_adjust += adjustment;
   this->*glob_size_to_adjust += adjustment;

   assert(this->*active_size_to_adjust >= 0);
   assert(this->*glob_size_to_adjust >= 0);

   for (DistributedTree* child_ : children) {
      auto* child = dynamic_cast<DistributedTreeCallbacks*>(child_);
      assert(child->children.empty());

      if (isInVector(rankMe, child->myProcs)) {
         child->*active_size_to_adjust += adjustment;
         child->*glob_size_to_adjust += adjustment;
      }
      assert(child->*active_size_to_adjust >= 0);
      assert(child->*glob_size_to_adjust >= 0);
   }
}

std::pair<int, int> DistributedTreeCallbacks::adjustSizesAfterSplit(const std::vector<unsigned int>& two_links_children_eq,
      const std::vector<unsigned int>& two_links_children_ineq) {
   assert(two_links_children_eq.size() == two_links_children_ineq.size());

   unsigned int my_two_links_eq{0};
   unsigned int my_two_links_ineq{0};
   unsigned int sum_two_links_children_eq{0};
   unsigned int sum_two_links_children_ineq{0};

   for (unsigned int i = 0; i < children.size(); ++i) {
      const auto& inner_leaf = dynamic_cast<const DistributedTreeCallbacks&>(*children[i]);
      if (isInVector(rankMe, inner_leaf.myProcs)) {
         my_two_links_eq += two_links_children_eq[i];
         my_two_links_ineq += two_links_children_ineq[i];
      }
      sum_two_links_children_eq += two_links_children_eq[i];
      sum_two_links_children_ineq += two_links_children_ineq[i];
   }
   const unsigned int not_my_two_links_eq = sum_two_links_children_eq - my_two_links_eq;
   const unsigned int not_my_two_links_ineq = sum_two_links_children_ineq - my_two_links_ineq;

   assert(sum_two_links_children_eq <= unsigned(myl_active));
   assert(sum_two_links_children_ineq <= unsigned(mzl_active));
   assert(not_my_two_links_eq <= MYL);
   assert(not_my_two_links_ineq <= MZL);

   myl_active -= static_cast<int>(sum_two_links_children_eq);
   mzl_active -= static_cast<int>(sum_two_links_children_ineq);
   MYL -= not_my_two_links_eq;
   MZL -= not_my_two_links_ineq;

   /* recompute sizes for the children */
   for (size_t i = 0; i < children.size(); ++i) {
      auto& inner_leaf = dynamic_cast<DistributedTreeCallbacks&>(*children[i]);
      auto& sub_root = dynamic_cast<DistributedTreeCallbacks&>(*inner_leaf.sub_root);
      assert(inner_leaf.is_hierarchical_inner_leaf);

      if (isInVector(rankMe, inner_leaf.myProcs)) {
         assert(isInVector(rankMe, sub_root.myProcs));

         inner_leaf.MYL = myl_active + two_links_children_eq[i];
         inner_leaf.myl_active = myl_active;
         if (MYL >= 0) {
            sub_root.adjust_active_size_by(&DistributedTreeCallbacks::myl_active, &DistributedTree::MYL,
               -myl_active - static_cast<int>(sum_two_links_children_eq - two_links_children_eq[i]));
         }

         inner_leaf.MZL = mzl_active + two_links_children_ineq[i];
         inner_leaf.mzl_active = mzl_active;
         if (MZL >= 0) {
            sub_root.adjust_active_size_by(&DistributedTreeCallbacks::mzl_active, &DistributedTree::MZL,
               -mzl_active - static_cast<int>(sum_two_links_children_ineq - two_links_children_ineq[i]));
         }

         inner_leaf.np = nx_active;

         inner_leaf.N = sub_root.N;
         inner_leaf.nx_active = static_cast<int>(sub_root.N);

         inner_leaf.MY = sub_root.MY;
         inner_leaf.my_active = static_cast<int>(sub_root.MY);

         inner_leaf.MZ = sub_root.MZ;
         inner_leaf.mz_active = static_cast<int>(sub_root.MZ);
      }
      else {
         assert(!isInVector(rankMe, sub_root.myProcs));

         // we need to store the info about how many two links we pushed to that child (and forgot about)
         inner_leaf.MYL = myl_active + two_links_children_eq[i];
         inner_leaf.MZL = mzl_active + two_links_children_ineq[i];
         sub_root.MYL = two_links_children_eq[i];
         sub_root.MZL = two_links_children_ineq[i];

         inner_leaf.N = inner_leaf.MY = inner_leaf.MZ = 0;
         inner_leaf.nx_active = inner_leaf.my_active = inner_leaf.mz_active = inner_leaf.myl_active = inner_leaf.mzl_active = 0;

#ifndef NDEBUG
         assert(sub_root.N == 0 && sub_root.MY == 0 && sub_root.MZ == 0);
         assert(
               sub_root.nx_active == 0 && sub_root.my_active == 0 && sub_root.mz_active == 0 && sub_root.myl_active == 0 && sub_root.mzl_active == 0);

         for (const auto& child_ : sub_root.children) {
            const auto& child = dynamic_cast<const DistributedTreeCallbacks*>(child_);
            assert(child->N == 0 && child->MY == 0 && child->MZ == 0 && child->MYL == 0 && child->MZL == 0);
            assert(child->nx_active == 0 && child->my_active == 0 && child->mz_active == 0 && child->myl_active == 0 && child->mzl_active == 0);
         }
#endif
      }
   }
   return std::make_pair(not_my_two_links_eq, not_my_two_links_ineq);
}

std::pair<int, int> DistributedTreeCallbacks::splitTree(int n_layers, DistributedQP* data_to_split) {
   if (n_layers == 1 || commWrkrs == MPI_COMM_NULL)
      return std::make_pair(0, 0);

   const std::vector<int>& twoLinksStartBlockA = data_to_split->getTwoLinksStartBlockA();
   const std::vector<int>& twoLinksStartBlockC = data_to_split->getTwoLinksStartBlockC();

   assert(!is_hierarchical_root);

#ifndef NDEBUG
   const size_t n_old_leafs = children.size();
#endif

   createSubcommunicatorsAndChildren(n_layers, map_node_sub_root);

   if (n_layers == 1) {
      if (PIPS_MPIgetRank(commWrkrs) == 0)
         std::cout << "No split applied!\n";
      return std::make_pair(0, 0);
   }

   assert(map_node_sub_root.size() == n_old_leafs);

   std::vector<unsigned int> two_links_children_eq, two_links_children_ineq;
   unsigned int two_links_root_eq{0};
   unsigned int two_links_root_ineq{0};

   countTwoLinksForChildTrees(twoLinksStartBlockA, twoLinksStartBlockC, two_links_children_eq, two_links_children_ineq, two_links_root_eq,
         two_links_root_ineq);

   const unsigned int two_links_children_eq_sum = std::accumulate(two_links_children_eq.begin(), two_links_children_eq.end(), unsigned(0));
   const unsigned int two_links_children_ineq_sum = std::accumulate(two_links_children_ineq.begin(), two_links_children_ineq.end(), unsigned(0));

   if (rankMe == 0 && !pipsipmpp_options::get_bool_parameter("SILENT")) {
      std::cout << "Splitting node into " << children.size() << " subroots\n";
      std::cout << "Splitting " << two_links_children_eq_sum + two_links_root_eq << " equality two-links into " << two_links_root_eq << " root and "
                << two_links_children_eq_sum << " child links\n";
      std::cout << "Splitting " << two_links_children_ineq_sum + two_links_root_ineq << " inequality two-links into " << two_links_root_ineq
                << " root and " << two_links_children_ineq_sum << " child links\n";
   }

   std::pair<int, int> deleted_myl_mzl = adjustSizesAfterSplit(two_links_children_eq, two_links_children_ineq);
   data_to_split->splitDataAccordingToTree();
   assertTreeStructureCorrect();

   assert(children.size() == data_to_split->children.size());

   std::pair<int, int> deleted_children(0, 0);
   for (size_t i = 0; i < children.size(); ++i) {
      if (n_layers - 1 > 1) {
         if (children[i]->sub_root->getCommWorkers() != MPI_COMM_NULL) {
            std::pair<int, int> child_deleted = children[i]->sub_root->splitTree(n_layers - 1, data_to_split->children[i]);
            children[i]->MYL -= child_deleted.first;
            children[i]->MZL -= child_deleted.second;

            deleted_children = std::make_pair(deleted_children.first + child_deleted.first, deleted_children.second + child_deleted.second);
            data_to_split->children[i]->splitStringMatricesAccordingToSubtreeStructure();
         }
      }
   }
   MYL -= deleted_children.first;
   MZL -= deleted_children.second;

   assertTreeStructureCorrect();
   data_to_split->recomputeSize();
   return std::make_pair<int, int>(deleted_myl_mzl.first + deleted_children.first, deleted_myl_mzl.second + deleted_children.second);
}

DistributedTree* DistributedTreeCallbacks::switchToHierarchicalTree(DistributedQP*& data_to_split) {
   assert(data_to_split->exploitingLinkStructure());

   const int n_layers = pipsipmpp_options::get_int_parameter("HIERARCHICAL_APPROACH_N_LAYERS");

   assert(!is_hierarchical_root);
   assert(np == -1);
   assertTreeStructureCorrect();

   /* distributed preconditioner must be deactivated */
   assert(!distributedPreconditionerActive());

   if (n_layers >= 1) {
      if (PIPS_MPIgetRank() == 0) {
         std::cout << "Building hierarchical data_to_split\n";
         std::cout << "Adding " << n_layers << " layers to hierarchical data_to_split\n";
      }

      if (n_layers > 1 && pipsipmpp_options::get_bool_parameter("HIERARCHICAL_APPLY_SPLIT")) {
         splitTree(n_layers, data_to_split);
         assertTreeStructureCorrect();
         printProcessTree();

         if (map_node_sub_root.empty()) {
            if (PIPS_MPIgetRank() == 0)
               std::cout << "Not a single split has been applied - hierarchical approach will have no additional sparse layers\n";
            pipsipmpp_options::set_int_parameter("HIERARCHICAL_APPROACH_N_LAYERS", 1);
         }
      }

      const int nx_to_shave = data_to_split->getNGlobalVars();
      const int myl_to_shave = data_to_split->getNGlobalEQConss();
      const int mzl_to_shave = data_to_split->getNGlobalINEQConss();

      assert(nx_to_shave >= 0);
      assert(myl_to_shave >= 0);
      assert(mzl_to_shave >= 0);

      assert(nx_to_shave <= nx_active);
      assert(myl_to_shave <= myl_active);
      assert(mzl_to_shave <= mzl_active);

      auto* top_layer = dynamic_cast<DistributedTreeCallbacks*>(shaveDenseBorder(nx_to_shave, myl_to_shave, mzl_to_shave));
      data_to_split = data_to_split->shaveDenseBorder(top_layer);

      if (PIPS_MPIgetRank() == 0)
         std::cout << "Hierarchical data_to_split built\n";

      return top_layer;
   }
   else
      return this;
}

std::vector<MPI_Comm> DistributedTreeCallbacks::getChildComms() const {
   std::vector<MPI_Comm> comms(children.size(), MPI_COMM_NULL);

   for (unsigned int i = 0; i < children.size(); ++i)
      comms[i] = children[i]->commWrkrs;

   return comms;
}
