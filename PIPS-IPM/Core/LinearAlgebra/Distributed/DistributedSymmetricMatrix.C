#include "DistributedSymmetricMatrix.h"
#include "DistributedVector.h"
#include "DoubleMatrixTypes.h"
#include "BorderedSymmetricMatrix.h"
#include "StripMatrix.h"
#include <cassert>

DistributedSymmetricMatrix::DistributedSymmetricMatrix(SymmetricMatrix* diag_, SparseMatrix* border_, MPI_Comm mpiComm_) : diag{diag_}, border{border_}, mpiComm{mpiComm_},
      iAmDistrib(PIPS_MPIgetDistributed(mpiComm)) {
   assert(diag);
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
DistributedSymmetricMatrix::DistributedSymmetricMatrix(long long global_n, int local_n, int local_nnz, MPI_Comm mpiComm_) : n(global_n), mpiComm(mpiComm_),
      iAmDistrib(PIPS_MPIgetDistributed(mpiComm)) {
   diag = new SparseSymmetricMatrix(local_n, local_nnz);
   // the cross Hessian is nullptr for the root node; it may be also nullptr for
   // children in the case when the Hessian does not have cross terms and
   // the children are created with this constructor. The border will be
   // set up to correct sizes later for this case.
}

void DistributedSymmetricMatrix::AddChild(DistributedSymmetricMatrix* child) {
   child->parent = this;
   assert(this->border == nullptr);

   if (!child->border)
      child->border = new SparseMatrix(static_cast<int>(child->diag->size()), static_cast<int>(this->diag->size()), 0);

   children.push_back(child);
}

DistributedSymmetricMatrix::~DistributedSymmetricMatrix() {
   for (auto & it : children)
      delete it;

   delete diag;
   delete border;
}

SymmetricMatrix* DistributedSymmetricMatrix::clone() const {
   SymmetricMatrix* diag_clone = diag->clone();
   SparseMatrix* border_clone = border ? dynamic_cast<SparseMatrix*>(border->cloneFull()) : nullptr;

   auto* clone = new DistributedSymmetricMatrix(diag_clone, border_clone, mpiComm);

   for (auto it : children) {
      auto* child = dynamic_cast<DistributedSymmetricMatrix*>(it->clone());
      clone->AddChild(child);
      clone->n += child->n;
   }

   assert(size() == clone->size());
   return clone;
}

void DistributedSymmetricMatrix::recomputeSize() {
   assert(diag);
   n = 0;
   for (auto& child : children) {
      if ( !child->is_a(kStochSymDummyMatrix)) {
         child->recomputeSize();
      }

      n += child->size();
   }

   if (diag->is_a(kStochSymMatrix) && !diag->is_a(kStochSymDummyMatrix)){
      dynamic_cast<DistributedSymmetricMatrix*>(diag)->recomputeSize();
   }

   n += diag->size();
}

int DistributedSymmetricMatrix::is_a(int type) const {
   return type == kStochSymMatrix || type == kSymMatrix;
}

void DistributedSymmetricMatrix::getSize(long long& m_, long long& n_) const {
   m_ = n;
   n_ = n;
}

void DistributedSymmetricMatrix::getSize(int& m_, int& n_) const {
   m_ = static_cast<int>(n);
   n_ = static_cast<int>(n);
}

long long DistributedSymmetricMatrix::size() const {
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
void DistributedSymmetricMatrix::mult(double beta, Vector<double>& y_, double alpha, const Vector<double>& x_) const {
   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_);

   if (0.0 == alpha) {
      y.first->scale(beta);
      return;
   }

   if (!parent) {
      if (PIPS_MPIiAmSpecial(iAmDistrib, mpiComm))
         diag->mult(beta, *y.first, alpha, *x.first);
      else
         y.first->setToZero();
   }
   else
      diag->mult(beta, *y.first, alpha, *x.first);

   // y0 = y0 + alpha * border^T * xi
   if (border) {
      assert(parent);
      assert(x.parent);
      assert(y.parent);

      border->transMult(1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.first);
      border->mult(1.0, *y.first, alpha, *x.getLinkingVecNotHierarchicalTop());
   }

   assert(y.children.size() == children.size());
   assert(x.children.size() == children.size());

   // recursively multiply the children
   for (size_t it = 0; it < children.size(); it++)
      children[it]->mult(beta, *(y.children[it]), alpha, *(x.children[it]));

   if (iAmDistrib && !parent)
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector<double>*>(y.first)->elements(), y.first->length(), mpiComm);
}

/** y = beta * y + alpha * this^T * x */
void DistributedSymmetricMatrix::transMult(double beta, Vector<double>& y_, double alpha, const Vector<double>& x_) const {
   // We are symmetric, this^T = this, therefore call 'mult' method
   this->mult(beta, y_, alpha, x_);
}

/** the magnitude of the element in this matrix with largest absolute value.
   */
double DistributedSymmetricMatrix::inf_norm() const {
   double maxNorm = 0.0;

   for (auto it : children)
      maxNorm = std::max(maxNorm, it->inf_norm());

   if (iAmDistrib)
      PIPS_MPIgetMaxInPlace(maxNorm, mpiComm);

   maxNorm = std::max(maxNorm, diag->inf_norm());
   if (border)
      maxNorm = std::max(maxNorm, border->inf_norm());
   return maxNorm;
}

double DistributedSymmetricMatrix::abminnormNonZero(double tol) const {
   double min = std::numeric_limits<double>::infinity();

   for (auto it : children)
      min = std::min(min, it->abminnormNonZero(tol));

   if (iAmDistrib)
      PIPS_MPIgetMinInPlace(min, mpiComm);

   min = std::min(min, diag->abminnormNonZero(tol));
   if (border)
      min = std::min(min, diag->abminnormNonZero(tol));
   return min;
}

void DistributedSymmetricMatrix::writeToStreamDense(std::ostream& out) const {
   const int rank = PIPS_MPIgetRank(mpiComm);
   const int world_size = PIPS_MPIgetSize(mpiComm);

   int m_loc, n_loc;
   int offset = 0;
   std::stringstream sout;
   MPI_Status status;
   int l;

   /* this is at the root node - thus there is no border */
   assert(this->border == nullptr);

   if (iAmDistrib)
      MPI_Barrier(mpiComm);

   if (iAmDistrib && rank > 0)  // receive offset from previous process
      MPI_Recv(&offset, 1, MPI_INT, (rank - 1), 0, mpiComm, MPI_STATUS_IGNORE);
   else  //  !iAmDistrib || (iAmDistrib && rank == 0)
      this->diag->writeToStreamDense(sout);

   for (auto it : children) {
      it->writeToStreamDenseChild(sout, offset);
      it->diag->getSize(m_loc, n_loc);
      offset += n_loc;
   }

   if (iAmDistrib && rank > 0) {
      std::string str = sout.str();
      // send string to rank ZERO to print it there:
      MPI_Ssend(str.c_str(), static_cast<int>(str.length()), MPI_CHAR, 0, rank, mpiComm);
      // send offset to next process:
      if (rank < world_size - 1)
         MPI_Ssend(&offset, 1, MPI_INT, rank + 1, 0, mpiComm);
   }
   else if (!iAmDistrib)
      out << sout.str();
   else if (iAmDistrib && rank == 0) {
      out << sout.str();
      MPI_Ssend(&offset, 1, MPI_INT, rank + 1, 0, mpiComm);

      for (int p = 1; p < world_size; p++) {
         MPI_Probe(p, p, mpiComm, &status);
         MPI_Get_count(&status, MPI_CHAR, &l);
         char* buf = new char[l];
         MPI_Recv(buf, l, MPI_CHAR, p, p, mpiComm, &status);
         std::string rowPartFromP(buf, l);
         out << rowPartFromP;
         delete[] buf;
      }
   }

   if (iAmDistrib)
      MPI_Barrier(mpiComm);
   std::cout << " done\n";
}

void DistributedSymmetricMatrix::writeToStreamDenseChild(std::stringstream& out, int offset) const {
   if (diag->is_a(kSparseSymMatrix)) {
      for (int r = 0; r < diag->size(); r++) {
         if (border)
            border->writeToStreamDenseRow(out, r);

         for (int i = 0; i < offset; i++)
            out << '\t';

         dynamic_cast<const SparseSymmetricMatrix&>(*diag).writeToStreamDenseRow(out, r);
         out << "\n";
      }
   }
   else
      dynamic_cast<const DistributedSymmetricMatrix&>(*diag).writeToStreamDense(out);
}


void DistributedSymmetricMatrix::getDiagonal(Vector<double>& vec_) const {
   auto& vec = dynamic_cast<DistributedVector<double>&>(vec_);
   assert(children.size() == vec.children.size());

   diag->getDiagonal(*vec.first);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->getDiagonal(*vec.children[it]);
}

void DistributedSymmetricMatrix::setToDiagonal(const Vector<double>& vec_) {
   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_);
   assert(children.size() == vec.children.size());

   diag->setToDiagonal(*vec.first);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->setToDiagonal(*vec.children[it]);
}

void DistributedSymmetricMatrix::atPutDiagonal(int idiag, const Vector<double>& v_) {
   const auto& v = dynamic_cast<const DistributedVector<double>&>(v_);

   //check the tree compatibility
   int nChildren = static_cast<int>(children.size());
   assert(v.children.size() - nChildren == 0);

   //check the node size compatibility
   assert(this->diag->size() == v.first->length());

   diag->atPutDiagonal(idiag, *v.first);

   for (int it = 0; it < nChildren; it++)
      children[it]->atPutDiagonal(idiag, *v.children[it]);
}

void DistributedSymmetricMatrix::atAddDiagonal(int idiag, const Vector<double>& v_) {
   const auto& v = dynamic_cast<const DistributedVector<double>&>(v_);

   //check the tree compatibility
   int nChildren = static_cast<int>(children.size());
   assert(v.children.size() - nChildren == 0);

   //check the node size compatibility
   assert(this->diag->size() == v.first->length());

   diag->atAddDiagonal(idiag, *v.first);

   for (int it = 0; it < nChildren; it++)
      children[it]->atAddDiagonal(idiag, *v.children[it]);
}
void DistributedSymmetricMatrix::fromGetDiagonal(int, Vector<double>& x_) const {
   assert("The value of the parameter is not supported!");

   auto& x = dynamic_cast<DistributedVector<double>&>(x_);
   assert(x.children.size() == children.size());

   diag->getDiagonal(*x.first);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->getDiagonal(*x.children[it]);
}

void DistributedSymmetricMatrix::symmetricScale(const Vector<double>& vec_) {
   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_);
   assert(children.size() == vec.children.size());

   diag->symmetricScale(*vec.first);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->symmetricScale(*vec.children[it]);
}

void DistributedSymmetricMatrix::columnScale(const Vector<double>& vec_) {
   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_);
   assert(children.size() == vec.children.size());

   diag->columnScale(*vec.first);

   if (border)
      border->columnScale(*vec.getLinkingVecNotHierarchicalTop());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->columnScale(*vec.children[it]);
}

void DistributedSymmetricMatrix::rowScale(const Vector<double>& vec_) {
   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_);
   assert(children.size() == vec.children.size());

   diag->rowScale(*vec.first);
   if (border)
      border->rowScale(*vec.first);

   for (size_t it = 0; it < children.size(); it++)
      children[it]->rowScale(*vec.children[it]);
}


void DistributedSymmetricMatrix::scalarMult(double num) {
   diag->scalarMult(num);
   if (border)
      border->scalarMult(num);
   for (auto & it : children)
      it->scalarMult(num);
}

void DistributedSymmetricMatrix::deleteEmptyRowsCols(const Vector<int>& nnzVec, const Vector<int>* linkParent) {
   const auto& nnzVecStoch = dynamic_cast<const DistributedVector<int>&>(nnzVec);
   assert(children.size() == nnzVecStoch.children.size());

   const auto* vec = dynamic_cast<const SimpleVector<int>*>( nnzVecStoch.first );
   assert(vec);

   const long long n_old = n;
   for (size_t it = 0; it < children.size(); it++)
      children[it]->deleteEmptyRowsCols(*nnzVecStoch.children[it], vec);

   if (linkParent != nullptr) {
      assert(border);
      assert(children.empty());
      // adapt border
      assert(dynamic_cast<const SimpleVector<int>*>(linkParent));
      border->deleteEmptyRowsCols(*vec, *linkParent);
   }
   else
      assert(!border);

   assert(diag->is_a(kSparseSymMatrix));
   dynamic_cast<SparseSymmetricMatrix&>(*diag).deleteEmptyRowsCols(*vec);

   const int nnzs = vec->getNnzs();
   const int length = vec->length();

   assert(nnzs <= length);
   n -= (length - nnzs);

   assert(n <= n_old);
   assert(n >= 0);

   const long long n_changes = n_old - n;

   if (linkParent != nullptr) {
      assert(parent != nullptr);
      parent->n -= n_changes;
   }
}


int StochSymDummyMatrix::is_a(int type) const {
   return type == kStochSymDummyMatrix || type == kStochSymMatrix;
}

void DistributedSymmetricMatrix::splitMatrix(const std::vector<unsigned int>& map_blocks_children, const std::vector<MPI_Comm>& child_comms) {
   const unsigned int n_curr_children = children.size();
   assert(n_curr_children == map_blocks_children.size());

   const unsigned int n_new_children = getNDistinctValues(map_blocks_children);
   std::vector<DistributedSymmetricMatrix*> new_children(n_new_children);

   unsigned int begin_curr_child_blocks{0};
   unsigned int end_curr_child_blocks{0};
   for (unsigned int i = 0; i < n_new_children; ++i) {
      while (end_curr_child_blocks != (n_curr_children - 1) &&
             map_blocks_children[end_curr_child_blocks] == map_blocks_children[end_curr_child_blocks + 1])
         ++end_curr_child_blocks;

      const unsigned int n_blocks_for_child = end_curr_child_blocks - begin_curr_child_blocks + 1;

      DistributedSymmetricMatrix* diag_new = (child_comms[i] == MPI_COMM_NULL) ? nullptr : new DistributedSymmetricMatrix(0, 0, 0, child_comms[i]);

      /* shave off empty two link part from respective children and add them to the new root/remove them from the old root */
      for (unsigned int j = 0; j < n_blocks_for_child; ++j) {
         DistributedSymmetricMatrix* child = children.front();
         children.erase(children.begin());

         if (child_comms[i] == MPI_COMM_NULL)
            assert(child->mpiComm == MPI_COMM_NULL);

         if (diag_new)
            diag_new->AddChild(child);
//         else
//            delete child;
      }
      if (diag_new)
         diag_new->recomputeSize();

      /* create child holding the new Bmat and it's Blmat part */
      new_children[i] = (child_comms[i] != MPI_COMM_NULL) ? new DistributedSymmetricMatrix(diag_new, nullptr, child_comms[i]) : new StochSymDummyMatrix();
      if (child_comms[i] != MPI_COMM_NULL)
         dynamic_cast<DistributedSymmetricMatrix*>(new_children[i]->diag)->parent = new_children[i];

      ++end_curr_child_blocks;
      begin_curr_child_blocks = end_curr_child_blocks;
   }
   assert(children.empty());

   /* exchange children and recompute sizes */
   children.insert(children.end(), new_children.begin(), new_children.end());

   for (auto& child : children) {
      child->recomputeSize();
      child->parent = this;
   }
   this->recomputeSize();
}

BorderedSymmetricMatrix* DistributedSymmetricMatrix::raiseBorder(int n_vars) {
   assert(parent == nullptr);
   assert(border == nullptr);

   StripMatrix* border_vertical = shaveBorder(n_vars);

   BorderedSymmetricMatrix* const border_layer = new BorderedSymmetricMatrix(this, border_vertical, new SparseSymmetricMatrix(n_vars, 0, false), mpiComm);

   DistributedSymmetricMatrix* me = this;
   IotrAddRef(&me);

   assert(n >= 0);

   return border_layer;
}

StripMatrix* DistributedSymmetricMatrix::shaveBorder(int n_vars) {
   assert(diag->is_a(kSparseSymMatrix));

   SparseMatrix* const border_top_left = parent ? new SparseMatrix(0, n_vars, 0) : dynamic_cast<SparseSymmetricMatrix*>(diag)->shaveSymLeftBottom(
         n_vars);
   auto* const border_vertical = new StripMatrix(true, border_top_left, nullptr, mpiComm);

   for (auto & it : children)
      border_vertical->addChild(it->shaveBorder2(n_vars));

   n -= n_vars;
   border_vertical->recomputeNonzeros();
   return border_vertical;
}

StripMatrix* DistributedSymmetricMatrix::shaveBorder2(int n_vars) {
   n -= n_vars;

   if (border) {
      SparseMatrix* const border_block = border->shaveLeft(n_vars);
      return new StripMatrix(true, border_block, nullptr, mpiComm);
   }
   else {
      assert(children.empty());
      assert(diag->is_a(kStochSymMatrix));

      auto& diags = dynamic_cast<DistributedSymmetricMatrix&>(*diag);
      return new StripMatrix(true, diags.shaveBorder(n_vars), nullptr, mpiComm);
   }
}