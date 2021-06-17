#include "DistributedMatrix.h"
#include "DoubleMatrixTypes.h"
#include "DistributedVector.h"
#include "SimpleVector.hpp"
#include "pipsdef.h"
#include "PIPSIPMppOptions.h"
#include <limits>
#include <algorithm>
#include <numeric>
#include <memory>

DistributedMatrix::DistributedMatrix(std::unique_ptr<GeneralMatrix> Amat_in, std::unique_ptr<GeneralMatrix> Bmat_in,
   std::unique_ptr<GeneralMatrix> Blmat_in, MPI_Comm mpiComm_, bool inner_leaf, bool inner_root) : Amat{std::move(Amat_in)},
   Bmat{std::move(Bmat_in)}, Blmat{std::move(Blmat_in)}, mpiComm{mpiComm_}, iAmDistrib{PIPS_MPIgetDistributed(mpiComm)},
   inner_leaf{inner_leaf},
   inner_root{inner_root} {
   assert(Amat);
   assert(Bmat);
   assert(Blmat);

#ifndef NDEBUG
   const auto[mA, nA] = Amat->n_rows_columns();
   const auto[mB, nB] = Bmat->n_rows_columns();
   const long long nBl = Blmat->n_columns();

   if (nA == 0 && mA == 0)
      assert(nBl == nB);
   else {
      assert(mA == mB);
      assert(nB == nBl);
   }
#endif

   recomputeSize();
}

DistributedMatrix::DistributedMatrix(long long global_m, long long global_n, int A_m, int A_n, int A_nnz, int B_m,
   int B_n, int B_nnz, MPI_Comm mpiComm_)
   : m(global_m), n(global_n), mpiComm(mpiComm_), iAmDistrib(PIPS_MPIgetDistributed(mpiComm)) {
   Amat = std::make_unique<SparseMatrix>(A_m, A_n, A_nnz);
   Bmat = std::make_unique<SparseMatrix>(B_m, B_n, B_nnz);
   Blmat = std::make_unique<SparseMatrix>(0, B_n, 0);
}

DistributedMatrix::DistributedMatrix(long long global_m, long long global_n, int A_m, int A_n, int A_nnz, int B_m,
   int B_n, int B_nnz, int Bl_m, int Bl_n,
   int Bl_nnz, MPI_Comm mpiComm_) : m(global_m), n(global_n), mpiComm(mpiComm_),
   iAmDistrib(PIPS_MPIgetDistributed(mpiComm)) {
   Amat = std::make_unique<SparseMatrix>(A_m, A_n, A_nnz);
   Bmat = std::make_unique<SparseMatrix>(B_m, B_n, B_nnz);
   Blmat = std::make_unique<SparseMatrix>(Bl_m, Bl_n, Bl_nnz);
}

bool DistributedMatrix::amatEmpty() const {
   const auto[mA, nA] = Amat->n_rows_columns();
   return mA <= 0 || nA <= 0;
}

bool DistributedMatrix::hasSparseMatrices() const {
   return Amat->is_a(kSparseGenMatrix) && Bmat->is_a(kSparseGenMatrix) && Blmat->is_a(kSparseGenMatrix);
}

std::unique_ptr<GeneralMatrix> DistributedMatrix::cloneFull(bool switchToDynamicStorage) const {
//   auto clone = std::make_unique<DistributedMatrix>(m, n, mpiComm);
   assert(hasSparseMatrices());

   // clone submatrices
   auto Amat_clone = Amat->cloneFull(switchToDynamicStorage);
   auto Bmat_clone = Bmat->cloneFull(switchToDynamicStorage);
   auto Blmat_clone = Blmat->cloneFull(switchToDynamicStorage);

   auto clone = std::make_unique<DistributedMatrix>(std::move(Amat_clone), std::move(Bmat_clone), std::move(Blmat_clone), mpiComm);
   for (const auto& it : children)
   {
      std::shared_ptr<DistributedMatrix> child_clone{dynamic_cast<DistributedMatrix*>(it->cloneFull(switchToDynamicStorage).release())};
      clone->children.push_back(std::move(child_clone));
   }
   return clone;
}

/* creates an empty copy of the matrix with n = 0 for all submatrices and m (cols) as before */
std::unique_ptr<GeneralMatrix> DistributedMatrix::cloneEmptyRows(bool switchToDynamicStorage) const {
   assert(hasSparseMatrices());

   // clone submatrices
   auto Amat_clone = Amat->cloneEmptyRows(switchToDynamicStorage);
   auto Bmat_clone = Bmat->cloneEmptyRows(switchToDynamicStorage);
   auto Blmat_clone = Blmat->cloneEmptyRows(switchToDynamicStorage);

   auto clone = std::make_unique<DistributedMatrix>(std::move(Amat_clone), std::move(Bmat_clone), std::move(Blmat_clone), mpiComm);
   for (const auto& it : children){
      std::shared_ptr<DistributedMatrix> child_clone{dynamic_cast<DistributedMatrix*>(it->cloneEmptyRows(switchToDynamicStorage).release())};
      clone->children.push_back(std::move(child_clone));
   }

   return clone;
}

void DistributedMatrix::AddChild(const std::shared_ptr<DistributedMatrix>& child) {
   children.push_back(child);
}

int DistributedMatrix::is_a(int type) const {
   return type == kDistributedMatrix || type == kGenMatrix;
}

int StochGenDummyMatrix::is_a(int type) const {
   return type == kStochGenDummyMatrix || type == kGenMatrix || type == kDistributedMatrix;
}

std::pair<long long, long long> DistributedMatrix::n_rows_columns() const {
   return {m, n};
}

long long DistributedMatrix::n_rows() const {
   return m;
}

long long DistributedMatrix::n_columns() const {
   return n;
}

void DistributedMatrix::columnScale2(const Vector<double>& vec) {
   const auto& scalevec = dynamic_cast<const DistributedVector<double>&>(vec);
   assert(scalevec.children.empty() && children.empty());

   Bmat->columnScale(*scalevec.first);

   if (!amatEmpty())
      Amat->columnScale(*scalevec.getLinkingVecNotHierarchicalTop());

   Blmat->columnScale(*scalevec.first);
}

void DistributedMatrix::columnScale(const Vector<double>& vec) {
   const auto& scalevec = dynamic_cast<const DistributedVector<double>&>(vec);

   assert(amatEmpty());

   Bmat->columnScale(*scalevec.getLinkingVecNotHierarchicalTop());
   Blmat->columnScale(*scalevec.getLinkingVecNotHierarchicalTop());

   assert(children.size() == scalevec.children.size());
   for (size_t it = 0; it < children.size(); it++)
      children[it]->columnScale2(*(scalevec.children[it]));
}

void DistributedMatrix::rowScale2(const Vector<double>& vec, const Vector<double>* linkingvec) {
   const auto& scalevec = dynamic_cast<const DistributedVector<double>&>(vec);

   assert(scalevec.children.empty() && children.empty());

   if (!amatEmpty())
      Amat->rowScale(*scalevec.first);

   Bmat->rowScale(*scalevec.first);

   if (linkingvec)
      Blmat->rowScale(*linkingvec);
}

void DistributedMatrix::rowScale(const Vector<double>& vec) {
   const auto& scalevec = dynamic_cast<const DistributedVector<double>&>(vec);

   Bmat->rowScale(*scalevec.first);
   if (scalevec.last)
      Blmat->rowScale(*scalevec.last);

   assert(children.size() == scalevec.children.size());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->rowScale2(*(scalevec.children[it]), scalevec.last.get());
}

void DistributedMatrix::scalarMult(double num) {
   Amat->scalarMult(num);
   Bmat->scalarMult(num);
   Blmat->scalarMult(num);

   for (auto& it : children)
      it->scalarMult(num);
}

void DistributedMatrix::getDiagonal(Vector<double>& vec_) const {
   auto& vec = dynamic_cast<DistributedVector<double>&>(vec_);

   Bmat->getDiagonal(*vec.first);

   assert(children.size() == vec.children.size());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->getDiagonal(*vec.children[it]);
}

void DistributedMatrix::setToDiagonal(const Vector<double>& vec_) {
   const auto& vec = dynamic_cast<const DistributedVector<double>&>(vec_);

   Bmat->setToDiagonal(*vec.first);

   assert(children.size() == vec.children.size());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->setToDiagonal(*vec.children[it]);
}

/* y = beta * y + alpha * this * x */
void DistributedMatrix::mult(double beta, Vector<double>& y_, double alpha, const Vector<double>& x_) const {
   if (0.0 == alpha) {
      y_.scale(beta);
      return;
   }

   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_);

   assert(amatEmpty());
   Bmat->mult(beta, *y.first, alpha, *x.getLinkingVecNotHierarchicalTop());

   if (y.last) {
      if (iAmSpecial(iAmDistrib, mpiComm))
         Blmat->mult(beta, *y.last, alpha, *x.getLinkingVecNotHierarchicalTop());
      else
         y.last->setToZero();
   }

   assert(y.children.size() == children.size());
   assert(x.children.size() == children.size());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->mult2(beta, *y.children[it], alpha, *x.children[it], y.last.get());

   if (iAmDistrib && y.last)
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector<double>&>(*y.last).elements(), y.last->length(), mpiComm);
}


/* mult method for children; needed only for linking constraints */
void
DistributedMatrix::mult2(double beta, DistributedVector<double>& y, double alpha, const DistributedVector<double>& x,
   Vector<double>* yparentl_) const {
   assert(alpha != 0.0);
   assert(children.empty());
   assert(y.children.size() == children.size());
   assert(x.children.size() == children.size());
   assert(x.first);
   assert(y.first);

   Bmat->mult(beta, *y.first, alpha, *x.first);

   if (yparentl_) {
      Blmat->mult(1.0, *yparentl_, alpha, *x.first);
      if (!iAmSpecial(iAmDistrib, mpiComm))
         yparentl_->setToZero();
   }

   if (!amatEmpty()) {
      const Vector<double>* link_vec = x.getLinkingVecNotHierarchicalTop();
      assert(link_vec != x.first.get());
      Amat->mult(1.0, *y.first, alpha, *link_vec);
   }
}


void DistributedMatrix::transMult(double beta, Vector<double>& y_, double alpha, const Vector<double>& x_) const {
   if (0.0 == alpha) {
      y_.scale(beta);
      return;
   }

   const auto& x = dynamic_cast<const DistributedVector<double>&>(x_);
   auto& y = dynamic_cast<DistributedVector<double>&>(y_);

   const bool at_root = y.first.get() == y.getLinkingVecNotHierarchicalTop();
   assert(y.first);
   assert(x.first);

   if (iAmSpecial(iAmDistrib, mpiComm)) {
      Bmat->transMult(at_root ? beta : 1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.first);

      if (x.last)
         Blmat->transMult(1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.last);
   } else if (at_root)
      y.first->setToZero();

   assert(y.children.size() == children.size());
   assert(x.children.size() == children.size());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->transMult2(beta, *y.children[it], alpha, *x.children[it], x.last.get());

   if (iAmDistrib && y.first.get() == y.getLinkingVecNotHierarchicalTop())
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector<double>&>(*y.first).elements(), y.first->length(), mpiComm);
}

void
DistributedMatrix::transMult2(double beta, DistributedVector<double>& y, double alpha,
   const DistributedVector<double>& x, const Vector<double>* xvecl) const {
   assert(alpha != 0.0);
   assert(x.first);
   assert(y.first);
   assert(y.children.size() - children.size() == 0);
   assert(x.children.size() - children.size() == 0);
   assert(children.empty());

   if (!amatEmpty())
      Amat->transMult(1.0, *y.getLinkingVecNotHierarchicalTop(), alpha, *x.first);

   Bmat->transMult(beta, *y.first, alpha, *x.first);
   if (xvecl)
      Blmat->transMult(1.0, *y.first, alpha, *xvecl);
}

double DistributedMatrix::inf_norm() const {
   double nrm = 0.0;

   for (const auto& it : children)
      nrm = std::max(nrm, it->inf_norm());

   if (iAmDistrib)
      PIPS_MPIgetMaxInPlace(nrm, mpiComm);

   nrm = std::max(nrm, std::max(Amat->inf_norm(), Bmat->inf_norm()));
   nrm = std::max(nrm, Blmat->inf_norm());

   return nrm;
}

double DistributedMatrix::abminnormNonZero(double tol) const {
   double nrm = std::numeric_limits<double>::infinity();

   for (const auto& it : children)
      nrm = std::min(nrm, it->abminnormNonZero(tol));

   if (iAmDistrib)
      PIPS_MPIgetMinInPlace(nrm, mpiComm);

   nrm = std::min(nrm, std::min(Amat->abminnormNonZero(tol), Bmat->abminnormNonZero(tol)));
   nrm = std::min(nrm, Blmat->abminnormNonZero(tol));

   return nrm;
}

void DistributedMatrix::getLinkVarsNnz(std::vector<int>& vec) const {
   assert(hasSparseMatrices());

   for (const auto& it : children)
      it->getLinkVarsNnzChild(vec);

   if (iAmDistrib) {
      int* buffer = new int[vec.size()];
      MPI_Allreduce(&vec[0], buffer, static_cast<int>(vec.size()), MPI_INT, MPI_SUM, mpiComm);

      std::memcpy(&vec[0], buffer, vec.size() * sizeof(int));

      delete[] buffer;
   }
}

void DistributedMatrix::getLinkVarsNnzChild(std::vector<int>& vec) const {
   assert(children.empty());
   assert(hasSparseMatrices());

   dynamic_cast<const SparseMatrix&>(*Amat).getLinkVarsNnz(vec);
}

void DistributedMatrix::write_to_streamDenseBorderedChild(const StripMatrix& border_left, std::ostream& out,
   int offset) const {
   assert(border_left.children.size() == this->children.size());

   if (Bmat->is_a(kDistributedMatrix)) {
      assert(border_left.first->is_a(kStripMatrix));
      dynamic_cast<DistributedMatrix&>(*Bmat).write_to_streamDenseBordered(dynamic_cast<StripMatrix&>(*border_left.first),
         out, offset);
   }
      /// Border.mat | Amat | offset | Bmat ///
   else {
      assert(border_left.first->is_a(kSparseGenMatrix));
      assert(hasSparseMatrices());
      assert(PIPS_MPIgetRank(mpiComm) == 0);

      const auto mB = Bmat->n_rows();

      assert(mB == Amat->n_rows() && Amat->n_rows() == border_left.first->n_rows());

      for (int row = 0; row < mB; ++row) {
         border_left.first->write_to_streamDenseRow(out, row);
         out << "|\t";
         Amat->write_to_streamDenseRow(out, row);

         for (int i = 0; i < offset; ++i)
            out << "\t";
         Bmat->write_to_streamDenseRow(out, row);
         out << "\n";
      }
   }
}


void
DistributedMatrix::write_to_streamDenseBordered(const StripMatrix& border_left, std::ostream& out, int offset) const {
   assert(border_left.children.size() == this->children.size());
   assert(hasSparseMatrices());
   assert(!children.empty());

   const int original_offset = offset;
   const int my_rank = PIPS_MPIgetRank(mpiComm);

   if (iAmDistrib)
      MPI_Barrier(mpiComm);

   const auto mBmat = this->Bmat->n_rows();
   assert(mBmat == border_left.first->n_rows());

   assert(Bmat->is_a(kSparseGenMatrix));
   assert(border_left.first->is_a(kSparseGenMatrix));

   /// Border.mat | Bmat ///
   if (my_rank == 0) {
      for (int i = 0; i < mBmat; ++i) {
         dynamic_cast<const SparseMatrix&>(*border_left.first).write_to_streamDenseRow(out, i);
         out << "|\t";
         dynamic_cast<const SparseMatrix&>(*this->Bmat).write_to_streamDenseRow(out, i);
         out << "\n";
      }
   }

   /// write children ///
   /// BorderChild | offset | Child ///
   std::ostringstream child_stream{};
   if (!amatEmpty())
      assert(children.empty());

   for (size_t it = 0; it < children.size(); it++) {
      MPI_Barrier(mpiComm);
      const StripMatrix& border_child = *border_left.children[it];

      children[it]->write_to_streamDenseBorderedChild(border_child, child_stream, offset);
      offset += static_cast<int>(PIPS_MPIgetSum(children[it]->n_columns(), mpiComm));

      MPI_Barrier(mpiComm);
   }

   const std::string children_string = child_stream.str();
   const std::string all_children = PIPS_MPIallgatherString(children_string, mpiComm);
   if (my_rank == 0)
      out << all_children;

   /// border.bl_mat | Blmat | offset | children ///
   const auto mlink = this->Blmat->n_rows();
   if (mlink > 0) {
      assert(border_left.last);
      assert(border_left.last->n_rows() == mlink);

      // for each row r do:
      for (int r = 0; r < mlink; r++) {
         MPI_Barrier(mpiComm);
         std::ostringstream link_row_stream;

         if (my_rank == 0) {
            dynamic_cast<const SparseMatrix&>(*border_left.last).write_to_streamDenseRow(link_row_stream, r);
            link_row_stream << "|\t";
            dynamic_cast<const SparseMatrix&>(*this->Blmat).write_to_streamDenseRow(link_row_stream, r);
            for (int i = 0; i < original_offset; ++i)
               link_row_stream << "\t";
         }

         for (const auto& it : children)
            it->Blmat->write_to_streamDenseRow(link_row_stream, r);

         const std::string children_link_row = link_row_stream.str();
         const std::string link_row = PIPS_MPIallgatherString(children_link_row, mpiComm);

         if (my_rank == 0) {
            out << link_row;
            out << "\n";
         }
         MPI_Barrier(mpiComm);
      }

      std::ostringstream dashed_line_stream;
      if (my_rank == 0) {
         dynamic_cast<const SparseMatrix&>(*border_left.last).writeDashedLineToStream(dashed_line_stream);
         dashed_line_stream << "|\t";
      }
      writeDashedLineToStream(dashed_line_stream, original_offset);

      if (my_rank == 0) {
         const std::string dashed_line = dashed_line_stream.str();
         out << dashed_line << "\n";
      }
   }
   if (iAmDistrib)
      MPI_Barrier(mpiComm);
}

void DistributedMatrix::writeDashedLineToStream(std::ostream& out, int offset) const {
   const int my_rank = PIPS_MPIgetRank(mpiComm);

   MPI_Barrier(mpiComm);
   std::ostringstream link_row_stream;

   if (my_rank == 0) {
      dynamic_cast<const SparseMatrix&>(*this->Blmat).writeDashedLineToStream(link_row_stream);
      for (int i = 0; i < offset; ++i)
         link_row_stream << "-\t";
   }

   for (const auto& it : children)
      it->Blmat->writeDashedLineToStream(link_row_stream);

   const std::string children_link_row = link_row_stream.str();
   const std::string link_row = PIPS_MPIallgatherString(children_link_row, mpiComm);

   if (my_rank == 0) {
      out << link_row;
      out << "\n";
   }
   MPI_Barrier(mpiComm);

}

void DistributedMatrix::write_to_streamDense(std::ostream& out, int offset) const {
   assert(hasSparseMatrices());
   assert(!children.empty());

   const int original_offset = offset;
   const int my_rank = PIPS_MPIgetRank(mpiComm);

   if (iAmDistrib)
      MPI_Barrier(mpiComm);

   assert(Bmat->is_a(kSparseGenMatrix));

   /// Bmat ///
   if (my_rank == 0) {
      for (int i = 0; i < Bmat->n_rows(); ++i) {
         for (int j = 0; j < offset; ++j)
            out << "\t";
         dynamic_cast<const SparseMatrix&>(*this->Bmat).write_to_streamDenseRow(out, i);
         out << "\n";
      }
   }

   /// write children ///
   /// offset | Child ///
   std::ostringstream child_stream{};
   if (!amatEmpty())
      assert(children.empty());

   for (const auto& it : children) {
      MPI_Barrier(mpiComm);
      it->write_to_streamDenseChild(child_stream, offset);
      offset += static_cast<int>(PIPS_MPIgetSum(it->n_rows(), mpiComm));

      MPI_Barrier(mpiComm);
   }

   const std::string children_string = child_stream.str();
   const std::string all_children = PIPS_MPIallgatherString(children_string, mpiComm);
   if (my_rank == 0)
      out << all_children;

   /// Blmat | offset | children ///
   const auto mlink = this->Blmat->n_rows();
   if (mlink > 0) {
      // for each row r do:
      for (int r = 0; r < mlink; r++) {
         MPI_Barrier(mpiComm);
         std::ostringstream link_row_stream;

         if (my_rank == 0) {
            dynamic_cast<const SparseMatrix&>(*this->Blmat).write_to_streamDenseRow(link_row_stream, r);
            for (int i = 0; i < original_offset; ++i)
               link_row_stream << "\t";
         }

         for (const auto& it : children)
            it->Blmat->write_to_streamDenseRow(link_row_stream, r);

         const std::string children_link_row = link_row_stream.str();
         const std::string link_row = PIPS_MPIallgatherString(children_link_row, mpiComm);

         if (my_rank == 0) {
            out << link_row;
            out << "\n";
         }
         MPI_Barrier(mpiComm);
      }
   }
   if (iAmDistrib)
      MPI_Barrier(mpiComm);
}

/** writes child matrix blocks, offset indicates the offset between A and B block. */
void DistributedMatrix::write_to_streamDenseChild(std::ostream& out, int offset) const {
   if (Bmat->is_a(kDistributedMatrix))
      dynamic_cast<DistributedMatrix&>(*Bmat).write_to_streamDense(out, offset);
      /// Border.mat | Amat | offset | Bmat ///
   else {
      assert(hasSparseMatrices());
      assert(PIPS_MPIgetRank(mpiComm) == 0);
      assert(Bmat->n_rows() == Amat->n_rows());

      for (int row = 0; row < Bmat->n_rows(); ++row) {
         Amat->write_to_streamDenseRow(out, row);
         for (int i = 0; i < offset; ++i)
            out << "\t";
         Bmat->write_to_streamDenseRow(out, row);
         out << "\n";
      }
   }
}

void DistributedMatrix::writeMPSformatRows(std::ostream& out, int rowType, const Vector<double>* irhs) const {
   assert(hasSparseMatrices());

   const int myRank = PIPS_MPIgetRank(mpiComm);
   std::string rt;
   if (rowType == 0)
      rt = "E";
   else if (rowType == 1)
      rt = "L";
   else if (rowType == 2)
      rt = "G";
   else
      assert(0);

   const auto* irhsStoch = dynamic_cast<const DistributedVector<double>*>(irhs);

   if (myRank == 0) {
      // A_0 block:
      const auto mB = this->Bmat->n_rows();
      for (int i = 0; i < mB; i++) {
         if (!irhs || (irhsStoch && dynamic_cast<const SimpleVector<double>&>(*irhsStoch->first)[i] != 0.0))
            out << " " << rt << " row_" << rt << "_" << "R" << "_" << i << "\n";
      }
      // linking rows:
      if (Blmat) {
         const auto mBl = this->Blmat->n_rows();
         for (int i = 0; i < mBl; i++) {
            if (!irhs || (irhsStoch && dynamic_cast<const SimpleVector<double>&>(*irhsStoch->last)[i] != 0.0))
               out << " " << rt << " row_" << rt << "_" << "L" << "_" << i << "\n";
            if (!irhs || (irhsStoch && dynamic_cast<const SimpleVector<double>&>(*irhsStoch->last)[i] != 0.0))
               out << " " << rt << " row_" << rt << "_" << "L" << "_" << i << "\n";
         }
      }
   }
   for (size_t it = 0; it < children.size(); it++) {
      const auto mA = children[it]->Amat->n_rows();
      for (int i = 0; i < mA; i++) {
         if (!irhs ||
            (irhsStoch && dynamic_cast<const SimpleVector<double>&>(*irhsStoch->children[it]->first)[i] != 0.0))
            out << " " << rt << " row_" << rt << "_" << it << "_" << i << "\n";
      }
   }
}

void DistributedMatrix::initTransposed(bool dynamic) const {
   assert(hasSparseMatrices());
   dynamic_cast<const SparseMatrix&>(*Bmat).initTransposed(dynamic);
   dynamic_cast<const SparseMatrix&>(*Blmat).initTransposed(dynamic);

   for (const auto& it : children)
      it->initTransposedChild(dynamic);
}

void DistributedMatrix::deleteTransposed() const {
   assert(hasSparseMatrices());

   dynamic_cast<const SparseMatrix&>(*Amat).deleteTransposed();
   dynamic_cast<const SparseMatrix&>(*Bmat).deleteTransposed();
   dynamic_cast<const SparseMatrix&>(*Blmat).deleteTransposed();

   for (const auto& it : children)
      it->deleteTransposed();
}

void DistributedMatrix::initTransposedChild(bool dynamic) const {
   assert(hasSparseMatrices());
   dynamic_cast<SparseMatrix&>(*Amat).initTransposed(dynamic);
   dynamic_cast<SparseMatrix&>(*Bmat).initTransposed(dynamic);

   if (Blmat != nullptr)
      dynamic_cast<SparseMatrix&>(*Blmat).initTransposed(dynamic);
}

int DistributedMatrix::numberOfNonZeros() const {
   assert(hasSparseMatrices());
   int nnz = 0;

   for (const auto& it : children)
      nnz += it->numberOfNonZeros();

   if (iAmDistrib)
      PIPS_MPIgetSumInPlace(nnz, mpiComm);

   nnz += Amat->numberOfNonZeros() + Bmat->numberOfNonZeros() + Blmat->numberOfNonZeros();

   return nnz;
}

void DistributedMatrix::getNnzPerRow(Vector<int>& nnzVec, Vector<int>* linkParent) const {
   assert(hasSparseMatrices());
   auto& nnzVecStoch = dynamic_cast<DistributedVector<int>&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   SimpleVector<int>* nnzvecl = nullptr;

   dynamic_cast<SparseMatrix&>(*Bmat).addNnzPerRow(*(nnzVecStoch.first));

   if (linkParent != nullptr)
      dynamic_cast<SparseMatrix&>(*Amat).addNnzPerRow(*(nnzVecStoch.first));

   /* with linking constraints? */
   if (nnzVecStoch.last || linkParent) {
      assert(nnzVecStoch.last == nullptr || linkParent == nullptr);

      if (linkParent)
         nnzvecl = dynamic_cast<SimpleVector<int>*>(linkParent);
      else
         nnzvecl = dynamic_cast<SimpleVector<int>*>(nnzVecStoch.last.get());

      if (linkParent != nullptr || iAmSpecial(iAmDistrib, mpiComm))
         dynamic_cast<SparseMatrix&>(*Blmat).addNnzPerRow(*nnzvecl);
   }


   for (size_t it = 0; it < children.size(); it++)
      children[it]->getNnzPerRow(*(nnzVecStoch.children[it]), nnzvecl);

   // distributed, with linking constraints, and at root?
   if (iAmDistrib && nnzVecStoch.last != nullptr && linkParent == nullptr) {
      PIPS_MPIsumArrayInPlace(nnzvecl->elements(), nnzvecl->length(), mpiComm);
   }
}

void DistributedMatrix::sum_transform_rows(Vector<double>& result_, const std::function<double(const double&)>& transform) const {
   const bool at_root = !children.empty();

   auto& result = dynamic_cast<DistributedVector<double>&>(result_);

   if (at_root)
   {
      assert(amatEmpty());
      assert(result.children.size() == children.size());
   }

   const bool has_linking = at_root ? result.last != nullptr : result.parent->last != nullptr;

   Bmat->sum_transform_rows(*result.first, transform);

   if (!amatEmpty()) {
      assert(!Bmat->is_a(kDistributedMatrix));
      Amat->sum_transform_rows(*result.first, transform);
   }

   if (has_linking) {
      if (at_root && iAmSpecial(iAmDistrib, mpiComm)) {
         Blmat->sum_transform_rows(*result.last, transform);
      } else {
         Blmat->sum_transform_rows(*result.parent->last, transform);
      }
   }

   for (size_t it = 0; it < children.size(); it++){
      children[it]->sum_transform_rows(*(result.children[it]), transform);
   }

   if (at_root && iAmDistrib) {
      PIPS_MPIsumArrayInPlace(dynamic_cast<SimpleVector<double>&>(*result.last).elements(), result.last->length(), mpiComm);
   }
}

void DistributedMatrix::getNnzPerCol(Vector<int>& nnzVec, Vector<int>* linkParent) const {
   assert(hasSparseMatrices());
   auto& nnzVecStoch = dynamic_cast<DistributedVector<int>&>(nnzVec);

   // assert tree compatibility
   assert(nnzVecStoch.children.size() == children.size());

   auto* vec = dynamic_cast<SimpleVector<int>*>(nnzVecStoch.first.get());

   if (iAmSpecial(iAmDistrib, mpiComm) || linkParent != nullptr) {
      dynamic_cast<SparseMatrix&>(*Bmat).addNnzPerCol(*(vec));

      /* with linking constraints? */
      if (Blmat->n_rows() > 0)
         dynamic_cast<SparseMatrix&>(*Blmat).addNnzPerCol(*vec);
   }

   // not at root?
   if (linkParent != nullptr)
      dynamic_cast<SparseMatrix&>(*Amat).addNnzPerCol(*linkParent);
   else {
      for (size_t it = 0; it < children.size(); it++)
         children[it]->getNnzPerCol(*(nnzVecStoch.children[it]), vec);
   }

   // distributed and at root?
   if (iAmDistrib && linkParent == nullptr) {
      PIPS_MPIsumArrayInPlace(vec->elements(), vec->length(), mpiComm);
   }
}

void DistributedMatrix::getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* col_scale_,
   Vector<double>& minmax_) const {
   assert(amatEmpty());

   auto& minmax = dynamic_cast<DistributedVector<double>&>(minmax_);

   const bool scale = col_scale_;
   const bool has_linking = minmax.last != nullptr;

   const DistributedVector<double>* const col_scale = scale ? dynamic_cast<const DistributedVector<double>*>(col_scale_)
      : nullptr;
   const Vector<double>* const col_scale_vec = scale ? col_scale->getLinkingVecNotHierarchicalTop() : nullptr;

   Bmat->getRowMinMaxVec(getMin, initializeVec, col_scale_vec, *minmax.first);

   if (has_linking) {
      if (initializeVec) {
         if (getMin)
            minmax.last->setToConstant(std::numeric_limits<double>::max());
         else
            minmax.last->setToZero();
      }

      if (iAmSpecial(iAmDistrib, mpiComm))
         Blmat->getRowMinMaxVec(getMin, false, col_scale_vec, *minmax.last);
   }

   assert(minmax.children.size() == children.size());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->getRowMinMaxVecChild(getMin, initializeVec, scale ? col_scale->children[it].get() : nullptr,
         *(minmax.children[it]), minmax.last.get());

   if (iAmDistrib) {
      if (getMin)
         PIPS_MPIminArrayInPlace(dynamic_cast<SimpleVector<double>&>(*minmax.last).elements(), minmax.last->length(),
            mpiComm);
      else
         PIPS_MPImaxArrayInPlace(dynamic_cast<SimpleVector<double>&>(*minmax.last).elements(), minmax.last->length(),
            mpiComm);
   }
}

void DistributedMatrix::getRowMinMaxVecChild(bool getMin, bool initializeVec, const Vector<double>* col_scale_,
   Vector<double>& minmax_,
   Vector<double>* minmax_linking_cons) const {
   assert(children.empty());
   auto& minmax = dynamic_cast<DistributedVector<double>&>(minmax_);

   const bool scale = col_scale_;
   const bool has_linking = minmax_linking_cons;

   const auto* const col_scale = dynamic_cast<const DistributedVector<double>*>(col_scale_);

   const Vector<double>* const col_scale_vec = scale ? col_scale->first.get() : nullptr;
   const Vector<double>* const col_scale_linkingvar_vec = scale ? col_scale->getLinkingVecNotHierarchicalTop()
      : nullptr;

   Bmat->getRowMinMaxVec(getMin, initializeVec, col_scale_vec, *(minmax.first));

   if (!amatEmpty()) {
      assert(!Bmat->is_a(kDistributedMatrix));
      Amat->getRowMinMaxVec(getMin, false, col_scale_linkingvar_vec, *(minmax.first));
   }

   if (has_linking)
      Blmat->getRowMinMaxVec(getMin, false, col_scale_vec, *minmax_linking_cons);
}

void DistributedMatrix::getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec_,
   Vector<double>& minmaxVec_) const {
   assert(amatEmpty());
   auto& minmaxVec = dynamic_cast<DistributedVector<double>&>(minmaxVec_);
   const auto* rowScaleVec = dynamic_cast<const DistributedVector<double>*>(rowScaleVec_);

   const bool scale = rowScaleVec;
   const bool has_linking = Blmat->n_rows() > 0;

   const Vector<double>* row_scale_vec = scale ? rowScaleVec->first.get() : nullptr;
   const Vector<double>* row_scale_link = scale ? rowScaleVec->last.get() : nullptr;

   if (minmaxVec.first.get() == minmaxVec.getLinkingVecNotHierarchicalTop())
      Bmat->getColMinMaxVec(getMin, initializeVec, row_scale_vec, *minmaxVec.getLinkingVecNotHierarchicalTop());
   else
      Bmat->getColMinMaxVec(getMin, false, row_scale_vec, *minmaxVec.getLinkingVecNotHierarchicalTop());

   if (has_linking)
      Blmat->getColMinMaxVec(getMin, false, row_scale_link, *minmaxVec.getLinkingVecNotHierarchicalTop());

   assert(minmaxVec.children.size() == children.size());

   for (size_t it = 0; it < children.size(); it++)
      children[it]->getColMinMaxVecChild(getMin, initializeVec, scale ? rowScaleVec->children[it].get() : nullptr,
         row_scale_link,
         *(minmaxVec.children[it]));

   if (iAmDistrib) {
      auto& mvec = dynamic_cast<SimpleVector<double>&>(*minmaxVec.first);
      if (getMin)
         PIPS_MPIminArrayInPlace(mvec.elements(), mvec.length(), mpiComm);
      else
         PIPS_MPImaxArrayInPlace(mvec.elements(), mvec.length(), mpiComm);
   }
}

void DistributedMatrix::getColMinMaxVecChild(bool getMin, bool initializeVec, const Vector<double>* rowScale_,
   const Vector<double>* rowScaleParent,
   Vector<double>& minmaxVec_) const {
   assert(children.empty());
   auto& minmaxVec = dynamic_cast<DistributedVector<double>&>(minmaxVec_);

   const bool scale = rowScale_;
   const bool has_linking = Blmat->n_rows() > 0;

   const auto* rowScale = dynamic_cast<const DistributedVector<double>*>(rowScale_);
   const Vector<double>* row_scale_vec = scale ? rowScale->first.get() : nullptr;

   Bmat->getColMinMaxVec(getMin, initializeVec, row_scale_vec, *minmaxVec.first);

   if (!amatEmpty()) {
      assert(!Bmat->is_a(kDistributedMatrix));
      Amat->getColMinMaxVec(getMin, false, row_scale_vec, *minmaxVec.getLinkingVecNotHierarchicalTop());
   }

   if (has_linking)
      Blmat->getColMinMaxVec(getMin, false, rowScaleParent, *minmaxVec.first);
}

void DistributedMatrix::addRowSums(Vector<double>& sumVec, Vector<double>* linkParent) const {
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      assert(false && "TODO : hierarchical version");
   assert(hasSparseMatrices());

   auto& sumVecStoch = dynamic_cast<DistributedVector<double>&>(sumVec);
   SimpleVector<double>* mvecl = nullptr;

   // assert tree compatibility
   assert(sumVecStoch.children.size() == children.size());

   Bmat->addRowSums(*sumVecStoch.first);

   // not at root?
   if (linkParent != nullptr)
      Amat->addRowSums(*sumVecStoch.first);

   /* with linking constraints? */
   if (sumVecStoch.last || linkParent) {
      assert(sumVecStoch.last == nullptr || linkParent == nullptr);

      // at root?
      if (linkParent == nullptr)
         mvecl = dynamic_cast<SimpleVector<double>*>(sumVecStoch.last.get());
      else
         mvecl = dynamic_cast<SimpleVector<double>*>(linkParent);

      if (linkParent != nullptr || iAmSpecial(iAmDistrib, mpiComm))
         Blmat->addRowSums(*mvecl);
   }

   for (size_t it = 0; it < children.size(); it++)
      children[it]->addRowSums(*(sumVecStoch.children[it]), mvecl);

   // distributed, with linking constraints, and at root?
   if (iAmDistrib && sumVecStoch.last != nullptr && linkParent == nullptr) {
      assert(mvecl != nullptr);

      // sum up linking constraints vectors
      const int locn = mvecl->length();
      auto* buffer = new double[locn];

      MPI_Allreduce(mvecl->elements(), buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      mvecl->copyFromArray(buffer);

      delete[] buffer;
   }
}

void DistributedMatrix::addColSums(Vector<double>& sumVec, Vector<double>* linkParent) const {
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      assert(false && "TODO : hierarchical version");
   assert(hasSparseMatrices());

   auto& sumVecStoch = dynamic_cast<DistributedVector<double>&>(sumVec);

   // assert tree compatibility
   assert(sumVecStoch.children.size() == children.size());

   auto* const mvec = dynamic_cast<SimpleVector<double>*>(sumVecStoch.first.get());

   if (iAmSpecial(iAmDistrib, mpiComm) || linkParent != nullptr)
      Bmat->addColSums(*mvec);

   /* with linking constraints? */
   if (Blmat->n_rows() > 0 && (iAmSpecial(iAmDistrib, mpiComm) || linkParent != nullptr))
      Blmat->addColSums(*mvec);

   // not at root?
   if (linkParent != nullptr)
      Amat->addColSums(*linkParent);
   else {
      for (size_t it = 0; it < children.size(); it++)
         children[it]->addColSums(*(sumVecStoch.children[it]), mvec);
   }

   // distributed and at root?
   if (iAmDistrib && linkParent == nullptr) {
      const int locn = mvec->length();
      auto* const entries = mvec->elements();
      auto* buffer = new double[locn];

      MPI_Allreduce(entries, buffer, locn, MPI_DOUBLE, MPI_SUM, mpiComm);

      mvec->copyFromArray(buffer);

      delete[] buffer;
   }
}

void DistributedMatrix::initStaticStorageFromDynamic(const Vector<int>& rowNnzVec, const Vector<int>& colNnzVec,
   const Vector<int>* rowLinkVec,
   const Vector<int>* colParentVec) {
   assert(hasSparseMatrices());

   const auto& rowNnzVecStoch = dynamic_cast<const DistributedVector<int>&>(rowNnzVec);
   const auto& colNnzVecStoch = dynamic_cast<const DistributedVector<int>&>(colNnzVec);

   assert(rowNnzVecStoch.children.size() == colNnzVecStoch.children.size());

   const auto* const rowvec = dynamic_cast<const SimpleVector<int>*>(rowNnzVecStoch.first.get());
   const auto* const colvec = dynamic_cast<const SimpleVector<int>*>(colNnzVecStoch.first.get());

   const auto* const rowlink = dynamic_cast<const SimpleVector<int>*>(rowNnzVecStoch.last.get());
   assert(rowvec);
   assert(colvec);

   dynamic_cast<SparseMatrix&>(*Amat).initStaticStorageFromDynamic(*rowvec,
      colParentVec); // initialized with colVec == nullptr for parent
   dynamic_cast<SparseMatrix&>(*Bmat).initStaticStorageFromDynamic(*rowvec, colvec);

   // at root?
   if (colParentVec == nullptr) {
      assert(rowLinkVec == nullptr);

      if (rowlink != nullptr)
         dynamic_cast<SparseMatrix&>(*Blmat).initStaticStorageFromDynamic(*rowlink, colvec);

      for (size_t it = 0; it < children.size(); it++)
         children[it]->initStaticStorageFromDynamic(*(rowNnzVecStoch.children[it]), *(colNnzVecStoch.children[it]),
            rowlink, colvec);
   } else {
      assert(children.empty());
      if (rowLinkVec != nullptr)
         dynamic_cast<SparseMatrix&>(*Blmat).initStaticStorageFromDynamic(*rowLinkVec, colvec);
   }

}

void DistributedMatrix::freeDynamicStorage() {
   assert(hasSparseMatrices());

   dynamic_cast<SparseMatrix&>(*Amat).freeDynamicStorage();
   dynamic_cast<SparseMatrix&>(*Bmat).freeDynamicStorage();
   dynamic_cast<SparseMatrix&>(*Blmat).freeDynamicStorage();

   for (auto& it : children)
      it->freeDynamicStorage();
}

void DistributedMatrix::recomputeSize(DistributedMatrix* parent) {
   m = 0;
   n = 0;

   if (Bmat->is_a(kDistributedMatrix)) {
      assert(children.empty());
      assert(inner_leaf);

      dynamic_cast<DistributedMatrix&>(*Bmat).recomputeSize();
   }

   if (!inner_root)
      std::tie(m, n) = Bmat->n_rows_columns();

   assert(m >= 0);
   assert(n >= 0);

   for (auto& it : children) {
      it->recomputeSize(this);

      const auto[m_child, n_child] = it->n_rows_columns();

      m += m_child;
      n += n_child;
   }

   if (!parent) {
      m += Blmat->n_rows();
   }
}

void DistributedMatrix::updateKLinkConsCount(std::vector<int>& linkCount) const {
   assert(hasSparseMatrices());

   if (!Blmat)
      return;

   const auto m_Blmat = Blmat->n_rows();
   assert(m_Blmat > 0);
   assert(linkCount.size() == size_t(m_Blmat));

   for (const auto& it : children) {
      if (!(it->is_a(kStochGenDummyMatrix))) {
         assert(it->Blmat);
         dynamic_cast<SparseMatrix&>(*it->Blmat).updateNonEmptyRowsCount(linkCount);
      }
   }

   if (iAmDistrib)
      MPI_Allreduce(MPI_IN_PLACE, &linkCount[0], m_Blmat, MPI_INT, MPI_SUM, mpiComm);
}

void DistributedMatrix::updateKLinkVarsCount(std::vector<int>& link_block_count) const {
   assert(hasSparseMatrices());
   const auto n_Bmat = Bmat->n_columns();

   if (n_Bmat == 0)
      return;

   assert(link_block_count.size() == size_t(n_Bmat));

   for (const auto& it : children) {
      if (!(it->is_a(kStochGenDummyMatrix))) {
         dynamic_cast<SparseMatrix&>(*it->Amat).getTranspose().updateNonEmptyRowsCount(link_block_count);
         dynamic_cast<SparseMatrix&>(*it->Amat).deleteTransposed();
      }
   }

   if (iAmDistrib)
      PIPS_MPIsumArrayInPlace(link_block_count, mpiComm);
}

void
DistributedMatrix::get2LinkStartBlocksAndCountsNew(std::vector<int>& block_start, std::vector<int>& block_count) const {
   assert(hasSparseMatrices());
   block_start.clear();
   block_count.clear();

   if (Blmat == nullptr)
      return;

   const auto m_Blmat = Blmat->n_rows();
   if (m_Blmat == 0)
      return;
   assert(m_Blmat > 0);

   const int n_blocks = static_cast<int>(children.size());
   block_count.resize(m_Blmat);
   std::fill(block_count.begin(), block_count.end(), 0);

   /* init with max + 1 and max - 1 for allreduce later */
   block_start.resize(m_Blmat);
   std::fill(block_start.begin(), block_start.end(), n_blocks);
   std::vector<int> block_end(m_Blmat, -1);

   std::vector<bool> is_2_link(m_Blmat, false);

   for (size_t it = 0; it < children.size(); it++)
      if (!(children[it]->is_a(kStochGenDummyMatrix))) {
         assert(children[it]->Blmat);
         dynamic_cast<const SparseMatrix&>(*children[it]->Blmat).updateNonEmptyRowsCountNew(static_cast<int>(it),
            block_count, block_start, block_end);
      }

   if (iAmDistrib) {
      // TODO : one can filter the non-candidates locally first on all processes
      PIPS_MPIminArrayInPlace(block_start, mpiComm);
      PIPS_MPImaxArrayInPlace(block_end, mpiComm);
      PIPS_MPIsumArrayInPlace(block_count, mpiComm);
   }

   for (int it = 0; it < m_Blmat; ++it) {
      const int start = block_start[it];
      const int end = block_end[it];

      assert(start == n_blocks || (0 <= start && start < n_blocks));
      if (start == n_blocks)
         assert(end == -1);
      else
         assert(start <= end && end < n_blocks);

      if (end == start + 1)
         assert(block_count[it] == 2);
      else {
         /* not a consecutive 2 link */
         block_start[it] = -1;
      }
   }
}

std::vector<int> DistributedMatrix::get2LinkStartBlocks() const {
   assert(hasSparseMatrices());
   if (Blmat == nullptr)
      return std::vector<int>();

   const auto m_loc = Blmat->n_rows();

   if (m_loc == 0)
      return std::vector<int>();
   assert(m_loc > 0);

   std::vector<int> linkBlockCount(m_loc, 0);
   std::vector<int> linkBlockStart(m_loc, -1);
   std::vector<int> linkBlockEnd(m_loc, -1);
   std::vector<bool> is2link(m_loc, false);

   for (size_t it = 0; it < children.size(); it++)
      if (!(children[it]->is_a(kStochGenDummyMatrix))) {
         assert(children[it]->Blmat);
         dynamic_cast<const SparseMatrix&>(*children[it]->Blmat).updateNonEmptyRowsCount(static_cast<int>(it),
            linkBlockCount, linkBlockStart, linkBlockEnd);
      }

   if (iAmDistrib)
      PIPS_MPIsumArrayInPlace(linkBlockCount, mpiComm);

   /* filter out process local two links already */
   for (int i = 0; i < m_loc; i++) {
      assert(linkBlockEnd[i] == -1 || linkBlockStart[i] <= linkBlockEnd[i]);

      if (linkBlockCount[i] == 2 && (linkBlockEnd[i] - linkBlockStart[i]) == 1) {
         assert(linkBlockStart[i] >= 0 && linkBlockEnd[i] >= 0);
         is2link[i] = true;
      }
   }


   if (iAmDistrib) {
      // find 2-links between different processes

      const int size = PIPS_MPIgetSize(mpiComm);
      assert(size > 0);

      // 1. allgather number of local 2-link candidates
      std::vector<int> localCandsIdx;
      std::vector<int> localCandsBlock;
      std::vector<int> candsPerProc(size, -1);

      /* a local candidate is a linking row that appears in exactly two blocks, starts on this process but where the second block is is not stored on this process */
      for (int i = 0; i < m_loc; i++)
         if (linkBlockCount[i] == 2 && linkBlockStart[i] >= 0 && linkBlockEnd[i] == -1) {
            assert(!is2link[i]);

            localCandsIdx.push_back(i);
            assert(unsigned(linkBlockStart[i]) < children.size());

            localCandsBlock.push_back(linkBlockStart[i]);
         }

      const int localcount = static_cast<int>(localCandsIdx.size());

      PIPS_MPIallgather(&localcount, 1, &candsPerProc[0], 1, mpiComm);

#ifndef NDEBUG
      for (int i : candsPerProc)
         assert(i >= 0);
#endif

      // 2. allgatherv 2-link candidates
      std::vector<int> displacements(size + 1, 0);
      for (int i = 1; i <= size; i++)
         displacements[i] = candsPerProc[i - 1] + displacements[i - 1];

      const int nAllCands = displacements[size];

      std::vector<int> allCandsRow(nAllCands, -1);
      std::vector<int> allCandsBlock(nAllCands, -1);

      MPI_Allgatherv(&localCandsIdx[0], localcount, MPI_INT, &allCandsRow[0], &candsPerProc[0], &displacements[0],
         MPI_INT, mpiComm);

      MPI_Allgatherv(&localCandsBlock[0], localcount, MPI_INT, &allCandsBlock[0], &candsPerProc[0], &displacements[0],
         MPI_INT, mpiComm);

#ifndef NDEBUG
      for (size_t i = 0; i < allCandsRow.size(); i++)
         assert(allCandsRow[i] >= 0 && allCandsRow[i] < m_loc && allCandsBlock[i] >= 0 &&
            allCandsBlock[i] < static_cast<int>(children.size()));
#endif


      // 3. check which candidates are indeed 2-links
      std::vector<int> blocksHash(m_loc, -1);
      for (int i = 0; i < size - 1; i++) {
         // hash
         for (int j = displacements[i]; j < displacements[i + 1]; j++)
            blocksHash[allCandsRow[j]] = allCandsBlock[j];

         // compare with next
         for (int j = displacements[i + 1]; j < displacements[i + 2]; j++) {
            assert(allCandsBlock[j] > 0);
            const int candRow = allCandsRow[j];
            const int candBlock = allCandsBlock[j];
            if (blocksHash[candRow] >= 0) {
               assert(blocksHash[candRow] != candBlock);

               const int startBlock = std::min(blocksHash[candRow], candBlock);
               const int endBlock = std::max(blocksHash[candRow], candBlock);

               assert(startBlock >= 0 && endBlock >= 0);

               if (endBlock - startBlock != 1)
                  continue;

               assert(!is2link[candRow]);
               is2link[candRow] = true;

               // start block owned by this MPI process?
               if (!children[startBlock]->is_a(kStochGenDummyMatrix)) {
                  assert(children[endBlock]->is_a(kStochGenDummyMatrix));
                  linkBlockStart[candRow] = startBlock;
               } else {
                  assert(children[startBlock]->is_a(kStochGenDummyMatrix));
                  linkBlockStart[candRow] = -1;
               }
            }
         }

         // un-hash
         for (int j = displacements[i]; j < displacements[i + 1]; j++)
            blocksHash[allCandsRow[j]] = -1;
      }
   }

   // correct block identifier
   for (int i = 0; i < m_loc; i++)
      if (!is2link[i])
         linkBlockStart[i] = -1;

   if (iAmDistrib)
      PIPS_MPImaxArrayInPlace(linkBlockStart, mpiComm);

   return linkBlockStart;
}


void DistributedMatrix::permuteLinkingVars(const std::vector<unsigned int>& permvec) {
   assert(hasSparseMatrices());
   if (Blmat)
      dynamic_cast<SparseMatrix&>(*Blmat).permuteCols(permvec);

   dynamic_cast<SparseMatrix&>(*Bmat).permuteCols(permvec);

   for (auto& it : children)
      it->permuteLinkingVarsChild(permvec);
}

void DistributedMatrix::permuteLinkingVarsChild(const std::vector<unsigned int>& permvec) {
   assert(hasSparseMatrices());

   dynamic_cast<SparseMatrix&>(*Amat).permuteCols(permvec);

   assert(children.empty());
}

void DistributedMatrix::permuteLinkingCons(const std::vector<unsigned int>& permvec) {
   assert(hasSparseMatrices());
   if (Blmat)
      dynamic_cast<SparseMatrix&>(*Blmat).permuteRows(permvec);

   for (auto& it : children)
      it->permuteLinkingCons(permvec);
}


void DistributedMatrix::updateTransposed() const {
   assert(hasSparseMatrices());
   dynamic_cast<const SparseMatrix&>(*Amat).updateTransposed();
   dynamic_cast<const SparseMatrix&>(*Bmat).updateTransposed();
   dynamic_cast<const SparseMatrix&>(*Blmat).updateTransposed();

   for (const auto& it : children)
      it->updateTransposed();
}

/* check whether root node date is same in all processes
 *
 * todo: check this
 * todo: make better use of std::vector and iterators
 * root node data is Amat (empty), Bmat and Blmat of root node. Children not checked.
 */
bool DistributedMatrix::isRootNodeInSync() const {
   assert(hasSparseMatrices());
   bool in_sync = true;

   assert(Amat);
   assert(Bmat);
   assert(Blmat);

   /* no need to check if not distributed or not in root node */
   if (!iAmDistrib || children.empty())
      return in_sync;

   const int my_rank = PIPS_MPIgetRank(mpiComm);
   const auto& bmat_sp = dynamic_cast<const SparseMatrix&>(*Bmat);
   const auto& blmat_sp = dynamic_cast<const SparseMatrix&>(*Blmat);

   /* since we are in root node Amat should be empty */
   assert(Amat->numberOfNonZeros() == 0);

   /* static storage */
   const int lenght_entries_bmat = bmat_sp.getStorage().len;
   const int length_columns_bmat = bmat_sp.getStorage().len;
   const int lenght_rowoffest_bmat = bmat_sp.getStorage().m + 1;

   const int lenght_entries_blmat = blmat_sp.getStorage().len;
   const int length_columns_blmat = blmat_sp.getStorage().len;
   const int lenght_rowoffest_blmat = blmat_sp.getStorage().m + 1;

   const long long count_row_cols =
      length_columns_bmat + lenght_rowoffest_bmat + length_columns_blmat + lenght_rowoffest_blmat;
   const long long count_entries = lenght_entries_bmat + lenght_entries_blmat;

   assert(count_row_cols < std::numeric_limits<int>::max());
   assert(count_entries < std::numeric_limits<int>::max());

   std::vector<double> sendbuf_entries(count_entries, 0.0);
   std::vector<double> recvbuf_entries(count_entries, 0.0);

   std::vector<int> sendbuf_row_col(count_row_cols, 0);
   std::vector<int> recvbuf_row_col(count_row_cols, 0);

   /* fill Bmat into send buffers */
   const double* M = bmat_sp.getStorage().M;
   const int* krowM = bmat_sp.getStorage().krowM;
   const int* jColM = bmat_sp.getStorage().jcolM;

   std::copy(M, M + lenght_entries_bmat, sendbuf_entries.begin());

   std::copy(krowM, krowM + lenght_rowoffest_bmat, sendbuf_row_col.begin());
   std::copy(jColM, jColM + lenght_entries_bmat, sendbuf_row_col.begin() + lenght_rowoffest_bmat);

   /* fill Blmat into send buffers */
   const double* Ml = blmat_sp.getStorage().M;
   const int* krowMl = blmat_sp.getStorage().krowM;
   const int* jColMl = blmat_sp.getStorage().jcolM;

   std::copy(Ml, Ml + lenght_entries_blmat, sendbuf_entries.begin() + lenght_entries_bmat);
   std::copy(krowMl, krowMl + lenght_rowoffest_blmat,
      sendbuf_row_col.begin() + lenght_rowoffest_bmat + lenght_entries_bmat);
   std::copy(jColMl, jColMl + lenght_entries_blmat,
      sendbuf_row_col.begin() + lenght_rowoffest_bmat + lenght_entries_bmat + lenght_rowoffest_blmat);

   /* Reduce Bmat and Blmat buffers */
   MPI_Allreduce(&sendbuf_entries[0], &recvbuf_entries[0], static_cast<int>(count_entries), MPI_DOUBLE, MPI_MAX,
      mpiComm);

   MPI_Allreduce(&sendbuf_row_col[0], &recvbuf_row_col[0], static_cast<int>(count_row_cols), MPI_INT, MPI_MAX,
      mpiComm);

   /* check recvbuf_entries */
   for (int i = 0; i < count_entries; ++i) {
      if (!PIPSisEQ(sendbuf_entries[i], recvbuf_entries[i])) {
         /* someone else had a higher value here */
         if (my_rank == 0)
            std::cout << "matrix entries out of sync\n";
         in_sync = false;
         break;
      }
   }

   for (int i = 0; i < count_row_cols; ++i) {
      if (!PIPSisEQ(sendbuf_row_col[i], recvbuf_row_col[i])) {
         /* someone else had a higher value here */
         if (my_rank == 0)
            std::cout << "matrix indices (col or row) out of sync\n";
         in_sync = false;
      }
   }

   /* if stoch mat has dynamic storage also check that */
   if (bmat_sp.hasDynamicStorage() || blmat_sp.hasDynamicStorage()) {
      assert(bmat_sp.hasDynamicStorage());
      assert(blmat_sp.hasDynamicStorage());

      const SparseStorageDynamic& Bmat_dyn = bmat_sp.getStorageDynamic();
      const SparseStorageDynamic& Blmat_dyn = blmat_sp.getStorageDynamic();

      /* dynamic storage */
      int bmat_dyn_len = 0;
      for (int i = 0; i < Bmat_dyn.n_rows(); ++i)
         bmat_dyn_len += (Bmat_dyn.getRowPtr(i).end - Bmat_dyn.getRowPtr(i).start);

      int blmat_dyn_len = 0;
      for (int i = 0; i < Blmat_dyn.n_rows(); ++i)
         blmat_dyn_len += (Blmat_dyn.getRowPtr(i).end - Blmat_dyn.getRowPtr(i).start);

      const int lenght_entries_bmat_dynamic = bmat_dyn_len;
      const int length_columns_bmat_dynamic = bmat_dyn_len;
      const int lenght_rowoffest_bmat_dynamic = Bmat_dyn.n_rows() + 1;

      const int lenght_entries_blmat_dynamic = blmat_dyn_len;
      const int length_columns_blmat_dynamic = blmat_dyn_len;
      const int lenght_rowoffest_blmat_dynamic = Blmat_dyn.n_rows() + 1;

      const long long count_row_cols_dyn =
         length_columns_bmat_dynamic + 2 * lenght_rowoffest_bmat_dynamic + length_columns_blmat_dynamic +
            2 * lenght_rowoffest_blmat_dynamic;
      const long long count_entries_dyn = lenght_entries_bmat_dynamic + lenght_entries_blmat_dynamic;

      assert(count_row_cols_dyn < std::numeric_limits<int>::max());
      assert(count_entries_dyn < std::numeric_limits<int>::max());

      std::vector<double> sendbuf_entries_dynamic(count_entries_dyn, 0.0);
      std::vector<double> recvbuf_entries_dynamic(count_entries_dyn, 0.0);

      std::vector<int> sendbuf_row_col_dynamic(count_row_cols_dyn, 0);
      std::vector<int> recvbuf_row_coldynamic(count_row_cols_dyn, 0);;

      /* fill Bmat into send buffers */
      const double* M_dyn = Bmat_dyn.getMat();
      const int* jColM_dyn = Bmat_dyn.getJcolM();

      int count_entries_2 = 0;
      int count_row_col = 0;

      /* entries Bmat into double array */
      for (int i = 0; i < Bmat_dyn.n_rows(); ++i) {
         for (int j = Bmat_dyn.getRowPtr(i).start; j < Bmat_dyn.getRowPtr(i).end; ++j) {
            sendbuf_entries_dynamic[count_entries_2] = M_dyn[j];
            count_entries_2++;
         }
      }
      assert(count_entries_2 == lenght_entries_bmat_dynamic);

      /* row pointers Bmat into int array */
      for (int i = 0; i < lenght_rowoffest_bmat_dynamic; ++i) {
         sendbuf_row_col_dynamic[count_row_col] = Bmat_dyn.getRowPtr()->start;
         sendbuf_row_col_dynamic[count_row_col + 1] = Bmat_dyn.getRowPtr()->end;
         count_row_col += 2;
      }
      assert(count_row_col == 2 * lenght_rowoffest_bmat_dynamic);

      /* col indices of Bmat into int array */
      for (int i = 0; i < Bmat_dyn.n_rows(); ++i) {
         for (int j = Bmat_dyn.getRowPtr(i).start; j < Bmat_dyn.getRowPtr(i).end; ++j) {
            sendbuf_row_col_dynamic[count_row_col] = jColM_dyn[j];
            count_row_col++;
         }
      }
      assert(count_row_col == 2 * lenght_rowoffest_bmat_dynamic + length_columns_bmat_dynamic);

      /* fill Blmat into send buffers */
      const double* Ml_dyn = Blmat_dyn.getMat();
      const int* jColMl_dyn = Blmat_dyn.getJcolM();

      /* entries Blmat into double array */
      for (int i = 0; i < Blmat_dyn.n_rows(); ++i) {
         for (int j = Blmat_dyn.getRowPtr(i).start; j < Blmat_dyn.getRowPtr(i).end; ++j) {
            sendbuf_entries_dynamic[count_entries_2] = Ml_dyn[j];
            count_entries_2++;
         }
      }
      assert(count_entries_2 == lenght_entries_bmat_dynamic + lenght_entries_blmat_dynamic);

      /* row pointers Blmat into int array */
      for (int i = 0; i < lenght_rowoffest_blmat_dynamic; ++i) {
         assert(2 * lenght_rowoffest_bmat_dynamic + lenght_entries_bmat_dynamic + 2 * i + 1 < count_row_cols_dyn);
         sendbuf_row_col_dynamic[count_row_col] = Blmat_dyn.getRowPtr()->start;
         sendbuf_row_col_dynamic[count_row_col + 1] = Blmat_dyn.getRowPtr()->end;
         count_row_col += 2;
      }
      assert(count_row_col ==
         2 * lenght_rowoffest_bmat_dynamic + length_columns_bmat_dynamic + 2 * lenght_rowoffest_blmat_dynamic);

      /* col indices of Bmat into int array */
      for (int i = 0; i < Blmat_dyn.n_rows(); ++i) {
         for (int j = Blmat_dyn.getRowPtr(i).start; j < Blmat_dyn.getRowPtr(i).end; ++j) {
            sendbuf_row_col_dynamic[count_row_col] = jColMl_dyn[j];
            count_row_col++;
         }
      }
      assert(count_row_col ==
         2 * lenght_rowoffest_bmat_dynamic + length_columns_bmat_dynamic + 2 * lenght_rowoffest_blmat_dynamic +
            length_columns_blmat_dynamic);

      /* Reduce Bmat and Blmat buffers */
      MPI_Allreduce(&sendbuf_entries_dynamic[0], &recvbuf_entries_dynamic[0], static_cast<int>(count_entries_dyn),
         MPI_DOUBLE, MPI_MAX, mpiComm);

      MPI_Allreduce(&sendbuf_row_col_dynamic[0], &recvbuf_row_coldynamic[0], static_cast<int>(count_row_cols_dyn),
         MPI_INT, MPI_MAX, mpiComm);

      /* check recvbuf_entries */
      for (int i = 0; i < count_entries_dyn; ++i) {
         if (!PIPSisEQ(sendbuf_entries_dynamic[i], recvbuf_entries_dynamic[i])) {
            /* someone else had a higher value here */
            if (my_rank == 0)
               std::cout << "matrix entries in dynamic storage out of sync\n";
            in_sync = false;
         }
      }
      for (int i = 0; i < count_row_cols_dyn; ++i) {
         if (!PIPSisEQ(sendbuf_row_col_dynamic[i], recvbuf_row_coldynamic[i])) {
            /* someone else had a higher value here */
            if (my_rank == 0)
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
int DistributedMatrix::appendRow(const DistributedMatrix& matrix_row, int child, int row, bool linking) {
   assert(hasSparseMatrices());
   // todo: check that matrix is in correct format
   assert(matrix_row.children.size() == children.size());
   assert(!children.empty());
   assert(-1 <= child && child <= (int) children.size());

   int index_row;

   // append row to all matrices necessary
   // todo maybe this can be done nicer - maybe we can just recursively call some method also on the dummies
   if (linking) {
      index_row = dynamic_cast<SparseMatrix&>(*Blmat).appendRow(dynamic_cast<const SparseMatrix&>(*matrix_row.Blmat),
         row);

      for (unsigned int i = 0; i < children.size(); ++i) {
         if (!children[i]->is_a(kStochGenDummyMatrix)) {
            assert(!matrix_row.children[i]->is_a(kStochGenDummyMatrix));
            dynamic_cast<SparseMatrix&>(*children[i]->Blmat).appendRow(
               dynamic_cast<const SparseMatrix&>(*matrix_row.children[i]->Blmat), row);
         }
      }
   } else {
      if (child != -1) {
         index_row = dynamic_cast<SparseMatrix&>(*children[child]->Amat).appendRow(
            dynamic_cast<const SparseMatrix&>(*matrix_row.children[child]->Amat), row);
#ifndef NDEBUG
         const int index_row1 = dynamic_cast<SparseMatrix&>(*children[child]->Bmat).appendRow(
            dynamic_cast<const SparseMatrix&>(*matrix_row.children[child]->Bmat), row);
#else
         dynamic_cast<SparseMatrix&>(*children[child]->Bmat).appendRow( dynamic_cast<SparseMatrix&>(*matrix_row.children[child]->Bmat), row );
#endif
         assert(index_row1 == index_row);
      } else
         index_row = dynamic_cast<SparseMatrix&>(*Bmat).appendRow(dynamic_cast<const SparseMatrix&>(*matrix_row.Bmat),
            row);
   }

   return index_row;
};

/* y += alpha RowAt(child, row, linking) */
void
DistributedMatrix::axpyWithRowAt(double alpha, DistributedVector<double>* y, SimpleVector<double>* y_linking, int child,
   int row, bool linking) const {
   assert(hasSparseMatrices());
   assert(y);
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(y->children.size() == children.size());

   /* go through all available children and calculate y += alpha * rowAt(row) */
   if (linking) {
      assert(Blmat);
      if (y_linking)
         dynamic_cast<const SparseMatrix&>(*Blmat).axpyWithRowAt(alpha, *y_linking, row);
      else {
         assert(y->first);
         dynamic_cast<const SparseMatrix&>(*Blmat).axpyWithRowAt(alpha, dynamic_cast<SimpleVector<double>&>(*y->first),
            row);
      }

      for (unsigned int i = 0; i < children.size(); ++i) {
         if (!children[i]->is_a(kStochGenDummyMatrix)) {
            assert(children[i]->Blmat);
            assert(y->children[i]->first);
            dynamic_cast<const SparseMatrix&>(*children[i]->Blmat).axpyWithRowAt(alpha,
               dynamic_cast<SimpleVector<double>&>(*y->children[i]->first), row);
         }
      }
   } else {
      if (child == -1) {
         assert(Bmat);
         if (y_linking)
            dynamic_cast<const SparseMatrix&>(*Bmat).axpyWithRowAt(alpha, *y_linking, row);
         else {
            assert(y->first);
            dynamic_cast<const SparseMatrix&>(*Bmat).axpyWithRowAt(alpha,
               dynamic_cast<SimpleVector<double>&>(*y->first), row);
         }
      } else {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);

         assert(y->children[child]->first);
         dynamic_cast<const SparseMatrix&>(*children[child]->Bmat).axpyWithRowAt(alpha,
            dynamic_cast<SimpleVector<double>&>(*y->children[child]->first), row);

         if (y_linking)
            dynamic_cast<const SparseMatrix&>(*children[child]->Amat).axpyWithRowAt(alpha, *y_linking, row);
         else {
            assert(y->first);
            dynamic_cast<const SparseMatrix&>(*children[child]->Amat).axpyWithRowAt(alpha,
               dynamic_cast<SimpleVector<double>&>(*y->first), row);
         }
      }
   }
}

void
DistributedMatrix::axpyWithRowAtPosNeg(double alpha, DistributedVector<double>* y_pos, SimpleVector<double>* y_link_pos,
   DistributedVector<double>* y_neg, SimpleVector<double>* y_link_neg, int child, int row, bool linking) const {
   assert(hasSparseMatrices());
   assert(y_pos && y_neg);
   assert((y_link_neg && y_link_pos) || (!y_link_neg && !y_link_pos));
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(y_neg->children.size() == children.size());
   assert(y_pos->children.size() == children.size());

   /* go through all available children and calculate y += alpha * rowAt(row) */
   if (linking) {
      assert(Blmat);
      if (y_link_pos)
         dynamic_cast<const SparseMatrix&>(*Blmat).axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
      else {
         assert(y_pos->first);
         assert(y_neg->first);
         dynamic_cast<const SparseMatrix&>(*Blmat).axpyWithRowAtPosNeg(alpha,
            dynamic_cast<SimpleVector<double>&>(*y_pos->first),
            dynamic_cast<SimpleVector<double>&>(*y_neg->first), row);
      }

      for (unsigned int i = 0; i < children.size(); ++i) {
         if (!children[i]->is_a(kStochGenDummyMatrix)) {
            assert(children[i]->Blmat);
            assert(y_pos->children[i]->first);
            assert(y_neg->children[i]->first);
            dynamic_cast<const SparseMatrix&>(*children[i]->Blmat).axpyWithRowAtPosNeg(alpha,
               dynamic_cast<SimpleVector<double>&>(*y_pos->children[i]->first),
               dynamic_cast<SimpleVector<double>&>(*y_neg->children[i]->first),
               row);
         }
      }
   } else {
      if (child == -1) {
         assert(Bmat);
         if (y_link_pos)
            dynamic_cast<const SparseMatrix&>(*Bmat).axpyWithRowAtPosNeg(alpha, *y_link_pos, *y_link_neg, row);
         else {
            assert(y_pos->first);
            assert(y_neg->first);
            dynamic_cast<const SparseMatrix&>(*Bmat).axpyWithRowAtPosNeg(alpha,
               dynamic_cast<SimpleVector<double>&>(*y_pos->first),
               dynamic_cast<SimpleVector<double>&>(*y_neg->first), row);
         }
      } else {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);

         assert(y_pos->children[child]->first);
         assert(y_neg->children[child]->first);
         dynamic_cast<const SparseMatrix&>(*children[child]->Bmat).axpyWithRowAtPosNeg(alpha,
            dynamic_cast<SimpleVector<double>&>(*y_pos->children[child]->first),
            dynamic_cast<SimpleVector<double>&>(*y_neg->children[child]->first), row);

         if (y_link_pos)
            dynamic_cast<const SparseMatrix&>(*children[child]->Amat).axpyWithRowAtPosNeg(alpha, *y_link_pos,
               *y_link_neg, row);
         else {
            assert(y_pos->first);
            assert(y_neg->first);
            dynamic_cast<const SparseMatrix&>(*children[child]->Amat).axpyWithRowAtPosNeg(alpha,
               dynamic_cast<SimpleVector<double>&>(*y_pos->first), dynamic_cast<SimpleVector<double>&>(*y_neg->first),
               row);
         }
      }
   }
}

double
DistributedMatrix::localRowTimesVec(const DistributedVector<double>& vec, int child, int row, bool linking) const {
   assert(hasSparseMatrices());
   assert(-1 <= child && child < static_cast<int>(children.size()));
   assert(vec.children.size() == children.size());

   double res = 0.0;

   /* go through all available children and multiply the first times row in submatrix */
   if (linking) {
      assert(Blmat);
      assert(vec.first);
      res += dynamic_cast<const SparseMatrix&>(*Blmat).localRowTimesVec(
         dynamic_cast<const SimpleVector<double>&>(*vec.first), row);

      for (unsigned int i = 0; i < children.size(); ++i) {
         if (!children[i]->is_a(kStochGenDummyMatrix)) {
            assert(children[i]->Blmat);
            assert(vec.children[i]->first);
            res += dynamic_cast<const SparseMatrix&>(*children[i]->Blmat).localRowTimesVec(
               dynamic_cast<const SimpleVector<double>&>(*vec.children[i]->first), row);
         }
      }
   } else {
      if (child == -1) {
         assert(Bmat);
         assert(vec.first);
         res += dynamic_cast<const SparseMatrix&>(*Bmat).localRowTimesVec(
            dynamic_cast<const SimpleVector<double>&>(*vec.first), row);
      } else {
         assert(children[child]->Amat);
         assert(children[child]->Bmat);
         assert(vec.first);
         assert(vec.children[child]->first);
         res += dynamic_cast<const SparseMatrix&>(*children[child]->Amat).localRowTimesVec(
            dynamic_cast<const SimpleVector<double>&>(*vec.first),
            row);
         res += dynamic_cast<const SparseMatrix&>(*children[child]->Bmat).localRowTimesVec(
            dynamic_cast<const SimpleVector<double>&>(*vec.children[child]->first), row);
      }
   }

   return res;
}

// TODO specify border and left from DistributedQP...
std::unique_ptr<BorderedMatrix> DistributedMatrix::raiseBorder(int m_conss, int n_vars) {
#ifndef NDEBUG
   const auto[m_link, n_link] = Blmat->n_rows_columns();
   assert(m_conss <= m_link && n_vars <= n_link);
#endif

   std::unique_ptr<SparseMatrix> A_left{dynamic_cast<SparseMatrix&>(*Bmat).shaveLeft(n_vars)};

   std::unique_ptr<SparseMatrix> Bl_left_top{dynamic_cast<SparseMatrix&>(*Blmat).shaveLeft(n_vars)};
   std::unique_ptr<GeneralMatrix> bottom_left_block = Bl_left_top->shaveBottom(m_conss);

   std::unique_ptr<GeneralMatrix> Bl_right_bottom = Blmat->shaveBottom(m_conss);

   auto border_bottom = std::make_unique<StripMatrix>(false, std::move(Bl_right_bottom), nullptr, mpiComm);
   auto border_left = std::make_unique<StripMatrix>(true, std::move(A_left), std::move(Bl_left_top), mpiComm);

   for (auto& it : children)
      it->shaveBorder(m_conss, n_vars, border_left.get(), border_bottom.get());
   border_left->recomputeNonzeros();
   border_bottom->recomputeNonzeros();

   m -= m_conss;
   n -= n_vars;

   auto bordered_matrix = std::make_unique<BorderedMatrix>(std::dynamic_pointer_cast<DistributedMatrix>(shared_from_this()), std::move(border_left),
      std::move(border_bottom), std::move(bottom_left_block), mpiComm);

   assert(m >= 0 && n >= 0);

   return bordered_matrix;
}

void DistributedMatrix::shaveBorder(int m_conss, int n_vars, StripMatrix* border_left, StripMatrix* border_bottom) {
   if (Bmat->is_a(kDistributedMatrix)) {
      assert(amatEmpty());
      border_left->addChild(std::make_unique<StripMatrix>(true, dynamic_cast<DistributedMatrix&>(*Bmat).shaveLeftBorder(n_vars), nullptr,
         mpiComm));
      border_bottom->addChild(std::make_unique<StripMatrix>(false, Blmat->shaveBottom(m_conss), nullptr, mpiComm));
   } else {
      assert(hasSparseMatrices());
      assert(children.empty());
      assert(PIPS_MPIgetSize(mpiComm) == 1);

      std::unique_ptr<SparseMatrix> border_a_mat = dynamic_cast<SparseMatrix&>(*Amat).shaveLeft(n_vars);
      std::unique_ptr<GeneralMatrix> border_bl_mat = dynamic_cast<SparseMatrix&>(*Blmat).shaveBottom(m_conss);

      std::unique_ptr<StripMatrix> border_left_child = std::make_unique<StripMatrix>(true, std::move(border_a_mat), nullptr,
         mpiComm);
      std::unique_ptr<StripMatrix> border_bottom_child = std::make_unique<StripMatrix>(false, std::move(border_bl_mat), nullptr,
         mpiComm);

      border_left->addChild(std::move(border_left_child));
      border_bottom->addChild(std::move(border_bottom_child));
   }
}

std::unique_ptr<StripMatrix> DistributedMatrix::shaveLeftBorder(int n_vars) {
   assert(!children.empty());
   assert(hasSparseMatrices());
   assert(amatEmpty());

   std::unique_ptr<SparseMatrix> border_b_mat{dynamic_cast<SparseMatrix&>(*Bmat).shaveLeft(n_vars)};
   std::unique_ptr<SparseMatrix> border_bl_mat{dynamic_cast<SparseMatrix&>(*Blmat).shaveLeft(n_vars)};

   auto border = std::make_unique<StripMatrix>(true, std::move(border_b_mat), std::move(border_bl_mat), mpiComm);

   for (auto& child : children) {
      border->addChild(child->shaveLeftBorderChild(n_vars));
   }

   border->recomputeNonzeros();
   return border;
}

std::unique_ptr<StripMatrix> DistributedMatrix::shaveLeftBorderChild(int n_vars) {
   assert(children.empty());

   if (Bmat->is_a(kDistributedMatrix)) {
      assert(amatEmpty());
      return std::make_unique<StripMatrix>(true,
         std::unique_ptr<StripMatrix>(dynamic_cast<DistributedMatrix&>(*Bmat).shaveLeftBorder(n_vars)), nullptr,
         mpiComm);
   } else
      return std::make_unique<StripMatrix>(true, dynamic_cast<SparseMatrix&>(*Amat).shaveLeft(n_vars), nullptr, mpiComm);
}

std::unique_ptr<StripMatrix> DistributedMatrix::shaveLinkingConstraints(unsigned int n_conss) {
   assert(hasSparseMatrices());

   std::unique_ptr<GeneralMatrix> border_bl_mat = Blmat->shaveBottom(static_cast<int>(n_conss));
   auto border = std::make_unique<StripMatrix>(false, std::move(border_bl_mat), nullptr, mpiComm);

   if (children.empty())
      assert(PIPS_MPIgetSize(mpiComm) == 1);
   for (auto& child : children) {
      std::unique_ptr<StripMatrix> border_child{child->shaveLinkingConstraints(n_conss)};
      border->addChild(std::move(border_child));
   }

   border->recomputeNonzeros();
   return border;
}

void DistributedMatrix::splitMatrix(const std::vector<int>& twolinks_start_in_block,
   const std::vector<unsigned int>& map_blocks_children,
   unsigned int n_links_in_root, const std::vector<MPI_Comm>& child_comms) {
   const unsigned int n_curr_children = children.size();

   assert(hasSparseMatrices());
   assert(n_curr_children == map_blocks_children.size());
   assert(n_curr_children == twolinks_start_in_block.size());
   assert(twolinks_start_in_block.back() == 0);

   auto[m_links_left, nBl] = Blmat->n_rows_columns();
   assert(std::accumulate(twolinks_start_in_block.begin(), twolinks_start_in_block.end(), 0) <= m_links_left);

   const unsigned int n_new_children = getNDistinctValues(map_blocks_children);
   std::vector<std::shared_ptr<DistributedMatrix>> new_children(n_new_children);

   std::unique_ptr<StripMatrix> Blmat_new{shaveLinkingConstraints(n_links_in_root)};

   Blmat_new->combineChildrenInNewChildren(map_blocks_children, child_comms);

   std::unique_ptr<GeneralMatrix> Blmat_leftover = std::move(Blmat);

#ifndef NDEBUG
   int n_child_links_sum{0};
   const unsigned int n_links_orig = m_links_left;
#endif
   m_links_left -= static_cast<int>(n_links_in_root);
   /* for each future new child collect its children and add them to the new child
    * then shave off the linking constraints that stay at the new child's level
    */
   unsigned int m_links_so_far{0};
   unsigned int begin_curr_child_blocks{0};
   unsigned int end_curr_child_blocks{0};
   for (unsigned int i = 0; i < n_new_children; ++i) {
      while (end_curr_child_blocks != (n_curr_children - 1) &&
         map_blocks_children[end_curr_child_blocks] == map_blocks_children[end_curr_child_blocks + 1])
         ++end_curr_child_blocks;

      const int n_links_for_child = std::accumulate(twolinks_start_in_block.begin() + begin_curr_child_blocks,
         twolinks_start_in_block.begin() + end_curr_child_blocks, 0);
      const unsigned int n_blocks_for_child = end_curr_child_blocks - begin_curr_child_blocks + 1;

#ifndef NDEBUG
      n_child_links_sum += n_links_for_child;
#endif
      /* combine children in new DistributedMatrix Bmat_loc */
      /* create root node with only Blmat */
      std::unique_ptr<GeneralMatrix> new_Blmat_leftover = Blmat_leftover->shaveBottom(m_links_left - n_links_for_child);

      if (child_comms[i] == MPI_COMM_NULL) {
         Blmat_leftover = nullptr;
      }

      std::unique_ptr<DistributedMatrix> Bmat_loc{
         (child_comms[i] == MPI_COMM_NULL) ? nullptr : new DistributedMatrix(std::make_unique<SparseMatrix>(0, 0, 0),
            std::make_unique<SparseMatrix>(0, nBl, 0), std::move(Blmat_leftover), child_comms[i], false, true)};

      Blmat_leftover = std::move(new_Blmat_leftover);

      /* shave off empty two link part from respective children and add them to the new root/remove them from the old root */
      for (unsigned int j = 0; j < n_blocks_for_child; ++j) {
         assert(m_links_left >= n_links_for_child);

         std::shared_ptr<DistributedMatrix> child = children.front();
         children.erase(children.begin());

         if (child_comms[i] == MPI_COMM_NULL)
            assert(child->mpiComm == MPI_COMM_NULL);

         if (child->mpiComm != MPI_COMM_NULL) {
            dynamic_cast<SparseMatrix&>(*child->Blmat).dropNEmptyRowsTop(static_cast<int>(m_links_so_far));
            dynamic_cast<SparseMatrix&>(*child->Blmat).dropNEmptyRowsBottom(m_links_left - n_links_for_child);
         }

#ifndef NDEBUG
         if (child->mpiComm != MPI_COMM_NULL) {
            assert(child->Blmat->n_rows() == n_links_for_child);
         }
#endif
         if (Bmat_loc)
            Bmat_loc->AddChild(child);
      }
      if (Bmat_loc)
         Bmat_loc->recomputeSize();

      /* create child holding the new Bmat_loc and it's Blmat part */
      if (child_comms[i] == MPI_COMM_NULL) {
         assert(Blmat_new->children[i]->is_a(kStringGenDummyMatrix));
         Blmat_new->children[i] = nullptr;
      } else {
         assert(Blmat_new->children[i]->is_a(kStripMatrix));
         assert(!Blmat_new->children[i]->is_a(kStringGenDummyMatrix));
      }

      new_children[i].reset(
         (child_comms[i] != MPI_COMM_NULL) ? new DistributedMatrix(std::make_unique<SparseMatrix>(0, 0, 0),
            std::move(Bmat_loc), std::move(Blmat_new->children[i]), child_comms[i], true, false) : new StochGenDummyMatrix());

      ++end_curr_child_blocks;
      begin_curr_child_blocks = end_curr_child_blocks;
      m_links_left -= n_links_for_child;
      m_links_so_far += n_links_for_child;
   }

   assert(n_child_links_sum + n_links_in_root == n_links_orig);
   assert(children.empty());
   assert(m_links_left == 0);

   /* exchange children and recompute sizes */
   children.insert(children.end(), new_children.begin(), new_children.end());

   Blmat_new->children.clear();
   Blmat = std::move(Blmat_new->first);
   Blmat_new->first = nullptr;

   recomputeSize();
}
