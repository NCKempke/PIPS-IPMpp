#ifndef DistributedMatrix_H
#define DistributedMatrix_H

#include "Vector.hpp"
#include "AbstractMatrix.h"
#include "SparseMatrix.h"
#include "BorderedMatrix.h"
#include "StripMatrix.h"
#include "DistributedTree.h"

#include "mpi.h"

#include <vector>

class DistributedMatrix : public GeneralMatrix {
protected:
   DistributedMatrix() = default;
public:

   DistributedMatrix(std::unique_ptr<GeneralMatrix> Amat, std::unique_ptr<GeneralMatrix> Bmat, std::unique_ptr<GeneralMatrix> Blmat, MPI_Comm mpiComm_, bool inner_leaf = false, bool inner_root = false);

   /** Constructs a matrix having local A and B blocks having the sizes and number of nz specified by
    *  A_m, A_n, A_nnz and B_m, B_n, B_nnz.
    *  Also sets the global sizes to 'global_m' and 'global_n'.
    *  The matrix that will be created  has no children, just local data.
    */
   DistributedMatrix(long long global_m, long long global_n, int A_m, int A_n, int A_nnz, int B_m, int B_n, int B_nnz, MPI_Comm mpiComm_);

   /** Constructs a matrix with local A, B, and Bl (linking constraints) blocks having the sizes and number of nz specified by
       A_m, A_n, A_nnz, B_m, B_n, B_nnz, and Bl_m, Bl_n, Bl_nnz. Otherwise, identical to the above constructor */
   DistributedMatrix(long long global_m, long long global_n, int A_m, int A_n, int A_nnz, int B_m, int B_n, int B_nnz, int Bl_m, int Bl_n, int Bl_nnz,
         MPI_Comm mpiComm_);

   // constructor for combining scenarios
   ~DistributedMatrix() override = default;

   using GeneralMatrix::cloneFull;
   using GeneralMatrix::cloneEmptyRows;
   [[nodiscard]] std::unique_ptr<GeneralMatrix> cloneEmptyRows(bool switchToDynamicStorage) const override;
   [[nodiscard]] std::unique_ptr<GeneralMatrix> cloneFull(bool switchToDynamicStorage) const override;

   virtual void AddChild(const std::shared_ptr<DistributedMatrix>& child);

   std::vector<std::shared_ptr<DistributedMatrix>> children;
   std::unique_ptr<GeneralMatrix> Amat{};
   std::unique_ptr<GeneralMatrix> Bmat{};
   std::unique_ptr<GeneralMatrix> Blmat{};

   long long m{-1};
   long long n{-1};
   MPI_Comm mpiComm{MPI_COMM_NULL};
   int iAmDistrib{false};

   /* is this matrix an inner matrix of the matrix hierarchy - if not, then its children hold the local Amat, Bmat and Blmat */
   const bool inner_leaf{false};
   const bool inner_root{false};
private:
   [[nodiscard]] virtual bool hasSparseMatrices() const;

   /** trans mult method for children with linking constraints */
   virtual void transMult2(double beta, DistributedVector<double>& y, double alpha, const DistributedVector<double>& x, const Vector<double>* xvecl) const;

   virtual void mult2(double beta, DistributedVector<double>& y, double alpha, const DistributedVector<double>& x, Vector<double>* yparentl_) const;

   /** column scale method for children */
   virtual void columnScale2(const Vector<double>& vec);

   /** row scale method for children */
   virtual void rowScale2(const Vector<double>& vec, const Vector<double>* linkingvec);

   virtual void getNnzPerRow(Vector<int>& nnzVec, Vector<int>* linkParent) const;

   virtual void getNnzPerCol(Vector<int>& nnzVec, Vector<int>* linkParent) const;

   void sum_transform_columns(Vector<double>& result, const std::function<double(const double&)>& transform, Vector<double>* link_parent) const;

   virtual void addRowSums(Vector<double>& sumVec, Vector<double>* linkParent) const;
   virtual void addColSums(Vector<double>& sumVec, Vector<double>* linkParent) const;

   virtual void initTransposedChild(bool dynamic) const;
   virtual void initStaticStorageFromDynamic(const Vector<int>& rowNnzVec, const Vector<int>& colNnzVec, const Vector<int>* rowLinkVec,
         const Vector<int>* colParentVec);

   virtual void permuteLinkingVarsChild(const std::vector<unsigned int>& permvec);
   virtual void getLinkVarsNnzChild(std::vector<int>& vec) const;

public:
   virtual void updateTransposed() const;

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override;
   [[nodiscard]] long long n_rows() const override;
   [[nodiscard]] long long n_columns() const override;

   /** The actual number of structural non-zero elements in this sparse
    *  matrix. This includes so-called "accidental" zeros, elements that
    *  are treated as non-zero even though their value happens to be zero.
    */
   [[nodiscard]] int numberOfNonZeros() const override;

   [[nodiscard]] int is_a(int matType) const override;

   void atPutDense(int, int, const double*, int, int, int) override { assert("Not implemented" && 0); };
   void fromGetDense(int, int, double*, int, int, int) const override { assert("Not implemented" && 0); };

   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void symmetricScale(const Vector<double>&) override { assert("Not implemented" && 0); };
   void scalarMult(double num) override;

   void fromGetSpRow(int, int, double[], int, int[], int&, int, int&) const override { assert("Not implemented" && 0); };
   void atPutSubmatrix(int, int, const AbstractMatrix&, int, int, int, int) override { assert("Not implemented" && 0); };
   void atPutSpRow(int, const double[], int, const int[], int&) override { assert("Not implemented" && 0); };
   void putSparseTriple(const int[], int, const int[], const double[], int&) override { assert("Not implemented" && 0); };

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;

   /** y = beta * y + alpha * this * x */
   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;

   virtual void getLinkVarsNnz(std::vector<int>& vec) const;

   void write_to_stream(std::ostream&) const override { assert("Not implemented" && 0); };
   virtual void write_to_streamDense(std::ostream& out, int offset) const;
   virtual void write_to_streamDenseBordered(const StripMatrix& border, std::ostream& out, int offset) const;
   void write_to_streamDense(std::ostream& out) const override { write_to_streamDense(out, 0); };
   void writeDashedLineToStream(std::ostream& out) const override { writeDashedLineToStream(out, 0); };
   virtual void writeDashedLineToStream(std::ostream& out, int offset) const;

   void writeMPSformatRows(std::ostream& out, int rowType, const Vector<double>* irhs) const override;

   /** initialize (dynamic) transposed matrices for A, B, Bl */
   virtual void initTransposed() const {
      initTransposed(false);
   }
   virtual void initTransposed(bool dynamic) const;
   virtual void deleteTransposed() const;

   void atPutDiagonal(int, const Vector<double>&) override { assert("Not implemented" && 0); };
   void atAddDiagonal(int, const Vector<double>&) override { assert("Not implemented" && 0); };
   void fromGetDiagonal(int, Vector<double>&) const override { assert("Not implemented" && 0); };
   void matTransDMultMat(const Vector<double>&, SymmetricMatrix**) const override { assert("Not implemented" && 0); };
   void matTransDinvMultMat(const Vector<double>&, SymmetricMatrix**) const override { assert("Not implemented" && 0); };

   void getNnzPerRow(Vector<int>& nnzVec) const override {
      getNnzPerRow(nnzVec, nullptr);
   };

   void getNnzPerCol(Vector<int>& nnzVec) const override {
      getNnzPerCol(nnzVec, nullptr);
   };

   /** fill vector with absolute minimum/maximum value of each row */
   void getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) const override;
   /** fill vector with absolute minimum/maximum value of each column */
   void getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) const override;

   void sum_transform_rows(Vector<double>& result, const std::function<double(const double&)>& transform) const override;
   void sum_transform_columns(Vector<double>& result, const std::function<double(const double&)>& transform) const override {
      sum_transform_columns(result, transform, nullptr);
   };

   void addRowSums(Vector<double>& sumVec) const override { addRowSums(sumVec, nullptr); };
   void addColSums(Vector<double>& sumVec) const override { addColSums(sumVec, nullptr); };

   virtual void initStaticStorageFromDynamic(const Vector<int>& rowNnzVec, const Vector<int>& colNnzVec) {
      initStaticStorageFromDynamic(rowNnzVec, colNnzVec, nullptr, nullptr);
   };
   virtual void freeDynamicStorage();
   void recomputeSize(DistributedMatrix* parent = nullptr);

   /** returns Simple Vector indicating which linking rows have entries in exactly two blocks (indicated by 1.0 versus 0.0)*/
   virtual void get2LinkStartBlocksAndCountsNew(std::vector<int>& block_start, std::vector<int>& block_count) const;
   [[nodiscard]] virtual std::vector<int> get2LinkStartBlocks() const;

   virtual void updateKLinkVarsCount(std::vector<int>& linkCount) const;
   virtual void updateKLinkConsCount(std::vector<int>& linkCount) const;

   virtual void permuteLinkingVars(const std::vector<unsigned int>& permvec);
   virtual void permuteLinkingCons(const std::vector<unsigned int>& permvec);

   [[nodiscard]] virtual bool isRootNodeInSync() const;

   virtual int appendRow(const DistributedMatrix& matrix_row, int child, int row, bool linking);

   /** calculate first^T * row where row is linking or not in child child and with row index row
    *  for a linking row only the available blocks will be multiplied - currently only possible for dynamic storage! (since
    *  this was its foremost usecase)
    */
   [[nodiscard]] virtual double localRowTimesVec(const DistributedVector<double>& vec, int child, int row, bool linking) const;

   /* y += alpha * RowAt(child, row, linking) */
   virtual void axpyWithRowAt(double alpha, DistributedVector<double>* y, SimpleVector<double>* y_linking, int child, int row, bool linking) const;
   virtual void
   axpyWithRowAtPosNeg(double alpha, DistributedVector<double>* y_pos, SimpleVector<double>* y_link_pos, DistributedVector<double>* y_neg,
         SimpleVector<double>* y_link_neg, int child, int row, bool linking) const;

   [[nodiscard]] virtual std::unique_ptr<BorderedMatrix> raiseBorder(int m_conss, int n_vars);

   [[nodiscard]] virtual std::unique_ptr<StripMatrix> shaveLinkingConstraints(unsigned int n_conss);
   virtual void
   splitMatrix(const std::vector<int>& twolinks_start_in_block, const std::vector<unsigned int>& map_blocks_children, unsigned int n_links_in_root,
         const std::vector<MPI_Comm>& child_comms);


protected:
   virtual void write_to_streamDenseChild(std::ostream& out, int offset) const;
   virtual void write_to_streamDenseBorderedChild(const StripMatrix& border_left, std::ostream& out, int offset) const;

   /* internal methods for linking cons and hierarchical structure */
   virtual void getRowMinMaxVecChild(bool getMin, bool initializeVec, const Vector<double>* colScaleVec_, Vector<double>& minmaxVec_,
         Vector<double>* minmax_link_parent) const;
   virtual void getColMinMaxVecChild(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, const Vector<double>* rowScaleParent,
         Vector<double>& minmaxVec) const;

   [[nodiscard]] virtual bool amatEmpty() const;
   virtual void shaveBorder(int m_conss, int n_vars, StripMatrix* border_left, StripMatrix* border_bottom);
   [[nodiscard]] virtual std::unique_ptr<StripMatrix> shaveLeftBorder(int n_vars);
   [[nodiscard]] virtual std::unique_ptr<StripMatrix> shaveLeftBorderChild(int n_vars);
};


/**
 * Dummy Class 
 */

class StochGenDummyMatrix : public DistributedMatrix {

protected:

public:

   StochGenDummyMatrix() : DistributedMatrix(0, 0, 0, 0, 0, 0, 0, 0, MPI_COMM_NULL) {};
   ~StochGenDummyMatrix() override = default;
   void AddChild(const std::shared_ptr<DistributedMatrix>& ) override {};

public:
   void updateTransposed() const override {};

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override { return {0,0}; };
   [[nodiscard]] long long n_rows() const override { return 0; };
   [[nodiscard]] long long n_columns() const override { return 0; };

   using GeneralMatrix::cloneFull;
   using GeneralMatrix::cloneEmptyRows;
   [[nodiscard]] std::unique_ptr<GeneralMatrix> cloneEmptyRows(bool) const override { return std::make_unique<StochGenDummyMatrix>(); };
   [[nodiscard]] std::unique_ptr<GeneralMatrix> cloneFull(bool) const override { return std::make_unique<StochGenDummyMatrix>(); };


   /** The actual number of structural non-zero elements in this sparse
    *  matrix. This includes so-called "accidental" zeros, elements that
    *  are treated as non-zero even though their value happens to be zero.
    */
   [[nodiscard]] int numberOfNonZeros() const override { return 0; };

   [[nodiscard]] int is_a(int matType) const override;

   void columnScale(const Vector<double>&) override {};
   void rowScale(const Vector<double>&) override {};
   void scalarMult(double) override {};

   void getDiagonal(Vector<double>&) const override {};
   void setToDiagonal(const Vector<double>&) override {};

   void mult(double, Vector<double>&, double, const Vector<double>&) const override {};
   void mult2(double, DistributedVector<double>&, double, const DistributedVector<double>&, Vector<double>*) const override {};

   void transMult(double, Vector<double>&, double, const Vector<double>&) const override {};
   void transMult2(double, DistributedVector<double>&, double, const DistributedVector<double>&, const Vector<double>*) const override {};

   [[nodiscard]] double inf_norm() const override { return 0.0; };
   [[nodiscard]] double abminnormNonZero(double) const override { return std::numeric_limits<double>::infinity(); };

   void permuteLinkingVarsChild(const std::vector<unsigned int>&) override {};
   void getLinkVarsNnzChild(std::vector<int>&) const override {};

   void getLinkVarsNnz(std::vector<int>&) const override {};

   void write_to_streamDense(std::ostream&) const override {};
   void write_to_streamDense(std::ostream&, int) const override {};
   void write_to_streamDenseBordered(const StripMatrix&, std::ostream&, int) const override {};
   void writeDashedLineToStream(std::ostream&) const override {};
   void writeDashedLineToStream(std::ostream&, int) const override {};

   void writeMPSformatRows(std::ostream&, int, const Vector<double>*) const override {};

protected:
   void write_to_streamDenseChild(std::ostream&, int) const override {};
   void write_to_streamDenseBorderedChild(const StripMatrix&, std::ostream&, int) const override {};

public:
   void initTransposedChild(bool) const override {};

   void columnScale2(const Vector<double>&) override {};
   void rowScale2(const Vector<double>&, const Vector<double>*) override {};

   void initTransposed() const override {};
   void initTransposed(bool) const override {};
   void deleteTransposed() const override {};

   void sum_transform_rows(Vector<double>&, const std::function<double(const double&)>&) const override {};

   void getNnzPerRow(Vector<int>&, Vector<int>*) const override {};
   void getNnzPerCol(Vector<int>&, Vector<int>*) const override {};
   void getNnzPerRow(Vector<int>&) const override {};
   void getNnzPerCol(Vector<int>&) const override {};

   void getRowMinMaxVec(bool, bool, const Vector<double>*, Vector<double>&) const override {};
   void getColMinMaxVec(bool, bool, const Vector<double>*, Vector<double>&) const override {};

   void addRowSums(Vector<double>&, Vector<double>*) const override {};
   void addColSums(Vector<double>&, Vector<double>*) const override {};
   void addRowSums(Vector<double>&) const override {};
   void addColSums(Vector<double>&) const override {};

   void freeDynamicStorage() override {};

   void initStaticStorageFromDynamic(const Vector<int>&, const Vector<int>&) override {};
   void initStaticStorageFromDynamic(const Vector<int>&, const Vector<int>&, const Vector<int>*, const Vector<int>*) override {};

   [[nodiscard]] std::vector<int> get2LinkStartBlocks() const override { return std::vector<int>(); };

   void updateKLinkVarsCount(std::vector<int>&) const override {};
   void updateKLinkConsCount(std::vector<int>&) const override {};

   void permuteLinkingVars(const std::vector<unsigned int>&) override {};
   void permuteLinkingCons(const std::vector<unsigned int>&) override {};

   [[nodiscard]] bool isRootNodeInSync() const override { return true; };

   int appendRow(const DistributedMatrix&, int, int, bool) override {
      assert(0 && "CANNOT APPEND ROW TO DUMMY MATRIX");
      return -1;
   };

   [[nodiscard]] double localRowTimesVec(const DistributedVector<double>&, int, int, bool) const override {
      assert(0 && "CANNOT MULTIPLY ROW WITH DUMMY MATRIX");
      return -1;
   };

   void axpyWithRowAt(double, DistributedVector<double>*, SimpleVector<double>*, int, int, bool) const override {};
   void axpyWithRowAtPosNeg(double, DistributedVector<double>*, SimpleVector<double>*, DistributedVector<double>*, SimpleVector<double>*, int, int,
         bool) const override {};

   [[nodiscard]] std::unique_ptr<BorderedMatrix> raiseBorder(int, int) override {
      assert(0 && "CANNOT SHAVE BORDER OFF OF A DUMMY MATRIX");
      return nullptr;
   };
   std::unique_ptr<StripMatrix> shaveLinkingConstraints(unsigned int) override { return std::make_unique<StringGenDummyMatrix>(); };
   void splitMatrix(const std::vector<int>&, const std::vector<unsigned int>&, unsigned int, const std::vector<MPI_Comm>&) override {
      assert(0 && "CANNOT SHAVE BORDER OFF OF A DUMMY MATRIX");
   };

protected:
   void getRowMinMaxVecChild(bool, bool, const Vector<double>*, Vector<double>&, Vector<double>*) const override {};
   void getColMinMaxVecChild(bool, bool, const Vector<double>*, const Vector<double>*, Vector<double>&) const override {};

   void shaveBorder(int, int, StripMatrix* border_left, StripMatrix* border_bottom) override {
      border_left->addChild(std::make_unique<StringGenDummyMatrix>());
      border_bottom->addChild(std::make_unique<StringGenDummyMatrix>());
   };

   std::unique_ptr<StripMatrix> shaveLeftBorder(int) override { return std::make_unique<StringGenDummyMatrix>(); };
   std::unique_ptr<StripMatrix> shaveLeftBorderChild(int) override { return std::make_unique<StringGenDummyMatrix>(); };

};

#endif
