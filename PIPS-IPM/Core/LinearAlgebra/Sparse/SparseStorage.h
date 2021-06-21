/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSESTORAGE_H
#define SPARSESTORAGE_H

#include "../Abstract/AbstractMatrix.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

/** A class for managing the matrix elements used by sparse matrices.
 *  @ingroup SparseLinearAlgebra
 */
class SparseStorage : public DoubleStorage {
private:
   /** store absolute non-zero minimum entry of row i and first[i] in first[i]; empty rows get value 0.0  */
   void getRowMinVec(const double* colScaleVec, double* vec) const;

   /** store absolute non-zero maximum entry of row i and first[i] in first[i]; empty rows get value 0.0  */
   void getRowMaxVec(const double* colScaleVec, double* vec) const;

   class index_sort {
   private:
      const int* indices;
      const int maxsize;

   public:
      index_sort(const int* indices, int maxsize) : indices(indices), maxsize(maxsize) {}

      bool operator()(int i, int j) const {
         assert(i < maxsize && j < maxsize);
         return (indices[i] < indices[j]);
      }
   };

protected:
   int neverDeleteElts;

public:
   static int instances;

   int m{};
   int n{};
   int len{};
   int* jcolM{};
   int* krowM{};
   double* M{};

   SparseStorage(int m_, int n_, int len_);
   SparseStorage(int m_, int n_, int len_, int* krowM_, int* jcolM_, double* M_, int deleteElts = 0);

   void copyFrom(int* krowM_, int* jcolM_, double* M_) const;

   void shiftRows(int row, int shift, int& info);

   [[nodiscard]] std::pair<int,int> n_rows_columns() const override;
   [[nodiscard]] int n_rows() const override;
   [[nodiscard]] int n_columns() const override;

   [[nodiscard]] bool isValid() const;
   [[nodiscard]] bool isSorted() const;


   [[nodiscard]] int length() const { return len; };
   [[nodiscard]] int numberOfNonZeros() const;

   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const override;
   void atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) override;

   void putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) override;

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;

   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void symmetricScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   void atPutSpRow(int col, const double A[], int lenA, const int irowA[], int& info) override;
   void fromGetSpRow(int row, int col, double A[], int lenA, int irowA[], int& nnz, int rowExtent, int& info) const override;

   virtual void clear();

   virtual void mult(double beta, double y[], double alpha, const double x[]) const;
   virtual void mult_transform(double beta, double y[], double alpha, const double x[], const std::function<double(const double&)>& transform) const;
   virtual void multSym(double beta, double y[], double alpha, const double x[]) const;

   virtual void transMult(double beta, double y[], double alpha, const double x[] ) const;
   virtual void transpose_mult_transform(double beta, double y[], double alpha, const double x[], const std::function<double(const double&)>& transform) const;
   virtual void transMultD(double beta, double y[], double alpha, const double x[], const double d[]) const;

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int idiag, Vector<double>& v) const override;

   void atPutDiagonal(int idiag, const double x[], int extent);
   void atAddDiagonal(int idiag, const double x[], int extent);

   void diagonal_set_to_constant_from(int from, int length, double value) override;
   void diagonal_add_constant_from(int from, int length, double value) override;

   virtual void write_to_stream(std::ostream& out) const;
   void writeNNZpatternToStreamDense(std::ostream& out) const;
   virtual void write_to_streamDense(std::ostream& out) const;
   virtual void write_to_streamDenseRow(std::ostream& out, int rowidx) const;

   virtual void symmetrize(int& info);
   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;

   /** Computes the sparsity pattern of MtM = M^T * D * M
    *  where D=diag(d) is a diagonal matrix and M=this.
    *
    *  Find the nonzero pattern of the product matrix, allocate it and returns it.
    *
    *  Also allocates, builds and returns this^T since it is needed later for
    *   numerical multiplication.
    */
   void matTransDSymbMultMat(const double* /*d*/, const int* krowMt, const int* jcolMt, const double* /*dMt*/, int** krowMtM, int** jcolMtM, double** dMtM) const;


   /** Numerical multiplication MtM = M^T * D * M  where
    *  D=diag(d) is a diagonal matrix and M=this.
    *  M^T and MtM buffers should be allocated before calling this method by calling
    *  method matTransDSymbMultMat.
    */
   void matTransDMultMat(const double* d, const int* krowMt, const int* jcolMt, const double* dMt, int* krowMtM, int* jcolMtM, double* dMtM) const;
   void matTransDinvMultMat(const double* d, const int* krowMt, const int* jcolMt, const double* dMt, int* krowMtM, int* jcolMtM, double* dMtM) const;

   /** Builds the transpose: Mt = this^T */
   void transpose(int* krowMt, int* jcolMt, double* dMt) const;

   void reduceToLower();

   void multMatSymUpper(double beta, SparseStorage& y, double alpha, const double x[], int yrow, int ycolstart) const;

   void fromGetColBlock(int col, double* A, int lda, int colExtent, int* colSparsity, bool& allzero) const;

   void fromGetRowsBlock(double* rows_array_dense, size_t row_start, size_t n_rows, size_t array_line_size, size_t array_line_offest,
         int* row_sparsity) const;
   void fromGetRowsBlock(const int* rowIndices, int nRows, int arrayLineSize, int arrayLineOffset, double* rowsArrayDense,
         int* rowSparsity = nullptr) const;

   /** add nnz per row to given array (of size nRows) */
   void addNnzPerRow(int* vec) const { addNnzPerRow(vec, 0, m); };
   void addNnzPerRow(int* vec, int begin_rows, int end_rows) const;

   void getLinkVarsNnz(std::vector<int>& vec) const;

   void sum_transform_rows(Vector<double>& result, const std::function<double(const double&)>& transform) const override;

   /** add abs. sum per row to given array (of size nRows) */
   void addRowSums(double* vec) const;

   /** store absolute non-zero minimum/maximum entry of row i and first[i] in first[i];
    *  empty rows get value 0.0 for maximization and <double>::max() for minimization  */
   void getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const;

   void permuteRows(const std::vector<unsigned int>& permvec);
   void permuteCols(const std::vector<unsigned int>& permvec);

   void sortCols();

   void dump(const std::string& filename) const;

   void deleteEmptyRowsCols(const int* nnzRowVec, const int* nnzColVec);

   void getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const;

   void getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const;

   void deleteEmptyRows(int*& orgIndex);

   // should be used with care! other methods might nor work correctly todo: add flag to check in other methods
   void c2fortran();

   void fortran2c();

   [[nodiscard]] bool fortranIndexed() const;

   void set2FortranIndexed();

   void deleteZeroRowsColsSym(int*& new2orgIdx);


   /*
    * computes the full sparse matrix representation from a upper triangular symmetric sparse representation
    *
    * Must be square, the storage for the full representation will be allocated within the matrix and must be released later
    */
   void fullMatrixFromUpperTriangular(int*& rowPtrFull, int*& colIdxFull, double*& valuesFull) const;

   std::unique_ptr<SparseStorage> shaveLeft(int n_cols);
   std::unique_ptr<SparseStorage> shaveSymLeftBottom(int n);
   std::unique_ptr<SparseStorage> shaveBottom(int n_rows);
   void dropNEmptyRowsBottom(int n_rows);
   void dropNEmptyRowsTop(int n_rows);

   ~SparseStorage() override;

private:
   bool isFortranIndexed{false};
};

#endif
