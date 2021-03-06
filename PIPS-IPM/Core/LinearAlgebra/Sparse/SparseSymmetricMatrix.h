/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef SPARSESYMMATRIX_H
#define SPARSESYMMATRIX_H

#include "../Abstract/AbstractMatrix.h"
#include "./SparseStorage.h"

class SparseMatrix;

/** Represents sparse symmetric matrices stored in
 *  row-major Harwell-Boeing format.
 *  @ingroup SparseLinearAlgebra
 */
class SparseSymmetricMatrix : public SymmetricMatrix {
private:
   std::unique_ptr<SparseStorage> mStorage;

   // is lower part of matrix stored? (otherwise upper part is stored)
   const bool isLower;

public:
   SparseSymmetricMatrix();
   SparseSymmetricMatrix(const SparseSymmetricMatrix& mat);

   SparseSymmetricMatrix(int size, int nnz, bool isLower = true);
   SparseSymmetricMatrix(int size, int nnz, int krowM[], int jcolM[], double M[], int deleteElts = 0, bool isLower = true);
   SparseSymmetricMatrix(std::unique_ptr<SparseStorage> m_storage, bool is_lower_);

   SparseStorage& getStorage() { return *mStorage; }
   [[nodiscard]] const SparseStorage& getStorage() const { return *mStorage; }

   int* krowM() { return mStorage->krowM; }
   int* jcolM() { return mStorage->jcolM; }
   double* M() { return mStorage->M; }

   [[nodiscard]] const int* krowM() const { return mStorage->krowM; }
   [[nodiscard]] const int* jcolM() const { return mStorage->jcolM; }
   [[nodiscard]] const double* M() const { return mStorage->M; }

   [[nodiscard]] int is_a(int type) const override;

   void putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) override;
   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const override;
   void fromGetSpRow(int row, int col, double A[], int lenA, int jcolA[], int& nnz, int colExtent, int& info) const override;

   void symmetricScale(const Vector<double>& vec) override;
   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;

   void symAtPutSpRow(int col, const double A[], int lenA, const int jcolA[], int& info) override;

   virtual void symPutZeroes();

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override;
   [[nodiscard]] long long n_rows() const override;
   [[nodiscard]] long long n_columns() const override;
   [[nodiscard]] long long size() const override;

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;
   void diagonal_add_constant_from(int from, int length, double value) override;
   void diagonal_set_to_constant_from(int from, int length, double value) override;

   void symAtPutSubmatrix(int this_start_row, int this_start_col, const AbstractMatrix& matix, int matrix_start_row, int matrix_start_col, int n_rows, int n_col) override;


   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   void transpose_mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;

   void write_to_stream(std::ostream& out) const override;
   void writeNNZpatternToStreamDense(std::ostream& out) const;
   void write_to_streamDense(std::ostream& out) const override;
   void write_to_streamDenseRow(std::ostream& out, int row) const override;

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;

   void fromGetDiagonal(int idiag, Vector<double>& v) const override;

   /** The actual number of structural non-zero elements in this sparse
    *  matrix. This includes so-called "accidental" zeros, elements that
    *  are treated as non-zero even though their value happens to be zero.
    */
   [[nodiscard]] int numberOfNonZeros() const { return mStorage->numberOfNonZeros(); }

   /** Reduce the matrix to lower triangular */
   void reduceToLower();

   [[nodiscard]] bool is_lower() const { return isLower; };

   void deleteEmptyRowsCols(const Vector<int>& nnzVec);

   void deleteZeroRowsCols(int*& new2orgIdx);

   void getSparseTriplet_c2fortran(int*& irn, int*& jcn, double*& val) const;

   void getSparseTriplet_fortran2fortran(int*& irn, int*& jcn, double*& val) const;

   std::unique_ptr<SparseMatrix> shaveSymLeftBottom(int n_vars);

   ~SparseSymmetricMatrix() override = default;

   [[nodiscard]] std::unique_ptr<SymmetricMatrix> clone() const override;

   void append_empty_diagonal(int n_values);
};

#endif
