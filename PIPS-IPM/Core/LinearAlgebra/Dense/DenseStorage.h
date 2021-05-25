/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DENSEDOUBLEMATRIX_H
#define DENSEDOUBLEMATRIX_H

#include "AbstractMatrix.h"
#include "Vector.hpp"

class SparseStorage;

extern int DenseStorageInstances;

/** A class for manupulating the storage of dense matrices.
 *  @ingroup DenseLinearAlgebra
 */
class DenseStorage : public DoubleStorage {

protected:
   int neverDeleteElts;
public:
   int m;
   int n;
   double** M;

   DenseStorage(int m, int n);
   DenseStorage(double A[], int m, int n);

   ~DenseStorage() override;

   [[nodiscard]] std::pair<int,int> n_rows_columns() const override;
   [[nodiscard]] int n_rows() const override;
   [[nodiscard]] int n_columns() const override;

   void getDiagonal(Vector<double>& vec) const override;
   void setToDiagonal(const Vector<double>& vec) override;

   void atPutDense(int row, int col, const double* A, int lda, int rowExtent, int colExtent) override;

   void fromGetDense(int row, int col, double* A, int lda, int rowExtent, int colExtent) const override;

   void atPutZeros(int row, int col, int rowExtent, int colExtent);

   void putZeros();

   void atAddOuterProductOf(int row, int col, double alpha, double* x, int incx, int nx);

   void addToDiagonalAt(double alpha, double x[], int incx, int idiag, int extent);
   void fromGetSpRow(int row, int col, double A[], int lenA, int irowA[], int& nnz, int rowExtent, int& info) const override;

   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;
   void symmetricScale(const Vector<double>& vec) override;
   void scalarMult(double num) override;
   [[nodiscard]] double inf_norm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;
   [[nodiscard]] int non_zeros() const;

   void atPutSpRow(int col, const double A[], int lenA, const int irowA[], int& info) override;
   void putSparseTriple(const int irow[], int len, const int jcol[], const double A[], int& info) override;

   void atPutDiagonal(int idiag, const Vector<double>& v) override;
   void atAddDiagonal(int idiag, const Vector<double>& v) override;
   void fromGetDiagonal(int idiag, Vector<double>& v) const override;
   void atPutDiagonal(int idiag, const double x[], int incx, int extent);
   void atAddDiagonal(int idiag, const double x[], int incx, int extent);

   void diagonal_add_constant_from(int from, int length, double value) override;
   void diagonal_set_to_constant_from(int from, int length, double value) override;

   void fill_from_sparse(const SparseStorage& other);
   void fill_from_dense(const DenseStorage& other);
};

#endif
