/*
 * BorderedSymMatrix.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#ifndef BORDEREDSYMMETRICMATRIX_H
#define BORDEREDSYMMETRICMATRIX_H

#include "AbstractMatrix.h"
#include "SparseSymmetricMatrix.h"
#include "DistributedSymmetricMatrix.h"
#include "StripMatrix.h"
#include "DistributedVector.h"
#include <vector>
#include "mpi.h"

class StripMatrix;

class DistributedSymmetricMatrix;

/* representing a matrix of the type
 *
 * [  D   B ]
 * [ B^T  Q ]
 *
 */

class BorderedSymmetricMatrix : public SymmetricMatrix {
public:
   BorderedSymmetricMatrix(std::shared_ptr<DistributedSymmetricMatrix> inner_matrix, std::unique_ptr<StripMatrix> border_vertical,
      std::unique_ptr<SymmetricMatrix> top_left_block, MPI_Comm mpiComm_);

   ~BorderedSymmetricMatrix() override = default;

   [[nodiscard]] int is_a(int type) const override;

   /** y = beta * y + alpha * this * x */
   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   /** y = beta * y + alpha * this^T * x */
   void transpose_mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   void fromGetDiagonal(int idiag, Vector<double>& x) const override;

   [[nodiscard]] double inf_norm() const override;
   void scalarMult(double num) override;

   [[nodiscard]] std::pair<long long, long long> n_rows_columns() const override;
   [[nodiscard]] long long n_rows() const override;
   [[nodiscard]] long long n_columns() const override;
   [[nodiscard]] long long size() const override;

   // ignored methods
   [[nodiscard]] double abminnormNonZero(double) const override {
      assert(false && "TODO: implement");
      return 0.0;
   };
   void write_to_streamDense(std::ostream&) const override { assert("Not implemented" && 0); };
   void symAtPutSpRow(int, const double[], int, const int[], int&) override { assert("Not implemented" && 0); };
   void symAtPutSubmatrix(int, int, const AbstractMatrix&, int, int, int, int) override { assert("Not implemented" && 0); };
   void atPutDiagonal(int, const Vector<double>&) override { assert("Not implemented" && 0); };
   void atAddDiagonal(int, const Vector<double>&) override { assert("Not implemented" && 0); };
   void fromGetDense(int, int, double*, int, int, int) const override { assert("Not implemented" && 0); };
   void fromGetSpRow(int, int, double[], int, int[], int&, int, int&) const override { assert("Not implemented" && 0); };
   void getDiagonal(Vector<double>&) const override { assert("Not implemented" && 0); };
   void setToDiagonal(const Vector<double>&) override { assert("Not implemented" && 0); };
   void symmetricScale(const Vector<double>&) override { assert("Not implemented" && 0); };
   void columnScale(const Vector<double>&) override { assert("Not implemented" && 0); };
   void rowScale(const Vector<double>&) override { assert("Not implemented" && 0); };
   void write_to_stream(std::ostream&) const override { assert("Not implemented" && 0); };
   void putSparseTriple(const int[], int, const int[], const double[], int&) override { assert("Not implemented" && 0); };

   // TODO could be more general..
   std::shared_ptr<DistributedSymmetricMatrix> inner_matrix;
   std::unique_ptr<StripMatrix> border_vertical;

   std::unique_ptr<SymmetricMatrix> top_left_block;

protected:

   MPI_Comm mpiComm;
   const int iAmDistrib;

   long long n;
};

#endif /* BORDEREDSYMMETRICMATRIX_H */
