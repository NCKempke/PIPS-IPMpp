/*
 * BorderedSymMatrix.C
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDSYMMATRIX_H_
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDSYMMATRIX_H_

#include "DoubleMatrix.h"
#include "SparseSymMatrix.h"
#include "StochSymMatrix.h"
#include "StringGenMatrix.h"
#include "DistributedVector.h"
#include <vector>
#include "mpi.h"

class StringGenMatrix;

class StochSymMatrix;

/* representing a matrix of the type
 *
 * [  D   B ]
 * [ B^T  Q ]
 *
 */

class BorderedSymMatrix : public SymMatrix {
public:
   BorderedSymMatrix(StochSymMatrix* inner_matrix, StringGenMatrix* border_vertical, SymMatrix* bottom_block, MPI_Comm mpiComm_);

   virtual ~BorderedSymMatrix();

   int isKindOf(int type) const override;

   /** y = beta * y + alpha * this * x */
   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   /** y = beta * y + alpha * this^T * x */
   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   void fromGetDiagonal(int idiag, Vector<double>& x) override;

   double abmaxnorm() const override;
   void scalarMult(double num) override;
   void getSize(long long& m, long long& n) const override;
   void getSize(int& m, int& n) const override;
   long long size() const override;

   // ignored methods
   double abminnormNonZero(double) const override {
      assert(false && "TODO: implement");
      return 0.0;
   };
   void writeToStreamDense(std::ostream&) const override { assert("Not implemented" && 0); };
   void symAtPutSpRow(int, double[], int, int[], int&) override { assert("Not implemented" && 0); };
   void symAtPutSubmatrix(int, int, DoubleMatrix&, int, int, int, int) override { assert("Not implemented" && 0); };
   void atPutDiagonal(int, const Vector<double>&) override { assert("Not implemented" && 0); };
   void atAddDiagonal(int, const Vector<double>&) override { assert("Not implemented" && 0); };
   void fromGetDense(int, int, double*, int, int, int) override { assert("Not implemented" && 0); };
   void fromGetSpRow(int, int, double[], int, int[], int&, int, int&) override { assert("Not implemented" && 0); };
   void getDiagonal(Vector<double>&) override { assert("Not implemented" && 0); };
   void setToDiagonal(const Vector<double>&) override { assert("Not implemented" && 0); };
   void symmetricScale(const Vector<double>&) override { assert("Not implemented" && 0); };
   void columnScale(const Vector<double>&) override { assert("Not implemented" && 0); };
   void rowScale(const Vector<double>&) override { assert("Not implemented" && 0); };
   void writeToStream(std::ostream&) const override { assert("Not implemented" && 0); };
   void putSparseTriple(int[], int, int[], double[], int&) override { assert("Not implemented" && 0); };

   // TODO could be more general..
   StochSymMatrix* inner_matrix;
   StringGenMatrix* border_vertical;

   SymMatrix* top_left_block;

protected:

   MPI_Comm mpiComm;
   const int iAmDistrib;

   long long n;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDSYMMATRIX_H_ */
