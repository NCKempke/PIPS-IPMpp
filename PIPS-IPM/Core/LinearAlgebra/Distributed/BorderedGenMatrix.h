/*
 * BorderedGenMatrix.h
 *
 *  Created on: Sep 14, 2020
 *      Author: bzfkempk
 */

#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDGENMATRIX_H_
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDGENMATRIX_H_

#include "DoubleMatrix.h"

#include <vector>
#include "mpi.h"

class StringGenMatrix;

class StochGenMatrix;

class SparseGenMatrix;


/*
 * representing a matrix of type
 *
 * [ B  K ]
 * [ C  B']
 *
 * Where K is a StochGenMatrix
 */

// TODO : make more general ? K B B' and C can be any matrices in theory..
class BorderedGenMatrix : public GenMatrix {
public:
   StochGenMatrix* const inner_matrix{};
   StringGenMatrix* const border_left{};
   StringGenMatrix* const border_bottom{};

   // TODO: is SparseGenMatrix appropriate? What does this block look like -> it has parts of the diagonals in it for inequality linking constraints and nothing else?
   SparseGenMatrix* const bottom_left_block{};

protected:

   MPI_Comm mpi_comm{MPI_COMM_NULL};
   const bool distributed{false};
   const int rank{-1};

   long long m{0};
   long long n{0};

public:
   BorderedGenMatrix(StochGenMatrix* inner_matrix, StringGenMatrix* border_left, StringGenMatrix* border_bottom, SparseGenMatrix* bottom_left_block,
         MPI_Comm mpi_comm_);
   virtual ~BorderedGenMatrix();

   int isKindOf(int matrixType) const override;

   /** y = beta * y + alpha * this * x */
   void mult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;
   /** y = beta * y + alpha * this^T * x */
   void transMult(double beta, Vector<double>& y, double alpha, const Vector<double>& x) const override;

   double abmaxnorm() const override;
   void columnScale(const Vector<double>& vec) override;
   void rowScale(const Vector<double>& vec) override;

   void scalarMult(double num) override;

   void getSize(long long& m_, long long& n_) const override;
   void getSize(int& m_, int& n_) const override;

   void getRowMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* colScaleVec, Vector<double>& minmaxVec) override;
   void getColMinMaxVec(bool getMin, bool initializeVec, const Vector<double>* rowScaleVec, Vector<double>& minmaxVec) override;

   void writeToStreamDense(std::ostream& out) const override;

   void addRowSums(Vector<double>& vec) const override;
   void addColSums(Vector<double>& vec) const override;

   /* methods not needed for Hierarchical approach */
   double abminnormNonZero(double) const override {
      assert(false && "TODO: implement");
      return 0.0;
   };
   void writeToStream(std::ostream&) const override { assert(0 && "not implemented"); }; // TODO : implement maybe?
   void getDiagonal(Vector<double>&) override {
      assert(0 && "not implemented");
   }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
   void setToDiagonal(const Vector<double>&) override {
      assert(0 && "not implemented");
   }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
   void atPutDiagonal(int, const Vector<double>&) override {
      assert(0 && "not implemented");
   }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
   void atAddDiagonal(int, const Vector<double>&) override {
      assert(0 && "not implemented");
   }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?
   void fromGetDiagonal(int, Vector<double>&) override {
      assert(0 && "not implemented");
   }; // TODO : not sure - maybe we want this to get forwarded to the underlying matrix?

   void matTransDMultMat(Vector<double>&, SymMatrix**) override { assert(0 && "not implemented"); }; // TODO : needed?
   void matTransDinvMultMat(Vector<double>&, SymMatrix**) override { assert(0 && "not implemented"); }; // TODO : needed?

   void getNnzPerRow(Vector<int>&) override { assert(0 && "not implemented"); };
   void getNnzPerCol(Vector<int>&) override { assert(0 && "not implemented"); };
   void fromGetDense(int, int, double*, int, int, int) override { assert(0 && "not implemented"); };
   void fromGetSpRow(int, int, double*, int, int*, int&, int, int&) override { assert(0 && "not implemented"); };
   void putSparseTriple(int*, int, int*, double*, int&) override { assert(0 && "not implemented"); };
   void symmetricScale(const Vector<double>&) override { assert(0 && "not implemented"); };
   void atPutSubmatrix(int, int, DoubleMatrix&, int, int, int, int) override { assert(0 && "not implemented"); };
   void atPutDense(int, int, double*, int, int, int) override { assert(0 && "not implemented"); };
   void atPutSpRow(int, double*, int, int*, int&) override { assert(0 && "not implemented"); };
private:

   template<typename T>
   bool hasVecStructureForBorderedMat(const Vector<T>& vec, bool row_vec) const;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_BORDEREDGENMATRIX_H_ */
