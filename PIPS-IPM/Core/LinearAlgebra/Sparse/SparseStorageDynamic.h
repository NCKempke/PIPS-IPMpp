/*
 * SparseStorageDynamic.h
 *
 *  Created on: 07.02.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_
#define PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_

#include "pipsport.h"
#include "DoubleMatrix.h"
#include "SparseStorage.h"
#include <vector>

typedef struct {
   int start;
   int end;
} ROWPTRS;

struct first_is_smaller {
   bool operator()(const std::pair<int, double>& x, const std::pair<int, double>& y) const {
      return x.first < y.first;
   }
};

/** A class for managing the matrix elements used by sparse matrices.
 *  @ingroup SparseLinearAlgebra
 */
class SparseStorageDynamic : public DoubleStorage {

private:
   const double spareRatio;

   int m;      // rows
   int m_len;  // length row array
   int n;      // cols
   int len;    // length col/value array
   int len_free;

   ROWPTRS* rowptr;
   int* jcolM;
   double* M;

   /* doubles the size of rowptr */
   void extendStorageRows();

   /* compresses storage and doubles size of the col entry storage */
   void extendStorageValues();

   /* shifts all rows such that every row is again row + spareRatio length */
   void rebuildSpareStructure(int guaranteed_spare = 0);

public:
   static int instances;

   [[nodiscard]] int getM() const { return m; };
   [[nodiscard]] int getN() const { return n; };
   [[nodiscard]] int getNVals() const { return len - len_free; };

   [[nodiscard]] const ROWPTRS* getRowPtr() const { return rowptr; };
   [[nodiscard]] ROWPTRS getRowPtr(int i) const;

   [[nodiscard]] const int* getJcolM() const { return jcolM; };
   [[nodiscard]] int getJcolM(int i) const;

   [[nodiscard]] const double* getMat() const { return M; };
   [[nodiscard]] double getMat(int i) const;
   void setMat(int i, double val);

   explicit SparseStorageDynamic(const SparseStorage& storage, double spareRatio = 0.2);
   SparseStorageDynamic(int m, int n, int len, double spareRatio = 0.2);
   SparseStorageDynamic(const SparseStorageDynamic& dynamicStorage);

   ~SparseStorageDynamic() override;

   void atPutDense(int, int, const double*, int, int, int) override { assert(0 && "not implemented here"); };
   void putSparseTriple(const int[], int, const int[], const double[], int&) override { assert(false && "not implemented here"); };
   void fromGetDense(int, int, double*, int, int, int) const override { assert(0 && "not implemented here"); };
   void atPutSpRow(int, const double*, int, const int*, int&) override { assert(0 && "not implemented here"); };
   void fromGetSpRow(int, int, double*, int, int*, int&, int, int&) const override { assert(0 && "not implemented here"); };
   void getDiagonal(Vector<double>&) const override { assert(0 && "not implemented here"); };
   void setToDiagonal(const Vector<double>&) override { assert(0 && "not implemented here"); };
   void atPutDiagonal(int, const Vector<double>&) override { assert(0 && "not implemented here"); };
   void atAddDiagonal(int, const Vector<double>&) override { assert(0 && "not implemented here"); };
   void fromGetDiagonal(int, Vector<double>&) const override { assert(0 && "not implemented here"); };
   void symmetricScale(const Vector<double>&) override { assert(0 && "not implemented here"); };
   void columnScale(const Vector<double>&) override { assert(0 && "not implemented here"); };
   void rowScale(const Vector<double>&) override { assert(0 && "not implemented here"); };
   void scalarMult(double) override { assert(0 && "not implemented here"); };

   void getSize(int& m, int& n) const override;

   void removeEntryAtIndex(int row, int col_idx);
   void removeEntryAtRowCol(int row, int col);

   bool addColToRow(double coeff, int col, int row);

   void clearRow(int row);
   void clearCol(int col);

   void appendRow(const SparseStorageDynamic& storage, int row);

   double rowTimesVec(const double* vec, int length, int row) const;
   void axpyWithRowAt(double alpha, double* y, int length, int row) const;
   void axpyWithRowAtPosNeg(double alpha, double* y_pos, double* y_neg, int length, int row) const;

   void scaleRow(int row, double factor);

   void addNnzPerRow(int* vec) const { addNnzPerRow(vec, 0, m); };
   void addNnzPerRow(int* vec, int begin_rows, int end_rows) const;

   void writeToStreamDense(std::ostream& out) const;
   void writeToStreamDenseRow(std::ostream& out, int rowidx) const;

   void restoreOrder();

   [[nodiscard]] double abmaxnorm() const override;
   [[nodiscard]] double abminnormNonZero(double tol) const override;

   [[nodiscard]] SparseStorage* getStaticStorage(const int* rowNnz, const int* colNnz) const;
   [[nodiscard]] SparseStorageDynamic* getTranspose() const;

   void getRowMaxVec(const double* colScaleVec, double* vec) const;
   void getRowMinMaxVec(bool getMin, const double* colScaleVec, double* vec) const;
   void getRowMinVec(const double* colScaleVec, double* vec) const;

};

typedef SmartPointer<SparseStorageDynamic> SparseStorageDynamicHandle;

#endif /* PIPS_IPM_CORE_SPARSELINEARALGEBRA_SPARSESTORAGEDYNAMIC_H_ */
