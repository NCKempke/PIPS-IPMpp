/*
 * StochRowStorage.h
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHROWSTORAGE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHROWSTORAGE_H_

#include "DistributedMatrix.h"
#include "DistributedVector.h"
#include "SystemType.h"

#include <memory>

class StochRowStorage {
public:
   StochRowStorage(const DistributedMatrix& system_matrix);
   ~StochRowStorage() = default;

   int storeRow(const INDEX& row, const DistributedMatrix& matrix_row);

   /** y = beta * y + alpha * stored row */
   void axpyAtRow(double beta, DistributedVector<double>* y, DenseVector<double>* y_linking, double alpha, const INDEX& row) const;
   void axpyAtRowPosNeg(double beta, DistributedVector<double>* y_pos, DenseVector<double>* y_link_pos, DistributedVector<double>* y_neg,
         DenseVector<double>* y_link_neg, double alpha, const INDEX& row) const;

   double multRowTimesVec(const INDEX& row, const DistributedVector<double>& vec) const;
   double getRowCoefficientAtColumn(const INDEX& row, const INDEX& col) const;

   // todo : deleteRowFromStorage
private:

   double multLinkingRowTimesVecWithoutBl0(int row, const DistributedVector<double>& vec) const;
   std::unique_ptr<DistributedMatrix> row_storage{};

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHROWSTORAGE_H_ */
