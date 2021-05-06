/*
 * StochColumnStorage.h
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_

#include "DistributedMatrix.h"
#include "DistributedVector.h"
#include "SystemType.h"

#include <memory>

class StochColumnStorage {
public:
   StochColumnStorage(const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part);
   ~StochColumnStorage() = default;

   int storeCol(const INDEX& col, const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part);

   /** y = beta * first + alpha * stored - either of the vectors can be nullptr and will then not be considered */
   /** there is the option of specifying a special vector to store the linking row results in - if no vector is given the corresponding link vector will be ignored */
   void axpyAtCol(double beta, DistributedVector<double>* eq_vec, DistributedVector<double>* ineq_vec, SimpleVector<double>* eq_link,
         SimpleVector<double>* ineq_link, double alpha, const INDEX& col) const;

   double multColTimesVec(const INDEX& col, const DistributedVector<double>& vec_eq, const DistributedVector<double>& vec_ineq) const;

   // todo: delete Column from storage
private:
   // todo : assert that transposed is initalized
   /// all columns get stored as rows
   /// for a linking column the additional block B0 is stored in the additional associated SparseGenMatrix
   std::unique_ptr<SparseMatrix> B0_eq{};
   std::unique_ptr<DistributedMatrix> stored_cols_eq{};

   std::unique_ptr<SparseMatrix> B0_ineq{};
   std::unique_ptr<DistributedMatrix> stored_cols_ineq{};

   const unsigned int nChildren;

   int storeLinkingCol(int col, const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part);
   int storeLocalCol(const INDEX& col, const DistributedMatrix& matrix_eq_part, const DistributedMatrix& matrix_ineq_part);

   double multiplyLocalColTimesVec(const INDEX& col, const DistributedVector<double>& vec_eq, const DistributedVector<double>& vec_ineq) const;
   double multiplyLinkingColTimesVec(int col, const DistributedVector<double>& vec_eq, const DistributedVector<double>& vec_ineq) const;
   double
   multiplyLinkingColTimesVecWithoutRootNode(int col, const DistributedVector<double>& vec_eq, const DistributedVector<double>& vec_ineq) const;

   void createStorageMatrix(SystemType system_type, const DistributedMatrix& sys_matrix);

   /* calculate first + alpha * stored for system */
   void axpyAtCol(DistributedVector<double>& vec, SimpleVector<double>* vec_link, double alpha, const INDEX& col, SystemType system_type) const;
};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_ */
