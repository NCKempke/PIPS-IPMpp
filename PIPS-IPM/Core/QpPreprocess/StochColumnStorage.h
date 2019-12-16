/*
 * StochColumnStorage.h
 *
 *  Created on: 12.12.2019
 *      Author: Nils-Christian Kempke
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_

#include "StochGenMatrix.h"
#include "StochVector.h"

class StochColumnStorage
{
public:
   StochColumnStorage( const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);
   ~StochColumnStorage();

   int storeCol( int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);

   double multColTimesVec( int node, int col, const StochVector& vec_eq, const StochVector& vec_ineq ) const;
   double getColCoefficientAtRow( int node, int col, int row) const;

   // todo: delete Column from storage
private:
   // todo : assert that transposed is initalized
   /// all columns get stored as rows
   /// for a linking column the additional block B0 is stored in the additional associated SparseGenMatrix
   SparseGenMatrixHandle b0_block_linking_cols_eq;
   StochGenMatrixHandle stored_cols_eq;

   StochGenMatrixHandle stored_cols_ineq;
   SparseGenMatrixHandle b0_block_linking_cols_ineq;

   int storeLinkingCol(int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);
   int storeLocalCol(int node, int col, const StochGenMatrix& matrix_eq_part, const StochGenMatrix& matrix_ineq_part);
   void createStorageMatrix(SparseGenMatrix* b0_block_storage, StochGenMatrix* col_storage, const StochGenMatrix& sys_matrix);
};




#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHCOLUMNSTORAGE_H_ */
