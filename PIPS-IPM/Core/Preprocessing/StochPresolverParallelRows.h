/*
 * StochPresolverParallelRows.h
 *
 *  Created on: 02.05.2018
 *      Author: bzfuslus
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERPARALLELROWS_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERPARALLELROWS_H_

#include "StochPresolverBase.h"

#include <boost/unordered_set.hpp>

namespace rowlib {
   static const double offset_hash_double = 0.127;

   struct rowWithColInd {
      const int id;
      const int offset_nA;
      const int lengthA;
      const int* const colIndicesA;
      const double* const norm_entriesA;
      const int lengthB;
      const int* const colIndicesB;
      const double* const norm_entriesB;

      rowWithColInd(int id, int offset, int lenA, const int* const colA, const double* const entA, int lenB, const int* const colB,
            const double* const entB) : id(id), offset_nA(offset), lengthA(lenA), colIndicesA(colA), norm_entriesA(entA), lengthB(lenB),
            colIndicesB(colB), norm_entriesB(entB) {}

      friend std::ostream& operator<<(std::ostream& out, const rowWithColInd& row) {
         out << "ID: " << row.id << ",offset_nA: " << row.offset_nA << ", A_col[";
         for (int i = 0; i < row.lengthA; ++i) {
            if (i != 0)
               out << ",";
            out << row.colIndicesA[i];
         }
         out << "], A_entr[";
         for (int i = 0; i < row.lengthA; ++i) {
            if (i != 0)
               out << ",";
            out << row.norm_entriesA[i];
         }
         out << "], B_col[";
         for (int i = 0; i < row.lengthB; ++i) {
            if (i != 0)
               out << ",";
            out << row.colIndicesB[i];
         }
         out << "], B_entr[";
         for (int i = 0; i < row.lengthB; ++i) {
            if (i != 0)
               out << ",";
            out << row.norm_entriesB[i];
         }
         out << "]";
         return out;
      }


   };

   bool operator==(rowWithColInd const& a, rowWithColInd const& b);
   std::size_t hash_value(rowWithColInd const& b);

   struct rowWithEntries {
      const int id;
      const int offset_nA;
      const int lengthA;
      const int* const colIndicesA;
      const double* const norm_entriesA;
      const int lengthB;
      const int* const colIndicesB;
      const double* const norm_entriesB;

      rowWithEntries(int id, int offset, int lenA, const int* const colA, const double* const entA, int lenB, const int* const colB,
            const double* const entB) : id(id), offset_nA(offset), lengthA(lenA), colIndicesA(colA), norm_entriesA(entA), lengthB(lenB),
            colIndicesB(colB), norm_entriesB(entB) {}
   };

   bool operator==(rowWithEntries const& a, rowWithEntries const& b);
   std::size_t hash_value(rowWithEntries const& b);
}

class StochPresolverParallelRows : public StochPresolverBase {
public:
   StochPresolverParallelRows(PresolveData& presolve_data, const DistributedProblem& origProb);

   ~StochPresolverParallelRows();

   // remove parallel rows
   bool applyPresolving() override;

private:

   /** tolerance for comparing two double values in two different rows and for them being considered equal */
   const double limit_tol_compare_entries;
   int n_rows_removed{0};

   /// extension to the pointer set from StochPresolverBase to point to C and A at the same moment rather than
   /// distinguishing between EQUALITY and INEQUALITY constraints
   const SparseStorageDynamic* currCmat{};
   const SparseStorageDynamic* currCmatTrans{};
   const SparseStorageDynamic* currDmat{};
   const SparseStorageDynamic* currDmatTrans{};
   const DenseVector<int>* currNnzRowC{};

   // pointers to the normalized and copied matrix blocks
   std::unique_ptr<SparseStorageDynamic> norm_Amat{};
   std::unique_ptr<SparseStorageDynamic> norm_Bmat{};
   std::unique_ptr<SparseStorageDynamic> norm_Cmat{};
   std::unique_ptr<SparseStorageDynamic> norm_Dmat{};
   std::unique_ptr<DenseVector<double>> norm_b{};
   std::unique_ptr<DenseVector<double>> norm_clow{};
   std::unique_ptr<DenseVector<double>> norm_cupp{};
   std::unique_ptr<DenseVector<double>> norm_iclow{};
   std::unique_ptr<DenseVector<double>> norm_icupp{};
   std::unique_ptr<DenseVector<double>> norm_factorC{};
   std::unique_ptr<DenseVector<double>> norm_factorA{};

   // data for the nearly parallel row case
   std::unique_ptr<DenseVector<int>> rowContainsSingletonVariableA{};
   std::unique_ptr<DenseVector<int>> rowContainsSingletonVariableC{};
   std::unique_ptr<DenseVector<double>> singletonCoeffsColParent{};
   std::unique_ptr<DenseVector<double>> singletonCoeffsColChild{};
   std::unique_ptr<DenseVector<int>> normNnzRowA{};
   std::unique_ptr<DenseVector<int>> normNnzRowC{};
   std::unique_ptr<DenseVector<int>> normNnzColParent{};
   std::unique_ptr<DenseVector<int>> normNnzColChild{};
   std::unique_ptr<SparseStorageDynamic> norm_AmatTrans{};
   std::unique_ptr<SparseStorageDynamic> norm_BmatTrans{};
   std::unique_ptr<SparseStorageDynamic> norm_CmatTrans{};
   std::unique_ptr<SparseStorageDynamic> norm_DmatTrans{};

   // number of rows of the A or B block
   int mA{0};
   // number of columns of the A or C block
   int nA{0};

   // unordered set?
   boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> > row_support_hashtable;
   boost::unordered_set<rowlib::rowWithEntries, boost::hash<rowlib::rowWithEntries> > row_coefficients_hashtable;

   void setNormalizedPointers(int node);
   void setNormalizedPointersMatrices(int node);
   void setNormalizedPointersMatrixBounds(int node);
   void setNormalizedNormFactors(int node);
   void setNormalizedSingletonFlags(int node);
   void setNormalizedReductionPointers(int node);
   void updateExtendedPointersForCurrentNode(int node);

   void removeSingletonVars();
   void removeEntry(int colIdx, DenseVector<int>& rowContainsSingletonVar, SparseStorageDynamic& matrix, SparseStorageDynamic& matrixTrans,
         DenseVector<int>& nnzRow, DenseVector<int>& nnzCol, bool parent);
   double removeEntryInDynamicStorage(SparseStorageDynamic& storage, int row, int col) const;

   void normalizeBlocksRowwise(SystemType system_type, SparseStorageDynamic* a_mat, SparseStorageDynamic* b_mat, DenseVector<double>* cupp,
         DenseVector<double>* clow, DenseVector<double>* icupp, DenseVector<double>* iclow) const;
   void
   insertRowsIntoHashtable(boost::unordered_set<rowlib::rowWithColInd, boost::hash<rowlib::rowWithColInd> >& rows, const SparseStorageDynamic* Ablock,
         const SparseStorageDynamic* Bblock, SystemType system_type, const DenseVector<int>* nnz_row_norm, const DenseVector<int>* nnz_row_orig);
   void compareRowsInCoeffHashTable(int& nRowElims, int it);
   bool checkRowsAreParallel(const rowlib::rowWithEntries& row1, const rowlib::rowWithEntries& row2) const;

   void tightenOriginalBoundsOfRow1(const INDEX& row1, const INDEX& row2) const;

   double getSingletonCoefficient(const INDEX& col) const;

   INDEX getRowSingletonVariable(const INDEX& row) const;
   bool rowContainsSingletonVariable(const INDEX& row) const;

   void tightenBoundsForSingleVar(int singleColIdx, double newxlow, double newxupp);

   bool twoParallelEqualityRows(const INDEX& row1, const INDEX& row2) const;
   bool parallelEqualityAndInequalityRow(const INDEX& row_eq, const INDEX& row_ineq) const;
   bool twoParallelInequalityRows(const INDEX& row1, const INDEX& row2) const;

   bool twoNearlyParallelEqualityRows(const INDEX& row1, const INDEX& row2) const;
   bool nearlyParallelEqualityAndInequalityRow(const INDEX& row_eq, const INDEX& row_ineq) const;
   bool twoNearlyParallelInequalityRows(const INDEX& row1, const INDEX& row2) const;

   void tightenLinkingVarsBounds();

};

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERPARALLELROWS_H_ */
