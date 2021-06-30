/*
 * StochPresolverBase.h
 *
 *  Created on: 06.04.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_

#include "DistributedVector.h"
#include "DistributedMatrix.h"
#include "PresolveData.h"
#include "DistributedProblem.hpp"
#include "SystemType.h"
#include "StochPostsolver.h"
#include <vector>

class StochPresolverBase {
public:
   StochPresolverBase(PresolveData& presolve_data, const DistributedProblem& origProb);
   virtual ~StochPresolverBase() = default;

   virtual bool applyPresolving() = 0;
   void countRowsCols(); // theoretically const but sets pointers

protected:
   void updatePointersForCurrentNode(int node, SystemType system_type);

private:
   void countRowsBlock(int& n_rows_total, int& n_rows_empty, int& n_rows_onesided, int& n_rows_boxed, int& n_rows_fixed, int& n_rows_singleton,
         SystemType system_type, BlockType block_type) const;
   void countBoxedColumns(int& n_cols_total, int& n_cols_empty, int& n_cols_free, int& n_cols_onesided, int& n_cols_boxed, int& n_cols_singleton,
         int& n_cols_orig_free, int& n_cols_orig_free_removed, const SimpleVector<double>& ixlow_orig, const SimpleVector<double>& ixupp_orig,
         bool at_root_node) const;

   void setPointersMatrices(const GeneralMatrix& mat, int node);
   void setPointersMatrixBoundsActivities(SystemType system_type, int node);
   void setPointersVarBounds(int node);
   void setPointersObjective(int node);
   void setReductionPointers(SystemType system_type, int node);
   void setPointersToNull();

protected:
   const int my_rank;
   const bool distributed;

   const int verbosity;

   const double INF_NEG;
   const double INF_POS;

   const int n_linking_vars;
   const int n_linking_rows_eq;
   const int n_linking_rows_ineq;

   /* not owned by the class itself - given from the outside */
   PresolveData& presolve_data;

   // pointers to the currently needed matrices and vectors for presolving
   const DistributedProblem& origProb;

   const SparseStorageDynamic* currAmat;
   const SparseStorageDynamic* currAmatTrans;
   const SparseStorageDynamic* currBmat;
   const SparseStorageDynamic* currBmatTrans;
   const SparseStorageDynamic* currBlmat;
   const SparseStorageDynamic* currBlmatTrans;

   const SimpleVector<double>* currxlowParent;
   const SimpleVector<double>* currIxlowParent;
   const SimpleVector<double>* currxuppParent;
   const SimpleVector<double>* currIxuppParent;
   const SimpleVector<double>* currxlowChild;
   const SimpleVector<double>* currIxlowChild;
   const SimpleVector<double>* currxuppChild;
   const SimpleVector<double>* currIxuppChild;

   const SimpleVector<double>* currEqRhs;
   const SimpleVector<double>* currIneqLhs;
   const SimpleVector<double>* currIclow;
   const SimpleVector<double>* currIneqRhs;
   const SimpleVector<double>* currIcupp;
   const SimpleVector<double>* currEqRhsLink;
   const SimpleVector<double>* currIneqLhsLink;
   const SimpleVector<double>* currIclowLink;
   const SimpleVector<double>* currIneqRhsLink;
   const SimpleVector<double>* currIcuppLink;

   const SimpleVector<double>* currgParent;
   const SimpleVector<double>* currgChild;

   const SimpleVector<int>* currNnzRow;
   const SimpleVector<int>* currNnzRowLink;

   const SimpleVector<int>* currNnzColParent;
   const SimpleVector<int>* currNnzColChild;

   /** the number of children */
   int nChildren;
   /** number of entry eliminations on this process in the current elimination routine */
   int localNelims;

};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
