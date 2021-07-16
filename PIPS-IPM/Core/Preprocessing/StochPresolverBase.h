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
         int& n_cols_orig_free, int& n_cols_orig_free_removed, const DenseVector<double>& ixlow_orig, const DenseVector<double>& ixupp_orig,
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

   const DenseVector<double>* currxlowParent;
   const DenseVector<double>* currIxlowParent;
   const DenseVector<double>* currxuppParent;
   const DenseVector<double>* currIxuppParent;
   const DenseVector<double>* currxlowChild;
   const DenseVector<double>* currIxlowChild;
   const DenseVector<double>* currxuppChild;
   const DenseVector<double>* currIxuppChild;

   const DenseVector<double>* currEqRhs;
   const DenseVector<double>* currIneqLhs;
   const DenseVector<double>* currIclow;
   const DenseVector<double>* currIneqRhs;
   const DenseVector<double>* currIcupp;
   const DenseVector<double>* currEqRhsLink;
   const DenseVector<double>* currIneqLhsLink;
   const DenseVector<double>* currIclowLink;
   const DenseVector<double>* currIneqRhsLink;
   const DenseVector<double>* currIcuppLink;

   const DenseVector<double>* currgParent;
   const DenseVector<double>* currgChild;

   const DenseVector<int>* currNnzRow;
   const DenseVector<int>* currNnzRowLink;

   const DenseVector<int>* currNnzColParent;
   const DenseVector<int>* currNnzColChild;

   /** the number of children */
   int nChildren;
   /** number of entry eliminations on this process in the current elimination routine */
   int localNelims;

};


#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVERBASE_H_ */
