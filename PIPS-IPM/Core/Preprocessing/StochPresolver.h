/*
 * StochPresolver.h
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_
#define PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_

#include "QpPresolver.h"
#include "PresolveData.h"
#include <vector>
#include <memory>

class DistributedTree;

class StochPresolverBase;

/**  * @defgroup QpPreprocess
 *
 * QP presolver
 * @{
 */

/**
 * Derived class for QP presolvers.
 */
class StochPresolver : public QpPresolver {
private:
   const int my_rank = -1;
   /** limit for max rounds to apply all presolvers */
   const int limit_max_rounds = -1;
   /** should free variables' bounds be reset after presolve (given the row implying these bounds was not removed */
   const bool reset_free_variables_after_presolve = false;
   /** should the problem be written to std::cout before and after presolve */
   const bool print_problem = false;
   /** should the presolved problem be written out in MPS format */
   const bool write_presolved_problem = false;

   const int verbosity = -1;

   /* tree belonging to origData and preDistributedQP */
   DistributedTree* const tree;

   PresolveData preDistributedQP;
   std::vector<std::unique_ptr<StochPresolverBase>> presolvers;

   void resetFreeVariables();
public:

   StochPresolver(DistributedTree* tree, const Problem& prob, Postsolver* postsolver);
   ~StochPresolver() override = default;

   Problem* presolve() override;
};

//@}

#endif /* PIPS_IPM_CORE_QPPREPROCESS_STOCHPRESOLVER_H_ */
