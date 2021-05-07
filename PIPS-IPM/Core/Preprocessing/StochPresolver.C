/*
 * StochPresolver.C
 *
 *  Created on: 26.01.2018
 *      Author: bzfrehfe
 */

#include "StochPresolver.h"
#include <cassert>
#include <iostream>
#include <string>
#include <memory>

#include "DistributedQP.hpp"
#include "DistributedTreeCallbacks.h"

#include "DistributedOptions.h"
#include "DistributedMatrix.h"
#include "Vector.hpp"
#include "StochPostsolver.h"
#include "StochPresolverBoundStrengthening.h"
#include "StochPresolverModelCleanup.h"
#include "StochPresolverColumnFixation.h"
#include "StochPresolverSingletonRows.h"
#include "StochPresolverSingletonColumns.h"
#include "StochPresolverParallelRows.h"

StochPresolver::StochPresolver(DistributedTree* tree_, const Problem& prob, Postsolver* postsolver = nullptr) : QpPresolver(prob, postsolver),
      my_rank(PIPS_MPIgetRank(MPI_COMM_WORLD)), limit_max_rounds(pips_options::get_int_parameter("PRESOLVE_MAX_ROUNDS")),
      reset_free_variables_after_presolve(pips_options::get_bool_parameter("PRESOLVE_RESET_FREE_VARIABLES")),
      print_problem(pips_options::get_bool_parameter("PRESOLVE_PRINT_PROBLEM")),
      write_presolved_problem(pips_options::get_bool_parameter("PRESOLVE_WRITE_PRESOLVED_PROBLEM_MPS")),
      verbosity(pips_options::get_int_parameter("PRESOLVE_VERBOSITY")), tree(tree_),
      preDistributedQP(dynamic_cast<const DistributedQP&>(origprob), dynamic_cast<StochPostsolver*>(postsolver)) {
   const DistributedQP& sorigprob = dynamic_cast<const DistributedQP&>(origprob);

   if (pips_options::get_bool_parameter("PRESOLVE_SINGLETON_ROWS"))
      presolvers.emplace_back(std::make_unique<StochPresolverSingletonRows>(preDistributedQP, sorigprob));

   if (pips_options::get_bool_parameter("PRESOLVE_COLUMN_FIXATION"))
      presolvers.emplace_back(std::make_unique<StochPresolverColumnFixation>(preDistributedQP, sorigprob));

   if (pips_options::get_bool_parameter("PRESOLVE_PARALLEL_ROWS"))
      presolvers.emplace_back(std::make_unique<StochPresolverParallelRows>(preDistributedQP, sorigprob));

   if (pips_options::get_bool_parameter("PRESOLVE_SINGLETON_COLUMNS"))
      presolvers.emplace_back(std::make_unique<StochPresolverSingletonColumns>(preDistributedQP, sorigprob));

   if (pips_options::get_bool_parameter("PRESOLVE_BOUND_STRENGTHENING"))
      presolvers.emplace_back(std::make_unique<StochPresolverBoundStrengthening>(preDistributedQP, sorigprob));
}

Problem* StochPresolver::presolve() {
   if (my_rank == 0)
      std::cout << "start stoch presolving\n";
   preDistributedQP.printRowColStats();

   const DistributedQP& sorigprob = dynamic_cast<const DistributedQP&>(origprob);
   sorigprob.printRanges();

   assert(sorigprob.isRootNodeInSync());
   assert(preDistributedQP.getPresProb().isRootNodeInSync());

   if (print_problem)
      sorigprob.writeToStreamDense(std::cout);

   /* initialize model clean up (necessary presolver) */
   StochPresolverModelCleanup presolverCleanup(preDistributedQP, sorigprob);

   if (my_rank == 0 && verbosity > 1)
      std::cout << "--- Before Presolving:\n";
   presolverCleanup.countRowsCols();

   // some while iterating over the list over and over until either every presolver says I'm done or some iterlimit is reached?
   presolverCleanup.applyPresolving();

   for (int i = 0; i < limit_max_rounds; ++i) {
      bool success = false;
      for (unsigned int i = 0; i < presolvers.size(); ++i) {
         bool presolver_success = presolvers[i]->applyPresolving();
         success = success || presolver_success;
      }

      presolverCleanup.applyPresolving();
   }

   if (my_rank == 0 && verbosity > 1)
      std::cout << "--- After Presolving:\n";
   presolverCleanup.countRowsCols();
   if (my_rank == 0)
      std::cout << "Objective offset: " << preDistributedQP.getObjOffset() << "\n";
   assert(preDistributedQP.getPresProb().isRootNodeInSync());

   if (reset_free_variables_after_presolve)
      resetFreeVariables();

   /* finalize data and switch tree to new presolved data */
   DistributedQP* finalPreDistributedQP = preDistributedQP.finalize();

   assert(tree != nullptr);
   assert(tree == finalPreDistributedQP->stochNode);

   DistributedTreeCallbacks& callbackTree = dynamic_cast<DistributedTreeCallbacks&>(*tree);
   callbackTree.initPresolvedData(*finalPreDistributedQP);
   callbackTree.switchToPresolvedData();

   /* change original bounds and set ixlow ixupp */
//   finalPreDistributedQP->xlowerBound().setNotIndicatedEntriesToVal( -1e10, *finalPreDistributedQP->ixlow );
//   finalPreDistributedQP->xupperBound().setNotIndicatedEntriesToVal( 1e10, *finalPreDistributedQP->ixupp );
//
//   Vector<double>* ixupp_inv = finalPreDistributedQP->ixupp->clone();
//   ixupp_inv->setToZero();
//   ixupp_inv->setNotIndicatedEntriesToVal(1.0, *finalPreDistributedQP->ixupp);
//
//   Vector<double>* ixlow_inv = finalPreDistributedQP->ixlow->clone();
//   ixlow_inv->setToZero();
//   ixlow_inv->setNotIndicatedEntriesToVal(1.0, *finalPreDistributedQP->ixlow);
//
//   finalPreDistributedQP->ixlow->setToConstant(1);
//   finalPreDistributedQP->ixupp->setToConstant(1);

   assert(finalPreDistributedQP);
   assert(finalPreDistributedQP->isRootNodeInSync());

   if (print_problem)
      finalPreDistributedQP->writeToStreamDense(std::cout);

   if (write_presolved_problem) {
      if (PIPS_MPIgetSize() > 1) {
         if (my_rank == 0)
            std::cout << "MPS format writer only available using one process!\n";
      }
      else {
         std::ofstream of("presolved.mps");

         if (of.is_open()) {
            //finalPreDistributedQP->writeMPSformat(of);
         }
         else if (my_rank == 0)
            std::cout << "Could not open presolved.mps to write out presolved problem!!\n";
      }
   }

   if (my_rank == 0)
      std::cout << "end stoch presolving\n";

   preDistributedQP.printRowColStats();
   finalPreDistributedQP->printRanges();

   return finalPreDistributedQP;
}

void StochPresolver::resetFreeVariables() {
   if (my_rank == 0)
      std::cout << "Resetting bounds found in bound strengthening\n";

   const DistributedQP& sorigprob = dynamic_cast<const DistributedQP&>(origprob);

   preDistributedQP.resetOriginallyFreeVarsBounds(sorigprob);
}
