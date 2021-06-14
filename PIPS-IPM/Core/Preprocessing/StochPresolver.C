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
#include "PIPSIPMppOptions.h"
#include "DistributedMatrix.h"
#include "Vector.hpp"
#include "StochPostsolver.h"
#include "StochPresolverBoundStrengthening.h"
#include "StochPresolverModelCleanup.h"
#include "StochPresolverColumnFixation.h"
#include "StochPresolverSingletonRows.h"
#include "StochPresolverSingletonColumns.h"
#include "StochPresolverParallelRows.h"

StochPresolver::StochPresolver(DistributedTree* tree_, const Problem& prob, Postsolver* postsolver = nullptr) : Presolver(prob, postsolver),
      my_rank(PIPS_MPIgetRank(MPI_COMM_WORLD)), limit_max_rounds(pipsipmpp_options::get_int_parameter("PRESOLVE_MAX_ROUNDS")),
      reset_free_variables_after_presolve(pipsipmpp_options::get_bool_parameter("PRESOLVE_RESET_FREE_VARIABLES")),
      print_problem(pipsipmpp_options::get_bool_parameter("PRESOLVE_PRINT_PROBLEM")),
      write_presolved_problem(pipsipmpp_options::get_bool_parameter("PRESOLVE_WRITE_PRESOLVED_PROBLEM_MPS")),
      transform_inequalities_to_equalities{pipsipmpp_options::get_bool_parameter("PRESOLVE_TRANSFROM_INEQUALITIES_INTO_EQUALITIES")},
      verbosity(pipsipmpp_options::get_int_parameter("PRESOLVE_VERBOSITY")), tree(tree_), original_problem(dynamic_cast<const DistributedQP&>(prob)),
      presolve_data(dynamic_cast<const DistributedQP&>(original_problem), dynamic_cast<StochPostsolver*>(postsolver)) {
   const auto& sorigprob = dynamic_cast<const DistributedQP&>(original_problem);

   if (pipsipmpp_options::get_bool_parameter("PRESOLVE_SINGLETON_ROWS"))
      presolvers.emplace_back(std::make_unique<StochPresolverSingletonRows>(presolve_data, sorigprob));

   if (pipsipmpp_options::get_bool_parameter("PRESOLVE_COLUMN_FIXATION"))
      presolvers.emplace_back(std::make_unique<StochPresolverColumnFixation>(presolve_data, sorigprob));

   if (pipsipmpp_options::get_bool_parameter("PRESOLVE_PARALLEL_ROWS"))
      presolvers.emplace_back(std::make_unique<StochPresolverParallelRows>(presolve_data, sorigprob));

   if (pipsipmpp_options::get_bool_parameter("PRESOLVE_SINGLETON_COLUMNS"))
      presolvers.emplace_back(std::make_unique<StochPresolverSingletonColumns>(presolve_data, sorigprob));

   if (pipsipmpp_options::get_bool_parameter("PRESOLVE_BOUND_STRENGTHENING"))
      presolvers.emplace_back(std::make_unique<StochPresolverBoundStrengthening>(presolve_data, sorigprob));
}

Problem* StochPresolver::presolve() {
   if (my_rank == 0)
      std::cout << "starting distributed presolving\n";
   presolve_data.printRowColStats();
   original_problem.printRanges();

   assert(original_problem.isRootNodeInSync());
   assert(presolve_data.getPresProb().isRootNodeInSync());
   if (print_problem)
      original_problem.write_to_streamDense(std::cout);

   run_presolve_loop();

   if (reset_free_variables_after_presolve) {
      resetFreeVariables();
   }

   if (transform_inequalities_to_equalities) {
      if (my_rank == 0)
         std::cout << "Transforming all inequality rows to equalities + slack variable...";
      presolve_data.transfrom_ineqalities_to_equalities();
      if (my_rank == 0)
         std::cout << " done\n";
   }

   /* finalize data and switch tree to new presolved data */
   auto* finalPreDistributedQP = dynamic_cast<DistributedQP*>(presolve_data.finalize());

   assert(tree != nullptr);
   assert(tree == finalPreDistributedQP->stochNode);

   auto& callbackTree = dynamic_cast<DistributedTreeCallbacks&>(*tree);
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
      finalPreDistributedQP->write_to_streamDense(std::cout);

   if (write_presolved_problem) {
      this->write_presolved_problem_to_file();
   }

   if (my_rank == 0)
      std::cout << "finished distributed presolving\n";
   presolve_data.printRowColStats();
   finalPreDistributedQP->printRanges();

   return finalPreDistributedQP;
}

void StochPresolver::resetFreeVariables() {
   if (my_rank == 0)
      std::cout << "Resetting bounds found in bound strengthening\n";

   const auto& sorigprob = dynamic_cast<const DistributedQP&>(original_problem);

   presolve_data.resetOriginallyFreeVarsBounds(sorigprob);
}

void StochPresolver::write_presolved_problem_to_file() {
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

void StochPresolver::run_presolve_loop() {
   /* initialize model clean up (necessary presolver) */
   StochPresolverModelCleanup presolverCleanup(presolve_data, original_problem);

   if (my_rank == 0 && verbosity > 1)
      std::cout << "--- Before Presolving:\n";
   presolverCleanup.countRowsCols();

   // some while iterating over the list over and over until either every presolver says I'm done or some iterlimit is reached?
   presolverCleanup.applyPresolving();

   for (int i = 0; i < limit_max_rounds; ++i) {
      bool success = false;
      for (auto & presolver : presolvers) {
         bool presolver_success = presolver->applyPresolving();
         success = success || presolver_success;
      }

      presolverCleanup.applyPresolving();
   }

   if (my_rank == 0 && verbosity > 1)
      std::cout << "--- After Presolving:\n";
   presolverCleanup.countRowsCols();
   if (my_rank == 0)
      std::cout << "Objective offset: " << presolve_data.getObjOffset() << "\n";
   assert(presolve_data.getPresProb().isRootNodeInSync());
}