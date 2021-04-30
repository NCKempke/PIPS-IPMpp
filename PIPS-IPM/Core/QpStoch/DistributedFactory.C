/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "DistributedFactory.h"
#include "DistributedQP.hpp"
#include "sTreeCallbacks.h"
#include "StochInputTree.h"
#include "StochSymMatrix.h"
#include "StochGenMatrix.h"
#include "StochVector.h"
#include "DistributedVariables.h"
#include "DistributedResiduals.hpp"
#include "sLinsysRoot.h"
#include "sLinsysLeaf.h"
#include "StochOptions.h"
#include "mpi.h"
#include <stdio.h>
#include "sLinsysRootAug.h"
#include "sLinsysRootAugHierInner.h"
#include "sLinsysRootBordered.h"

class PardisoProjectSolver;
class PardisoMKLSolver;
class MumpsSolverLeaf;
class PardisoSchurSolver;
class Ma27Solver;
class Ma57Solver;

#ifdef WITH_MA57
#include "Ma57Solver.h"
#endif

#ifdef WITH_MA27
#include "Ma27Solver.h"
#endif

#ifdef WITH_PARDISO
#include "sLinsysLeafSchurSlv.h"
#include "PardisoProjectSolver.h"
#include "../LinearSolvers/PardisoSolver/PardisoSchurSolver/PardisoProjectSchurSolver.h"

#endif

#ifdef WITH_MKL_PARDISO
#include "sLinsysLeafSchurSlv.h"
#include "PardisoMKLSolver.h"
#include "../LinearSolvers/PardisoSolver/PardisoSchurSolver/PardisoMKLSchurSolver.h"
#endif

#ifdef WITH_MUMPS
#include "sLinsysLeafMumps.h"
#include "../LinearSolvers/MumpsSolver/MumpsSolverLeaf.h"
#endif


DistributedFactory::DistributedFactory(StochInputTree* inputTree, MPI_Comm comm) : tree(new sTreeCallbacks(inputTree)) {
   tree->assignProcesses(comm);

   tree->computeGlobalSizes();
   // now the sizes of the problem are available, set them for the parent class
   tree->getGlobalSizes(nx, my, mz);
}

DistributedFactory::~DistributedFactory() {
   if (tree)
      delete tree;
}

DoubleLinearSolver* DistributedFactory::make_leaf_solver(const DoubleMatrix* kkt_) {
   const SparseSymMatrix* kkt = dynamic_cast<const SparseSymMatrix*>(kkt_);
   assert(kkt);

   const SolverType leaf_solver = pips_options::getSolverLeaf();

   if (!pips_options::getBoolParameter("SC_COMPUTE_BLOCKWISE")) {
      if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new MumpsSolverLeaf(kkt);
#endif
      }
      else if (leaf_solver == SolverType::SOLVER_PARDISO) {
#ifdef WITH_PARDISO
         return new PardisoProjectSchurSolver(kkt);
#endif
      }
      else if (leaf_solver == SolverType::SOLVER_MKL_PARDISO) {
#ifdef WITH_MKL_PARDISO
         return new PardisoMKLSchurSolver(kkt);
#endif
      }

      PIPS_MPIabortIf(true, "No leaf solver for Schur Complement computation could be found - should not happen..");
   }
   else {
      if (leaf_solver == SolverType::SOLVER_PARDISO) {
#ifdef WITH_PARDISO
         return new PardisoProjectSolver(kkt);
#endif
      }
      else if (leaf_solver == SolverType::SOLVER_MKL_PARDISO) {
#ifdef WITH_MKL_PARDISO
         return new PardisoMKLSolver(kkt);
#endif
      }
      else if (leaf_solver == SolverType::SOLVER_MA57) {
#ifdef WITH_MA57
         return new Ma57Solver(kkt);
#endif
      }
      else if (leaf_solver == SolverType::SOLVER_MA27) {
#ifdef WITH_MA27
         return new Ma27Solver(kkt);
#endif
      }
      else if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new MumpsSolverLeaf(kkt);
#endif
      }

      PIPS_MPIabortIf(true, "No leaf solver for Blockwise Schur Complement computation could be found - should not happen..");
   }
   return nullptr;
}

sLinsysLeaf* DistributedFactory::make_linear_system_leaf(DistributedQP* problem, OoqpVector* primal_diagonal, OoqpVector* dq, OoqpVector* nomegaInv,
   OoqpVector* primal_regularization, OoqpVector* dual_equality_regularization, OoqpVector* dual_inequality_regularization, OoqpVector* rhs) {
   assert(problem);
   static bool printed = false;
   const SolverType leaf_solver = pips_options::getSolverLeaf();

   if (!pips_options::getBoolParameter("SC_COMPUTE_BLOCKWISE")) {
      if (PIPS_MPIgetRank() == 0 && !printed)
         std::cout << "Using " << leaf_solver << " for leaf schur complement computation - sFactory\n";
      printed = true;

      if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new sLinsysLeafMumps(this, problem, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
#endif
      }
      else if (leaf_solver == SolverType::SOLVER_PARDISO || leaf_solver == SolverType::SOLVER_MKL_PARDISO) {
#if defined(WITH_PARDISO) or defined(WITH_MKL_PARDISO)
         return new sLinsysLeafSchurSlv(this, problem, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
#endif
      }
      else {
         if (PIPS_MPIgetRank() == 0) {
            std::cout << "WARNING: did not specify SC_COMPUTE_BLOCKWISE but " << leaf_solver
                      << " which can only compute the Schur Complement blockwise\n";
            std::cout << "WARNING: checking for PARDISO or MUMPS instead...\n";
         }

         SolverType solver = SolverType::SOLVER_NONE;
         if (pips_options::isSolverAvailable(SolverType::SOLVER_PARDISO))
            solver = SolverType::SOLVER_PARDISO;
         else if (pips_options::isSolverAvailable(SolverType::SOLVER_MKL_PARDISO))
            solver = SolverType::SOLVER_MKL_PARDISO;
         else if (pips_options::isSolverAvailable(SolverType::SOLVER_MUMPS))
            solver = SolverType::SOLVER_MUMPS;

         if (solver != SolverType::SOLVER_NONE) {
            if (PIPS_MPIgetRank() == 0)
               std::cout << " Found solver " << solver << " - using that for leaf computations\n";
            pips_options::setIntParameter("LINEAR_LEAF_SOLVER", solver);
#if defined(WITH_PARDISO) or defined(WITH_MKL_PARDISO)
            return new sLinsysLeafSchurSlv(this, problem, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
#endif
         }

         PIPS_MPIabortIf(true, "Error: Could not find suitable solver - please specify SC_COMPUTE_BLOCKWISE");
         return nullptr;
      }
   }
   else {
      if (PIPS_MPIgetRank() == 0 && !printed)
         std::cout << "Using " << leaf_solver << " for blockwise Schur Complement computation - deactivating distributed preconditioner - sFactory\n";
      pips_options::setBoolParameter("PRECONDITION_DISTRIBUTED", false);
      printed = true;

      if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new sLinsysLeafMumps(this, problem, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
#endif
      }
      else
         return new sLinsysLeaf(this, problem, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
   }
   return nullptr;
}

Problem* DistributedFactory::make_problem() {
#ifdef TIMING
   double t2 = MPI_Wtime();
#endif

   StochGenMatrixHandle A(tree->createA());
   StochVectorHandle b(tree->createb());

   StochGenMatrixHandle C(tree->createC());
   StochVectorHandle clow(tree->createclow());
   StochVectorHandle iclow(tree->createiclow());
   StochVectorHandle cupp(tree->createcupp());
   StochVectorHandle icupp(tree->createicupp());

   StochSymMatrixHandle Q(tree->createQ());
   StochVectorHandle c(tree->createc());

   StochVectorHandle xlow(tree->createxlow());
   StochVectorHandle ixlow(tree->createixlow());
   StochVectorHandle xupp(tree->createxupp());
   StochVectorHandle ixupp(tree->createixupp());

#ifdef TIMING
   MPI_Barrier( tree->getCommWrkrs() );
   t2 = MPI_Wtime() - t2;
   if ( PIPS_MPIgetRank(tree->getCommWorkers() == 0) )
         std::cout << "IO second part took " << t2 << " sec\n";
#endif

   this->problem = new DistributedQP(tree, c, Q, xlow, ixlow, xupp, ixupp, A, b, C, clow, iclow, cupp, icupp);
   return this->problem;
}

// TODO adjust this for hierarchical approach
Variables* DistributedFactory::make_variables(Problem& problem) {
   OoqpVectorHandle x = OoqpVectorHandle(make_primal_vector());
   OoqpVectorHandle s = OoqpVectorHandle(make_inequalities_dual_vector());
   OoqpVectorHandle y = OoqpVectorHandle(make_equalities_dual_vector());
   OoqpVectorHandle z = OoqpVectorHandle(make_inequalities_dual_vector());
   OoqpVectorHandle v = OoqpVectorHandle(make_primal_vector());
   OoqpVectorHandle gamma = OoqpVectorHandle(make_primal_vector());
   OoqpVectorHandle w = OoqpVectorHandle(make_primal_vector());
   OoqpVectorHandle phi = OoqpVectorHandle(make_primal_vector());
   OoqpVectorHandle t = OoqpVectorHandle(make_inequalities_dual_vector());
   OoqpVectorHandle lambda = OoqpVectorHandle(make_inequalities_dual_vector());
   OoqpVectorHandle u = OoqpVectorHandle(make_inequalities_dual_vector());
   OoqpVectorHandle pi = OoqpVectorHandle(make_inequalities_dual_vector());

   DistributedVariables* variables = new DistributedVariables(tree, x, s, y, z, v, gamma, w, phi, t, lambda, u, pi, problem.ixlow, problem.ixlow->numberOfNonzeros(), problem.ixupp,
         problem.ixupp->numberOfNonzeros(), problem.iclow, problem.iclow->numberOfNonzeros(), problem.icupp, problem.icupp->numberOfNonzeros());
   registeredVars.push_back(variables);
   return variables;
}

Residuals* DistributedFactory::make_residuals(Problem& problem) {
   residuals = new DistributedResiduals(tree, problem.ixlow, problem.ixupp, problem.iclow, problem.icupp);
   return residuals;
}

AbstractLinearSystem* DistributedFactory::make_linear_system(Problem&) {
   if (pips_options::getBoolParameter("HIERARCHICAL"))
      linsys = newLinsysRootHierarchical();
   else
      linsys = make_linear_system_root();
   return linsys;
}

OoqpVector* DistributedFactory::make_primal_vector() const {
   assert(!la);
   return tree->newPrimalVector();
}

OoqpVector* DistributedFactory::make_equalities_dual_vector() const {
   assert(!la);
   return tree->newDualYVector();
}

OoqpVector* DistributedFactory::make_inequalities_dual_vector() const {
   assert(!la);
   return tree->newDualZVector();
}

OoqpVector* DistributedFactory::make_right_hand_side() const {
   assert(!la);
   return tree->newRhs();
}

void DistributedFactory::iterate_started() {
   iterTmMonitor.recIterateTm_start();
   tree->startMonitors();
}

void DistributedFactory::iterate_ended() {
   tree->stopMonitors();

   if (tree->balanceLoad()) {
      printf("Should not get here! OMG OMG OMG\n");
   }
   //logging and monitoring
   iterTmMonitor.recIterateTm_stop();
   m_tmTotal += iterTmMonitor.tmIterate;

   if (PIPS_MPIgetRank() == 0) {
#ifdef TIMING
      extern double g_iterNumber;
      printf("TIME %g SOFAR %g ITER %d\n", iterTmMonitor.tmIterate, m_tmTotal, (int)g_iterNumber);
      //#elseif STOCH_TESTING
      //printf("ITERATION WALLTIME: iter=%g  Total=%g\n", iterTmMonitor.tmIterate, m_tmTotal);
#endif
   }
}

sLinsysRoot* DistributedFactory::make_linear_system_root() {
   assert(problem);
   return new sLinsysRootAug(this, problem);
}

sLinsysRoot* DistributedFactory::make_linear_system_root(DistributedQP* prob, OoqpVector* primal_diagonal, OoqpVector* dq, OoqpVector* nomegaInv,
   OoqpVector* primal_regularization, OoqpVector* dual_equality_regularization, OoqpVector* dual_inequality_regularization, OoqpVector* rhs) {
   if (prob->isHierarchyInnerLeaf())
      return new sLinsysRootAugHierInner(this, prob, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
   else
      return new sLinsysRootAug(this, prob, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs, true);
}

DoubleLinearSolver* DistributedFactory::make_root_solver() { return nullptr; };

sLinsysRoot* DistributedFactory::newLinsysRootHierarchical() {
   return new sLinsysRootBordered(this, problem);
}

Problem* DistributedFactory::switchToHierarchicalData(Problem*) {
   hier_tree_swap.reset(tree->clone());

   tree = tree->switchToHierarchicalTree(problem);

   assert(tree->getChildren().size() == 1);
   assert(tree->isHierarchicalRoot());

   assert(problem->isHierarchyRoot());
   return problem;
}

void DistributedFactory::switchToOriginalTree() {
   assert(hier_tree_swap);

   sTree* tmp = tree;
   tree = hier_tree_swap.get();
   hier_tree_swap.release();
   hier_tree_swap.reset(tmp);
}
