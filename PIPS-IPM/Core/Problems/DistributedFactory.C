/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include <memory>
#include "DistributedFactory.h"
#include "DistributedQP.hpp"
#include "DistributedTreeCallbacks.h"
#include "DistributedInputTree.h"
#include "DistributedSymmetricMatrix.h"
#include "DistributedMatrix.h"
#include "DistributedVector.h"
#include "DistributedVariables.h"
#include "DistributedResiduals.hpp"
#include "DistributedRootLinearSystem.h"
#include "DistributedLeafLinearSystem.h"
#include "PIPSIPMppOptions.h"
#include "mpi.h"
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


DistributedFactory::DistributedFactory(DistributedInputTree* inputTree, MPI_Comm comm) : tree(
   new DistributedTreeCallbacks(inputTree)) {
   tree->assignProcesses(comm);
   tree->computeGlobalSizes();
   // now the sizes of the problem are available, set them for the parent class
   tree->getGlobalSizes(nx, my, mz);
}

// TODO make tree unique_ptr
DistributedFactory::~DistributedFactory() {
   delete tree;
}

DoubleLinearSolver* DistributedFactory::make_leaf_solver(const AbstractMatrix* kkt_) {
   const auto& kkt = dynamic_cast<const SparseSymmetricMatrix&>(*kkt_);

   const SolverType leaf_solver = pipsipmpp_options::get_solver_leaf();

   if (!pipsipmpp_options::get_bool_parameter("SC_COMPUTE_BLOCKWISE")) {
      if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new MumpsSolverLeaf(kkt);
#endif
      } else if (leaf_solver == SolverType::SOLVER_PARDISO) {
#ifdef WITH_PARDISO
         return new PardisoProjectSchurSolver(kkt);
#endif
      } else if (leaf_solver == SolverType::SOLVER_MKL_PARDISO) {
#ifdef WITH_MKL_PARDISO
         return new PardisoMKLSchurSolver(kkt);
#endif
      }

      PIPS_MPIabortIf(true, "No leaf solver for Schur Complement computation could be found - should not happen..");
   } else {
      if (leaf_solver == SolverType::SOLVER_PARDISO) {
#ifdef WITH_PARDISO
         return new PardisoProjectSolver(kkt);
#endif
      } else if (leaf_solver == SolverType::SOLVER_MKL_PARDISO) {
#ifdef WITH_MKL_PARDISO
         return new PardisoMKLSolver(kkt);
#endif
      } else if (leaf_solver == SolverType::SOLVER_MA57) {
#ifdef WITH_MA57
         return new Ma57Solver(kkt);
#endif
      } else if (leaf_solver == SolverType::SOLVER_MA27) {
#ifdef WITH_MA27
         return new Ma27Solver(kkt);
#endif
      } else if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new MumpsSolverLeaf(kkt);
#endif
      }

      PIPS_MPIabortIf(true,
         "No leaf solver for Blockwise Schur Complement computation could be found - should not happen..");
   }
   return nullptr;
}

DistributedLeafLinearSystem*
DistributedFactory::make_linear_system_leaf(DistributedQP* problem, std::shared_ptr<Vector<double>> primal_diagonal,
   std::shared_ptr<Vector<double>> dq,
   std::shared_ptr<Vector<double>> nomegaInv,
   std::shared_ptr<Vector<double>> primal_regularization, std::shared_ptr<Vector<double>> dual_equality_regularization,
   std::shared_ptr<Vector<double>> dual_inequality_regularization,
   std::shared_ptr<Vector<double>> rhs) {
   assert(problem);
   assert(primal_diagonal);
   static bool printed = false;
   const SolverType leaf_solver = pipsipmpp_options::get_solver_leaf();

   if (!pipsipmpp_options::get_bool_parameter("SC_COMPUTE_BLOCKWISE")) {
      if (PIPS_MPIgetRank() == 0 && !printed)
         std::cout << "Using " << leaf_solver << " for leaf schur complement computation - sFactory\n";
      printed = true;

      if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new sLinsysLeafMumps(this, problem, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
#endif
      } else if (leaf_solver == SolverType::SOLVER_PARDISO || leaf_solver == SolverType::SOLVER_MKL_PARDISO) {
#if defined(WITH_PARDISO) or defined(WITH_MKL_PARDISO)
         return new sLinsysLeafSchurSlv(this, problem, std::move(primal_diagonal), std::move(dq), std::move(nomegaInv),
            std::move(primal_regularization),
            std::move(dual_equality_regularization),
            std::move(dual_inequality_regularization), std::move(rhs));
#endif
      } else {
         if (PIPS_MPIgetRank() == 0) {
            std::cout << "WARNING: did not specify SC_COMPUTE_BLOCKWISE but " << leaf_solver
                      << " which can only compute the Schur Complement blockwise\n";
            std::cout << "WARNING: checking for PARDISO or MUMPS instead...\n";
         }

         SolverType solver = SolverType::SOLVER_NONE;
         if (pipsipmpp_options::is_solver_available(SolverType::SOLVER_PARDISO))
            solver = SolverType::SOLVER_PARDISO;
         else if (pipsipmpp_options::is_solver_available(SolverType::SOLVER_MKL_PARDISO))
            solver = SolverType::SOLVER_MKL_PARDISO;
         else if (pipsipmpp_options::is_solver_available(SolverType::SOLVER_MUMPS))
            solver = SolverType::SOLVER_MUMPS;

         if (solver != SolverType::SOLVER_NONE) {
            if (PIPS_MPIgetRank() == 0)
               std::cout << " Found solver " << solver << " - using that for leaf computations\n";
            pipsipmpp_options::set_int_parameter("LINEAR_LEAF_SOLVER", solver);
#if defined(WITH_PARDISO) or defined(WITH_MKL_PARDISO)
            return new sLinsysLeafSchurSlv(this, problem, std::move(primal_diagonal), std::move(dq),
               std::move(nomegaInv), std::move(primal_regularization),
               std::move(dual_equality_regularization),
               std::move(dual_inequality_regularization), std::move(rhs));
#endif
         }

         PIPS_MPIabortIf(true, "Error: Could not find suitable solver - please specify SC_COMPUTE_BLOCKWISE");
         return nullptr;
      }
   } else {
      if (PIPS_MPIgetRank() == 0 && !printed)
         std::cout << "Using " << leaf_solver
                   << " for blockwise Schur Complement computation - deactivating distributed preconditioner - sFactory\n";
      pipsipmpp_options::set_bool_parameter("PRECONDITION_DISTRIBUTED", false);
      printed = true;

      if (leaf_solver == SolverType::SOLVER_MUMPS) {
#ifdef WITH_MUMPS
         return new sLinsysLeafMumps(this, problem, primal_diagonal, dq, nomegaInv, primal_regularization, dual_equality_regularization, dual_inequality_regularization, rhs);
#endif
      } else
         return new DistributedLeafLinearSystem(this, problem, std::move(primal_diagonal), std::move(dq),
            std::move(nomegaInv), std::move(primal_regularization),
            std::move(dual_equality_regularization),
            std::move(dual_inequality_regularization), std::move(rhs));
   }
   return nullptr;
}

Problem* DistributedFactory::make_problem() const {
#ifdef TIMING
   double t2 = MPI_Wtime();
#endif

   std::shared_ptr<GeneralMatrix> A(tree->createA());
   std::shared_ptr<DistributedVector<double>> b(tree->createb());

   std::shared_ptr<GeneralMatrix> C(tree->createC());
   std::shared_ptr<DistributedVector<double>> clow(tree->createclow());
   std::shared_ptr<DistributedVector<double>> iclow(tree->createiclow());
   std::shared_ptr<DistributedVector<double>> cupp(tree->createcupp());
   std::shared_ptr<DistributedVector<double>> icupp(tree->createicupp());

   std::shared_ptr<SymmetricMatrix> Q(tree->createQ());
   std::shared_ptr<DistributedVector<double>> c(tree->createc());

   std::shared_ptr<DistributedVector<double>> xlow(tree->createxlow());
   std::shared_ptr<DistributedVector<double>> ixlow(tree->createixlow());
   std::shared_ptr<DistributedVector<double>> xupp(tree->createxupp());
   std::shared_ptr<DistributedVector<double>> ixupp(tree->createixupp());

#ifdef TIMING
   MPI_Barrier( tree->getCommWrkrs() );
   t2 = MPI_Wtime() - t2;
   if ( PIPS_MPIgetRank(tree->getCommWorkers() == 0) )
         std::cout << "IO second part took " << t2 << " sec\n";
#endif

   return new DistributedQP(tree, std::move(c), std::move(Q), std::move(xlow), std::move(ixlow), std::move(xupp),
      std::move(ixupp), std::move(A), std::move(b), std::move(C), std::move(clow), std::move(iclow), std::move(cupp),
      std::move(icupp));
}

Variables* DistributedFactory::make_variables(Problem& problem) const {

   std::unique_ptr<Vector<double>> x(make_primal_vector());
   std::unique_ptr<Vector<double>> s(make_inequalities_dual_vector());
   std::unique_ptr<Vector<double>> y(make_equalities_dual_vector());
   std::unique_ptr<Vector<double>> z(make_inequalities_dual_vector());
   std::unique_ptr<Vector<double>> v(make_primal_vector());
   std::unique_ptr<Vector<double>> gamma(make_primal_vector());
   std::unique_ptr<Vector<double>> w(make_primal_vector());
   std::unique_ptr<Vector<double>> phi(make_primal_vector());
   std::unique_ptr<Vector<double>> t(make_inequalities_dual_vector());
   std::unique_ptr<Vector<double>> lambda(make_inequalities_dual_vector());
   std::unique_ptr<Vector<double>> u(make_inequalities_dual_vector());
   std::unique_ptr<Vector<double>> pi(make_inequalities_dual_vector());

   assert(problem.ixlow && problem.ixupp && problem.iclow && problem.icupp);
   auto* variables = new DistributedVariables(tree, std::move(x), std::move(s), std::move(y), std::move(z),
      std::move(v), std::move(gamma), std::move(w), std::move(phi), std::move(t), std::move(lambda), std::move(u),
      std::move(pi), problem.ixlow, problem.ixlow->number_nonzeros(), problem.ixupp, problem.ixupp->number_nonzeros(),
      problem.iclow, problem.iclow->number_nonzeros(), problem.icupp, problem.icupp->number_nonzeros());
   return variables;
}

Residuals* DistributedFactory::make_residuals(Problem& problem) const {

   std::unique_ptr<Vector<double>> lagrangian_gradient{tree->new_primal_vector()};

   std::unique_ptr<Vector<double>> rA{tree->newDualYVector()};
   std::unique_ptr<Vector<double>> rC{tree->newDualZVector()};

   std::unique_ptr<Vector<double>> rz{tree->newDualZVector()};

   const bool mclow_empty = problem.mclow <= 0;
   std::unique_ptr<Vector<double>> rt{tree->newDualZVector(mclow_empty)};
   std::unique_ptr<Vector<double>> rlambda{tree->newDualZVector(mclow_empty)};

   const bool mcupp_empty = problem.mcupp <= 0;
   std::unique_ptr<Vector<double>> ru{tree->newDualZVector(mcupp_empty)};
   std::unique_ptr<Vector<double>> rpi{tree->newDualZVector(mcupp_empty)};

   const bool nxlow_empty = problem.nxlow <= 0;
   std::unique_ptr<Vector<double>> rv{tree->new_primal_vector(nxlow_empty)};
   std::unique_ptr<Vector<double>> rgamma{tree->new_primal_vector(nxlow_empty)};

   const bool nxupp_empty = problem.nxupp <= 0;
   std::unique_ptr<Vector<double>> rw{tree->new_primal_vector(nxupp_empty)};
   std::unique_ptr<Vector<double>> rphi{tree->new_primal_vector(nxupp_empty)};

   assert(problem.ixlow && problem.ixupp && problem.iclow && problem.icupp);
   return new DistributedResiduals(std::move(lagrangian_gradient), std::move(rA), std::move(rC), std::move(rz),
      std::move(rt), std::move(rlambda), std::move(ru), std::move(rpi), std::move(rv), std::move(rgamma), std::move(rw),
      std::move(rphi), problem.ixlow, problem.ixupp, problem.iclow, problem.icupp);
}

std::unique_ptr<AbstractLinearSystem> DistributedFactory::make_linear_system(Problem& problem_in) {
   auto* problem = dynamic_cast<DistributedQP*>(&problem_in);
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      return make_root_hierarchical_linear_system(problem);
   else
      return make_linear_system_root(problem);
}

Vector<double>* DistributedFactory::make_primal_vector() const {
   return tree->new_primal_vector();
}

Vector<double>* DistributedFactory::make_equalities_dual_vector() const {
   return tree->newDualYVector();
}

Vector<double>* DistributedFactory::make_inequalities_dual_vector() const {
   return tree->newDualZVector();
}

Vector<double>* DistributedFactory::make_right_hand_side() const {
   return tree->newRhs();
}

void DistributedFactory::iterate_started() {
   timer.start();
   tree->startMonitors();
}

void DistributedFactory::iterate_ended() {
   tree->stopMonitors();

   if (tree->balanceLoad()) {
      printf("Should not get here! OMG OMG OMG\n");
   }
   //logging and monitoring
   timer.stop();
   total_time += timer.end_time;

   if (PIPS_MPIgetRank() == 0) {
#ifdef TIMING
      extern double g_iterNumber;
      printf("TIME %g SOFAR %g ITER %d\n", iterTmMonitor.tmIterate, m_tmTotal, (int)g_iterNumber);
      //#elseif STOCH_TESTING
      //printf("ITERATION WALLTIME: iter=%g  Total=%g\n", iterTmMonitor.tmIterate, m_tmTotal);
#endif
   }
}

std::unique_ptr<DistributedRootLinearSystem> DistributedFactory::make_linear_system_root(DistributedQP* problem) {
   assert(problem);
   return std::make_unique<sLinsysRootAug>(this, problem);
}

DistributedRootLinearSystem*
DistributedFactory::make_linear_system_root(DistributedQP* prob, std::shared_ptr<Vector<double>> primal_diagonal,
   std::shared_ptr<Vector<double>> dq,
   std::shared_ptr<Vector<double>> nomegaInv,
   std::shared_ptr<Vector<double>> primal_regularization, std::shared_ptr<Vector<double>> dual_equality_regularization,
   std::shared_ptr<Vector<double>> dual_inequality_regularization,
   std::shared_ptr<Vector<double>> rhs) {
   if (prob->isHierarchyInnerLeaf())
      return new sLinsysRootAugHierInner(this, prob, std::move(primal_diagonal), std::move(dq), std::move(nomegaInv),
         std::move(primal_regularization),
         std::move(dual_equality_regularization),
         std::move(dual_inequality_regularization), std::move(rhs));
   else
      return new sLinsysRootAug(this, prob, std::move(primal_diagonal), std::move(dq), std::move(nomegaInv),
         std::move(primal_regularization),
         std::move(dual_equality_regularization),
         std::move(dual_inequality_regularization), rhs, true);
}

std::unique_ptr<DistributedRootLinearSystem>
DistributedFactory::make_root_hierarchical_linear_system(DistributedQP* problem) {
   return std::make_unique<sLinsysRootBordered>(this, problem);
}

DistributedQP* DistributedFactory::switchToHierarchicalData(DistributedQP* problem) {
   hier_tree_swap.reset(tree->clone());

   tree = tree->switchToHierarchicalTree(problem);

   assert(tree->getChildren().size() == 1);
   assert(tree->isHierarchicalRoot());

   assert(problem->isHierarchyRoot());
   return problem;
}

void DistributedFactory::switchToOriginalTree() {
   assert(hier_tree_swap);

   DistributedTree* tmp = tree;
   tree = hier_tree_swap.get();
   hier_tree_swap.release();
   hier_tree_swap.reset(tmp);
}
