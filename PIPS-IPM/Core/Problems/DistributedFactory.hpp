/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef DISTRIBUTEDFACTORY_H
#define DISTRIBUTEDFACTORY_H

#include "DistributedTree.h"
#include "ProblemFactory.h"

// save diagnostic state
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"

#include "mpi.h"
// turn the warnings back 
#pragma GCC diagnostic pop

#include <cassert>
#include <memory>

class QP;

class DistributedProblem;

class Variables;

class DistributedInputTree;

class DistributedSymmetricMatrix;

class DistributedResiduals;

class DistributedVariables;

class DistributedRootLinearSystem;

class DistributedLeafLinearSystem;

class DoubleLinearSolver;

class AbstractMatrix;

class Problem;

class Residuals;

class AbstractLinearSystem;

class RegularizationStrategy;

#include "StochResourcesMonitor.hpp"

class DistributedFactory : public ProblemFactory {
public:
   explicit DistributedFactory(DistributedInputTree* tree, MPI_Comm comm = MPI_COMM_WORLD);

   [[nodiscard]] std::unique_ptr<Problem> make_problem() const override;

   [[nodiscard]] std::unique_ptr<Residuals> make_residuals(const Problem& problem) const override;

   [[nodiscard]] std::unique_ptr<Variables> make_variables(const Problem& problem) const override;

   [[nodiscard]] std::unique_ptr<AbstractLinearSystem> make_linear_system(Problem& problem) const override;

   [[nodiscard]] std::unique_ptr<RegularizationStrategy> make_regularization_strategy(unsigned int positive_eigenvalues, unsigned int negative_eigenvalues) const override;

   /** create x shaped vector using tree */
   [[nodiscard]] std::unique_ptr<Vector<double>> make_primal_vector() const override;

   /** create dual A shaped vector using tree */
   [[nodiscard]] std::unique_ptr<Vector<double>> make_equalities_dual_vector() const override;

   /** create dual C shaped vector using tree */
   [[nodiscard]] std::unique_ptr<Vector<double>> make_inequalities_dual_vector() const override;

   /** create rhs for augmented system using tree */
   [[nodiscard]] std::unique_ptr<Vector<double>> make_right_hand_side() const override;

   [[nodiscard]] std::unique_ptr<DistributedRootLinearSystem> make_linear_system_root(DistributedProblem* problem) const;

   [[nodiscard]] std::unique_ptr<DistributedRootLinearSystem> make_root_hierarchical_linear_system(DistributedProblem* problem) const;

   [[nodiscard]] std::unique_ptr<DistributedRootLinearSystem>
   make_linear_system_root(DistributedProblem* problem, std::shared_ptr<Vector<double>> primal_diagonal,
      std::shared_ptr<Vector<double>> dq, std::shared_ptr<Vector<double>> nomegaInv,
      std::shared_ptr<Vector<double>> primal_regularization,
      std::shared_ptr<Vector<double>> dual_equality_regularization,
      std::shared_ptr<Vector<double>> dual_inequality_regularization,
      std::shared_ptr<Vector<double>> rhs) const;

   DistributedProblem* switchToHierarchicalData(DistributedProblem* problem);

   void switchToOriginalTree();

   [[nodiscard]] std::unique_ptr<DistributedLeafLinearSystem>
   make_linear_system_leaf(DistributedProblem* problem, std::shared_ptr<Vector<double>> primal_diagonal, std::shared_ptr<Vector<double>> dq,
      std::shared_ptr<Vector<double>> nomegaInv,
      std::shared_ptr<Vector<double>> primal_regularization, std::shared_ptr<Vector<double>> dual_equality_regularization,
      std::shared_ptr<Vector<double>> dual_inequality_regularization,
      std::shared_ptr<Vector<double>> rhs) const;

   static std::unique_ptr<DoubleLinearSolver> make_leaf_solver(const AbstractMatrix* kkt);

   std::unique_ptr<DistributedTree> tree{};

   void iterate_started();

   void iterate_ended();

   DistributedResiduals* residuals{};

   //DistributedRootLinearSystem* linear_system{};

   Timer timer;
   double total_time{0.0};

   ~DistributedFactory() override = default;

protected:
   /** number of elements in x */
   long long nx{0};

   /** number of rows in A and b including linking rows (sFactory..) */
   long long my{0};

   /** number of rows in C including linking rows */
   long long mz{0};

   std::unique_ptr<DistributedTree> hier_tree_swap{};

   DistributedFactory() = default;
};

#endif
