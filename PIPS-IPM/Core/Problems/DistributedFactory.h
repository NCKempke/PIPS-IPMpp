/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef DISTRIBUTEDFACTORY_H
#define DISTRIBUTEDFACTORY_H

#include "DistributedTree.h"

// save diagnostic state
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsuggest-override"

#include "mpi.h"
// turn the warnings back 
#pragma GCC diagnostic pop

#include <cassert>
#include <memory>

class QP;

class DistributedQP;

class Variables;

class StochInputTree;

class StochSymMatrix;

class DistributedResiduals;

class DistributedVariables;

class DistributedRootLinearSystem;

class DistributedLeafLinearSystem;

class DoubleLinearSolver;

class DoubleMatrix;

class Problem;

class Residuals;

class AbstractLinearSystem;

#include "StochResourcesMonitor.h"

class DistributedFactory {
public:
   DistributedFactory(StochInputTree* tree, MPI_Comm comm = MPI_COMM_WORLD);
   Problem* make_problem();

   Residuals* make_residuals(Problem& problem);

   Variables* make_variables(Problem& problem);

   AbstractLinearSystem* make_linear_system(Problem& problem);

   /** create x shaped vector using tree */
   Vector<double>* make_primal_vector() const;

   /** create dual A shaped vector using tree */
   Vector<double>* make_equalities_dual_vector() const;

   /** create dual C shaped vector using tree */
   Vector<double>* make_inequalities_dual_vector() const;

   /** create rhs for augmented system using tree */
   Vector<double>* make_right_hand_side() const;

   DistributedRootLinearSystem* make_linear_system_root();

   DistributedRootLinearSystem* make_root_hierarchical_linear_system();

   DistributedRootLinearSystem* make_linear_system_root(DistributedQP* problem, Vector<double>* primal_diagonal, Vector<double>* dq, Vector<double>* nomegaInv,
         Vector<double>* primal_regularization, Vector<double>* dual_equality_regularization, Vector<double>* dual_inequality_regularization,
         Vector<double>* rhs);

   Problem* switchToHierarchicalData(Problem* problem);

   void switchToOriginalTree();

   DistributedLeafLinearSystem*
   make_linear_system_leaf(DistributedQP* problem, Vector<double>* primal_diagonal, Vector<double>* dq, Vector<double>* nomegaInv,
         Vector<double>* primal_regularization, Vector<double>* dual_equality_regularization, Vector<double>* dual_inequality_regularization,
         Vector<double>* rhs);

   DoubleLinearSolver* make_leaf_solver(const DoubleMatrix* kkt);

   DistributedTree* tree{};
   DistributedQP* problem{};

   void iterate_started();

   void iterate_ended();

   DistributedResiduals* residuals{};
   std::vector<DistributedVariables*> registered_variables;

   DistributedRootLinearSystem* linear_system{};

   StochIterateResourcesMonitor iterTmMonitor;
   double m_tmTotal{0.0};

   ~DistributedFactory();

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
