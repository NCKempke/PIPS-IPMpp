/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SPSTOCHFACTORY
#define SPSTOCHFACTORY

#include "ProblemFactory.h"
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

class DistributedLinearSystem;

class DistributedRootLinearSystem;

class DistributedLeafLinearSystem;

class DoubleLinearSolver;

class DoubleMatrix;

#include "StochResourcesMonitor.h"

class DistributedFactory : public ProblemFactory {
public:

   DistributedFactory(StochInputTree* tree, MPI_Comm comm = MPI_COMM_WORLD);
   virtual Problem* make_problem();

   Residuals* make_residuals(Problem& problem) override;

   Variables* make_variables(Problem& problem) override;

   AbstractLinearSystem* make_linear_system(Problem& problem) override;

   /** create x shaped vector using tree */
   Vector<double>* make_primal_vector() const override;

   /** create dual A shaped vector using tree */
   Vector<double>* make_equalities_dual_vector() const override;

   /** create dual C shaped vector using tree */
   Vector<double>* make_inequalities_dual_vector() const override;

   /** create rhs for augmented system using tree */
   Vector<double>* make_right_hand_side() const override;

   DoubleLinearSolver* make_root_solver();

   DistributedRootLinearSystem* make_linear_system_root();

   DistributedRootLinearSystem* newLinsysRootHierarchical();

   DistributedRootLinearSystem* make_linear_system_root(DistributedQP* problem, Vector<double>* primal_diagonal, Vector<double>* dq, Vector<double>* nomegaInv,
         Vector<double>* primal_regularization, Vector<double>* dual_equality_regularization, Vector<double>* dual_inequality_regularization,
         Vector<double>* rhs);

   Problem* switchToHierarchicalData(Problem* problem);

   void switchToOriginalTree();

   void join_right_hand_side(Vector<double>&, const Vector<double>&, const Vector<double>&, const Vector<double>&) const override {
      assert(0 && "not implemented here");
   };

   void separate_variables(Vector<double>&, Vector<double>&, Vector<double>&, const Vector<double>&) const override {
      assert(0 && "not implemented here");
   };

   virtual DistributedLeafLinearSystem*
   make_linear_system_leaf(DistributedQP* problem, Vector<double>* primal_diagonal, Vector<double>* dq, Vector<double>* nomegaInv,
         Vector<double>* primal_regularization, Vector<double>* dual_equality_regularization, Vector<double>* dual_inequality_regularization,
         Vector<double>* rhs);

   virtual DoubleLinearSolver* make_leaf_solver(const DoubleMatrix* kkt);

   DistributedTree* tree{};
   DistributedQP* problem{};

   virtual void iterate_started();

   virtual void iterate_ended();

   DistributedResiduals* residuals{};
   std::vector<DistributedVariables*> registeredVars;

   DistributedRootLinearSystem* linsys{};

   StochIterateResourcesMonitor iterTmMonitor;
   double m_tmTotal{0.0};

   ~DistributedFactory() override;

protected:
   std::unique_ptr<DistributedTree> hier_tree_swap{};
   DistributedFactory() = default;
};

#endif
