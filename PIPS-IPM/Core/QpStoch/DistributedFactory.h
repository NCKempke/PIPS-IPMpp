/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SPSTOCHFACTORY
#define SPSTOCHFACTORY

#include "ProblemFactory.h"
#include "sTree.h"

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

class sLinsysRoot;

class sLinsysLeaf;

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
   OoqpVector* make_primal_vector() const override;

   /** create dual A shaped vector using tree */
   OoqpVector* make_equalities_dual_vector() const override;

   /** create dual C shaped vector using tree */
   OoqpVector* make_inequalities_dual_vector() const override;

   /** create rhs for augmented system using tree */
   OoqpVector* make_right_hand_side() const override;

   DoubleLinearSolver* make_root_solver();

   sLinsysRoot* make_linear_system_root();

   sLinsysRoot* newLinsysRootHierarchical();

   sLinsysRoot* make_linear_system_root(DistributedQP* problem, OoqpVector* dd, OoqpVector* dq, OoqpVector* nomegaInv, OoqpVector* regP, OoqpVector* regDy, OoqpVector* regDz, OoqpVector* rhs);

   Problem* switchToHierarchicalData(Problem* problem);

   void switchToOriginalTree();

   void join_right_hand_side(OoqpVector&, const OoqpVector&, const OoqpVector&, const OoqpVector&) const override { assert(0 && "not implemented here"); };

   void separate_variables(OoqpVector&, OoqpVector&, OoqpVector&, const OoqpVector&) const override { assert(0 && "not implemented here"); };

   virtual sLinsysLeaf* make_linear_system_leaf(DistributedQP* problem, OoqpVector* dd, OoqpVector* dq, OoqpVector* nomegaInv, OoqpVector* regP, OoqpVector* regDy, OoqpVector* regDz, OoqpVector* rhs);

   virtual DoubleLinearSolver* make_leaf_solver(const DoubleMatrix* kkt);

   sTree* tree{};
   DistributedQP* problem{};

   virtual void iterate_started();

   virtual void iterate_ended();

   DistributedResiduals* residuals{};
   std::vector<DistributedVariables*> registeredVars;

   sLinsysRoot* linsys{};

   StochIterateResourcesMonitor iterTmMonitor;
   double m_tmTotal{0.0};

   ~DistributedFactory() override;

protected:
   std::unique_ptr<sTree> hier_tree_swap{};
   DistributedFactory() = default;
};

#endif
