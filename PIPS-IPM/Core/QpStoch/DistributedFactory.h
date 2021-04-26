/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SPSTOCHFACTORY
#define SPSTOCHFACTORY

#include "ProblemFormulation.h"
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

class sVars;

class sLinsys;

class sLinsysRoot;

class sLinsysLeaf;

class DoubleLinearSolver;

class DoubleMatrix;

#include "StochResourcesMonitor.h"

class DistributedFactory : public ProblemFormulation {
public:

   DistributedFactory(StochInputTree*, MPI_Comm comm = MPI_COMM_WORLD);

protected:
   DistributedFactory() = default;

   ~DistributedFactory() override;

public:

   virtual Problem* create_problem();

   Residuals* makeResiduals(Problem* prob_in) override;

   Variables* makeVariables(Problem* prob_in) override;

   LinearSystem* makeLinsys(Problem* prob_in) override;

   /** create x shaped vector using tree */
   OoqpVector* makePrimalVector() const override;

   /** create dual A shaped vector using tree */
   OoqpVector* makeDualYVector() const override;

   /** create dual C shaped vector using tree */
   OoqpVector* makeDualZVector() const override;

   /** create rhs for augmented system using tree */
   OoqpVector* makeRhs() const override;

   virtual sLinsysRoot* newLinsysRootHierarchical() {
      assert(0 && "not implemented here");
      return nullptr;
   }

   virtual Problem* switchToHierarchicalData(Problem* /*prob_in*/ ) {
      assert(0 && "not implemented here");
      return nullptr;
   }

   virtual void switchToOriginalTree() { assert(0 && "not implemented here"); }

   void joinRHS(OoqpVector&, const OoqpVector&, const OoqpVector&, const OoqpVector&) const override { assert(0 && "not implemented here"); };

   void separateVars(OoqpVector&, OoqpVector&, OoqpVector&, const OoqpVector&) const override { assert(0 && "not implemented here"); };

   virtual sLinsysRoot* newLinsysRoot() = 0;

   virtual sLinsysRoot* newLinsysRoot(DistributedQP* prob, OoqpVector* dd, OoqpVector* dq, OoqpVector* nomegaInv, OoqpVector* rhs) = 0;

   virtual sLinsysLeaf* newLinsysLeaf(DistributedQP* prob, OoqpVector* dd, OoqpVector* dq, OoqpVector* nomegaInv, OoqpVector* rhs);

   virtual DoubleLinearSolver* newRootSolver() = 0;

   virtual DoubleLinearSolver* newLeafSolver(const DoubleMatrix* kkt);

   sTree* tree{};
   DistributedQP* data{};

   virtual void iterateStarted();

   virtual void iterateEnded();

   DistributedResiduals* resid{};
   std::vector<sVars*> registeredVars;

   sLinsysRoot* linsys{};

   StochIterateResourcesMonitor iterTmMonitor;
   double m_tmTotal{0.0};

protected:
   std::unique_ptr<sTree> hier_tree_swap{};
};

#endif
