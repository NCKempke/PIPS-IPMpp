/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SPSTOCHFACTORY
#define SPSTOCHFACTORY

#include "QpGen.h"
#include "mpi.h"
#include <cassert>
#include <memory>

class QpGenData;
class sData;

class QpGenVars;
class StochInputTree;
class sTree;
class StochSymMatrix;
class sResiduals;
class sVars;
class sLinsys;
class sLinsysRoot;
class sLinsysLeaf;

#include "StochResourcesMonitor.h"

class sFactory : public QpGen
{
   public:

      sFactory( StochInputTree*, MPI_Comm comm = MPI_COMM_WORLD );

   protected:
      sFactory() = default;
      ~sFactory() override;

 public:

  virtual Data* makeData();
  virtual Residuals * makeResiduals( Data * prob_in );
  virtual Variables * makeVariables( Data * prob_in );
  virtual LinearSystem* makeLinsys( Data * prob_in );

  /** create x shaped vector using tree */
  OoqpVector* makePrimalVector() const override;
  /** create dual A shaped vector using tree */
  OoqpVector* makeDualYVector() const override;
  /** create dual C shaped vector using tree */
  OoqpVector* makeDualZVector() const override;
  /** create rhs for augmented system using tree */
  OoqpVector* makeRhs() const override;


  virtual sLinsysRoot* newLinsysRootHierarchical() { assert( 0 && "not implemented here" ); return nullptr; }
  virtual Data* switchToHierarchicalData( Data* /*prob_in*/ ) { assert( 0 && "not implemented here" ); return nullptr; }

  virtual void switchToOriginalTree() { assert( 0 && "not implemented here" ); }


  void joinRHS( OoqpVector&, const OoqpVector&, const OoqpVector&, const OoqpVector&) const override
  { assert(0 && "not implemented here"); };


  void separateVars( OoqpVector&, OoqpVector&, OoqpVector&, const OoqpVector& ) const override
  { assert(0 && "not implemented here"); };

  virtual sLinsysRoot* newLinsysRoot() = 0;
  virtual sLinsysRoot* newLinsysRoot(sData* prob, OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs) = 0;
  virtual sLinsysLeaf* newLinsysLeaf(sData* prob, OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);



  sTree * tree{};
  sData * data{};

  virtual void iterateStarted();
  virtual void iterateEnded();

  sResiduals *resid{};
  std::vector<sVars*> registeredVars;

  sLinsysRoot* linsys{};

  StochIterateResourcesMonitor iterTmMonitor;
  double m_tmTotal{0.0};

 protected:
  std::unique_ptr<sTree> hier_tree_swap{};
};

#endif
