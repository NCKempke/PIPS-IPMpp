/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SPSTOCHFACTORY
#define SPSTOCHFACTORY

#include "QpGen.h"

// save diagnostic state
#pragma GCC diagnostic push 
#pragma GCC diagnostic ignored "-Wsuggest-override"
#include "mpi.h"
// turn the warnings back 
#pragma GCC diagnostic pop

#include <cassert>

class QpGenData;
class sData;

class QpGenVars;
class StochInputTree;
//class stochasticInput;
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
  Residuals * makeResiduals( Data * prob_in ) override;
  Variables * makeVariables( Data * prob_in ) override;
  LinearSystem* makeLinsys( Data * prob_in ) override;

  virtual sLinsysRoot* newLinsysRootHierarchical() { assert( 0 && "not implemented here" ); return nullptr; }
  virtual Data* switchToHierarchicalData( Data* /*prob_in*/ ) { assert( 0 && "not implemented here" ); return nullptr; }

  virtual void collapseHierarchicalTree() { assert( 0 && "not implemented here" ); }

  void joinRHS( OoqpVector&, const OoqpVector&, const OoqpVector&, const OoqpVector&) const override
  { assert(0 && "not implemented here"); };

  void separateVars( OoqpVector&, OoqpVector&, OoqpVector&, const OoqpVector& ) const override
  { assert(0 && "not implemented here"); };

  virtual sLinsysRoot* newLinsysRoot() = 0;
  virtual sLinsysRoot* newLinsysRoot(sData* prob, sTree* tree_,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs) = 0;

  virtual sLinsysLeaf* newLinsysLeaf(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
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
};

#endif
