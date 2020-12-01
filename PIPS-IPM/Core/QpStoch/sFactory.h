/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SPSTOCHFACTORY
#define SPSTOCHFACTORY

#include "QpGen.h"
#include "mpi.h"
#include <cassert>

class QpGenData;
class sData;

class QpGenVars;
class StochInputTree;
class stochasticInput;
class sTree;
class StochSymMatrix;
class sResiduals;
class sVars;
class sLinsys;
class sLinsysRoot;
class sLinsysLeaf;

#include "StochResourcesMonitor.h"

class sFactory : public QpGen {
 public:
  sFactory( stochasticInput&, MPI_Comm comm = MPI_COMM_WORLD );

  /** This is a obsolete constructor since it uses sTreeCallbacks to create
   *   data objects
   */
  sFactory( StochInputTree*, MPI_Comm comm = MPI_COMM_WORLD );

 protected:
  sFactory();

 public:

  virtual ~sFactory();

  virtual Data* makeData();
  virtual Residuals * makeResiduals( Data * prob_in );
  virtual Variables * makeVariables( Data * prob_in );
  virtual LinearSystem* makeLinsys( Data * prob_in );

  virtual sLinsysRoot* newLinsysRootHierarchical() { assert( 0 && "not implemented here" ); return nullptr; }
  virtual Data* switchToHierarchicalData( Data* prob_in ) { assert( 0 && "not implemented here" ); return nullptr; }

  virtual void collapseHierarchicalTree() { assert( 0 && "not implemented here" ); }


  void joinRHS( OoqpVector& rhs_in, const OoqpVector& rhs1_in,
			const OoqpVector& rhs2_in, const OoqpVector& rhs3_in ) const override;

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
	      OoqpVector& z_in, const OoqpVector& vars_in ) const override;

  virtual sLinsysRoot* newLinsysRoot() = 0;
  virtual sLinsysRoot* newLinsysRoot(sData* prob, sTree* tree_,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs) = 0;

  virtual sLinsysLeaf* newLinsysLeaf(sData* prob,
				     OoqpVector* dd,OoqpVector* dq,
				     OoqpVector* nomegaInv, OoqpVector* rhs);


  sTree * tree;
  sData * data;

  virtual void iterateStarted();
  virtual void iterateEnded();
  void writeProblemToStream(ostream& out, bool printRhs) const;

  sResiduals *resid;
  vector<sVars*> registeredVars;

  sLinsysRoot* linsys;

  StochIterateResourcesMonitor iterTmMonitor;
  double m_tmTotal;
};

#endif
