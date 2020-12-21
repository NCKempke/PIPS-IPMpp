#ifndef SPSTOCHQPGENFACTORY
#define SPSTOCHQPGENFACTORY

#include "QpGen.h"
#include "mpi.h"

class QpGenData;
class QpGenStochData;

class QpGenVars;
class StochInputTree;
class StochTree;
class StochSymMatrix;
class QpGenResiduals2;
class QpGenStochVars;
class QpGenStochLinsys;
class QpGenStochLinsysRoot;
class QpGenStochLinsysLeaf;
#include "StochResourcesMonitor.h"

class QpGenStoch : public QpGen {
 protected:
  int m_blocks;

  int nnzQ, nnzA, nnzC;

 public:
  QpGenStoch( StochInputTree* );
 protected:
  QpGenStoch( int nx_, int my_, int mz_, int nnzQ_, int nnzA_, int nnzC_ );
  QpGenStoch();

 public:

  virtual ~QpGenStoch();

  virtual Data  * makeData();

  virtual Residuals     * makeResiduals( Data * prob_in );
  virtual Variables     * makeVariables( Data * prob_in );

  virtual LinearSystem* makeLinsys( Data * prob_in );


  void joinRHS( OoqpVector& rhs_in, const  OoqpVector& rhs1_in,
			const OoqpVector& rhs2_in, const OoqpVector& rhs3_in ) const;

  void separateVars( OoqpVector& x_in, OoqpVector& y_in,
			     OoqpVector& z_in, const OoqpVector& vars_in ) const;

  virtual QpGenStochLinsysRoot* newLinsysRoot() = 0;
  virtual QpGenStochLinsysRoot* newLinsysRoot(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs) = 0;

  virtual QpGenStochLinsysLeaf* newLinsysLeaf(QpGenStochData* prob,
					      OoqpVector* dd,OoqpVector* dq,
					      OoqpVector* nomegaInv, OoqpVector* rhs);


  StochTree* tree;
  QpGenStochData * data;
  //  Variables

  virtual void iterateStarted();
  virtual void iterateEnded();

  QpGenResiduals2 *resid;
  std::vector<QpGenStochVars*> registeredVars;

  QpGenStochLinsys* linsys;

  StochIterateResourcesMonitor iterTmMonitor;
  double m_tmTotal;
};

#endif
