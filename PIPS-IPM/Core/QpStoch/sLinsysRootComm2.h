/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SROOTLINSYSCOMM2
#define SROOTLINSYSCOMM2

#include "sLinsysRoot.h"

class sFactory;
class sData;

// DEBUG only
//#include "ScaDenSymMatrix.h"

/**
 * ROOT (= NON-leaf) linear system
 */
class sLinsysRootComm2 : public sLinsysRoot {
 public:
  sLinsysRootComm2(sFactory * factory_, sData * prob_);
  sLinsysRootComm2(sFactory* factory,
         sTree* tree_,
	      sData* prob_,
	      OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_,
	      OoqpVector* rhs_);

  SymMatrix* createKKT (sData* prob) override = 0;

  void factor2(sData *prob, Variables *vars) override;
  /* Atoms methods of FACTOR2 for a non-leaf linear system */
  using sLinsysRoot::reduceKKT;
  void reduceKKT() override;

  using sLinsysRoot::factorizeKKT;
  void factorizeKKT() override;

  void finalizeKKT(sData* prob, Variables* vars) override = 0;

  void Dsolve ( sData *prob, OoqpVector& x) override;
  void Lsolve ( sData *prob, OoqpVector& x ) override;
  void Ltsolve( sData *prob, OoqpVector& x ) override;
  void solveReduced( sData *prob, SimpleVector& b) override = 0;
 public:
  virtual ~sLinsysRootComm2();

  void submatrixReduce(DenseSymMatrix* A,
			  int row, int col, int drow, int dcol,
			  MPI_Comm comm);
 protected:
  void assembleLocalKKT( sData* prob ) override;
};

#endif
