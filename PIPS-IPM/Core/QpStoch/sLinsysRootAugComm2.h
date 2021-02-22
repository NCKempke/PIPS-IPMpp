/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYSAUGCOMM2
#define SAUGLINSYSAUGCOMM2

#include "sLinsysRootComm2.h"
class sData;
/** 
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAugComm2 : public sLinsysRootComm2 {
 protected:
  sLinsysRootAugComm2() {};

  SymMatrix* createKKT (sData* prob) override;
  DoubleLinearSolver* createSolver (sData* prob, SymMatrix* kktmat) override;
 public:

  sLinsysRootAugComm2(sFactory * factory_, sData * prob_);
  sLinsysRootAugComm2(sFactory* factory,
		 sData* prob_,				    
		 OoqpVector* dd_, OoqpVector* dq_, 
		 OoqpVector* nomegaInv_,
		 OoqpVector* rhs_);
  virtual ~sLinsysRootAugComm2();

  void finalizeKKT( sData* prob, Variables* vars) override;
 protected:
  void solveReduced( sData *prob, SimpleVector& b) override;
  void solveReducedLinkCons( sData *prob, SimpleVector& b) override;
  void solveWithIterRef( sData *prob, SimpleVector& b);
  void solveWithBiCGStab( sData *prob, SimpleVector& b);

  /** y = beta*y - alpha* SC * x */
  void SCmult ( double beta, SimpleVector& y, double alpha, SimpleVector& x, sData* prob);

  SymMatrix* CtDC;
  SimpleVector* redRhs;
};

#endif

