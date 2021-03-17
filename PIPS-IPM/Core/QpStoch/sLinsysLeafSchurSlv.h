/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS_SCHURSLV
#define STOCHLEAFLINSYS_SCHURSLV

#include "sLinsysLeaf.h"

class StochTree;
class sFactory;
class sData;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeafSchurSlv : public sLinsysLeaf
{
 public:
  sLinsysLeafSchurSlv(sFactory* factory,
		      sData* prob_,				    
		      OoqpVector* dd_, OoqpVector* dq_, 
		      OoqpVector* nomegaInv_,
		      OoqpVector* rhs_) : sLinsysLeaf(factory, prob_, dd_, dq_, nomegaInv_, rhs_) {};

  void factor2(sData *prob, Variables *vars) override;
  void addTermToDenseSchurCompl(sData *prob, 
				DenseSymMatrix& SC) override;
  void addTermToSparseSchurCompl(sData *prob,
            SparseSymMatrix& SC) override;

 private:
  bool switchedToSafeSlv{false};

}; 

#endif
