/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS_SCHURSLV
#define STOCHLEAFLINSYS_SCHURSLV

#include "sLinsysLeaf.h"

class StochTree;
class sFactory;
class DistributedQP;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeafSchurSlv : public sLinsysLeaf
{
 public:
  sLinsysLeafSchurSlv(sFactory* factory,
		      DistributedQP* prob_,				    
		      OoqpVector* dd_, OoqpVector* dq_, 
		      OoqpVector* nomegaInv_,
		      OoqpVector* primal_reg_,
		      OoqpVector* dual_y_reg_,
		      OoqpVector* dual_z_reg_,
		      OoqpVector* rhs_
         ) : sLinsysLeaf(factory, prob_, dd_, dq_, nomegaInv_, primal_reg_, dual_y_reg_, dual_z_reg_, rhs_) {};

  void factor2(DistributedQP *prob, Variables *vars) override;
  void addTermToDenseSchurCompl(DistributedQP *prob, 
				DenseSymMatrix& SC) override;
  void addTermToSparseSchurCompl(DistributedQP *prob,
            SparseSymMatrix& SC) override;

 private:
  bool switchedToSafeSlv{false};

}; 


#endif
