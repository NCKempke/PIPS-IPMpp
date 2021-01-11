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
  template<class LINSOLVER>
  sLinsysLeafSchurSlv(sFactory* factory,
		      sData* prob_,				    
		      OoqpVector* dd_, OoqpVector* dq_, 
		      OoqpVector* nomegaInv_,
		      OoqpVector* rhs_,
            OoqpVector* primal_reg,
            OoqpVector* dual_y_reg,
            OoqpVector* dual_z_reg,
            LINSOLVER* solver
         );

  void factor2(sData *prob, Variables *vars);
  void addTermToDenseSchurCompl(sData *prob, 
				DenseSymMatrix& SC);
  void addTermToSparseSchurCompl(sData *prob,
            SparseSymMatrix& SC);

 private:
  bool switchedToSafeSlv;

}; 
template<class LINSOLVER>
sLinsysLeafSchurSlv::sLinsysLeafSchurSlv(sFactory* factory,
					 sData* prob,
					 OoqpVector* dd_, 
					 OoqpVector* dq_, 
					 OoqpVector* nomegaInv_,
					 OoqpVector* rhs_,
					 OoqpVector* primal_reg,
					 OoqpVector* dual_y_reg,
					 OoqpVector* dual_z_reg,
					 LINSOLVER* s)
: sLinsysLeaf(factory, prob, dd_, dq_, nomegaInv_, rhs_, primal_reg, dual_y_reg, dual_z_reg, s),
  switchedToSafeSlv(false)
{}


#endif
