/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS_SCHURSLV
#define STOCHLEAFLINSYS_SCHURSLV

#include "DistributedLeafLinearSystem.h"

class StochTree;

class DistributedFactory;

class DistributedQP;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeafSchurSlv : public DistributedLeafLinearSystem {
public:
   sLinsysLeafSchurSlv(DistributedFactory* factory, DistributedQP* prob_, Vector<double>* dd_, Vector<double>* dq_, Vector<double>* nomegaInv_,
         Vector<double>* primal_reg_, Vector<double>* dual_y_reg_, Vector<double>* dual_z_reg_, Vector<double>* rhs_) : DistributedLeafLinearSystem(factory, prob_,
         dd_, dq_, nomegaInv_, primal_reg_, dual_y_reg_, dual_z_reg_, rhs_) {};

   void factor2(DistributedQP* prob, Variables* vars) override;
   void addTermToDenseSchurCompl(DistributedQP* prob, DenseSymMatrix& SC) override;
   void addTermToSparseSchurCompl(DistributedQP* prob, SparseSymMatrix& SC) override;

private:
   bool switchedToSafeSlv{false};

};


#endif
