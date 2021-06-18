/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS_SCHURSLV
#define STOCHLEAFLINSYS_SCHURSLV

#include "DistributedLeafLinearSystem.h"

class DistributedFactory;

class DistributedQP;

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to QpGenSparseLinsys.
 */
class sLinsysLeafSchurSlv : public DistributedLeafLinearSystem {
public:
   sLinsysLeafSchurSlv(const DistributedFactory& factory, DistributedQP* prob_, std::shared_ptr<Vector<double>> dd_,
      std::shared_ptr<Vector<double>> dq_, std::shared_ptr<Vector<double>> nomegaInv_,
      std::shared_ptr<Vector<double>> primal_reg_, std::shared_ptr<Vector<double>> dual_y_reg_,
      std::shared_ptr<Vector<double>> dual_z_reg_, std::shared_ptr<Vector<double>> rhs_) : DistributedLeafLinearSystem(
      factory, prob_,
      std::move(dd_), std::move(dq_), std::move(nomegaInv_), std::move(primal_reg_), std::move(dual_y_reg_), std::move(dual_z_reg_), std::move(rhs_)) {};

   void factor2() override;

   void addTermToDenseSchurCompl(DenseSymmetricMatrix& SC) override;

   void addTermToSparseSchurCompl(SparseSymmetricMatrix& SC) override;

private:
   bool switchedToSafeSlv{false};

};


#endif
