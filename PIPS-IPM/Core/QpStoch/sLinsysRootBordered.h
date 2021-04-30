/*
 * sLinsysRootBordered.h
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */
#ifndef PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_
#define PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_

#include "sLinsysRoot.h"
#include <memory>

class sLinsysRootBordered : public sLinsysRoot {
public:
   sLinsysRootBordered(DistributedFactory* factory_, DistributedQP* prob_);

   ~sLinsysRootBordered() override = default;

      void add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization, double dual_inequality_regularization) override;

      void finalizeKKT(DistributedQP* prob, Variables*) override;

   void Lsolve(DistributedQP*, OoqpVector& x) override;
   void Dsolve(DistributedQP*, OoqpVector& x) override;
   void Ltsolve(DistributedQP*, OoqpVector& v) override;

   void computeInnerSystemRightHandSide(DistributedVector<double>& rhs_inner, const SimpleVector<double>& x0, bool) override;
protected:
   SymMatrix* createKKT(DistributedQP*);
   void assembleLocalKKT(DistributedQP* prob) override;
   void reduceKKT(DistributedQP*) override;

   DoubleLinearSolver* createSolver(DistributedQP*, const SymMatrix* kktmat);
private:
   void computeSchurCompRightHandSide(const DistributedVector<double>& rhs_inner, SimpleVector<double>& b0);

};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_ */
