/*
 * sLinsysRootBordered.h
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */
#ifndef PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_
#define PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_

#include "DistributedRootLinearSystem.h"
#include <memory>

class sLinsysRootBordered : public DistributedRootLinearSystem {
public:
   sLinsysRootBordered(DistributedFactory* factory_, DistributedQP* prob_);

   ~sLinsysRootBordered() override = default;

   void
   add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization, double dual_inequality_regularization) override;

   void reset_regularization_local_kkt() override;

   void finalizeKKT() override;

   void Lsolve(Vector<double>& x) override;
   void Dsolve(Vector<double>& x) override;
   void Ltsolve(Vector<double>& v) override;

   void computeInnerSystemRightHandSide(DistributedVector<double>& rhs_inner, const SimpleVector<double>& x0, bool) override;
protected:
   SymmetricMatrix* createKKT();
   void assembleLocalKKT() override;
   void reduceKKT() override;

   DoubleLinearSolver* createSolver(const SymmetricMatrix* kktmat);
private:
   void computeSchurCompRightHandSide(const DistributedVector<double>& rhs_inner, SimpleVector<double>& b0);

};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_ */
