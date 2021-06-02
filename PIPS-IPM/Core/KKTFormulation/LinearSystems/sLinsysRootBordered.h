/*
 * sLinsysRootBordered.h
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */
#ifndef PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_
#define PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_

#include "sLinsysRootAug.h"
#include <memory>

class sLinsysRootBordered : public sLinsysRootAug {
public:
   sLinsysRootBordered(DistributedFactory* factory_, DistributedQP* prob_);

   ~sLinsysRootBordered() override = default;

   void Lsolve(Vector<double>& x) override;
   void Dsolve(Vector<double>& x) override;
   void Ltsolve(Vector<double>& v) override;

   void computeInnerSystemRightHandSide(DistributedVector<double>& rhs_inner, const SimpleVector<double>& x0, bool) override;
protected:
   void assembleLocalKKT() override;

   DoubleLinearSolver* createSolver(const SymmetricMatrix* kktmat);
private:
   void computeSchurCompRightHandSide(const DistributedVector<double>& rhs_inner, SimpleVector<double>& b0);

};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_ */
