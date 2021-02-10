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

class sLinsysRootBordered : public sLinsysRoot
{
   public:
      sLinsysRootBordered(sFactory * factory_, sData * prob_);

      ~sLinsysRootBordered() override = default;

      void finalizeKKT(sData* prob, Variables* ) override;

      void Lsolve(sData*, OoqpVector& x) override;
      void Dsolve(sData*, OoqpVector& x) override;
      void Ltsolve(sData*, OoqpVector& v) override;

      void computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& x0 ) override;
   protected:
      SymMatrix* createKKT(sData*);
      void assembleLocalKKT(sData* prob) override;
      void reduceKKT(sData*) override;

      DoubleLinearSolver* createSolver(sData*, SymMatrix* kktmat);
   private:
      void computeSchurCompRightHandSide( const StochVector& rhs_inner, SimpleVector& b0 );

};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_ */
