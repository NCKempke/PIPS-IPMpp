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

      void finalizeKKT(sData* prob, Variables* vars) override;
      void solveReduced( sData *prob, SimpleVector& b) override
      { assert ( false && "should not end up here" ); };

      void Lsolve(sData *prob, OoqpVector& x) override;
      void Dsolve(sData *prob, OoqpVector& x) override;
      void Ltsolve(sData *prob, OoqpVector& v) override;

   protected:
      SymMatrix* createKKT(sData* prob) override;
      void assembleLocalKKT(sData* prob) override;
      void reduceKKT(sData* prob) override;

      DoubleLinearSolver* createSolver(sData* prob, SymMatrix* kktmat) override;

   private:
      void computeSchurCompRightHandSide( const StochVector& rhs_inner, SimpleVector& b0 );
      void computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& x0 );

      std::unique_ptr<StochVector> sol_inner{};
};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_ */
