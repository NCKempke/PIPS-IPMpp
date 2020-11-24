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

class sLinsysRootBordered : public sLinsysRootAug
{
   public:
      sLinsysRootBordered(sFactory * factory_, sData * prob_);

      virtual ~sLinsysRootBordered() = default;

      void solveReduced( sData *prob, SimpleVector& b) override
      { assert ( false && "should not end up here" ); };

      void Lsolve(sData *prob, OoqpVector& x) override;
      void Dsolve(sData *prob, OoqpVector& x) override;
      void Ltsolve(sData *prob, OoqpVector& v) override;

   protected:
      void assembleLocalKKT(sData* prob) override;
      void reduceKKT(sData* prob) override;

   private:
      void computeSchurCompRightHandSide( const StochVector& rhs_inner, SimpleVector& b0 );
      void computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& x0 );

      std::unique_ptr<StochVector> sol_inner{};
};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_ */
