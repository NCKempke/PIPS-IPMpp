/*
 * sLinsysRootBordered.h
 *
 *  Created on: Sep 22, 2020
 *      Author: bzfkempk
 */
#ifndef PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_
#define PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_

#include "sLinsysRoot.h"

class sLinsysRootBordered : public sLinsysRoot
{
   public:
      sLinsysRootBordered(sFactory * factory_, sData * prob_);
      sLinsysRootBordered(sFactory* factory, sData* prob_, OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_, OoqpVector* rhs_)
         { assert("should not be called" && 0); };

      virtual ~sLinsysRootBordered();

      void finalizeKKT(sData* prob, Variables* vars) override;
      void solveReduced( sData *prob, SimpleVector& b) override;

   protected:

      SymMatrix* createKKT(sData* prob) override;
      DoubleLinearSolver* createSolver(sData* prob, SymMatrix* kktmat) override;


   private:

};

#endif /* PIPS_IPM_CORE_QPSTOCH_SLINSYSROOTBORDERED_H_ */
