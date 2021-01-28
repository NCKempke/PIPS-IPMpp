/*
 * sLinsysRootAugHierInner.h
 *
 *  Created on: 27.01.2021
 *      Author: bzfkempk
 */
#ifndef PIPS_IPM_CORE_STOCHLINEARALGEBRA_SLINSYSROOTAUGHIERINNER
#define PIPS_IPM_CORE_STOCHLINEARALGEBRA_SLINSYSROOTAUGHIERINNER

#include "sLinsysRootAug.h"

class sLinsysRootAugHierInner : public sLinsysRootAug
{
   public:
      sLinsysRootAugHierInner(sFactory* factory, sData* prob_, OoqpVector* dd_,
            OoqpVector* dq_, OoqpVector* nomegaInv_, OoqpVector* rhs_);

      ~sLinsysRootAugHierInner() override = default;

      void addTermToSchurComplBlocked(sData* prob, bool sparseSC, SymMatrix& SC, bool linking_only ) override;

      void assembleLocalKKT( sData* prob ) override;
      void putXDiagonal( OoqpVector& xdiag ) override;
      void putZDiagonal( OoqpVector& zdiag ) override;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_SLINSYSROOTAUGHIERINNER */
