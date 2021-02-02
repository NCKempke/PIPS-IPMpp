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

      void assembleLocalKKT( sData* prob ) override;

      void finalizeZ0Hierarchical( DenseGenMatrix& buffer, BorderLinsys& ) override;

      void addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, bool use_local_RAC_mat ) override;

      void addTermToSchurComplBlocked( sData* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC_mat ) override;

      void putXDiagonal( OoqpVector& xdiag ) override;
      void putZDiagonal( OoqpVector& zdiag ) override;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_SLINSYSROOTAUGHIERINNER */
