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

      void finalizeZ0Hierarchical( DenseGenMatrix& buffer, BorderLinsys&, std::vector<BorderMod>& Br_mod_border ) override;

      void addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border ) override;

      void addTermToSchurComplBlocked( sData* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC ) override;
      void LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
            int parent_nx, int parent_my, int parent_mz, bool sparse_res, bool sym_res ) override;

      void putXDiagonal( OoqpVector& xdiag ) override;
      void putZDiagonal( OoqpVector& zdiag ) override;
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_SLINSYSROOTAUGHIERINNER */
