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

      void Lsolve(sData *prob, OoqpVector& x) override;
      void Ltsolve( sData *prob, OoqpVector& x ) override;
      void Ltsolve2(sData*, StochVector& x, SimpleVector& x0, bool use_local_RAC ) override;

      using sLinsys::LsolveHierarchyBorder;
      void LsolveHierarchyBorder( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool two_link_border ) override;

      using sLinsys::LtsolveHierarchyBorder;
      void LtsolveHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
            std::vector<BorderMod>& br_mod_border, bool sym_res, bool sparse_res, bool two_link_border ) override;

      void computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& b0, bool use_local_RAC ) override;

      void addLniziLinkCons( sData *prob, OoqpVector& z0, OoqpVector& zi, bool use_local_RAC ) override;
      void addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border ) override;
      void addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border ) override;

      void addInnerBorderKiInvBrToRes( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool use_local_RAC ) override;

      void addTermToSchurComplBlocked( sData* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC ) override;
      void LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
            bool sparse_res, bool sym_res, bool use_local_RAC ) override;

      void putXDiagonal( OoqpVector& xdiag ) override;
      void putZDiagonal( OoqpVector& zdiag ) override;

   private:
    void createSolversAndKKts(sData* prob);
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_SLINSYSROOTAUGHIERINNER */
