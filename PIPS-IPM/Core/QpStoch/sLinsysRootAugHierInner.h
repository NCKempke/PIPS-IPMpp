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
      sLinsysRootAugHierInner(sFactory* factory, DistributedQP* prob_, OoqpVector* dd_,
            OoqpVector* dq_, OoqpVector* nomegaInv_, OoqpVector* regP_, OoqpVector* regDy_, OoqpVector* regDz_, OoqpVector* rhs_);
      ~sLinsysRootAugHierInner() override = default;

      void assembleLocalKKT( DistributedQP* prob ) override;

      void Lsolve(DistributedQP *prob, OoqpVector& x) override;
      void Ltsolve( DistributedQP *prob, OoqpVector& x ) override;
      void Ltsolve2(DistributedQP*, StochVector& x, SimpleVector& x0, bool use_local_RAC ) override;

      using sLinsys::LsolveHierarchyBorder;
      void LsolveHierarchyBorder( DenseGenMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool two_link_border, int begin_cols, int end_cols ) override;

      using sLinsys::LtsolveHierarchyBorder;
      void LtsolveHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
            std::vector<BorderMod>& br_mod_border, bool sym_res, bool sparse_res, int begin_cols, int end_cols ) override;

      void computeInnerSystemRightHandSide( StochVector& rhs_inner, const SimpleVector& b0, bool use_local_RAC ) override;

      void addLniziLinkCons( DistributedQP *prob, OoqpVector& z0, OoqpVector& zi, bool use_local_RAC ) override;
      void addBorderTimesRhsToB0( StochVector& rhs, SimpleVector& b0, BorderLinsys& border ) override;
      void addBorderX0ToRhs( StochVector& rhs, const SimpleVector& x0, BorderLinsys& border ) override;

      void addInnerBorderKiInvBrToRes( DoubleMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool use_local_RAC, bool sparse_res, bool sym_res, int begin_cols, int end_cols,
            int n_empty_rows_inner_border ) override;

      void addTermToSchurComplBlocked( DistributedQP* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC, int n_empty_rows_inner_border ) override;
      void LniTransMultHierarchyBorder( DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
            bool sparse_res, bool sym_res, bool use_local_RAC, int begin_cols, int end_cols, int n_empty_rows_inner_border ) override;

      void putXDiagonal( const OoqpVector& xdiag ) override;
      void putZDiagonal( const OoqpVector& zdiag ) override;

   private:
    void createSolversAndKKts(DistributedQP* prob);
};

#endif /* PIPS_IPM_CORE_STOCHLINEARALGEBRA_SLINSYSROOTAUGHIERINNER */
