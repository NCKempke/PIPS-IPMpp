/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef STOCHLEAFLINSYS
#define STOCHLEAFLINSYS

#include "DistributedLinearSystem.h"
#include "sTree.h"
#include "DistributedFactory.h"
#include "DistributedQP.hpp"
#include "SparseSymMatrix.h"
#include "SparseGenMatrix.h"

#include "pipsport.h"

#include "omp.h"

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to SparseLinearSystem.
 */
class sLinsysLeaf : public DistributedLinearSystem {
public:
   sLinsysLeaf(DistributedFactory* factory, DistributedQP* prob_, OoqpVector* dd_, OoqpVector* dq_, OoqpVector* nomegaInv_, OoqpVector* primal_reg_,
         OoqpVector* dual_y_reg_, OoqpVector* dual_z_reg_, OoqpVector* rhs_);

   ~sLinsysLeaf() override = default;

   void factor2(DistributedQP* prob, Variables* vars) override;
   void assembleKKT(DistributedQP*, Variables*) override {};
   void allreduceAndFactorKKT(DistributedQP* prob, Variables* vars) override { factor2(prob, vars); };

   void Lsolve(DistributedQP*, OoqpVector&) override {};
   void Dsolve(DistributedQP*, OoqpVector& x) override;
   void Ltsolve(DistributedQP*, OoqpVector&) override {};

   //void Lsolve2 ( OoqpVector& x ) override;
   //void Dsolve2 ( OoqpVector& x ) override;
   void Ltsolve2(DistributedQP* prob, DistributedVector<double>& x, SimpleVector<double>& xp, bool) override;

   void putXDiagonal(const OoqpVector& xdiag_) override;
   void putZDiagonal(const OoqpVector& zdiag_) override;

   void addRegularization(OoqpVector& regP_, OoqpVector& regDy_, OoqpVector& regDz_) const override;
   void addRegularizationsToKKTs(const OoqpVector& regP_, const OoqpVector& regDy_, const OoqpVector& regDz_) override;

   //void Ltsolve_internal(  DistributedQP *prob, DistributedVector<double>& x, SimpleVector<double>& xp);
   void deleteChildren() override;

   void addTermToSchurComplBlocked(DistributedQP* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC, int) override;

   void addLniziLinkCons(DistributedQP* prob, OoqpVector& z0_, OoqpVector& zi_, bool) override;

   void addInnerBorderKiInvBrToRes(DoubleMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool, bool sparse_res, bool sym_res,
         int begin_cols, int end_cols, int) override;
   void
   LniTransMultHierarchyBorder(DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
         bool sparse_res, bool sym_res, bool, int begin_cols, int end_cols, int n_empty_rows_inner_border) override;

protected:

   static void mySymAtPutSubmatrix(SymMatrix& kkt, GenMatrix& B, GenMatrix&, int locnx, int locmy, int);

   void addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0, BorderLinsys& border) override;
   void addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0, BorderLinsys& border) override;
private:
   void addBorderTimesRhsToB0(SimpleVector<double>& rhs, SimpleVector<double>& b0, BorderBiBlock& border);
   void addBorderX0ToRhs(SimpleVector<double>& rhs, const SimpleVector<double>& x0, BorderBiBlock& border);

   /* compute result += B_inner^T K^-1 Br */
   void addInnerBorderKiInvBrToRes(DenseGenMatrix& result, BorderLinsys& Br, int begin_cols, int end_cols);

   /* compute result += B_inner^T K^-1 ( Br - Br_mod_border ) */
   void addLeftBorderKiInvBrToRes(DoubleMatrix& result, BorderBiBlock& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sparse_res,
         bool sym_res, int begin_cols_br, int end_cols_br, int begin_cols_res, int end_cols_res);
};

#endif
