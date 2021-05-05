/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef DISTRIBUTEDLEAFLINEARSYSTEM_H
#define DISTRIBUTEDLEAFLINEARSYSTEM_H

#include "DistributedLinearSystem.h"
#include "DistributedTree.h"
#include "DistributedFactory.h"
#include "DistributedQP.hpp"
#include "SparseSymmetricMatrix.h"
#include "SparseMatrix.h"

#include "pipsport.h"

#include "omp.h"

/** This class solves the linear system corresponding to a leaf node.
 *  It just redirects the call to SparseLinearSystem.
 */
class DistributedLeafLinearSystem : public DistributedLinearSystem {
public:
   DistributedLeafLinearSystem(DistributedFactory* factory, DistributedQP* prob_, Vector<double>* dd_, Vector<double>* dq_, Vector<double>* nomegaInv_,
         Vector<double>* primal_reg_, Vector<double>* dual_y_reg_, Vector<double>* dual_z_reg_, Vector<double>* rhs_);

   ~DistributedLeafLinearSystem() override = default;

   void factor2(DistributedQP* prob, Variables* vars) override;
   void assembleKKT(DistributedQP*, Variables*) override {};
   void allreduceAndFactorKKT(DistributedQP* prob, Variables* vars) override { factor2(prob, vars); };

   void Lsolve(DistributedQP*, Vector<double>&) override {};
   void Dsolve(DistributedQP*, Vector<double>& x) override;
   void Ltsolve(DistributedQP*, Vector<double>&) override {};

   //void Lsolve2 ( Vector<double>& x ) override;
   //void Dsolve2 ( Vector<double>& x ) override;
   void Ltsolve2(DistributedQP* prob, DistributedVector<double>& x, SimpleVector<double>& xp, bool) override;

   void put_primal_diagonal() override;
   void put_dual_inequalites_diagonal() override;
   void put_barrier_parameter(double barrier) override;
   void clear_dual_equality_diagonal() override;

   void
   add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization, double dual_inequality_regularization) override;

   //void Ltsolve_internal(  DistributedQP *prob, DistributedVector<double>& x, SimpleVector<double>& xp);
   void deleteChildren() override;

   void addTermToSchurComplBlocked(DistributedQP* prob, bool sparseSC, SymmetricMatrix& SC, bool use_local_RAC, int) override;

   void addLniziLinkCons(DistributedQP* prob, Vector<double>& z0_, Vector<double>& zi_, bool) override;

   void addInnerBorderKiInvBrToRes(AbstractMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool, bool sparse_res, bool sym_res,
         int begin_cols, int end_cols, int) override;
   void
   LniTransMultHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
         bool sparse_res, bool sym_res, bool, int begin_cols, int end_cols, int n_empty_rows_inner_border) override;

protected:
   void add_regularization_diagonal(int offset, double regularization, Vector<double>& regularization_vector);

   static void mySymAtPutSubmatrix(SymmetricMatrix& kkt, GeneralMatrix& B, GeneralMatrix&, int locnx, int locmy, int);

   void addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0, BorderLinsys& border) override;
   void addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0, BorderLinsys& border) override;
private:
   void addBorderTimesRhsToB0(SimpleVector<double>& rhs, SimpleVector<double>& b0, BorderBiBlock& border);
   void addBorderX0ToRhs(SimpleVector<double>& rhs, const SimpleVector<double>& x0, BorderBiBlock& border);

   /* compute result += B_inner^T K^-1 Br */
   void addInnerBorderKiInvBrToRes(DenseMatrix& result, BorderLinsys& Br, int begin_cols, int end_cols);

   /* compute result += B_inner^T K^-1 ( Br - Br_mod_border ) */
   void addLeftBorderKiInvBrToRes(AbstractMatrix& result, BorderBiBlock& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sparse_res,
         bool sym_res, int begin_cols_br, int end_cols_br, int begin_cols_res, int end_cols_res);
};

#endif
