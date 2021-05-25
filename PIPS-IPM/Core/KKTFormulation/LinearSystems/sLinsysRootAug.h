/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYS
#define SAUGLINSYS

#ifdef WITH_PARDISO

#include "PardisoProjectIndefSolver.h"

#endif

#ifdef WITH_MKL_PARDISO
#include "PardisoMKLIndefSolver.h"
#endif

#ifdef WITH_MA57

#include "Ma57SolverRoot.h"

#endif

#ifdef WITH_MA27

#include "Ma27SolverRoot.h"

#endif

#ifdef WITH_MUMPS
#include "MumpsSolverRoot.h"
#endif

#include "DistributedRootLinearSystem.h"
#include <memory>

class DistributedQP;

/**
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAug : public DistributedRootLinearSystem {

public:
   sLinsysRootAug(DistributedFactory* factory_, DistributedQP* prob_);
   sLinsysRootAug(DistributedFactory* factory, DistributedQP* prob_, Vector<double>* dd_, Vector<double>* dq_, Vector<double>* nomegaInv_,
         Vector<double>* regP, Vector<double>* regDy, Vector<double>* regDz, Vector<double>* rhs_, bool creat_solvers);
   ~sLinsysRootAug() override = default;

   void finalizeKKT(DistributedQP* prob, Variables* vars) override;
   void finalizeKKTdist(DistributedQP* prob) override;

   void Lsolve(DistributedQP* prob, Vector<double>& x) override;
   void Dsolve(DistributedQP* prob, Vector<double>& x) override;
   void Ltsolve(DistributedQP* prob, Vector<double>& x) override;

   using DistributedLinearSystem::LsolveHierarchyBorder;
   void LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool two_link_border, int begin_cols,
         int end_cols) override;

   using DistributedLinearSystem::LtsolveHierarchyBorder;
   void LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& br_mod_border,
         bool sym_res, bool sparse_res, int begin_cols, int end_cols) override;

   void
   add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization, double dual_inequality_regularization) override;

protected:
   SymmetricMatrix* createKKT(DistributedQP* prob) const;
   void createSolversSparse(SolverType solver);
   void createSolversDense();

   void assembleLocalKKT(DistributedQP* prob) override;
   void solveReducedLinkCons(DistributedQP* prob, SimpleVector<double>& b);
   void solveReducedLinkConsBlocked(DistributedQP* data, DenseMatrix& rhs_mat_transp, int rhs_start, int n_rhs);
   void addBlTKiInvBrToRes(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sym_res,
         bool sparse_res) override;
   void addBlTKiInvBrToResBlockwise(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sym_res,
         bool sparse_res, DenseMatrix& buffer_b0, int begin_cols, int end_cols);

private:
   void createSolversAndKKts(DistributedQP* prob);
   void finalizeKKTdense(DistributedQP* prob, Variables* vars);
   void finalizeKKTsparse(DistributedQP* prob, Variables* vars);
   void solveWithIterRef(DistributedQP* prob, SimpleVector<double>& b);
   void solveWithBiCGStab(DistributedQP* prob, SimpleVector<double>& b);

   void DsolveHierarchyBorder(DenseMatrix& b, int n_cols) override;

   void add_CtDC_to_sparse_schur_complement(const SymmetricMatrix& CtDC);
   void add_CtDC_to_dense_schur_complement(const SymmetricMatrix& CtDC);

   /** computes CtDC and stores it  in CtDC_loc + adds it to the Schur complement. If CtDC == nullptr it will also allocate CtDC */
   void compute_CtDC_and_add_to_Schur_complement(SymmetricMatrix*& CtDC_loc, const Vector<double>& diagonal);

   void clear_CtDC_from_sparse_schur_complement(const SymmetricMatrix& CtDC);
   void clear_CtDC_from_dense_schur_complement(const SymmetricMatrix& CtDC);

   /** computes CtDC and stores it  in CtDC_loc + adds it to the Schur complement. If CtDC == nullptr it will also allocate CtDC */
   void clear_CtDC_from_schur_complement(const SymmetricMatrix& CtDC_loc);

   // add specified columns of given matrix Ht (either Ft or Gt) to Schur complement
   void addLinkConsBlock0Matrix(DistributedQP* prob, SparseMatrix& Ht, int nHtOffsetCols, int nKktOffsetCols, int startCol, int endCol);

   /** y = beta*y - alpha* SC * x */
   void SCmult(double beta, SimpleVector<double>& y, double alpha, SimpleVector<double>& x, DistributedQP* prob);

   /** stores C^T (Omega + Reg)^-1 C from the KKT matrix */
   std::unique_ptr<SymmetricMatrix> CtDC;
   /** stores Omeag + Reg for computation of CtDC */
   std::unique_ptr<SimpleVector<double>> dual_inequality_non_link_diagonal_regularized;

   std::vector<double> reduced_rhss_blocked;
   std::unique_ptr<SimpleVector<double>> redRhs;
};

#endif

