/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYS
#define SAUGLINSYS

#include "DistributedRootLinearSystem.h"
#include <memory>

class DistributedQP;

/**
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAug : public DistributedRootLinearSystem {

public:
   sLinsysRootAug(DistributedFactory* factory_, DistributedQP* prob_, bool is_hierarchy_root = false);
   sLinsysRootAug(DistributedFactory* factory, DistributedQP* prob_, std::shared_ptr<Vector<double>> dd_, std::shared_ptr<Vector<double>> dq_, std::shared_ptr<Vector<double>> nomegaInv_,
         std::shared_ptr<Vector<double>> regP, std::shared_ptr<Vector<double>> regDy, std::shared_ptr<Vector<double>> regDz, std::shared_ptr<Vector<double>> rhs_, bool create_sub_root_solver);
   ~sLinsysRootAug() override = default;

   void finalizeKKT() override;
   void finalizeKKTdist() override;

   void Lsolve(Vector<double>& x) override;
   void Dsolve(Vector<double>& x) override;
   void Ltsolve(Vector<double>& x) override;

   using DistributedLinearSystem::LsolveHierarchyBorder;
   void LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool two_link_border, int begin_cols,
         int end_cols) override;

   using DistributedLinearSystem::LtsolveHierarchyBorder;
   void LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& br_mod_border,
         bool sym_res, bool sparse_res, int begin_cols, int end_cols) override;

   void
   add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization, double dual_inequality_regularization) override;

protected:
   void assembleLocalKKT() override;
   void solveReducedLinkCons(SimpleVector<double>& b);
   void solveReducedLinkConsBlocked(DenseMatrix& rhs_mat_transp, int rhs_start, int n_rhs);
   void addBlTKiInvBrToRes(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sym_res,
         bool sparse_res) override;
   void addBlTKiInvBrToResBlockwise(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sym_res,
         bool sparse_res, DenseMatrix& buffer_b0, int begin_cols, int end_cols);

private:
   void finalizeKKTdense();
   void finalizeKKTsparse();
   void solveWithIterRef(SimpleVector<double>& b);
   void solveWithBiCGStab(SimpleVector<double>& b);

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
   void addLinkConsBlock0Matrix(const SparseMatrix& Ht, int nHtOffsetCols, int nKktOffsetCols, int startCol, int endCol);

   void addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0, BorderLinsys& border) override;

   void addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0, BorderLinsys& border) override;

   /** y = beta*y - alpha* SC * x */
   void SCmult(double beta, SimpleVector<double>& y, double alpha, SimpleVector<double>& x);

   /** stores C^T (Omega + Reg)^-1 C from the KKT matrix */
   std::unique_ptr<SymmetricMatrix> CtDC;
   /** stores Omeag + Reg for computation of CtDC */
   std::unique_ptr<SimpleVector<double>> dual_inequality_non_link_diagonal_regularized;

   std::vector<double> reduced_rhss_blocked;
   std::unique_ptr<SimpleVector<double>> redRhs;
};

#endif

