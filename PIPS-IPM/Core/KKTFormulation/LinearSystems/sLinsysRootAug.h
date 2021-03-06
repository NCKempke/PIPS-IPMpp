/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#ifndef SAUGLINSYS
#define SAUGLINSYS

#include "DistributedRootLinearSystem.h"
#include <memory>

class DistributedProblem;

/**
 * ROOT (= NON-leaf) linear system in reduced augmented form
 */
class sLinsysRootAug : public DistributedRootLinearSystem {

public:
   sLinsysRootAug(const DistributedFactory& factory_, DistributedProblem* prob_, bool is_hierarchy_root = false);
   sLinsysRootAug(const DistributedFactory& factory, DistributedProblem* prob_, std::shared_ptr<Vector<double>> dd_, std::shared_ptr<Vector<double>> dq_, std::shared_ptr<Vector<double>> nomegaInv_,
         std::shared_ptr<Vector<double>> regP, std::shared_ptr<Vector<double>> regDy, std::shared_ptr<Vector<double>> regDz, std::shared_ptr<Vector<double>> rhs_, bool creat_solvers);
   ~sLinsysRootAug() override = default;

   void finalizeKKT() override;
   void finalizeKKTdist() override;

   void Lsolve(Vector<double>& x) override;
   void Dsolve(Vector<double>& x) override;
   void Ltsolve(Vector<double>& x) override;

   using DistributedLinearSystem::LsolveHierarchyBorder;
   void LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, int begin_cols,
         int end_cols) override;

   using DistributedLinearSystem::LtsolveHierarchyBorder;
   void LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& br_mod_border,
         bool sym_res, bool sparse_res, int begin_cols, int end_cols) override;

   void
   add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization, double dual_inequality_regularization) override;
   void put_dual_inequalites_diagonal() override;

protected:
   void assembleLocalKKT() override;
   void solveReducedLinkCons(DenseVector<double>& b);
   void solveReducedLinkConsBlocked(DenseMatrix& rhs_mat_transp, int rhs_start, int n_rhs);
   void addBlTKiInvBrToRes(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sym_res,
         bool sparse_res) override;
   void addBlTKiInvBrToResBlockwise(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool sym_res,
         bool sparse_res, DenseMatrix& buffer_b0, int begin_cols, int end_cols);

private:
   void finalizeKKTdense();
   void finalizeKKTsparse();
   void solveWithIterRef(DenseVector<double>& b);
   void solveWithBiCGStab(DenseVector<double>& b);

   void DsolveHierarchyBorder(DenseMatrix& b, int n_cols) override;

   void add_CtDC_to_sparse_schur_complement(const SymmetricMatrix& CtDC);
   void add_CtDC_to_dense_schur_complement(const SymmetricMatrix& CtDC);

   /** computes CtDC and stores it in CtDC_loc + adds it to the Schur complement. If CtDC == nullptr it will also allocate CtDC */
   void compute_CtDC_and_add_to_Schur_complement(SymmetricMatrix*& CtDC_loc, const Vector<double>& diagonal);

   void clear_CtDC_from_sparse_schur_complement(const SymmetricMatrix& CtDC);
   void clear_CtDC_from_dense_schur_complement(const SymmetricMatrix& CtDC);

   /** computes CtDC and stores it  in CtDC_loc + adds it to the Schur complement. If CtDC == nullptr it will also allocate CtDC */
   void clear_CtDC_from_schur_complement(const SymmetricMatrix& CtDC_loc);

   // add specified columns of given matrix Ht (either Ft or Gt) to Schur complement
   void addLinkConsBlock0Matrix(const SparseMatrix& Ht, int nHtOffsetCols, int nKktOffsetCols, int startCol, int endCol);

   void addBorderX0ToRhs(DistributedVector<double>& rhs, const DenseVector<double>& x0, BorderLinsys& border) override;

   void addBorderTimesRhsToB0(DistributedVector<double>& rhs, DenseVector<double>& b0, BorderLinsys& border) override;

   /** y = beta*y - alpha* SC * x */
   void SCmult(double beta, DenseVector<double>& y, double alpha, DenseVector<double>& x);

   void schur_complement_put_primal_block(const Vector<double>* diagonal, SymmetricMatrix& schur_complement) const;
   void schur_complement_add_CTDC_block(const Vector<double>* diagonal, SymmetricMatrix& schur_complement);
   void schur_complement_put_AT_block(SymmetricMatrix& schur_complement) const;

   // TODO : unused ..
   const bool eliminate_C_block_from_kkt{true};
   /** stores C^T (Omega + Reg)^-1 C from the KKT matrix */
   std::unique_ptr<SymmetricMatrix> CtDC;

   /** since we eliminate the C0 block we need a buffer for (D + reg)^-1 */
   std::unique_ptr<Vector<double>> dual_inequality_diagonal_regularized{};

   std::vector<double> reduced_rhss_blocked;
   std::unique_ptr<DenseVector<double>> redRhs;


};

#endif

