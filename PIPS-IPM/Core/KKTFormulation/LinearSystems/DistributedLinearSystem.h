/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SLINSYS
#define SLINSYS

#include "LinearSystem.h"
#include "DoubleLinearSolver.h"
#include "Vector.hpp"
#include "DenseSymmetricMatrix.h"
#include "SparseSymmetricMatrix.h"
#include "DenseMatrix.h"
#include "SimpleVector.hpp"
#include "DistributedVector.h"
#include "StripMatrix.h"
#include "PIPSIPMppOptions.h"

#include "RACFG_BLOCK.h"
#include "BorderMod_Block.h"

#include <vector>
#include <memory>
#include <tuple>

#include "mpi.h"

class DistributedTree;

class DistributedFactory;

class DistributedProblem;

class StochNodeResourcesMonitor;

class DistributedLinearSystem : public LinearSystem {

public:
   DistributedLinearSystem(const DistributedFactory& factory, DistributedProblem* prob, bool is_hierarchy_root = false);

   DistributedLinearSystem(const DistributedFactory& factory, DistributedProblem* prob, std::shared_ptr<Vector<double>> dd, std::shared_ptr<Vector<double>> dq, std::shared_ptr<Vector<double>> nomegaInv,
         std::shared_ptr<Vector<double>> primal_reg_, std::shared_ptr<Vector<double>> dual_y_reg_, std::shared_ptr<Vector<double>> dual_z_reg_, std::shared_ptr<Vector<double>> rhs, bool create_iter_ref_vecs);

   ~DistributedLinearSystem() override = default;

   void factorize(const Variables& variables) override;

   virtual void factor2() = 0;

   virtual void assembleKKT() = 0;

   virtual void allreduceAndFactorKKT() = 0;

   virtual void Lsolve(Vector<double>& x) = 0;

   virtual void Dsolve(Vector<double>& x) = 0;

   virtual void Ltsolve(Vector<double>& x) = 0;

   virtual void Ltsolve2(DistributedVector<double>& x, SimpleVector<double>& xp, bool use_local_RAC) = 0;

   void solveCompressed(Vector<double>& rhs) override;

   [[nodiscard]] virtual bool isDummy() const { return false; };

protected:
   int locnx{};
   int locmy{};
   int locmyl{};
   int locmz{};
   int locmzl{};

   const DistributedProblem* data{};

   int iAmDistrib;

   /* members for blockwise schur complement computation */
   bool computeBlockwiseSC{false};

   /* batches right hand sides are being processed int */
   int blocksizemax{0};
   std::vector<double> colsBlockDense;
   std::vector<int> colId;
   std::vector<int> colSparsity;

   /* is this linsys the overall root */
   const bool is_hierarchy_root{false};

   const int blocksize_hierarchical{20};
   const bool sc_compute_blockwise_hierarchical{false};
   std::unique_ptr<DenseMatrix> buffer_blocked_hierarchical{};

public:
   MPI_Comm mpiComm{MPI_COMM_NULL};
   const DistributedTree* distributed_tree{};
   StochNodeResourcesMonitor* resource_monitor{};
protected:
   /* depending on SC_HIERARCHICAL_COMPUTE_BLOCKWISE either allocated a full buffer of buffer_m rows or a smaller one - returns number of rows in buffer */
   int allocateAndZeroBlockedComputationsBuffer(int buffer_m, int buffer_n);

public:
   virtual void addLniziLinkCons(Vector<double>& /*z0*/, Vector<double>& /*zi*/, bool /*use_local_RAC*/) {
      assert(false && "not implemented here");
   };

   /* put BiT into res */
   virtual void putBiTBorder(DenseMatrix& res, const BorderBiBlock& BiT, int begin_rows, int end_rows) const;

   /* compute Bli^T X_i = Bli^T Ki^-1 (Bri - Bi_{inner} X0) and add it to SC */
   virtual void LniTransMultHierarchyBorder(AbstractMatrix& /*SC*/, const DenseMatrix& /*X0*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/,
         std::vector<BorderMod>& /*Br_mod_border*/, bool /*sparse_res*/, bool /*sym_res*/, bool /*use_local_RAC*/, int /*begin_cols*/,
         int /*end_cols*/, int /*n_empty_rows_inner_border*/) { assert(false && "not implemented here"); };

   /** y += alpha * Lni^T * x */
   virtual void LniTransMult(SimpleVector<double>& y, double alpha, SimpleVector<double>& x);

   /** Method(s) that use a memory-friendly mechanism for computing
    *  the terms from the Schur Complement
    */
   virtual void addTermToDenseSchurCompl(DenseSymmetricMatrix& SC);

   virtual void addTermToSchurComplBlocked(bool /*sparseSC*/, SymmetricMatrix& /*SC*/, bool /*use_local_RAC*/,
         int /*n_empty_rows_inner_border*/) { assert(0 && "not implemented here"); };

   virtual void
   computeInnerSystemRightHandSide(DistributedVector<double>& /*rhs_inner*/, const SimpleVector<double>& /*b0*/, bool /*use_local_RAC*/) {
      assert(false && "not implemented here");
   };

   /* add you part of the border times rhs to b0 */
   virtual void addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0, BorderLinsys& border) = 0;

   /* add you part of the border times rhs to b0 */
   virtual void addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0, BorderLinsys& border)  = 0;

   virtual void addBiTLeftKiBiRightToResBlockedParallelSolvers(bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
         /* const */ BorderBiBlock& border_right, AbstractMatrix& result, int begin_cols, int end_cols, int begin_rows_res, int end_rows_res);

   void addBiTLeftKiDenseToResBlockedParallelSolvers(bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
         /* const */ DenseMatrix& BT, AbstractMatrix& result, int begin_rows_res, int end_rows_res);

   virtual void addTermToSparseSchurCompl(SparseSymmetricMatrix& /*SC*/ ) { assert(0 && "not implemented here"); };

   /** Used in the iterative refinement for the dense Schur complement systems
    * Computes res += [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    */
   virtual void addTermToSchurResidual(SimpleVector<double>& res, SimpleVector<double>& x);

   // TODO only compute bottom left part for symmetric matrices
   /* compute result += Bl^T K^-1 Br where K is our own linear system */
   virtual void addBlTKiInvBrToRes(AbstractMatrix& /*result*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/,
         bool /*sym_res*/, bool /*sparse_res*/) { assert(false && "not implemented here"); }

   /* compute Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ) and add it up in result */
   virtual void
   LsolveHierarchyBorder(DenseMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/,
      int /*begin_cols*/, int /*end_cols*/) { assert(false && "not implemented here"); };

   virtual void
   LsolveHierarchyBorder(DenseMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*use_local_RAC*/,
         int /*begin_cols*/, int /*end_cols*/) { assert(false && "not implemented here"); };

   /* solve with SC and comput X_0 = SC^-1 B_0 */
   virtual void DsolveHierarchyBorder(DenseMatrix& /*buffer_b0*/, int /*n_cols*/) { assert(false && "not implemented here"); };

   /* compute RES += SUM_i Bli_^T X_i = Bli^T Ki^-1 ( ( Bri - sum_j Bmodij Xij ) - Bi_{inner} X0) */
   virtual void LtsolveHierarchyBorder(AbstractMatrix& /*res*/, const DenseMatrix& /*X0*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/,
         std::vector<BorderMod>& /*Br_mod_border*/, bool /*sym_res*/, bool /*sparse_res*/, int /*begin_cols*/, int /*end_cols*/) {
      assert(false && "not implemented here");
   };

   virtual void LtsolveHierarchyBorder(AbstractMatrix& /*res*/, const DenseMatrix& /*X0*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/,
         std::vector<BorderMod>& /*Br_mod_border*/, bool /*sym_res*/, bool /*sparse_res*/, bool /*use_local_RAC*/, int /*begin_cols*/,
         int /*end_cols*/) { assert(false && "not implemented here"); };

   /* compute Bi_{inner}^T Ki^{-1} ( Bri - sum_j Brmod_ij Xj )and add it to result */
   virtual void
   addInnerBorderKiInvBrToRes(AbstractMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*has_RAC*/,
         bool /*sparse_res*/, bool /*sym_res*/, int /*begin_cols*/, int /*end_cols*/, int /*n_empty_rows_inner_border*/) {
      assert(false && "not implemented here");
   };

protected:
   void addLeftBorderTimesDenseColsToResTransp(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col, int n_cols,
         bool sparse_res, bool sym_res, AbstractMatrix& res) const;

   static void addLeftBorderTimesDenseColsToResTranspSparse(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col, int n_cols,
         SparseSymmetricMatrix& res) ;

   void addLeftBorderTimesDenseColsToResTranspDense(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col, int n_cols,
         int n_cols_res, double** res) const;

   /* calculate res -= BT0 * X0 */
   void finalizeDenseBorderBlocked(BorderLinsys& B, const DenseMatrix& X, DenseMatrix& result, int begin_rows, int end_rows);

   /* calculate res -= X0 * BT */
   static void multRightDenseBorderBlocked(BorderBiBlock& BT, const DenseMatrix& X, DenseMatrix& result, int begin_rows, int end_rows);

   /* calculate res -= (sum_j X0jT * BjT ) */
   void multRightDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseMatrix& result, int begin_cols, int end_cols);

   /* calculate res -= (sum_j X0jT * B0JT ) */
   void finalizeDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseMatrix& result, int begin_rows, int end_rows);

};

#endif
