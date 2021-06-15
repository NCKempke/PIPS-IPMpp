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
#include "SimpleVector.h"
#include "DistributedVector.h"
#include "StripMatrix.h"

#include <vector>
#include <memory>
#include <tuple>

#include "mpi.h"

class DistributedTree;

class DistributedFactory;

class DistributedQP;

class DistributedLinearSystem : public LinearSystem {
public:
   template<typename T>
   struct RACFG_BLOCK {
   private:
      std::unique_ptr<T> dummy{new T()};

   public:
      const bool use_local_RAC{};
      const bool has_RAC{};
      /* represents a block like
       * [ R_i 0 F_i^T G_i^T ]             [ R_i^T A_i^T C_i^T ]
       * [ A_i 0   0     0   ] or possibly [   0     0     0   ]
       * [ C_i 0   0     0   ]             [  F_i    0     0   ]
       *                                   [  G_i    0     0   ]
       */
      const T& R;
      const T& A;
      const T& C;
      const T& F;
      const T& G;

      /* n_empty_rows gives the distance between RAC and F,G blocks */
      int n_empty_rows;

      [[nodiscard]] bool isEmpty() const;

      RACFG_BLOCK(const T& R, const T& A, const T& C, int n_empty_rows, const T& F, const T& G) : has_RAC{true}, R{R}, A{A}, C{C}, F{F}, G{G}, n_empty_rows{n_empty_rows} {
         assert(n_empty_rows >= 0);
      };

      RACFG_BLOCK(int n_empty_rows, const T& F, const T& G, bool use_local_RAC) : use_local_RAC{use_local_RAC}, has_RAC{false}, R{*dummy}, A{*dummy}, C{*dummy},
            F{F}, G{G}, n_empty_rows{n_empty_rows} { assert(n_empty_rows >= 0); };

      RACFG_BLOCK(const RACFG_BLOCK<T>& block) : use_local_RAC{block.use_local_RAC}, has_RAC{block.has_RAC}, R{block.R}, A{block.A}, C{block.C},
            F{block.F}, G{block.G}, n_empty_rows{block.n_empty_rows} { assert(n_empty_rows >= 0); };
   };

   using BorderLinsys = RACFG_BLOCK<StripMatrix>;
   using BorderBiBlock = RACFG_BLOCK<SparseMatrix>;

   static BorderLinsys getChild(BorderLinsys& border, unsigned int i) {
      const bool dummy = border.F.children[i]->is_a(kStringGenDummyMatrix);
      assert(i < border.F.children.size());
      if (border.has_RAC) {
         if (!dummy && border.F.children[i]->first->is_a(kStripMatrix))
            return BorderLinsys(dynamic_cast<StripMatrix&>(*border.R.children[i]->first),
                  dynamic_cast<StripMatrix&>(*border.A.children[i]->first), dynamic_cast<StripMatrix&>(*border.C.children[i]->first),
                  border.n_empty_rows, dynamic_cast<StripMatrix&>(*border.F.children[i]->first),
                  dynamic_cast<StripMatrix&>(*border.G.children[i]->first));
         else
            return BorderLinsys(*border.R.children[i], *border.A.children[i], *border.C.children[i], border.n_empty_rows, *border.F.children[i],
                  *border.G.children[i]);
      }
      else {
         if (!dummy && border.F.children[i]->first->is_a(kStripMatrix))
            return BorderLinsys(border.n_empty_rows, dynamic_cast<StripMatrix&>(*border.F.children[i]->first),
                  dynamic_cast<StripMatrix&>(*border.G.children[i]->first), border.use_local_RAC);
         else
            return BorderLinsys(border.n_empty_rows, *border.F.children[i], *border.G.children[i], border.use_local_RAC);
      }
   }

   template<typename T>
   struct BorderMod_Block {
   public:
      BorderLinsys border;
      const T& multiplier;

      BorderMod_Block(BorderLinsys& border_, const T& multiplier) : border{border_}, multiplier{multiplier} {};
   };

   using BorderMod = BorderMod_Block<DenseMatrix>;


   template<typename T>
   static BorderMod_Block<T> getChild(BorderMod_Block<T>& bordermod, unsigned int i) {
      BorderLinsys child = getChild(bordermod.border, i);
      return BorderMod_Block<T>(child, bordermod.multiplier);
   }

   DistributedLinearSystem(DistributedFactory* factory, DistributedQP* prob, bool is_hierarchy_root = false);

   DistributedLinearSystem(DistributedFactory* factory, DistributedQP* prob, std::shared_ptr<Vector<double>> dd, std::shared_ptr<Vector<double>> dq, std::shared_ptr<Vector<double>> nomegaInv,
         std::shared_ptr<Vector<double>> primal_reg_, std::shared_ptr<Vector<double>> dual_y_reg_, std::shared_ptr<Vector<double>> dual_z_reg_, std::shared_ptr<Vector<double>> rhs, bool create_iter_ref_vecs);

   ~DistributedLinearSystem() override = default;

   void factorize(Variables& variables) override;

   virtual void factor2() = 0;

   virtual void assembleKKT() = 0;

   virtual void allreduceAndFactorKKT() = 0;

   virtual void Lsolve(Vector<double>& x) = 0;

   virtual void Dsolve(Vector<double>& x) = 0;

   virtual void Ltsolve(Vector<double>& x) = 0;

   virtual void Ltsolve2(DistributedVector<double>& x, SimpleVector<double>& xp, bool use_local_RAC) = 0;

   void solveCompressed(Vector<double>& rhs) override;

   virtual void deleteChildren() = 0;

   [[nodiscard]] virtual bool isDummy() const { return false; };

protected:
   int locnx{};
   int locmy{};
   int locmyl{};
   int locmz{};
   int locmzl{};

   const DistributedQP* data{};

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
   DistributedTree* stochNode{};

protected:
   /* depending on SC_HIERARCHICAL_COMPUTE_BLOCKWISE either allocated a full buffer of buffer_m rows or a smaller one - returns number of rows in buffer */
   int allocateAndZeroBlockedComputationsBuffer(int buffer_m, int buffer_n);

public:
   virtual void addLnizi(Vector<double>& z0, Vector<double>& zi);

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

public:
   /* add you part of the border times rhs to b0 */
   virtual void addBorderTimesRhsToB0(DistributedVector<double>& /*rhs*/, SimpleVector<double>& /*b0*/, BorderLinsys& /*border*/ ) {
      assert(false && "not implemented here");
   };

   /* add you part of the border times rhs to b0 */
   virtual void addBorderX0ToRhs(DistributedVector<double>& /*rhs*/, const SimpleVector<double>& /*x0*/, BorderLinsys& /*border*/ ) {
      assert(false && "not implemented here");
   };

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
   LsolveHierarchyBorder(DenseMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*two_link_border*/,
         int /*begin_cols*/, int /*end_cols*/) { assert(false && "not implemented here"); };

   virtual void
   LsolveHierarchyBorder(DenseMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*use_local_RAC*/,
         bool /*two_link_border*/, int /*begin_cols*/, int /*end_cols*/) { assert(false && "not implemented here"); };

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

template<>
bool DistributedLinearSystem::RACFG_BLOCK<StripMatrix>::isEmpty() const;

template<>
bool DistributedLinearSystem::RACFG_BLOCK<SparseMatrix>::isEmpty() const;


#endif
