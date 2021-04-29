/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#ifndef SLINSYS
#define SLINSYS

#include "QpGenLinsys.h"
#include "DoubleLinearSolver.h"
#include "OoqpVectorHandle.h"
#include "DenseSymMatrix.h"
#include "SparseSymMatrix.h"
#include "DenseGenMatrix.h"
#include "SimpleVector.h"
#include "StochVector.h"
#include "StringGenMatrix.h"

#include <vector>
#include <memory>
#include <tuple>

#include "mpi.h"

class sTree;

class DistributedFactory;

class DistributedQP;

class StringSymMatrix;

class sLinsys : public QpGenLinsys {
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
      T& R;
      T& A;
      T& C;
      T& F;
      T& G;

      /* n_empty_rows gives the distance between RAC and F,G blocks */
      int n_empty_rows;

      bool isEmpty() const;

      RACFG_BLOCK(T& R, T& A, T& C, int n_empty_rows, T& F, T& G) : has_RAC{true}, R{R}, A{A}, C{C}, F{F}, G{G}, n_empty_rows{n_empty_rows} {
         assert(n_empty_rows >= 0);
      };

      RACFG_BLOCK(int n_empty_rows, T& F, T& G, bool use_local_RAC) : use_local_RAC{use_local_RAC}, has_RAC{false}, R{*dummy}, A{*dummy}, C{*dummy},
            F{F}, G{G}, n_empty_rows{n_empty_rows} { assert(n_empty_rows >= 0); };

      RACFG_BLOCK(const RACFG_BLOCK<T>& block) : use_local_RAC{block.use_local_RAC}, has_RAC{block.has_RAC}, R{block.R}, A{block.A}, C{block.C},
            F{block.F}, G{block.G}, n_empty_rows{block.n_empty_rows} { assert(n_empty_rows >= 0); };
   };

   using BorderLinsys = RACFG_BLOCK<StringGenMatrix>;
   using BorderBiBlock = RACFG_BLOCK<SparseGenMatrix>;

   static BorderLinsys getChild(BorderLinsys& border, unsigned int i) {
      const bool dummy = border.F.children[i]->isKindOf(kStringGenDummyMatrix);
      assert(i < border.F.children.size());
      if (border.has_RAC) {
         if (!dummy && border.F.children[i]->mat->isKindOf(kStringGenMatrix))
            return BorderLinsys(dynamic_cast<StringGenMatrix&>(*border.R.children[i]->mat),
                  dynamic_cast<StringGenMatrix&>(*border.A.children[i]->mat), dynamic_cast<StringGenMatrix&>(*border.C.children[i]->mat),
                  border.n_empty_rows, dynamic_cast<StringGenMatrix&>(*border.F.children[i]->mat),
                  dynamic_cast<StringGenMatrix&>(*border.G.children[i]->mat));
         else
            return BorderLinsys(*border.R.children[i], *border.A.children[i], *border.C.children[i], border.n_empty_rows, *border.F.children[i],
                  *border.G.children[i]);
      }
      else {
         if (!dummy && border.F.children[i]->mat->isKindOf(kStringGenMatrix))
            return BorderLinsys(border.n_empty_rows, dynamic_cast<StringGenMatrix&>(*border.F.children[i]->mat),
                  dynamic_cast<StringGenMatrix&>(*border.G.children[i]->mat), border.use_local_RAC);
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

   using BorderMod = BorderMod_Block<DenseGenMatrix>;
   using BorderModVector = BorderMod_Block<StochVector>;


   template<typename T>
   static BorderMod_Block<T> getChild(BorderMod_Block<T>& bordermod, unsigned int i) {
      BorderLinsys child = getChild(bordermod.border, i);
      return BorderMod_Block<T>(child, bordermod.multiplier);
   }

  sLinsys(DistributedFactory* factory, DistributedQP* prob, bool is_hierarchy_root = false);
  sLinsys(DistributedFactory* factory,
     DistributedQP* prob,
		   OoqpVector* dd,
		   OoqpVector* dq,
		   OoqpVector* nomegaInv,
		   OoqpVector* primal_reg_,
		   OoqpVector* dual_y_reg_,
		   OoqpVector* dual_z_reg_,
		   OoqpVector* rhs,
		   bool create_iter_ref_vecs
  );

   ~sLinsys() override = default;

   void factorize(Problem* prob, Variables* vars) override;

   virtual void factor2(DistributedQP* prob, Variables* vars) = 0;

   virtual void assembleKKT(DistributedQP* prob, Variables* vars) = 0;

   virtual void allreduceAndFactorKKT(DistributedQP* prob, Variables* vars) = 0;

   virtual void Lsolve(DistributedQP* prob, OoqpVector& x) = 0;

   virtual void Dsolve(DistributedQP* prob, OoqpVector& x) = 0;

   virtual void Ltsolve(DistributedQP* prob, OoqpVector& x) = 0;

   virtual void Ltsolve2(DistributedQP* prob, StochVector& x, SimpleVector& xp, bool use_local_RAC) = 0;

   void solveCompressed(OoqpVector& rhs) override;

   void joinRHS(OoqpVector& rhs_in, const OoqpVector& rhs1_in, const OoqpVector& rhs2_in, const OoqpVector& rhs3_in) const override;

   void separateVars(OoqpVector& x_in, OoqpVector& y_in, OoqpVector& z_in, const OoqpVector& vars_in) const override;

   virtual void deleteChildren() = 0;

   virtual bool isDummy() const { return false; };

protected:
   int locnx, locmy, locmyl, locmz, locmzl;
   DistributedQP* data{};

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

   /* symmetric Schur Complement / whole KKT system in lower triangular from */
   std::unique_ptr<SymMatrix> kkt{};
   std::unique_ptr<DoubleLinearSolver> solver{};

   const int blocksize_hierarchical{20};
   const bool sc_compute_blockwise_hierarchical{false};
   std::unique_ptr<DenseGenMatrix> buffer_blocked_hierarchical{};

public:
   MPI_Comm mpiComm{MPI_COMM_NULL};
   sTree* stochNode{};

protected:
   /* depending on SC_HIERARCHICAL_COMPUTE_BLOCKWISE either allocated a full buffer of buffer_m rows or a smaller one - returns number of rows in buffer */
   int allocateAndZeroBlockedComputationsBuffer(int buffer_m, int buffer_n);

public:
   virtual void addLnizi(DistributedQP* prob, OoqpVector& z0, OoqpVector& zi);

   virtual void addLniziLinkCons(DistributedQP*/*prob*/, OoqpVector& /*z0*/, OoqpVector& /*zi*/, bool /*use_local_RAC*/) {
      assert(false && "not implemented here");
   };


   /* put BiT into res */
   virtual void putBiTBorder(DenseGenMatrix& res, const BorderBiBlock& BiT, int begin_rows, int end_rows) const;

   /* compute Bli^T X_i = Bli^T Ki^-1 (Bri - Bi_{inner} X0) and add it to SC */
   virtual void LniTransMultHierarchyBorder(DoubleMatrix& /*SC*/, const DenseGenMatrix& /*X0*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/,
         std::vector<BorderMod>& /*Br_mod_border*/, bool /*sparse_res*/, bool /*sym_res*/, bool /*use_local_RAC*/, int /*begin_cols*/,
         int /*end_cols*/, int /*n_empty_rows_inner_border*/) { assert(false && "not implemented here"); };

   /** y += alpha * Lni^T * x */
   virtual void LniTransMult(DistributedQP* prob, SimpleVector& y, double alpha, SimpleVector& x);

   /** Methods that use dense matrices U and V to compute the
    *  terms from the Schur complement.
    */
   virtual void allocU(DenseGenMatrix** Ut, int np);

   virtual void allocV(DenseGenMatrix** V, int np);

   virtual void computeU_V(DistributedQP* prob, DenseGenMatrix* U, DenseGenMatrix* V);

   /** Method(s) that use a memory-friendly mechanism for computing
    *  the terms from the Schur Complement
    */
   virtual void addTermToDenseSchurCompl(DistributedQP* prob, DenseSymMatrix& SC);

   virtual void addTermToSchurComplBlocked(DistributedQP* /*prob*/, bool /*sparseSC*/, SymMatrix& /*SC*/, bool /*use_local_RAC*/,
         int /*n_empty_rows_inner_border*/) { assert(0 && "not implemented here"); };

   virtual void computeInnerSystemRightHandSide(StochVector& /*rhs_inner*/, const SimpleVector& /*b0*/, bool /*use_local_RAC*/) {
      assert(false && "not implemented here");
   };

public:

   /* add you part of the border times rhs to b0 */
   virtual void addBorderTimesRhsToB0(StochVector& /*rhs*/, SimpleVector& /*b0*/, BorderLinsys& /*border*/ ) {
      assert(false && "not implemented here");
   };

   /* add you part of the border times rhs to b0 */
   virtual void addBorderX0ToRhs(StochVector& /*rhs*/, const SimpleVector& /*x0*/, BorderLinsys& /*border*/ ) {
      assert(false && "not implemented here");
   };

   virtual void addBiTLeftKiBiRightToResBlockedParallelSolvers(bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
         /* const */ BorderBiBlock& border_right, DoubleMatrix& result, int begin_cols, int end_cols, int begin_rows_res, int end_rows_res);

   void addBiTLeftKiDenseToResBlockedParallelSolvers(bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
         /* const */ DenseGenMatrix& BT, DoubleMatrix& result, int begin_rows_res, int end_rows_res);

   virtual void addTermToSparseSchurCompl(DistributedQP* /*prob*/, SparseSymMatrix& /*SC*/ ) { assert(0 && "not implemented here"); };

   /** Used in the iterative refinement for the dense Schur complement systems
    * Computes res += [0 A^T C^T ]*inv(KKT)*[0;A;C] x
    */
   virtual void addTermToSchurResidual(DistributedQP* prob, SimpleVector& res, SimpleVector& x);

   // TODO only compute bottom left part for symmetric matrices
   /* compute result += Bl^T K^-1 Br where K is our own linear system */
   virtual void addBlTKiInvBrToRes(DoubleMatrix& /*result*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/,
         bool /*sym_res*/, bool /*sparse_res*/) { assert(false && "not implemented here"); }

   /* compute Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ) and add it up in result */
   virtual void
   LsolveHierarchyBorder(DenseGenMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*two_link_border*/,
         int /*begin_cols*/, int /*end_cols*/) { assert(false && "not implemented here"); };

   virtual void
   LsolveHierarchyBorder(DenseGenMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*use_local_RAC*/,
         bool /*two_link_border*/, int /*begin_cols*/, int /*end_cols*/) { assert(false && "not implemented here"); };

   /* solve with SC and comput X_0 = SC^-1 B_0 */
   virtual void DsolveHierarchyBorder(DenseGenMatrix& /*buffer_b0*/, int /*n_cols*/) { assert(false && "not implemented here"); };

   /* compute RES += SUM_i Bli_^T X_i = Bli^T Ki^-1 ( ( Bri - sum_j Bmodij Xij ) - Bi_{inner} X0) */
   virtual void LtsolveHierarchyBorder(DoubleMatrix& /*res*/, const DenseGenMatrix& /*X0*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/,
         std::vector<BorderMod>& /*Br_mod_border*/, bool /*sym_res*/, bool /*sparse_res*/, int /*begin_cols*/, int /*end_cols*/) {
      assert(false && "not implemented here");
   };

   virtual void LtsolveHierarchyBorder(DoubleMatrix& /*res*/, const DenseGenMatrix& /*X0*/, BorderLinsys& /*Bl*/, BorderLinsys& /*Br*/,
         std::vector<BorderMod>& /*Br_mod_border*/, bool /*sym_res*/, bool /*sparse_res*/, bool /*use_local_RAC*/, int /*begin_cols*/,
         int /*end_cols*/) { assert(false && "not implemented here"); };

   /* compute Bi_{inner}^T Ki^{-1} ( Bri - sum_j Brmod_ij Xj )and add it to result */
   virtual void
   addInnerBorderKiInvBrToRes(DoubleMatrix& /*result*/, BorderLinsys& /*Br*/, std::vector<BorderMod>& /*Br_mod_border*/, bool /*has_RAC*/,
         bool /*sparse_res*/, bool /*sym_res*/, int /*begin_cols*/, int /*end_cols*/, int /*n_empty_rows_inner_border*/) {
      assert(false && "not implemented here");
   };

protected:
   void addLeftBorderTimesDenseColsToResTransp(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col, int n_cols,
         bool sparse_res, bool sym_res, DoubleMatrix& res) const;

   void addLeftBorderTimesDenseColsToResTranspSparse(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col, int n_cols,
         SparseSymMatrix& res) const;

   void addLeftBorderTimesDenseColsToResTranspDense(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col, int n_cols,
         int n_cols_res, double** res) const;

   /* calculate res -= BT0 * X0 */
   void finalizeDenseBorderBlocked(BorderLinsys& B, const DenseGenMatrix& X, DenseGenMatrix& result, int begin_rows, int end_rows);

   /* calculate res -= X0 * BT */
   void multRightDenseBorderBlocked(BorderBiBlock& BT, const DenseGenMatrix& X, DenseGenMatrix& result, int begin_rows, int end_rows);

   /* calculate res -= (sum_j X0jT * BjT ) */
   void multRightDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseGenMatrix& result, int begin_cols, int end_cols);

   /* calculate res -= (sum_j X0jT * B0JT ) */
   void finalizeDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseGenMatrix& result, int begin_rows, int end_rows);

};

template<>
bool sLinsys::RACFG_BLOCK<StringGenMatrix>::isEmpty() const;

template<>
bool sLinsys::RACFG_BLOCK<SparseGenMatrix>::isEmpty() const;


#endif
