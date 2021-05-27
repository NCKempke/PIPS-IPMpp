/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include <memory>
#include <utility>

#include "DistributedLeafLinearSystem.h"

DistributedLeafLinearSystem::DistributedLeafLinearSystem(DistributedFactory* factory_, DistributedQP* prob,
   std::shared_ptr<Vector<double>> dd_, std::shared_ptr<Vector<double>> dq_, std::shared_ptr<Vector<double>> nomegaInv_,
   std::shared_ptr<Vector<double>> primal_reg_, std::shared_ptr<Vector<double>> dual_y_reg_,
   std::shared_ptr<Vector<double>> dual_z_reg_, std::shared_ptr<Vector<double>> rhs_) : DistributedLinearSystem(
   factory_,
   prob, std::move(dd_), std::move(dq_), std::move(nomegaInv_), std::move(primal_reg_), std::move(dual_y_reg_), std::move(dual_z_reg_), std::move(rhs_), false) {
#ifdef TIMING
   const int myRank = PIPS_MPIgetRank(mpiComm);
   const double t0 = MPI_Wtime();
#endif

   prob->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);
   const int n = locnx + locmy + locmz;

   int nnzQ, nnzB, nnzD;
   prob->getLocalNnz(nnzQ, nnzB, nnzD);

   if (apply_regularization) {
      regularization_strategy = std::make_unique<RegularizationStrategy>(static_cast<unsigned int>(locnx),
         static_cast<unsigned int>(locmy + locmz));
   }
#ifdef TIMING
   if( myRank == 0 )
      std::cout << "Rank 0: building local Schur matrix ..." << std::endl;
#endif

   /* allocate and copy lower triangular:
    *
    * [ Qq BiT DiT ]
    * [ Bi  0   0  ]
    * [ Di  0   0  ]
    *
    * where  Qq = Q + V^-1 Gamma + W^-1 Phi (so we estimate its nnzs as n_1 + nnzQ)
    */

   auto* kkt_sp = new SparseSymmetricMatrix(n, n + nnzQ + nnzB + nnzD);

   SimpleVector<double> v(n);
   v.setToZero();
   kkt_sp->setToDiagonal(v);

   kkt_sp->symAtPutSubmatrix(0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);

   // TODO this logic or is flawed - requires Bi to exist..
   if (locmz > 0) {
      kkt_sp->symAtPutSubmatrix(locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);
      kkt_sp->symAtPutSubmatrix(locnx + locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
   } else
      kkt_sp->symAtPutSubmatrix(locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);
//      mySymAtPutSubmatrix(*kkt_sp, data->getLocalB(), locnx, locmy, locmz);

#ifdef TIMING
   if( myRank == 0 ) std::cout << "Rank 0: finished " << std::endl;
#endif

   kkt.reset(kkt_sp);
   solver.reset(DistributedFactory::make_leaf_solver(kkt_sp));

#ifdef TIMING
   const double t1 = MPI_Wtime() - t0;
   if (myRank == 0) printf("Rank 0: new sLinsysLeaf took %f sec\n",t1);
#endif

   assert(this->primal_diagonal);
   mpiComm = (dynamic_cast<DistributedVector<double>&>(*this->primal_diagonal)).mpiComm;
}

void DistributedLeafLinearSystem::factor2() {
   // Diagonals were already updated, so
   // just trigger a local refactorization (if needed, depends on the type of lin solver).
   stochNode->resMon.recFactTmLocal_start();

   if (apply_regularization) {
      factorize_with_correct_inertia();
   } else {
      solver->matrixChanged();
   }

   stochNode->resMon.recFactTmLocal_stop();
}

void DistributedLeafLinearSystem::put_primal_diagonal() {
   assert(primal_diagonal);
   const auto& primal_diagonal_stoch = dynamic_cast<const DistributedVector<double>&>(*primal_diagonal);
   kkt->atPutDiagonal(0, *primal_diagonal_stoch.first);
}

void DistributedLeafLinearSystem::clear_dual_equality_diagonal() {
   kkt->diagonal_add_constant_from(locnx, locmy, 0.0);
};

void DistributedLeafLinearSystem::put_dual_inequalites_diagonal() {
   assert(nomegaInv);
   const auto& nomegaInv_stoch = dynamic_cast<const DistributedVector<double>&>(*nomegaInv);
   kkt->atPutDiagonal(locnx + locmy, *nomegaInv_stoch.first);
}

void DistributedLeafLinearSystem::put_barrier_parameter(double barrier) {
   this->barrier_parameter_current_iterate = barrier;
}

void DistributedLeafLinearSystem::add_regularization_diagonal(int offset, double regularization,
   Vector<double>& regularization_vector_) {
   assert(dynamic_cast<DistributedVector<double>&>(regularization_vector_).first);

   auto& regularization_vector = *dynamic_cast<DistributedVector<double>&>(regularization_vector_).first;

   regularization_vector.addConstant(regularization);
   kkt->diagonal_add_constant_from(offset, regularization_vector.length(), regularization);
}

/** adds regularization terms to primal, dualy and dualz vectors - these might depend on the level of linsys we are in */
void
DistributedLeafLinearSystem::add_regularization_local_kkt(double primal_regularization,
   double dual_equality_regularization, double dual_inequality_regularization) {
   assert(this->primal_regularization_diagonal);
   assert(this->dual_equality_regularization_diagonal);
   assert(this->dual_inequality_regularization_diagonal);

   std::cout << "regularizing leaf with " << primal_regularization << " " << dual_equality_regularization << " "
             << dual_inequality_regularization << std::endl;
   if (locnx > 0) {
      add_regularization_diagonal(0, primal_regularization, *this->primal_regularization_diagonal);
   }

   if (locmy > 0) {
      add_regularization_diagonal(locnx, -dual_equality_regularization, *this->dual_equality_regularization_diagonal);
   }

   if (locmz > 0) {
      add_regularization_diagonal(locnx + locmy, -dual_inequality_regularization,
         *this->dual_inequality_regularization_diagonal);
   }
}

void DistributedLeafLinearSystem::Dsolve(Vector<double>& x_in) {
   auto& x = dynamic_cast<DistributedVector<double>&>(x_in);
   assert(x.children.empty());
   stochNode->resMon.recDsolveTmChildren_start();
   solver->Dsolve(*x.first);
   stochNode->resMon.recDsolveTmChildren_stop();
}

void DistributedLeafLinearSystem::Ltsolve2(DistributedVector<double>& x, SimpleVector<double>& xp, bool) {
   auto& b = dynamic_cast<DistributedVector<double>&>(x);
   auto& bi = dynamic_cast<SimpleVector<double>&>(*b.first);
   assert(b.children.empty());

#ifdef TIMING
   stochNode->resMon.eLtsolve.clear();
   stochNode->resMon.recLtsolveTmLocal_start();
#endif

   //b_i -= Lni^T x0
   LniTransMult(bi, -1.0, xp);
   //  solver->Ltsolve(bi); -> empty
#ifdef TIMING
   stochNode->resMon.recLtsolveTmChildren_stop();
#endif
}

void DistributedLeafLinearSystem::deleteChildren() {}

/** sum up right hand side for (current) scenario i and add it to right hand side of scenario 0 */
void DistributedLeafLinearSystem::addLniziLinkCons(Vector<double>& z0_, Vector<double>& zi_, bool /*use_local_RAC*/) {
   auto& z0 = dynamic_cast<SimpleVector<double>&>(z0_);
   auto& zi = dynamic_cast<SimpleVector<double>&>(*dynamic_cast<DistributedVector<double>&>(zi_).first);

   solver->solve(zi);

   int nx0 = data->hasRAC() ? data->getLocalA().n_columns() : 0;

   SimpleVector<double> z01(&z0[0], nx0);
   SimpleVector<double> zi1(&zi[0], locnx);

   if (data->hasRAC()) {
      const SparseMatrix& A = data->getLocalA();
      const SparseMatrix& C = data->getLocalC();
      const SparseMatrix& R = data->getLocalCrossHessian();

      SimpleVector<double> zi2(&zi[locnx], locmy);
      SimpleVector<double> zi3(&zi[locnx + locmy], locmz);

      R.transMult(1.0, z01, -1.0, zi1);
      A.transMult(1.0, z01, -1.0, zi2);
      C.transMult(1.0, z01, -1.0, zi3);
   }

   if (locmyl > 0) {
      assert(locmyl >= 0);
      const int nxMyMz = z0.length() - locmyl - locmzl;

      SimpleVector<double> z0myl(&z0[nxMyMz], locmyl);
      const SparseMatrix& F = data->getLocalF();
      F.mult(1.0, z0myl, -1.0, zi1);
   }

   if (locmzl > 0) {
      assert(locmyl >= 0);
      const int nxMyMzMyl = z0.length() - locmzl;

      SimpleVector<double> z0mzl(&z0[nxMyMzMyl], locmzl);
      const SparseMatrix& G = data->getLocalG();
      G.mult(1.0, z0mzl, -1.0, zi1);
   }
}

void
DistributedLeafLinearSystem::addTermToSchurComplBlocked(bool sparseSC, SymmetricMatrix& SC, bool use_local_RAC, int) {
   const bool sc_is_sym = true;

   std::unique_ptr<BorderBiBlock> border_right{};
   std::unique_ptr<BorderBiBlock> border_left_transp{};

   if (use_local_RAC) {
      const int n_empty =
         SC.size() - data->getLocalCrossHessian().getStorage().n - data->getLocalF().getStorage().m -
            data->getLocalG().getStorage().m;

      assert(data->hasRAC());
      border_right = std::make_unique<BorderBiBlock>(
         data->getLocalCrossHessian(), data->getLocalA(), data->getLocalC(), n_empty, data->getLocalF().getTranspose(),
         data->getLocalG().getTranspose());

      border_left_transp = std::make_unique<BorderBiBlock>(
         data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(),
         data->getLocalC().getTranspose(),
         n_empty, data->getLocalF(), data->getLocalG());

   } else {
      assert(!data->hasRAC());
      const int n_empty = SC.size() - data->getLocalF().getStorage().m - data->getLocalG().getStorage().m;

      border_right = std::make_unique<BorderBiBlock>(n_empty, data->getLocalF().getTranspose(),
         data->getLocalG().getTranspose(), use_local_RAC);
      border_left_transp = std::make_unique<BorderBiBlock>(n_empty, data->getLocalF(), data->getLocalG(),
         use_local_RAC);
   }


   if (border_left_transp->isEmpty() || border_right->isEmpty())
      return;

   addBiTLeftKiBiRightToResBlockedParallelSolvers(sparseSC, sc_is_sym, *border_left_transp, *border_right, SC, 0,
      SC.size(), 0, SC.size());
}

void
DistributedLeafLinearSystem::mySymAtPutSubmatrix(SymmetricMatrix& kkt_, const GeneralMatrix& B_, int locnx, int locmy,
   int) {
   auto& kkt = dynamic_cast<SparseSymmetricMatrix&>(kkt_);
   auto& B = dynamic_cast<const SparseMatrix&>(B_);

   int* jcolK = kkt.jcolM();
   const int* jcolB = B.jcolM(); //int* jcolD = D.jcolM();
   int* krowK = kkt.krowM();
   const int* krowB = B.krowM(); //int* krowD =  D.krowM();
   double* MK = kkt.M();
   const double* MB = B.M();

   for (int i = 0; i < locmy; i++) {
      int itK = krowK[i + locnx];
      int j = krowB[i];

      for (; j < krowB[i + 1]; j++) {

         if (jcolB[j] < i + locnx) {
            jcolK[itK] = jcolB[j];
            MK[itK] = MB[j];
            itK++;
         }
      }
      jcolK[itK] = i + locnx;
      MK[itK] = 0.0;
      itK++;

      assert(j == krowB[i + 1]);

      krowK[i + locnx + 1] = itK;
   }
}

/* compute result += B_inner^T K^-1 Br */
void DistributedLeafLinearSystem::addInnerBorderKiInvBrToRes(DenseMatrix& result, BorderLinsys& Br, int begin_cols,
   int end_cols) {
   assert(Br.A.children.empty());

   /* empty dummy */
   std::vector<BorderMod> Br_mod_border;

   addInnerBorderKiInvBrToRes(result, Br, Br_mod_border, false, false, false, begin_cols, end_cols, 0);
}

/* compute result += [ Bl^T K^-1 ( Br - SUM_j Brmodj Xj ) ]^T = (Br^T - SUM_j Xj^T Brmodj^T) K^-1 Bl for cols begin_cols to end_cols in (Br - SUM_j Brmodj Xj) */
void DistributedLeafLinearSystem::addLeftBorderKiInvBrToRes(AbstractMatrix& result, BorderBiBlock& Bl, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border,
   bool sparse_res, bool sym_res, int begin_cols_br, int end_cols_br, int begin_cols_res, int end_cols_res) {
#ifndef NDEBUG
   const int n_cols = end_cols_br - begin_cols_br;
   assert(end_cols_res - begin_cols_res == n_cols);

   const auto mres = result.n_rows();
   assert(0 <= begin_cols_br && begin_cols_br <= end_cols_br);
   assert(0 <= begin_cols_res && begin_cols_res <= end_cols_res);
   assert(end_cols_res <= mres);
   assert(n_cols <= mres);

   int nx_border{0};
   if (Br.has_RAC)
      nx_border = Br.R.n_columns();
   else if (Br.use_local_RAC)
      nx_border = data->getLocalCrossHessian().n_columns();

   const int myl_border = Br.F.first->n_rows();
   const int mzl_border = Br.G.first->n_rows();

   if (sc_compute_blockwise_hierarchical)
      assert(end_cols_br <= nx_border + Br.n_empty_rows + myl_border + mzl_border);
   else {
      assert(end_cols_br == nx_border + Br.n_empty_rows + myl_border + mzl_border && begin_cols_br == 0);
   }
#endif

   if (Bl.isEmpty())
      return;

   /* buffer for (Bri - (sum_j Brmodj * Xmodj)_i - Bi_{inner} X0)^T = Bri^T - X0^T Bi_{inner}^T - (sum_j Xmodj^T Brmodj^T)_i */
   /* Bi buffer and X0 are in transposed form for memory alignment reasons when solving with K_i */
   const int m_buffer = allocateAndZeroBlockedComputationsBuffer(result.n_rows(), kkt->size());
   (void) m_buffer;
   assert(n_cols <= m_buffer);
   /* put cols from begin_cols to end_cold of (Bri)^T into buffer
    *
    *                [ RiT 0 AiT CiT ]
    * Bri^T        = [  Fi 0  0   0  ]
    *                [  Gi 0  0   0  ]
    */
   if (!Br_mod_border.empty()) {
      std::unique_ptr<BorderBiBlock> BriT{};

      if (Br.has_RAC)
         BriT = std::make_unique<BorderBiBlock>(
            dynamic_cast<SparseMatrix&>(*Br.R.first).getTranspose(),
            dynamic_cast<SparseMatrix&>(*Br.A.first).getTranspose(),
            dynamic_cast<SparseMatrix&>(*Br.C.first).getTranspose(), Br.n_empty_rows,
            dynamic_cast<SparseMatrix&>(*Br.F.first),
            dynamic_cast<SparseMatrix&>(*Br.G.first));
      else if (Br.use_local_RAC)
         BriT = std::make_unique<BorderBiBlock>(data->getLocalCrossHessian().getTranspose(),
            data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
            Br.n_empty_rows, dynamic_cast<SparseMatrix&>(*Br.F.first), dynamic_cast<SparseMatrix&>(*Br.G.first));
      else
         BriT = std::make_unique<BorderBiBlock>(Br.n_empty_rows, dynamic_cast<SparseMatrix&>(*Br.F.first),
            dynamic_cast<SparseMatrix&>(*Br.G.first), false);


      // TODO : return early if all Bordermods and Br were empty
      if (!BriT->isEmpty())
         putBiTBorder(*buffer_blocked_hierarchical, *BriT, begin_cols_br, end_cols_br);

      /* BiT_buffer += X_j^T Bmodj for all j */
      multRightDenseBorderModBlocked(Br_mod_border, *buffer_blocked_hierarchical, begin_cols_br, end_cols_br);

      /* compute B_{inner}^T Ki^-1 Bi_buffer = B_{inner}^T Ki^-1 (Bri^T - sumj Xj^T Bmodj^T) */
      addBiTLeftKiDenseToResBlockedParallelSolvers(sparse_res, sym_res, Bl, *buffer_blocked_hierarchical, result,
         begin_cols_res, end_cols_res);
   } else {
      std::unique_ptr<BorderBiBlock> BriT{};

      if (Br.has_RAC)
         BriT = std::make_unique<BorderBiBlock>(dynamic_cast<SparseMatrix&>(*Br.R.first),
            dynamic_cast<SparseMatrix&>(*Br.A.first),
            dynamic_cast<SparseMatrix&>(*Br.C.first), Br.n_empty_rows,
            dynamic_cast<SparseMatrix&>(*Br.F.first).getTranspose(),
            dynamic_cast<SparseMatrix&>(*Br.G.first).getTranspose());
      else if (Br.use_local_RAC)
         BriT = std::make_unique<BorderBiBlock>(data->getLocalCrossHessian(), data->getLocalA(), data->getLocalC(),
            Br.n_empty_rows,
            dynamic_cast<SparseMatrix&>(*Br.F.first).getTranspose(),
            dynamic_cast<SparseMatrix&>(*Br.G.first).getTranspose());
      else
         BriT = std::make_unique<BorderBiBlock>(Br.n_empty_rows,
            dynamic_cast<SparseMatrix&>(*Br.F.first).getTranspose(),
            dynamic_cast<SparseMatrix&>(*Br.G.first).getTranspose(), false);


      if (!Br.isEmpty())
         addBiTLeftKiBiRightToResBlockedParallelSolvers(sparse_res, sym_res, Bl, *BriT, result, begin_cols_br,
            end_cols_br, begin_cols_res,
            end_cols_res);
   }
}

/* compute result += [ B_{inner}^T K^-1 ( Br - SUM_j Brmodj Xj ) ]^T = (Br^T - SUM_j Xj^T Brmodj^T) K^-1 B_{inner} */
void DistributedLeafLinearSystem::addInnerBorderKiInvBrToRes(AbstractMatrix& result, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool, bool sparse_res,
   bool sym_res, int begin_cols, int end_cols, int) {
   const auto res_n = result.n_columns();
   const int n_empty = data->hasRAC() ? res_n - data->getLocalCrossHessian().getStorage().n -
      data->getLocalF().getStorage().m -
      data->getLocalG().getStorage().m : res_n - data->getLocalF().getStorage().m -
      data->getLocalG().getStorage().m;
   assert(n_empty >= 0);

   BorderBiBlock BiT_inner = data->hasRAC() ? BorderBiBlock(data->getLocalCrossHessian().getTranspose(),
      data->getLocalA().getTranspose(),
      data->getLocalC().getTranspose(), n_empty, data->getLocalF(), data->getLocalG()) : BorderBiBlock(n_empty,
      data->getLocalF(),
      data->getLocalG(), false);

   addLeftBorderKiInvBrToRes(result, BiT_inner, Br, Br_mod_border, sparse_res, sym_res, begin_cols, end_cols, 0,
      end_cols - begin_cols);
}

/* compute res += [Bli^T X_i]^T = [ Bli^T Ki^-1 (Bri - Br_mod_border - Bi_{inner} X0) ]^T = (Bri^T - SUM_i Xi^T Brmodi^T - X0^T Bi_{inner}^T) Ki^{-1} Bli and add it to res
 * begin_cols to end_cols is the position in res */
void
DistributedLeafLinearSystem::LniTransMultHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl,
   BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool sparse_res, bool sym_res, bool, int begin_cols, int end_cols,
   int n_empty_rows_inner_border) {
   std::unique_ptr<BorderBiBlock> BliT{};

   if (Bl.has_RAC)
      BliT = std::make_unique<BorderBiBlock>(dynamic_cast<SparseMatrix&>(*Bl.R.first).getTranspose(),
         dynamic_cast<SparseMatrix&>(*Bl.A.first).getTranspose(),
         dynamic_cast<SparseMatrix&>(*Bl.C.first).getTranspose(), Bl.n_empty_rows,
         dynamic_cast<SparseMatrix&>(*Bl.F.first),
         dynamic_cast<SparseMatrix&>(*Bl.G.first));
   else if (Bl.use_local_RAC)
      BliT = std::make_unique<BorderBiBlock>(data->getLocalCrossHessian().getTranspose(),
         data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
         Bl.n_empty_rows, dynamic_cast<SparseMatrix&>(*Bl.F.first), dynamic_cast<SparseMatrix&>(*Bl.G.first));
   else
      BliT = std::make_unique<BorderBiBlock>(Bl.n_empty_rows, dynamic_cast<SparseMatrix&>(*Bl.F.first),
         dynamic_cast<SparseMatrix&>(*Bl.G.first), false);

   /* constructed to be able to call addLeftBorderKiInvBrToRes.. */
   assert(!data->isHierarchyInnerLeaf());

   // TODO : const cast...
   std::unique_ptr<const StripMatrix> localF_view(
      std::make_unique<const StripMatrix>(true, std::unique_ptr<SparseMatrix>(const_cast<SparseMatrix*>(&data->getLocalF())), nullptr, mpiComm, true));
   std::unique_ptr<const StripMatrix> localG_view(
      std::make_unique<const StripMatrix>(true, std::unique_ptr<SparseMatrix>(const_cast<SparseMatrix*>(&data->getLocalG())), nullptr, mpiComm, true));

   assert(n_empty_rows_inner_border >= 0);

   BorderLinsys B_inner(n_empty_rows_inner_border, *localF_view, *localG_view, data->hasRAC());
   if (!B_inner.isEmpty()) {
      BorderMod inner_mod(B_inner, X0);
      Br_mod_border.push_back(inner_mod);
   }

   /* sym res is a schur complement while the other is a buffer */
   if (sym_res)
      addLeftBorderKiInvBrToRes(res, *BliT, Br, Br_mod_border, sparse_res, sym_res, begin_cols, end_cols, begin_cols,
         end_cols);
   else
      addLeftBorderKiInvBrToRes(res, *BliT, Br, Br_mod_border, sparse_res, sym_res, begin_cols, end_cols, 0,
         end_cols - begin_cols);

}

void DistributedLeafLinearSystem::addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0,
   BorderLinsys& border) {
   assert(border.F.children.empty());
   assert(rhs.children.empty());

   assert(rhs.first);
   if (border.has_RAC) {
      assert(border.R.first);
      assert(border.A.first);
      assert(border.C.first);
   }
   assert(border.F.first);
   assert(border.G.first);

   std::unique_ptr<BorderBiBlock> border_block{};

   if (border.has_RAC)
      border_block = std::make_unique<BorderBiBlock>(dynamic_cast<SparseMatrix&>(*border.R.first),
         dynamic_cast<SparseMatrix&>(*border.A.first),
         dynamic_cast<SparseMatrix&>(*border.C.first), border.n_empty_rows,
         dynamic_cast<SparseMatrix&>(*border.F.first),
         dynamic_cast<SparseMatrix&>(*border.G.first));
   else if (border.use_local_RAC)
      border_block = std::make_unique<BorderBiBlock>(data->getLocalCrossHessian(), data->getLocalA(), data->getLocalC(),
         border.n_empty_rows,
         dynamic_cast<SparseMatrix&>(*border.F.first), dynamic_cast<SparseMatrix&>(*border.G.first));
   else
      border_block = std::make_unique<BorderBiBlock>(
         border.n_empty_rows, dynamic_cast<SparseMatrix&>(*border.F.first),
         dynamic_cast<SparseMatrix&>(*border.G.first),
         false);

   if (border_block->isEmpty())
      return;

   addBorderTimesRhsToB0(dynamic_cast<SimpleVector<double>&>(*rhs.first), b0, *border_block);
}

void DistributedLeafLinearSystem::addBorderTimesRhsToB0(SimpleVector<double>& rhs, SimpleVector<double>& b0,
   BorderBiBlock& border) {
   if (border.isEmpty())
      return;

   const auto[mFi, nFi] = border.F.n_rows_columns();
   const auto mGi = border.G.n_rows();

   int mRi{0};
   int nRi{0};
   int mAi{0};
   int mCi{0};
   if (border.has_RAC) {
      std::tie(mRi, nRi) = border.R.n_rows_columns();
      mAi = border.A.n_rows();
      mCi = border.C.n_rows();

      assert(nFi == mRi);
      assert(rhs.length() == mRi + mAi + mCi);
   } else
      assert(rhs.length() >= nFi);

   assert(b0.length() >= nRi + mFi + mGi);
   const int nb0 = b0.length();

   auto& zi = dynamic_cast<SimpleVector<double>&>(rhs);
   SimpleVector<double> zi1(&zi[0], nFi);

   if (border.has_RAC) {
      SimpleVector<double> zi2(&zi[mRi], mAi);
      SimpleVector<double> zi3(&zi[mRi + mAi], mCi);

      SimpleVector<double> b1(&b0[0], nRi);

      border.R.transMult(1.0, b1, -1.0, zi1);
      border.A.transMult(1.0, b1, -1.0, zi2);
      border.C.transMult(1.0, b1, -1.0, zi3);
   }

   SimpleVector<double> b2(&b0[nb0 - mFi - mGi], mFi);
   SimpleVector<double> b3(&b0[nb0 - mGi], mGi);

   border.F.mult(1.0, b2, -1.0, zi1);
   border.G.mult(1.0, b3, -1.0, zi1);
}

void DistributedLeafLinearSystem::addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0,
   BorderLinsys& border) {
   assert(border.F.children.empty());
   assert(rhs.children.empty());

   assert(rhs.first);
   if (border.has_RAC) {
      assert(border.R.first);
      assert(border.A.first);
      assert(border.C.first);
   }
   assert(border.F.first);
   assert(border.G.first);

   std::unique_ptr<BorderBiBlock> border_block{};

   if (border.has_RAC)
      border_block = std::make_unique<BorderBiBlock>(dynamic_cast<SparseMatrix&>(*border.R.first),
         dynamic_cast<SparseMatrix&>(*border.A.first),
         dynamic_cast<SparseMatrix&>(*border.C.first), border.n_empty_rows,
         dynamic_cast<SparseMatrix&>(*border.F.first),
         dynamic_cast<SparseMatrix&>(*border.G.first));
   else if (border.use_local_RAC)
      border_block = std::make_unique<BorderBiBlock>(data->getLocalCrossHessian(), data->getLocalA(), data->getLocalC(),
         border.n_empty_rows,
         dynamic_cast<SparseMatrix&>(*border.F.first), dynamic_cast<SparseMatrix&>(*border.G.first));
   else
      border_block = std::make_unique<BorderBiBlock>(
         border.n_empty_rows, dynamic_cast<SparseMatrix&>(*border.F.first),
         dynamic_cast<SparseMatrix&>(*border.G.first),
         false);

   if (border_block->isEmpty())
      return;

   addBorderX0ToRhs(dynamic_cast<SimpleVector<double>&>(*rhs.first), x0, *border_block);
}

void DistributedLeafLinearSystem::addBorderX0ToRhs(SimpleVector<double>& rhs, const SimpleVector<double>& x0,
   BorderBiBlock& border) {
   const auto mFi = border.F.n_rows();
   const auto mGi = border.G.n_rows();
   const auto mRi = border.has_RAC ? border.R.n_rows() : 0;
   const auto mAi = border.has_RAC ? border.A.n_rows() : 0;

#ifndef NDEBUG
   const auto mCi = border.has_RAC ? border.C.n_rows() : 0;
   int nRi{0};
   if (border.has_RAC) {
      assert(border.F.n_columns() == mRi);
      assert(rhs.length() == mRi + mAi + mCi);
   } else
      assert(rhs.length() >= border.F.n_columns());

   assert(x0.length() >= nRi + mFi + mGi);
#endif
   if (border.isEmpty())
      return;

   const int nb0 = x0.length();

   double* rhsi1 = &rhs[0];

   if (border.has_RAC) {
      const double* x1 = &x0[0];

      double* rhsi2 = &rhs[mRi];
      double* rhsi3 = &rhs[mRi + mAi];
      border.R.mult(1.0, rhsi1, 1, -1.0, x1, 1);
      border.A.mult(1.0, rhsi2, 1, -1.0, x1, 1);
      border.C.mult(1.0, rhsi3, 1, -1.0, x1, 1);
   }

   const double* x2 = &x0[nb0 - mFi - mGi];
   const double* x3 = &x0[nb0 - mGi];

   border.F.transMult(1.0, rhsi1, 1, -1.0, x2, 1);
   border.G.transMult(1.0, rhsi1, 1, -1.0, x3, 1);
}
