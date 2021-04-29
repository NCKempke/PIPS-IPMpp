/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysLeaf.h"

sLinsysLeaf::sLinsysLeaf(DistributedFactory *factory_, DistributedQP* prob,
          OoqpVector* dd_,
          OoqpVector* dq_,
          OoqpVector* nomegaInv_,
          OoqpVector* primal_reg_,
          OoqpVector* dual_y_reg_,
          OoqpVector* dual_z_reg_,
          OoqpVector* rhs_
       )
  : DistributedLinearSystem(factory_, prob, dd_, dq_, nomegaInv_, primal_reg_, dual_y_reg_, dual_z_reg_, rhs_, false)
{
#ifdef TIMING
   const int myRank = PIPS_MPIgetRank(mpiComm);
   const double t0 = MPI_Wtime();
#endif

   prob->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);
   const int n = locnx + locmy + locmz;

   int nnzQ, nnzB, nnzD;
   prob->getLocalNnz(nnzQ, nnzB, nnzD);

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

   SparseSymMatrix* kkt_sp = new SparseSymMatrix(n, n + nnzQ + nnzB + nnzD);

   SimpleVector<double> v(n);
   v.setToZero();
   kkt_sp->setToDiagonal(v);

   kkt_sp->symAtPutSubmatrix(0, 0, prob->getLocalQ(), 0, 0, locnx, locnx);

   if (locmz > 0) {
      kkt_sp->symAtPutSubmatrix(locnx, 0, prob->getLocalB(), 0, 0, locmy, locnx);
      kkt_sp->symAtPutSubmatrix(locnx + locmy, 0, prob->getLocalD(), 0, 0, locmz, locnx);
   }
   else
      mySymAtPutSubmatrix(*kkt_sp, prob->getLocalB(), prob->getLocalD(), locnx, locmy, locmz);

#ifdef TIMING
   if( myRank == 0 ) std::cout << "Rank 0: finished " << std::endl;
#endif

   kkt.reset(kkt_sp);
   solver.reset(factory_->make_leaf_solver(kkt_sp));

#ifdef TIMING
   const double t1 = MPI_Wtime() - t0;
   if (myRank == 0) printf("Rank 0: new sLinsysLeaf took %f sec\n",t1);
#endif

   mpiComm = (dynamic_cast<StochVector*>(dd_))->mpiComm;
}

void sLinsysLeaf::factor2(DistributedQP*, Variables*) {
   // Diagonals were already updated, so
   // just trigger a local refactorization (if needed, depends on the type of lin solver).
   stochNode->resMon.recFactTmLocal_start();

   solver->matrixChanged();

   stochNode->resMon.recFactTmLocal_stop();
}

void sLinsysLeaf::putXDiagonal( const OoqpVector& xdiag_ )
{
  const StochVector& xdiag = dynamic_cast<const StochVector&>(xdiag_);
  kkt->atPutDiagonal( 0, *xdiag.first );
}

void sLinsysLeaf::putZDiagonal( const OoqpVector& zdiag_ )
{
  const StochVector& zdiag = dynamic_cast<const StochVector&>(zdiag_);
  kkt->atPutDiagonal( locnx + locmy, *zdiag.first );
}

/** adds regularization terms to primal, dualy and dualz vectors - these might depend on the level of linsys we are in */
void sLinsysLeaf::addRegularization( OoqpVector& regP_, OoqpVector& regDy_, OoqpVector& regDz_ ) const
{
   if( locnx > 0 )
      dynamic_cast<StochVector&>(regP_).first->setToConstant( primal_reg_val );

   if( locmy > 0 )
      dynamic_cast<StochVector&>(regDy_).first->setToConstant( dual_y_reg_val );

   if( locmz > 0 )
      dynamic_cast<StochVector&>(regDz_).first->setToConstant( dual_z_reg_val );
}

void sLinsysLeaf::addRegularizationsToKKTs( const OoqpVector& regP_, const OoqpVector& regDy_, const OoqpVector& regDz_ )
{
  const StochVector& regP = dynamic_cast<const StochVector&>(regP_);
  kkt->atAddDiagonal( 0, *regP.first );

  const StochVector& regDy = dynamic_cast<const StochVector&>(regDy_);
  kkt->atPutDiagonal( locnx, *regDy.first );

  const StochVector& regDz = dynamic_cast<const StochVector&>(regDz_);
  kkt->atAddDiagonal( locnx + locmy, *regDz.first );
}

void sLinsysLeaf::Dsolve( DistributedQP*, OoqpVector& x_in )
{
   StochVector& x = dynamic_cast<StochVector&>(x_in);
   assert(x.children.size() == 0);
   stochNode->resMon.recDsolveTmChildren_start();
   solver->Dsolve(*x.first);
   stochNode->resMon.recDsolveTmChildren_stop();
}

void sLinsysLeaf::Ltsolve2(DistributedQP* prob, StochVector& x, SimpleVector<double>& xp, bool) {
   StochVector& b = dynamic_cast<StochVector&>(x);
   SimpleVector<double>& bi = dynamic_cast<SimpleVector<double>&>(*b.first);
   assert(0 == b.children.size());

#ifdef TIMING
   stochNode->resMon.eLtsolve.clear();
   stochNode->resMon.recLtsolveTmLocal_start();
#endif

   //b_i -= Lni^T x0
   LniTransMult(prob, bi, -1.0, xp);
   //  solver->Ltsolve(bi); -> empty
#ifdef TIMING
   stochNode->resMon.recLtsolveTmChildren_stop();
#endif
}

void sLinsysLeaf::deleteChildren() {}

/** sum up right hand side for (current) scenario i and add it to right hand side of scenario 0 */
void sLinsysLeaf::addLniziLinkCons(DistributedQP* prob, OoqpVector& z0_, OoqpVector& zi_, bool /*use_local_RAC*/) {
   SimpleVector<double>& z0 = dynamic_cast<SimpleVector<double>&>(z0_);
   SimpleVector<double>& zi = dynamic_cast<SimpleVector<double>&>(*dynamic_cast<StochVector&>(zi_).first);

   solver->solve(zi);

   int dummy{0};
   int nx0{0};

   SparseGenMatrix& A = prob->getLocalA();
   if (data->hasRAC())
      A.getSize(dummy, nx0);

   SimpleVector<double> z01(&z0[0], nx0);
   SimpleVector<double> zi1(&zi[0], locnx);

   if (data->hasRAC()) {
      SparseGenMatrix& C = prob->getLocalC();
      SparseGenMatrix& R = prob->getLocalCrossHessian();

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
      SparseGenMatrix& F = prob->getLocalF();
      F.mult(1.0, z0myl, -1.0, zi1);
   }

   if (locmzl > 0) {
      assert(locmyl >= 0);
      const int nxMyMzMyl = z0.length() - locmzl;

      SimpleVector<double> z0mzl(&z0[nxMyMzMyl], locmzl);
      SparseGenMatrix& G = prob->getLocalG();
      G.mult(1.0, z0mzl, -1.0, zi1);
   }
}

void sLinsysLeaf::addTermToSchurComplBlocked(DistributedQP* prob, bool sparseSC, SymMatrix& SC, bool use_local_RAC, int) {
   assert(prob == data);

   const bool sc_is_sym = true;

   std::unique_ptr<BorderBiBlock> border_right{};
   std::unique_ptr<BorderBiBlock> border_left_transp{};

   if (use_local_RAC) {
      const int n_empty =
            SC.size() - prob->getLocalCrossHessian().getStorageRef().n - prob->getLocalF().getStorageRef().m - prob->getLocalG().getStorageRef().m;

      assert(prob->hasRAC());
      border_right.reset(
            new BorderBiBlock(prob->getLocalCrossHessian(), prob->getLocalA(), prob->getLocalC(), n_empty, prob->getLocalF().getTranspose(),
                  prob->getLocalG().getTranspose()));

      border_left_transp.reset(
            new BorderBiBlock(prob->getLocalCrossHessian().getTranspose(), prob->getLocalA().getTranspose(), prob->getLocalC().getTranspose(),
                  n_empty, prob->getLocalF(), prob->getLocalG()));

   }
   else {
      assert(!prob->hasRAC());
      const int n_empty = SC.size() - prob->getLocalF().getStorageRef().m - prob->getLocalG().getStorageRef().m;

      border_right.reset(new BorderBiBlock(n_empty, prob->getLocalF().getTranspose(), prob->getLocalG().getTranspose(), use_local_RAC));
      border_left_transp.reset(new BorderBiBlock(n_empty, prob->getLocalF(), prob->getLocalG(), use_local_RAC));
   }


   if (border_left_transp->isEmpty() || border_right->isEmpty())
      return;

   addBiTLeftKiBiRightToResBlockedParallelSolvers(sparseSC, sc_is_sym, *border_left_transp, *border_right, SC, 0, SC.size(), 0, SC.size());
}

void sLinsysLeaf::mySymAtPutSubmatrix(SymMatrix& kkt_, GenMatrix& B_, GenMatrix&, int locnx, int locmy, int) {
   SparseSymMatrix& kkt = dynamic_cast<SparseSymMatrix&>(kkt_);
   SparseGenMatrix& B = dynamic_cast<SparseGenMatrix&>(B_);
   //SparseGenMatrix& D   = reinterpret_cast<SparseGenMatrix&>(D_);

   int* jcolK = kkt.jcolM();
   int* jcolB = B.jcolM(); //int* jcolD = D.jcolM();
   int* krowK = kkt.krowM();
   int* krowB = B.krowM(); //int* krowD =  D.krowM();
   double* MK = kkt.M();
   double* MB = B.M();

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
void sLinsysLeaf::addInnerBorderKiInvBrToRes(DenseGenMatrix& result, BorderLinsys& Br, int begin_cols, int end_cols) {
   assert(Br.A.children.size() == 0);

   /* empty dummy */
   std::vector<BorderMod> Br_mod_border;

   addInnerBorderKiInvBrToRes(result, Br, Br_mod_border, false, false, false, begin_cols, end_cols, 0);
}

/* compute result += [ Bl^T K^-1 ( Br - SUM_j Brmodj Xj ) ]^T = (Br^T - SUM_j Xj^T Brmodj^T) K^-1 Bl for cols begin_cols to end_cols in (Br - SUM_j Brmodj Xj) */
void sLinsysLeaf::addLeftBorderKiInvBrToRes(DoubleMatrix& result, BorderBiBlock& Bl, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
      bool sparse_res, bool sym_res, int begin_cols_br, int end_cols_br, int begin_cols_res, int end_cols_res) {
   int dummy;
#ifndef NDEBUG
   const int n_cols = end_cols_br - begin_cols_br;
   assert(end_cols_res - begin_cols_res == n_cols);

   int mres, nres;
   result.getSize(mres, nres);
   assert(0 <= begin_cols_br && begin_cols_br <= end_cols_br);
   assert(0 <= begin_cols_res && begin_cols_res <= end_cols_res);
   assert(end_cols_res <= mres);
   assert(n_cols <= mres);

   int nx_border{0};
   if (Br.has_RAC)
      Br.R.getSize(dummy, nx_border);
   else if (Br.use_local_RAC)
      data->getLocalCrossHessian().getSize(dummy, nx_border);

   int myl_border, mzl_border;
   Br.F.mat->getSize(myl_border, dummy);
   Br.G.mat->getSize(mzl_border, dummy);

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
   const int n_buffer = kkt->size();
   int m_result;
   result.getSize(m_result, dummy);

#ifndef NDEBUG
   const int m_buffer = allocateAndZeroBlockedComputationsBuffer(m_result, n_buffer);
   assert(n_cols <= m_buffer);
#else
   allocateAndZeroBlockedComputationsBuffer(m_result, n_buffer);
#endif
   /* put cols from begin_cols to end_cold of (Bri)^T into buffer
    *
    *                [ RiT 0 AiT CiT ]
    * Bri^T        = [  Fi 0  0   0  ]
    *                [  Gi 0  0   0  ]
    */
   if (!Br_mod_border.empty()) {
      std::unique_ptr<BorderBiBlock> BriT{};

      if (Br.has_RAC)
         BriT.reset(
               new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*Br.R.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*Br.A.mat).getTranspose(),
                     dynamic_cast<SparseGenMatrix&>(*Br.C.mat).getTranspose(), Br.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Br.F.mat),
                     dynamic_cast<SparseGenMatrix&>(*Br.G.mat)));
      else if (Br.use_local_RAC)
         BriT.reset(new BorderBiBlock(data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
               Br.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat)));
      else
         BriT.reset(new BorderBiBlock(Br.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Br.F.mat), dynamic_cast<SparseGenMatrix&>(*Br.G.mat), false));


      // TODO : return early if all Bordermods and Br were empty
      if (!BriT->isEmpty())
         putBiTBorder(*buffer_blocked_hierarchical, *BriT, begin_cols_br, end_cols_br);

      /* BiT_buffer += X_j^T Bmodj for all j */
      multRightDenseBorderModBlocked(Br_mod_border, *buffer_blocked_hierarchical, begin_cols_br, end_cols_br);

      /* compute B_{inner}^T Ki^-1 Bi_buffer = B_{inner}^T Ki^-1 (Bri^T - sumj Xj^T Bmodj^T) */
      addBiTLeftKiDenseToResBlockedParallelSolvers(sparse_res, sym_res, Bl, *buffer_blocked_hierarchical, result, begin_cols_res, end_cols_res);
   }
   else {
      std::unique_ptr<BorderBiBlock> BriT{};

      if (Br.has_RAC)
         BriT.reset(new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*Br.R.mat), dynamic_cast<SparseGenMatrix&>(*Br.A.mat),
               dynamic_cast<SparseGenMatrix&>(*Br.C.mat), Br.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Br.F.mat).getTranspose(),
               dynamic_cast<SparseGenMatrix&>(*Br.G.mat).getTranspose()));
      else if (Br.use_local_RAC)
         BriT.reset(new BorderBiBlock(data->getLocalCrossHessian(), data->getLocalA(), data->getLocalC(), Br.n_empty_rows,
               dynamic_cast<SparseGenMatrix&>(*Br.F.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*Br.G.mat).getTranspose()));
      else
         BriT.reset(new BorderBiBlock(Br.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Br.F.mat).getTranspose(),
               dynamic_cast<SparseGenMatrix&>(*Br.G.mat).getTranspose(), false));


      if (!Br.isEmpty())
         addBiTLeftKiBiRightToResBlockedParallelSolvers(sparse_res, sym_res, Bl, *BriT, result, begin_cols_br, end_cols_br, begin_cols_res,
               end_cols_res);
   }
}

/* compute result += [ B_{inner}^T K^-1 ( Br - SUM_j Brmodj Xj ) ]^T = (Br^T - SUM_j Xj^T Brmodj^T) K^-1 B_{inner} */
void sLinsysLeaf::addInnerBorderKiInvBrToRes(DoubleMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border, bool, bool sparse_res,
      bool sym_res, int begin_cols, int end_cols, int) {
   int res_m, res_n;
   result.getSize(res_m, res_n);
   const int n_empty = data->hasRAC() ? res_n - data->getLocalCrossHessian().getStorageRef().n - data->getLocalF().getStorageRef().m -
                                        data->getLocalG().getStorageRef().m : res_n - data->getLocalF().getStorageRef().m -
                                                                              data->getLocalG().getStorageRef().m;
   assert(n_empty >= 0);

   BorderBiBlock BiT_inner = data->hasRAC() ? BorderBiBlock(data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(),
         data->getLocalC().getTranspose(), n_empty, data->getLocalF(), data->getLocalG()) : BorderBiBlock(n_empty, data->getLocalF(),
         data->getLocalG(), false);

   addLeftBorderKiInvBrToRes(result, BiT_inner, Br, Br_mod_border, sparse_res, sym_res, begin_cols, end_cols, 0, end_cols - begin_cols);
}

/* compute res += [Bli^T X_i]^T = [ Bli^T Ki^-1 (Bri - Br_mod_border - Bi_{inner} X0) ]^T = (Bri^T - SUM_i Xi^T Brmodi^T - X0^T Bi_{inner}^T) Ki^{-1} Bli and add it to res
 * begin_cols to end_cols is the position in res */
void sLinsysLeaf::LniTransMultHierarchyBorder(DoubleMatrix& res, const DenseGenMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
      std::vector<BorderMod>& Br_mod_border, bool sparse_res, bool sym_res, bool, int begin_cols, int end_cols, int n_empty_rows_inner_border) {
   std::unique_ptr<BorderBiBlock> BliT{};

   if (Bl.has_RAC)
      BliT.reset(new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*Bl.R.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*Bl.A.mat).getTranspose(),
            dynamic_cast<SparseGenMatrix&>(*Bl.C.mat).getTranspose(), Bl.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Bl.F.mat),
            dynamic_cast<SparseGenMatrix&>(*Bl.G.mat)));
   else if (Bl.use_local_RAC)
      BliT.reset(new BorderBiBlock(data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
            Bl.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Bl.F.mat), dynamic_cast<SparseGenMatrix&>(*Bl.G.mat)));
   else
      BliT.reset(new BorderBiBlock(Bl.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*Bl.F.mat), dynamic_cast<SparseGenMatrix&>(*Bl.G.mat), false));

   /* constructed to be able to call addLeftBorderKiInvBrToRes.. */
   std::unique_ptr<StringGenMatrix> localF_view(std::make_unique<StringGenMatrix>(true, &data->getLocalF(), nullptr, mpiComm, true));
   std::unique_ptr<StringGenMatrix> localG_view(std::make_unique<StringGenMatrix>(true, &data->getLocalG(), nullptr, mpiComm, true));

   assert(n_empty_rows_inner_border >= 0);

   BorderLinsys B_inner(n_empty_rows_inner_border, *localF_view, *localG_view, data->hasRAC());
   if (!B_inner.isEmpty()) {
      BorderMod inner_mod(B_inner, X0);
      Br_mod_border.push_back(inner_mod);
   }

   /* sym res is a schur complement while the other is a buffer */
   if (sym_res)
      addLeftBorderKiInvBrToRes(res, *BliT, Br, Br_mod_border, sparse_res, sym_res, begin_cols, end_cols, begin_cols, end_cols);
   else
      addLeftBorderKiInvBrToRes(res, *BliT, Br, Br_mod_border, sparse_res, sym_res, begin_cols, end_cols, 0, end_cols - begin_cols);

}

void sLinsysLeaf::addBorderTimesRhsToB0(StochVector& rhs, SimpleVector<double>& b0, BorderLinsys& border) {
   assert(border.F.children.size() == 0);
   assert(rhs.children.size() == 0);

   assert(rhs.first);
   if (border.has_RAC) {
      assert(border.R.mat);
      assert(border.A.mat);
      assert(border.C.mat);
   }
   assert(border.F.mat);
   assert(border.G.mat);

   std::unique_ptr<BorderBiBlock> border_block{};

   if (border.has_RAC)
      border_block.reset(new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*border.R.mat), dynamic_cast<SparseGenMatrix&>(*border.A.mat),
            dynamic_cast<SparseGenMatrix&>(*border.C.mat), border.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*border.F.mat),
            dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
   else if (border.use_local_RAC)
      border_block.reset(new BorderBiBlock(data->getLocalCrossHessian(), data->getLocalA(), data->getLocalC(), border.n_empty_rows,
            dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
   else
      border_block.reset(
            new BorderBiBlock(border.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat),
                  false));

   if (border_block->isEmpty())
      return;

   addBorderTimesRhsToB0(dynamic_cast<SimpleVector<double>&>(*rhs.first), b0, *border_block);
}

void sLinsysLeaf::addBorderTimesRhsToB0(SimpleVector<double>& rhs, SimpleVector<double>& b0, BorderBiBlock& border) {
   if (border.isEmpty())
      return;

   int mFi, nFi;
   border.F.getSize(mFi, nFi);
   int mGi, nGi;
   border.G.getSize(mGi, nGi);

   int mRi{0};
   int nRi{0};
   int mAi{0};
   int nAi{0};
   int mCi{0};
   int nCi{0};
   if (border.has_RAC) {
      border.R.getSize(mRi, nRi);
      border.A.getSize(mAi, nAi);
      border.C.getSize(mCi, nCi);

      assert(nFi == mRi);
      assert(rhs.length() == mRi + mAi + mCi);
   }
   else
      assert(rhs.length() >= nFi);

   assert(b0.length() >= nRi + mFi + mGi);
   const int nb0 = b0.length();

   SimpleVector<double>& zi = dynamic_cast<SimpleVector<double>&>(rhs);
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

void sLinsysLeaf::addBorderX0ToRhs(StochVector& rhs, const SimpleVector<double>& x0, BorderLinsys& border) {
   assert(border.F.children.size() == 0);
   assert(rhs.children.size() == 0);

   assert(rhs.first);
   if (border.has_RAC) {
      assert(border.R.mat);
      assert(border.A.mat);
      assert(border.C.mat);
   }
   assert(border.F.mat);
   assert(border.G.mat);

   std::unique_ptr<BorderBiBlock> border_block{};

   if (border.has_RAC)
      border_block.reset(new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*border.R.mat), dynamic_cast<SparseGenMatrix&>(*border.A.mat),
            dynamic_cast<SparseGenMatrix&>(*border.C.mat), border.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*border.F.mat),
            dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
   else if (border.use_local_RAC)
      border_block.reset(new BorderBiBlock(data->getLocalCrossHessian(), data->getLocalA(), data->getLocalC(), border.n_empty_rows,
            dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
   else
      border_block.reset(
            new BorderBiBlock(border.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat),
                  false));

   if (border_block->isEmpty())
      return;

   addBorderX0ToRhs(dynamic_cast<SimpleVector<double>&>(*rhs.first), x0, *border_block);
}

void sLinsysLeaf::addBorderX0ToRhs(SimpleVector<double>& rhs, const SimpleVector<double>& x0, BorderBiBlock& border) {
   int mFi, nFi;
   border.F.getSize(mFi, nFi);
   int mGi, nGi;
   border.G.getSize(mGi, nGi);
   int mRi{0};
   int nRi{0};
   int mAi{0};
   int nAi{0};
   int mCi{0};
   int nCi{0};

   if (border.has_RAC) {
      border.R.getSize(mRi, nRi);
      border.A.getSize(mAi, nAi);
      border.C.getSize(mCi, nCi);

      assert(nFi == mRi);
      assert(rhs.length() == mRi + mAi + mCi);
   }
   else
      assert(rhs.length() >= nFi);

   assert(x0.length() >= nRi + mFi + mGi);

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
