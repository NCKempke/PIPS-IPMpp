/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#include <memory>
#include <utility>

#include "PIPSIPMppOptions.h"
#include "BorderedSymmetricMatrix.h"
#include "DistributedLinearSystem.h"
#include "DistributedFactory.hpp"
#include "DistributedProblem.hpp"

DistributedLinearSystem::DistributedLinearSystem(const DistributedFactory& factory_, DistributedProblem* problem,
   bool is_hierarchy_root) : LinearSystem(
   factory_, *problem), data{problem},
   computeBlockwiseSC(pipsipmpp_options::get_bool_parameter("SC_COMPUTE_BLOCKWISE")),
   blocksizemax(pipsipmpp_options::get_int_parameter("SC_BLOCKWISE_BLOCKSIZE_MAX")),
   is_hierarchy_root(is_hierarchy_root),
   blocksize_hierarchical(pipsipmpp_options::get_int_parameter("SC_BLOCKSIZE_HIERARCHICAL")),
   sc_compute_blockwise_hierarchical{pipsipmpp_options::get_bool_parameter("SC_HIERARCHICAL_COMPUTE_BLOCKWISE")},
   stochNode{factory_.tree.get()} {
   if (sc_compute_blockwise_hierarchical && PIPS_MPIgetRank() == 0)
      std::cout << "Computing hierarchical Schur complements blockwise with buffersize " << blocksize_hierarchical
                << " (times # of available OMP threads)\n";

   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      assert(is_hierarchy_root);

   problem->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);

   //get the communicator from one of the vectors
   const auto& dds = dynamic_cast<const DistributedVector<double>&>(*primal_diagonal);
   this->mpiComm = dds.mpiComm;
   this->iAmDistrib = dds.iAmDistrib;
}

DistributedLinearSystem::DistributedLinearSystem(const DistributedFactory& factory_, DistributedProblem* problem,
   std::shared_ptr<Vector<double>> dd_, std::shared_ptr<Vector<double>> dq_,
   std::shared_ptr<Vector<double>> nomegaInv_, std::shared_ptr<Vector<double>> primal_reg_,
   std::shared_ptr<Vector<double>> dual_y_reg_, std::shared_ptr<Vector<double>> dual_z_reg_,
   std::shared_ptr<Vector<double>> rhs_,
   bool create_iter_ref_vecs) : LinearSystem(factory_, *problem, std::move(dd_), std::move(dq_), std::move(nomegaInv_),
   std::move(primal_reg_), std::move(dual_y_reg_), std::move(dual_z_reg_), std::move(rhs_), create_iter_ref_vecs), data{problem},
   computeBlockwiseSC(pipsipmpp_options::get_bool_parameter("SC_COMPUTE_BLOCKWISE")),
   blocksizemax(pipsipmpp_options::get_int_parameter("SC_BLOCKWISE_BLOCKSIZE_MAX")),
   blocksize_hierarchical(pipsipmpp_options::get_int_parameter("SC_BLOCKSIZE_HIERARCHICAL")),
   sc_compute_blockwise_hierarchical{pipsipmpp_options::get_bool_parameter("SC_HIERARCHICAL_COMPUTE_BLOCKWISE")},
   stochNode{factory_.tree.get()} {
   problem->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);

   if (primal_diagonal) {
      const auto& primal_diagonal_stoch = dynamic_cast<const DistributedVector<double>&>(*primal_diagonal);
      mpiComm = primal_diagonal_stoch.mpiComm;
      iAmDistrib = primal_diagonal_stoch.iAmDistrib;
   } else {
      mpiComm = MPI_COMM_NULL;
      iAmDistrib = false;
   }
}

void DistributedLinearSystem::factorize(const Variables& vars) {
#ifdef TIMING
   double tTot = MPI_Wtime();
#endif
   // the call to the the parent's method takes care of all necessary updates
   // to the KKT system (updating diagonals mainly). This is done recursively,
   // we don't have to worry about it anymore.
   LinearSystem::factorize(vars);

   // now DO THE LINEAR ALGEBRA!
   // in order to avoid a call to QpGenLinsys::factor, call factor2 method.
   factor2();

#ifdef TIMING
   tTot = MPI_Wtime() - tTot;
   MPI_Barrier(MPI_COMM_WORLD);
   const int rank == PIPS_MPIgetRANK(MPI_COMM_WORLD);
   //if( 128 * ( myRank / 128 ) == 0 )
   if( 0 == myRank )
       std::cout << "Outer fact. total time " << tTot << std::endl;
#endif
}

void
DistributedLinearSystem::finalizeDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseMatrix& result,
   int begin_rows, int end_rows) {
   /* compute BiT_buffer += X_j^T Bmodj for all j */
   for (auto& border_mod_block : border_mod) {
      if (border_mod_block.border.isEmpty())
         continue;
      finalizeDenseBorderBlocked(border_mod_block.border, border_mod_block.multiplier, result, begin_rows, end_rows);
   }
}

void
DistributedLinearSystem::multRightDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseMatrix& result,
   int begin_cols, int end_cols) {
   /* compute BiT_buffer += X_j^T Bmodj for all j */
   for (auto& border_mod_block : border_mod) {
      std::unique_ptr<BorderBiBlock> BiT_mod{};

      BorderLinsys& border = border_mod_block.border;

      if (border.isEmpty())
         continue;

      if (border.use_local_RAC)
         BiT_mod = std::make_unique<BorderBiBlock>(
            data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(),
            data->getLocalC().getTranspose(),
            border.n_empty_rows, dynamic_cast<SparseMatrix&>(*border.F.first),
            dynamic_cast<SparseMatrix&>(*border.G.first));
      else if (border.has_RAC)
         BiT_mod = std::make_unique<BorderBiBlock>(dynamic_cast<SparseMatrix&>(*border.R.first).getTranspose(),
            dynamic_cast<SparseMatrix&>(*border.A.first).getTranspose(),
            dynamic_cast<SparseMatrix&>(*border.C.first).getTranspose(),
            border.n_empty_rows, dynamic_cast<SparseMatrix&>(*border.F.first),
            dynamic_cast<SparseMatrix&>(*border.G.first));
      else
         BiT_mod = std::make_unique<BorderBiBlock>(
            border.n_empty_rows, dynamic_cast<SparseMatrix&>(*border.F.first),
            dynamic_cast<SparseMatrix&>(*border.G.first),
            false);

      multRightDenseBorderBlocked(*BiT_mod, border_mod_block.multiplier, result, begin_cols, end_cols);
   }
}

/* compute
 *              locnx locmy locmz locmyl locmzl
 * nx_border  [   0    A0T  C0T    F0VT G0VT ]
 * myl_border [  F0C    0    0      0    0   ]
 * mzl_border [  G0C    0    0      0    0   ]
 *
 *               [  0 F0C^T  G0C^T ]^T
 *               [ A0   0     0    ]
 *               [ C0   0     0    ]
 * result -= X * [ F0V  0     0    ]
 *               [ G0V  0     0    ]
 */
void
DistributedLinearSystem::finalizeDenseBorderBlocked(BorderLinsys& B, const DenseMatrix& X, DenseMatrix& result,
   int begin_rows, int end_rows) {
   const bool has_RAC = B.has_RAC;

   if (!B.has_RAC && !B.use_local_RAC)
      return;

   const SparseMatrix* F0cons_border = has_RAC ? dynamic_cast<SparseMatrix*>(B.F.first.get()) : nullptr;
   const SparseMatrix* G0cons_border = has_RAC ? dynamic_cast<SparseMatrix*>(B.G.first.get()) : nullptr;

   const SparseMatrix* A0_border = has_RAC ? dynamic_cast<SparseMatrix*>(B.A.first.get()) : nullptr;
   const SparseMatrix* C0_border = has_RAC ? dynamic_cast<SparseMatrix*>(B.C.first.get()) : nullptr;
   const SparseMatrix* F0vec_border = has_RAC ? dynamic_cast<SparseMatrix*>(B.A.last.get()) : &data->getLocalF();
   const SparseMatrix* G0vec_border = has_RAC ? dynamic_cast<SparseMatrix*>(B.C.last.get()) : &data->getLocalG();

   assert(F0vec_border);
   assert(G0vec_border);
   const int n_rows = end_rows - begin_rows;

   const auto mA0 = A0_border ? A0_border->n_rows() : 0;
   const auto mC0 = C0_border ? C0_border->n_rows() : 0;

   long long mF0C{0};
   long long nF0C{0};
   if (F0cons_border) {
      std::tie(mF0C, nF0C) = F0cons_border->n_rows_columns();
   }
   const auto[mF0V, nF0V] = F0vec_border->n_rows_columns();
   const auto nG0V = G0vec_border->n_columns();

   if (!has_RAC && nF0V == 0 && nG0V == 0) {
      return;
   }

#ifndef NDEBUG
   const auto nA0 = A0_border ? A0_border->n_columns() : 0;
   const auto nC0 = C0_border ? C0_border->n_columns() : 0;

   const auto mG0V = G0vec_border->n_rows();

   const auto[mX0, nX0] = X.n_rows_columns();
   const auto[mRes, nRes] = result.n_rows_columns();

   long long mG0C{0};
   long long nG0C{0};
   if (G0cons_border) {
      std::tie(mG0C, nG0C) = G0cons_border->n_rows_columns();
   }

   assert(n_rows <= mX0 && n_rows <= mRes);
   assert(nA0 == nC0);
   assert(nF0V == nG0V);

   if (has_RAC)
      assert(nA0 == nF0V);

   assert(nF0C + mA0 + mC0 + mF0V + mG0V == nRes);
   assert(nF0C == nG0C);

   if (has_RAC)
      assert(nX0 == nF0V + mF0C + mG0C);
   else
      assert(nX0 >= nF0V);
#endif

   if (has_RAC) {
      /*            [  0  ]
       * res -= X * [ F0C ]
       *            [ G0C ]
       */
      X.multMatAt(0, n_rows, nF0V, 1.0, 0, 0, result, -1.0, *F0cons_border);
      X.multMatAt(0, n_rows, nF0V + mF0C, 1.0, 0, 0, result, -1.0, *G0cons_border);

      /*            [ A0T ]
       * res -= X * [  0  ]
       *            [  0  ]
       */
      X.multMatAt(0, n_rows, 0, 1.0, 0, nF0C, result, -1.0, A0_border->getTranspose());

      /*            [ C0T ]
       * res -= X * [  0  ]
       *            [  0  ]
       */
      X.multMatAt(0, n_rows, 0, 1.0, 0, nF0C + mA0, result, -1.0, C0_border->getTranspose());
   }

   /*            [ F0VT ]
    * res -= X * [  0   ]
    *            [  0   ]
    */
   X.multMatAt(0, n_rows, 0, 1.0, 0, nF0C + mA0 + mC0, result, -1.0, F0vec_border->getTranspose());

   /*            [ G0VT ]
    * res -= X * [  0   ]
    *            [  0   ]
    */
   X.multMatAt(0, n_rows, 0, 1.0, 0, nF0C + mA0 + mC0 + mF0V, result, -1.0, G0vec_border->getTranspose());
}


/* calculate res -= X * BT */
void DistributedLinearSystem::multRightDenseBorderBlocked(BorderBiBlock& BT, const DenseMatrix& X, DenseMatrix& result,
   int begin_rows,
   int end_rows) {
   /*
    *        [  RiT   AiT   CiT ]
    * Bi^T = [   0     0     0  ]
    *        [   Fi    0     0  ]
    *        [   Gi    0     0  ]
    */
   const bool with_RAC = BT.has_RAC;

   const auto nX = X.n_columns();
   const auto nR = BT.R.n_columns();
   const auto nA = BT.A.n_columns();
   const auto mF = BT.F.n_rows();
   const auto mG = BT.G.n_rows();


#ifndef NDEBUG
   const auto mX = X.n_rows();
   const auto mR = BT.R.n_rows();
   const auto mA = BT.A.n_rows();
   const auto nF = BT.F.n_columns();
   const auto nG = BT.G.n_columns();

   const int n_empty_rows = BT.n_empty_rows;
   const auto[mRes, nRes] = result.n_rows_columns();

   assert(nF == nG);
   if (with_RAC) {
      const auto[mC, nC] = BT.C.n_rows_columns();
      assert(mR == mA);
      assert(mR == mC);
      assert(nR == nF);

      assert(nR + nA + nC == nRes);
      assert(mR + n_empty_rows + mF + mG == nX);
   } else {
      assert(nF <= nRes);
      assert(mF + mG == nX);
   }

   assert(0 <= begin_rows && begin_rows <= end_rows);
   assert(end_rows - begin_rows <= mX);
   assert(end_rows - begin_rows <= mRes);
#endif

   if (BT.isEmpty())
      return;

   // X from the right with each column of Bi^T todo add OMP to submethods

   const int n_rows = end_rows - begin_rows;

   /*            [ RiT ]
    *            [  0  ]
    * res -= X * [  0  ]
    *            [  Fi ]
    *            [  Gi ]
    */
   if (with_RAC)
      X.multMatAt(0, n_rows, 0, 1.0, 0, 0, result, -1.0, BT.R);

   if (mF > 0)
      X.multMatAt(0, n_rows, nX - mF - mG, 1.0, 0, 0, result, -1.0, BT.F);

   if (mG > 0)
      X.multMatAt(0, n_rows, nX - mG, 1.0, 0, 0, result, -1.0, BT.G);

   if (with_RAC) {
      /*            [ AiT ]
       *            [  0  ]
       * res -= X * [  0  ]
       *            [  0  ]
       *            [  0  ]
       */
      X.multMatAt(0, n_rows, 0, 1.0, 0, nR, result, -1.0, BT.A);

      /*            [ CiT ]
       *            [  0  ]
       * res -= X * [  0  ]
       *            [  0  ]
       *            [  0  ]
       */
      X.multMatAt(0, n_rows, 0, 1.0, 0, nR + nA, result, -1.0, BT.C);
   }
}

void
DistributedLinearSystem::putBiTBorder(DenseMatrix& res, const BorderBiBlock& BiT, int begin_rows, int end_rows) const {
   /* add (Bri)^T to res
    *
    *                [ RiT AiT CiT ]
    *                [  0   0   0  ]
    * Bri^T        = [  Fi  0   0  ]
    *                [  Gi  0   0  ]
    */

   const auto[mRt, nRt] = BiT.R.n_rows_columns();
   const auto nAt = BiT.A.n_columns();
   const auto mF = BiT.F.n_rows();

   const int n_empty_rows = BiT.n_empty_rows;

#ifndef NDEBUG
   const auto nF = BiT.F.n_columns();
   const auto[mG, nG] = BiT.G.n_rows_columns();

   const long long m_border = BiT.has_RAC ? mRt + n_empty_rows + mF + mG : n_empty_rows + mF + mG;

   const auto[mres, nres] = res.n_rows_columns();
   const auto nCt = BiT.C.n_columns();

   assert(mres >= end_rows - begin_rows);
   assert(nF == nG);
   if (BiT.has_RAC) {
      assert(nF == nRt);
      assert(nRt + nAt + nCt == nres);
   } else
      assert(nres >= nF);

   assert(0 <= begin_rows && begin_rows <= end_rows && end_rows <= m_border);
#endif
   if (BiT.isEmpty())
      return;

   const auto length_col = dynamic_cast<SparseSymmetricMatrix&>(*kkt).size();

   const int end_RAC = std::max(0ll, mRt);
   if (BiT.has_RAC && begin_rows < end_RAC) {
      const int begin_block_RAC = begin_rows;
      const int end_block_RAC = std::min(end_RAC, end_rows);
      const int n_rhs_block = end_block_RAC - begin_rows;

      BiT.R.fromGetRowsBlock(begin_block_RAC, n_rhs_block, length_col, 0, res[0]);
      BiT.A.fromGetRowsBlock(begin_block_RAC, n_rhs_block, length_col, nRt, res[0]);
      BiT.C.fromGetRowsBlock(begin_block_RAC, n_rhs_block, length_col, (nRt + nAt), res[0]);
   }

   const int begin_F = end_RAC + n_empty_rows;
   const int end_F = begin_F + mF;
   if (begin_rows < end_F && begin_F < end_rows) {
      const int begin_block_F = std::max(begin_rows, begin_F) - begin_F;
      const int end_block_F = std::min(end_rows, end_F) - begin_F;
      const int n_rhs_block = end_block_F - begin_block_F;

      const int res_start_row = std::max(begin_F - begin_rows, 0);
      assert(res_start_row + n_rhs_block <= mres);

      BiT.F.fromGetRowsBlock(begin_block_F, n_rhs_block, length_col, 0, res[res_start_row]);
   }

   const int begin_G = end_F;
   if (begin_G < end_rows) {
      const int begin_block_G = std::max(begin_rows, begin_G) - begin_G;
      const int end_block_G = end_rows - begin_G;
      const int n_rhs_block = end_block_G - begin_block_G;

      const int res_start_row = std::max(begin_G - begin_rows, 0);
      assert(res_start_row + n_rhs_block <= mres);

      BiT.G.fromGetRowsBlock(begin_block_G, n_rhs_block, length_col, 0, res[res_start_row]);
   }
}

void DistributedLinearSystem::solveCompressed(Vector<double>& rhs_) {
   auto& rhs = dynamic_cast<DistributedVector<double>&>(rhs_);
#ifdef TIMING
   //double tTot=MPI_Wtime();
#endif
   Lsolve(rhs);
   Dsolve(rhs);
   Ltsolve(rhs);
#ifdef TIMING
   //cout << "SolveCompressed took: " << (MPI_Wtime()-tTot) << endl;
#endif
}


/*
 *  y = alpha*Lni^T x + beta*y
 *
 *                       ( [ R 0 0 ]     )
 *  y = beta*y + Di\Li\ (  [ A 0 0 ] * x )
 *                      (  [ C 0 0 ]    )
 */
void DistributedLinearSystem::LniTransMult(SimpleVector<double>& y, double alpha, SimpleVector<double>& x) {
   const SparseMatrix& A = data->getLocalA();
   int N{0}, nx0{0};

   //get nx(parent) from the number of cols of A (or C). Use N as dummy
   if (data->hasRAC())
      std::tie(N, nx0) = A.n_rows_columns();
   // a mild assert
   assert(nx0 <= x.length());

   N = locnx + locmy + locmz;
   assert(y.length() == N);

   //!memopt
   SimpleVector<double> LniTx(N);

   SimpleVector<double> x1(&x[0], nx0);
   SimpleVector<double> LniTx1(&LniTx[0], locnx);

   LniTx1.setToZero();
   if (data->hasRAC()) {
      SimpleVector<double> LniTx2(&LniTx[locnx], locmy);
      SimpleVector<double> LniTx3(&LniTx[locnx + locmy], locmz);

      const SparseMatrix& C = data->getLocalC();
      const SparseMatrix& R = data->getLocalCrossHessian();
      R.mult(0.0, LniTx1, 1.0, x1);
      A.mult(0.0, LniTx2, 1.0, x1);
      C.mult(0.0, LniTx3, 1.0, x1);

   }

   if (locmyl > 0) {
      int nxMyMzP = x.length() - locmyl - locmzl;

      const SparseMatrix& F = data->getLocalF();
      SimpleVector<double> xlink(&x[nxMyMzP], locmyl);

      F.transpose_mult(1.0, LniTx1, 1.0, xlink);
   }

   if (locmzl > 0) {
      int nxMyMzMylP = x.length() - locmzl;

      const SparseMatrix& G = data->getLocalG();
      SimpleVector<double> xlink(&x[nxMyMzMylP], locmzl);

      G.transpose_mult(1.0, LniTx1, 1.0, xlink);
   }

//  solver->Lsolve(LniTx); -> empty
   solver->Dsolve(LniTx);
   y.axpy(alpha, LniTx);
}


/*
 * Computes res += [R^T A^T C^T ] * inv(KKT) * [R 0 F^T G^T ] x
 *                 [0           ]              [A           ]
 *                 [F           ]              [C           ]
 *                 [G           ]
 */

void DistributedLinearSystem::addTermToSchurResidual(SimpleVector<double>& res, SimpleVector<double>& x) {
   const SparseMatrix& A = data->getLocalA();
   const SparseMatrix& C = data->getLocalC();
   const SparseMatrix& F = data->getLocalF();
   const SparseMatrix& G = data->getLocalG();
   const SparseMatrix& R = data->getLocalCrossHessian();

#ifndef NDEBUG
   assert(A.n_rows() == locmy);
   assert(C.n_rows() == locmz);
   assert(F.n_rows() == locmyl);
   assert(G.n_rows() == locmzl);
   assert(R.n_rows() == locnx);

   assert(res.length() >= x.length());
   assert(x.length() >= A.n_columns());
#endif

   // res contains mz buffer part
   int N = locnx + locmy + locmz;
   SimpleVector<double> y(N);

   R.getStorage().mult(0.0, &y[0], 1.0, &x[0]);
   A.getStorage().mult(0.0, &y[locnx], 1.0, &x[0]);
   C.getStorage().mult(0.0, &y[locnx + locmy], 1.0, &x[0]);

   if (locmyl > 0) {
      assert(res.length() == x.length());
      F.getStorage().transMult(1.0, &y[0], 1.0, &x[x.length() - locmyl - locmzl]);
   }

   if (locmzl > 0) {
      assert(res.length() == x.length());
      G.getStorage().transMult(1.0, &y[0], 1.0, &x[x.length() - locmzl]);
   }

   //cout << "4 - y norm:" << y.twonorm() << endl;
   //printf("%g  %g  %g  %g\n", y[locnx+locmy+0], y[locnx+locmy+1], y[locnx+locmy+2], y[locnx+locmy+3]);
   solver->solve(y);

   R.getStorage().transMult(1.0, &res[0], 1.0, &y[0]);
   A.getStorage().transMult(1.0, &res[0], 1.0, &y[locnx]);
   C.getStorage().transMult(1.0, &res[0], 1.0, &y[locnx + locmy]);

   if (locmyl > 0)
      F.getStorage().mult(1.0, &res[res.length() - locmyl - locmzl], 1.0, &y[0]);

   if (locmzl > 0)
      G.getStorage().mult(1.0, &res[res.length() - locmzl], 1.0, &y[0]);
}

void DistributedLinearSystem::addTermToDenseSchurCompl(DenseSymmetricMatrix& SC) {
   const SparseMatrix& A = data->getLocalA();
   const SparseMatrix& C = data->getLocalC();
   const SparseMatrix& F = data->getLocalF();
   const SparseMatrix& G = data->getLocalG();
   const SparseMatrix& R = data->getLocalCrossHessian();

   const bool withR = (R.n_columns() != -1);
   const bool withA = (A.n_columns() != -1);
   const bool withC = (C.n_columns() != -1);


   assert(A.n_rows() == locmy);
   assert(locmyl >= 0);
   assert(locmzl >= 0);

   const int nxMyP = SC.size() - locmyl - locmzl;
   const int nxMyMzP = SC.size() - locmzl;

   auto nRAC = A.n_columns();

   if (nRAC == -1) {
      nRAC = C.n_columns();
   }
   if (nRAC == -1) {
      nRAC = SC.size();
   }

   assert(SC.size() >= nRAC);


   const int N = locnx + locmy + locmz;

   SimpleVector<double> col(N);
   SimpleVector<int> nnzPerColRAC(nRAC);

   if (withR)
      R.addNnzPerCol(nnzPerColRAC);

   if (withA)
      A.addNnzPerCol(nnzPerColRAC);

   if (withC)
      C.addNnzPerCol(nnzPerColRAC);

   const int withMyl = (locmyl > 0);
   const int withMzl = (locmzl > 0);

   for (int it = 0; it < nRAC; it++) {
      if (nnzPerColRAC[it] == 0)
         continue;

      double* const pcol = &col[0];

      for (int it1 = 0; it1 < locnx; it1++)
         pcol[it1] = 0.0;

      R.fromGetDense(0, it, &col[0], 1, locnx, 1);
      A.fromGetDense(0, it, &col[locnx], 1, locmy, 1);
      C.fromGetDense(0, it, &col[locnx + locmy], 1, locmz, 1);

      solver->solve(col);

      //here we have colGi = inv(H_i)* it-th col of Gi^t
      //now do colSC = Gi * inv(H_i)* it-th col of Gi^t

      // SC+=R*x
      R.getStorage().transMult(1.0, &SC[it][0], -1.0, &col[0]);

      // SC+=At*y
      A.getStorage().transMult(1.0, &SC[it][0], -1.0, &col[locnx]);

      // SC+=Ct*z
      C.getStorage().transMult(1.0, &SC[it][0], -1.0, &col[locnx + locmy]);

      // do we have linking equality constraints? If so, set SC+=F*x
      if (withMyl)
         F.getStorage().mult(1.0, &SC[it][nxMyP], -1.0, &col[0]);

      // do we have linking inequality constraints? If so, set SC+=G*x
      if (withMzl)
         G.getStorage().mult(1.0, &SC[it][nxMyMzP], -1.0, &col[0]);
   }

   // do we have linking equality constraints?
   if (withMyl) {
      SimpleVector<int> nnzPerColFt(locmyl);
      F.addNnzPerRow(nnzPerColFt);

      // do column-wise multiplication for columns containing Ft (F transposed)
      for (int it = 0; it < locmyl; it++) {

         if (nnzPerColFt[it] == 0)
            continue;

         double* pcol = &col[0];

         // get it'th column from Ft (i.e., it'th row from F)
         F.fromGetDense(it, 0, &col[0], 1, 1, locnx);

         for (int it1 = locnx; it1 < locnx + locmy + locmz; it1++)
            pcol[it1] = 0.0;

         solver->solve(col);

         R.getStorage().transMult(1.0, &SC[it + nxMyP][0], -1.0, &col[0]);
         A.getStorage().transMult(1.0, &SC[it + nxMyP][0], -1.0, &col[locnx]);
         C.getStorage().transMult(1.0, &SC[it + nxMyP][0], -1.0, &col[locnx + locmy]);

         // here we have colGi = inv(H_i)* (it + locnx + locmy)-th col of Gi^t
         // now do colSC = Gi * inv(H_i)* (it + locnx + locmy)-th col of Gi^t

         F.getStorage().mult(1.0, &SC[it + nxMyP][nxMyP], -1.0, &col[0]);

         if (withMzl)
            G.getStorage().mult(1.0, &SC[it + nxMyP][nxMyMzP], -1.0, &col[0]);
      }
   }

   // do we have linking inequality constraints?
   if (withMzl) {
      SimpleVector<int> nnzPerColGt(locmzl);
      G.addNnzPerRow(nnzPerColGt);

      // do column-wise multiplication for columns containing Gt (G transposed)
      for (int it = 0; it < locmzl; it++) {
         if (nnzPerColGt[it] == 0)
            continue;

         double* pcol = &col[0];

         // get it'th column from Gt (i.e., it'th row from G)
         G.fromGetDense(it, 0, &col[0], 1, 1, locnx);

         for (int it1 = locnx; it1 < locnx + locmy + locmz; it1++)
            pcol[it1] = 0.0;

         solver->solve(col);

         R.getStorage().transMult(1.0, &SC[it + nxMyMzP][0], -1.0, &col[0]);
         A.getStorage().transMult(1.0, &SC[it + nxMyMzP][0], -1.0, &col[locnx]);
         C.getStorage().transMult(1.0, &SC[it + nxMyMzP][0], -1.0, &col[locnx + locmy]);

         // here we have colGi = inv(H_i)* (it + locnx + locmy + locmyl)-th col of Gi^t
         // now do colSC = Gi * inv(H_i)* (it + locnx + locmy + locmyl)-th col of Gi^t

         if (withMyl)
            F.getStorage().mult(1.0, &SC[it + nxMyMzP][nxMyP], -1.0, &col[0]);

         G.getStorage().mult(1.0, &SC[it + nxMyMzP][nxMyMzP], -1.0, &col[0]);
      }
   }
}

/* res += [ Bl^T Ki^{-1} BT ]^T for cols begin_rows to end_rows in res */
void DistributedLinearSystem::addBiTLeftKiDenseToResBlockedParallelSolvers(bool sparse_res, bool sym_res,
   const BorderBiBlock& BlT,
   /* const */ DenseMatrix& BT, AbstractMatrix& result, int begin_rows_res, int end_rows_res) {
#ifndef NDEBUG
   const auto[m_res, n_res] = result.n_rows_columns();
   assert(m_res >= 0 && n_res >= 0);
   if (sym_res)
      assert(m_res == n_res);
#endif

   if (BlT.isEmpty())
      return;

   const auto nB = BT.n_columns();
   assert(0 <= begin_rows_res && begin_rows_res <= end_rows_res);
   const int n_cols = end_rows_res - begin_rows_res;
   assert(n_cols <= BT.n_rows());
   assert(end_rows_res <= m_res);

   const int chunk_length = blocksizemax * PIPSgetnOMPthreads();

   if (colId.empty() || colId.size() < static_cast<unsigned int>(chunk_length))
      colId.resize(chunk_length);

   /* note that B is in storage transposed for rhs access reasons */
#ifdef TIME_SCHUR
   const double t_start = omp_get_wtime();
#endif

   //
   //     res +=  Bl^T  K^-1 B
   //
   const int chunks = std::ceil(static_cast<double>(n_cols) / chunk_length);

   for (int i = 0; i < chunks; i++) {
      const int actual_blocksize = std::min((i + 1) * chunk_length, n_cols) - (i * chunk_length);

      double* colsBlockDense_loc = BT[i * chunk_length];
      solver->solve(actual_blocksize, colsBlockDense_loc, nullptr);

      for (int j = 0; j < actual_blocksize; ++j) {
         colId[j] = j + i * chunk_length + begin_rows_res;
         assert(colId[j] < m_res);
      }

      addLeftBorderTimesDenseColsToResTransp(BlT, colsBlockDense_loc, colId.data(), nB, actual_blocksize, sparse_res,
         sym_res, result);
   }

#ifdef TIME_SCHUR
   const double t_end = omp_get_wtime();
   std::cout << "t_end - t_start:" << (t_end - t_start) << std::endl;
   assert(0);
#endif
}

/* res is in transposed form here */
/*
 * calculate
 *
 * res += border_left_transposed K^-1 border_right only for begin_cols, end_cols
 *
 *                [ R 0 F G ]                          [ R A C ]
 * Border_right = [ A 0 0 0 ] border_left_transposed = [ 0 0 0 ]
 *                [ C 0 0 0 ]                          [ F 0 0 ]
 *                                                     [ G 0 0 ]
 */
void DistributedLinearSystem::addBiTLeftKiBiRightToResBlockedParallelSolvers(bool sparse_res, bool sym_res,
   const BorderBiBlock& border_left_transp,
   /* const */ BorderBiBlock& border_right, AbstractMatrix& result, int begin_cols, int end_cols, int begin_rows_res,
   int end_rows_res) {
   if (sparse_res)
      assert(sym_res);

   const auto nF_right = border_right.F.n_columns();
   const auto nG_right = border_right.G.n_columns();

   const auto[mR_r, nR_r] = border_right.R.n_rows_columns();
   const auto mA_r = border_right.A.n_rows();

   const bool with_RAC = border_right.has_RAC;
   const bool withF = (nF_right > 0);
   const bool withG = (nG_right > 0);

   const long long length_col = dynamic_cast<SparseSymmetricMatrix&>(*kkt).size();

#ifndef NDEBUG
   const auto mF_right = border_right.F.n_rows();
   const auto mG_right = border_right.G.n_rows();
   const auto nA_r = border_right.A.n_columns();

   const auto[n_res_tp, m_res_tp] = result.n_rows_columns();

   assert(0 <= begin_cols && begin_cols <= end_cols);
   assert(end_cols - begin_cols <= n_res_tp);

   const auto[mF_left, nF_left] = border_left_transp.F.n_rows_columns();
   const auto[mG_left, nG_left] = border_left_transp.G.n_rows_columns();

   if (border_left_transp.has_RAC) {
      const auto[mR_left, nR_left] = border_left_transp.R.n_rows_columns();
      const auto[mA_left, nA_left] = border_left_transp.A.n_rows_columns();
      const auto[mC_left, nC_left] = border_left_transp.C.n_rows_columns();

      assert(mR_left == mA_left);
      assert(mR_left == mC_left);
      assert(nR_left == nF_left);
      assert(nR_left == nG_left);

      assert(mR_left + border_left_transp.n_empty_rows + mF_left + mG_left == m_res_tp);
      assert(nR_left + nA_left + nC_left == length_col);
   } else {
      assert(nF_left == nG_left);
      assert(nF_left <= length_col);
      assert(border_left_transp.n_empty_rows + mF_left + mG_left == m_res_tp);
   }

   if (with_RAC) {
      const auto[mC_right, nC_right] = border_right.C.n_rows_columns();

      assert(nR_r == nA_r);
      assert(nR_r == nC_right);
      assert(mR_r == mF_right);
      assert(mR_r == mG_right);

      assert(end_cols <= nR_r + border_right.n_empty_rows + nF_right + nG_right);
      assert(mR_r + mA_r + mC_right == length_col);
   } else {
      assert(mF_right == mG_right);
      assert(mF_right <= length_col);

      /* less equal since there could be link0vars */
      assert(end_cols <= border_right.n_empty_rows + nF_right + nG_right);
   }
#endif

   if (border_left_transp.isEmpty() || border_right.isEmpty())
      return;

   const int chunk_length = blocksizemax * PIPSgetnOMPthreads();

   if (colsBlockDense.empty() || colsBlockDense.size() < static_cast<unsigned int>(chunk_length * length_col))
      colsBlockDense.resize(chunk_length * length_col);
   if (colId.empty() || colId.size() < static_cast<unsigned int>(chunk_length))
      colId.resize(chunk_length);

   // indicating whether a right hand side is zero - deactivated since problemlems in Pardiso
#if 0
   if( colSparsity.empty() )
      colSparsity.resize( length_col * chunk_length );
#endif

   int* colSparsity_ptr = nullptr;
   if (!colSparsity.empty())
      colSparsity_ptr = colSparsity.data();

#ifdef TIME_SCHUR
   const double t_start = omp_get_wtime();
#endif


   if (with_RAC) {
      //                       (R)
      //     SC +=  B^T  K^-1  (A)
      //                       (C)
      if (begin_cols < nR_r) {
         const int begin_block_RAC = begin_cols;
         const int end_block_RAC = std::min(end_cols, (int) nR_r);

         const int n_cols = end_block_RAC - begin_block_RAC;
         // TODO : add buffer for nonzeros and do not reallocate all the time
         SimpleVector<int> nnzPerColRAC(n_cols);

         border_right.R.addNnzPerCol(nnzPerColRAC, begin_block_RAC, end_block_RAC);
         border_right.A.addNnzPerCol(nnzPerColRAC, begin_block_RAC, end_block_RAC);
         border_right.C.addNnzPerCol(nnzPerColRAC, begin_block_RAC, end_block_RAC);

         const int chunks_RAC = std::ceil(static_cast<double>(n_cols) / chunk_length);

         for (int i = 0; i < chunks_RAC; i++) {
            assert(i * chunk_length + begin_block_RAC <= end_block_RAC);

            const int actual_blocksize =
               std::min((i + 1) * chunk_length + begin_block_RAC, end_block_RAC) - (i * chunk_length + begin_block_RAC);
            assert(0 <= actual_blocksize);
            assert(i * chunk_length + actual_blocksize <= n_cols);

            int nrhs = 0;

            for (int j = 0; j < actual_blocksize; ++j)
               if (nnzPerColRAC[j + i * chunk_length] != 0)
                  colId[nrhs++] = j + i * chunk_length + begin_block_RAC;

            if (nrhs == 0)
               continue;

            std::fill(colsBlockDense.begin(), colsBlockDense.end(), 0);

            border_right.R.fromGetColsBlock(colId.data(), nrhs, length_col, 0, colsBlockDense.data(), colSparsity_ptr);
            border_right.A.fromGetColsBlock(colId.data(), nrhs, length_col, mR_r, colsBlockDense.data(),
               colSparsity_ptr);
            border_right.C.fromGetColsBlock(colId.data(), nrhs, length_col, (mR_r + mA_r), colsBlockDense.data(),
               colSparsity_ptr);

            solver->solve(nrhs, colsBlockDense.data(), colSparsity_ptr);

            /* map indices back to buffer */
            for (int j = 0; j < nrhs; ++j) {
               colId[j] -= begin_block_RAC + begin_rows_res;
               assert(colId[j] < n_res_tp);
            }

            addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense.data(), colId.data(), length_col,
               nrhs, sparse_res, sym_res,
               result);
         }
      }
   }

   // do we have linking equality constraints?
   if (withF) {
      //                       (F^T)
      //     SC +=  B^T  K^-1  (0  )
      //                       (0  )
      const int begin_F = with_RAC ? nR_r + border_right.n_empty_rows : border_right.n_empty_rows;
      const int end_F = begin_F + nF_right;
      if (begin_F < end_cols && begin_cols < end_F) {
         const int begin_block_F = std::max(begin_F, begin_cols) - begin_F;
         const int end_block_F = std::min(end_cols, end_F) - begin_F;
         assert(0 <= begin_block_F && begin_block_F <= end_block_F && end_block_F <= nF_right);

         const int n_cols = end_block_F - begin_block_F;

         SimpleVector<int> nnzPerColFt(n_cols);
         border_right.F.addNnzPerCol(nnzPerColFt, begin_block_F, end_block_F);

         const int chunks_F = std::ceil(static_cast<double>(n_cols) / chunk_length);

         // do block-wise multiplication for columns of F^T part
         for (int i = 0; i < chunks_F; ++i) {
            const int actual_blocksize =
               std::min((i + 1) * chunk_length + begin_block_F, end_block_F) - (i * chunk_length + begin_block_F);

            int nrhs = 0;

            for (int j = 0; j < actual_blocksize; ++j)
               if (nnzPerColFt[j + i * chunk_length] != 0)
                  colId[nrhs++] = j + i * chunk_length + begin_block_F;

            if (nrhs == 0)
               continue;

            std::fill(colsBlockDense.begin(), colsBlockDense.end(), 0);

            // get column block from Ft (i.e., row block from F)
            border_right.F.fromGetColsBlock(colId.data(), nrhs, length_col, 0, colsBlockDense.data(), colSparsity_ptr);

            solver->solve(nrhs, colsBlockDense.data(), colSparsity_ptr);

            for (int j = 0; j < nrhs; ++j) {
               colId[j] += begin_F - begin_cols + begin_rows_res;
               assert(colId[j] < n_res_tp);
            }

            addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense.data(), colId.data(), length_col,
               nrhs, sparse_res, sym_res,
               result);
         }
      }
   }

   // do we have linking inequality constraints?
   if (withG) {
      //                       (G^T)
      //     SC +=  B^T  K^-1  (0  )
      //                       (0  )
      const int begin_G = with_RAC ? nR_r + border_right.n_empty_rows + nF_right : border_right.n_empty_rows + nF_right;
      const int end_G = begin_G + nG_right;

      if (begin_G < end_cols && begin_cols < end_G) {
         const int begin_block_G = std::max(begin_G, begin_cols) - begin_G;
         const int end_block_G = std::min(end_cols, end_G) - begin_G;
         assert(0 <= begin_block_G && begin_block_G <= end_block_G && end_block_G <= nG_right);

         const int n_cols = end_block_G - begin_block_G;

         SimpleVector<int> nnzPerColGt(n_cols);
         border_right.G.addNnzPerCol(nnzPerColGt, begin_block_G, end_block_G);

         const int chunks_G = std::ceil(static_cast<double>(n_cols) / chunk_length);

         // do block-wise multiplication for columns of G^T part
         for (int i = 0; i < chunks_G; ++i) {
            const int actual_blocksize =
               std::min((i + 1) * chunk_length + begin_block_G, end_block_G) - (i * chunk_length + begin_block_G);

            int nrhs = 0;

            for (int j = 0; j < actual_blocksize; ++j)
               if (nnzPerColGt[j + i * chunk_length] != 0)
                  colId[nrhs++] = j + i * chunk_length + begin_block_G;

            if (nrhs == 0)
               continue;

            std::fill(colsBlockDense.begin(), colsBlockDense.end(), 0);

            border_right.G.fromGetColsBlock(colId.data(), nrhs, length_col, 0, colsBlockDense.data(), colSparsity_ptr);

            solver->solve(nrhs, colsBlockDense.data(), colSparsity_ptr);

            for (int j = 0; j < nrhs; ++j) {
               colId[j] += begin_G - begin_cols + begin_rows_res;
               assert(colId[j] < n_res_tp);
            }

            addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense.data(), colId.data(), length_col,
               nrhs, sparse_res, sym_res,
               result);
         }
      }
   }

#if 0
   // debug stuff
   const int myrank = PIPS_MPIgetRank();

   if( myrank == 0 )
   {
      static int iteration = 0;
      ofstream myfile;
      char filename[50];
      sprintf(filename, "../blocked_%d_%d.txt", myrank, iteration);
      myfile.open(filename);
      iteration++;
      SC.write_to_stream(myfile); // todo write out in each iteration with global counter and MPI rank!
      myfile.close();

      assert(0);
   }
#endif


#ifdef TIME_SCHUR
   const double t_end = omp_get_wtime();
   std::cout << "t_end - t_start:" << (t_end - t_start) << std::endl;
   assert(0);
#endif
}

void
DistributedLinearSystem::addLeftBorderTimesDenseColsToResTranspSparse(const BorderBiBlock& Bl, const double* cols,
   const int* cols_id, int length_col,
   int n_cols, SparseSymmetricMatrix& res) {
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * cols = border_left * cols
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *  the size of the zero rows in border_left is determined by res and can be zero
    */
   const auto mF = Bl.F.n_rows();
   const auto mG = Bl.G.n_rows();
   const auto nRes = res.n_columns();

   const bool with_RAC = Bl.has_RAC;
   const bool with_F = mF > 0;
   const bool with_G = mG > 0;

#ifndef NDEBUG
   const auto nF = Bl.F.n_columns();
   const auto nG = Bl.G.n_columns();
   const auto mRes = res.n_rows();

   assert(mRes == nRes);
   if (with_RAC) {
      const auto[mR, nR] = Bl.R.n_rows_columns();
      const auto[mA, nA] = Bl.A.n_rows_columns();
      const auto[mC, nC] = Bl.C.n_rows_columns();

      assert(nF == nG && nF == nR);
      assert(length_col == nR + nA + nC);
      assert(nRes >= mR + mF + mG);
      assert(mR == mA && mA == mC);
   } else {
      /* >= since there could be 0linkvars */
      assert(nRes >= mF + mG);
   }
#endif

   if (Bl.isEmpty())
      return;

   // multiply each column with left_border and add if to res
#pragma omp parallel for schedule(dynamic, 10)
   for (int it_col = 0; it_col < n_cols; it_col++) {
      const double* const col = &cols[it_col * length_col];
      const int row_res = cols_id[it_col];

      assert(row_res < mRes);
      if (with_RAC) {
         const long long nR = Bl.R.n_columns();
         const long long nA = Bl.A.n_columns();

         Bl.R.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, 0);
         Bl.A.multMatSymUpper(1.0, res, -1.0, &col[nR], row_res, 0);
         Bl.C.multMatSymUpper(1.0, res, -1.0, &col[nR + nA], row_res, 0);
      }

      if (with_F)
         Bl.F.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, nRes - mF - mG);

      if (with_G)
         Bl.G.multMatSymUpper(1.0, res, -1.0, &col[0], row_res, nRes - mG);
   }
}

void
DistributedLinearSystem::addLeftBorderTimesDenseColsToResTranspDense(const BorderBiBlock& Bl, const double* cols,
   const int* cols_id, int length_col,
   int n_cols, int n_cols_res, double** res) const {
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * colsBlockDense = border_left * colsBlockDense
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *
    *  cols lie as rows in storage
    *  the size of the zero rows in border_left is determined by res and can be zero
    */
   const auto mF = Bl.F.n_rows();
   const auto mG = Bl.G.n_rows();

   const bool with_RAC = Bl.has_RAC;
   const bool with_F = mF > 0;
   const bool with_G = mG > 0;

#ifndef NDEBUG
   const auto nF = Bl.F.n_columns();
   const auto nG = Bl.G.n_columns();

   if (with_RAC) {
      const auto[mR, nR] = Bl.R.n_rows_columns();
      const auto[mA, nA] = Bl.A.n_rows_columns();
      const auto[mC, nC] = Bl.C.n_rows_columns();

      assert(mR == mA && mA == mC);
      assert(nF == nG && nF == nR);
      assert(length_col == nR + nA + nC);
      assert(n_cols_res >= mR + mF + mG);
   } else
      assert(n_cols_res >= mF + mG);

   assert(n_cols_res >= 1);
#endif

   if (Bl.isEmpty())
      return;

   // multiply each column with left factor of SC
#pragma omp parallel for schedule(dynamic, 10)
   for (int it_col = 0; it_col < n_cols; it_col++) {
      const double* const col = &cols[it_col * length_col];
      const int row_res = cols_id[it_col];

      if (with_RAC) {
         const auto nR = Bl.R.n_columns();
         const auto nA = Bl.A.n_columns();
         Bl.R.getStorage().mult(1.0, &res[row_res][0], -1.0, &col[0]);
         Bl.A.getStorage().mult(1.0, &res[row_res][0], -1.0, &col[nR]);
         Bl.C.getStorage().mult(1.0, &res[row_res][0], -1.0, &col[nR + nA]);
      }

      if (with_F)
         Bl.F.getStorage().mult(1.0, &res[row_res][n_cols_res - mF - mG], -1.0, &col[0]);
      if (with_G)
         Bl.G.getStorage().mult(1.0, &res[row_res][n_cols_res - mG], -1.0, &col[0]);
   }
}

void
DistributedLinearSystem::addLeftBorderTimesDenseColsToResTransp(const BorderBiBlock& border_left, const double* cols,
   const int* cols_id,
   int length_col, int blocksize, bool sparse_res, bool sym_res, AbstractMatrix& res) const {
   if (border_left.isEmpty())
      return;

   assert(cols_id && cols);

   if (sparse_res) {
      assert(sym_res);
      auto& res_sparse = dynamic_cast<SparseSymmetricMatrix&>(res);
      addLeftBorderTimesDenseColsToResTranspSparse(border_left, cols, cols_id, length_col, blocksize, res_sparse);
   } else {
      double** res_array;
      int res_ncols;

#ifndef NDEBUG
      int res_mrows;
      if (sym_res)
         res_mrows = dynamic_cast<DenseSymmetricMatrix&>(res).size();
      else
         res_mrows = dynamic_cast<DenseMatrix&>(res).mStorage->m;
      for (int i = 0; i < blocksize; ++i)
         assert(cols_id[i] < res_mrows);
#endif

      if (sym_res) {
         auto& res_dense = dynamic_cast<DenseSymmetricMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.size();
      } else {
         auto& res_dense = dynamic_cast<DenseMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.mStorage->n;
      }

      addLeftBorderTimesDenseColsToResTranspDense(border_left, cols, cols_id, length_col, blocksize, res_ncols,
         res_array);
   }
}

int DistributedLinearSystem::allocateAndZeroBlockedComputationsBuffer(int buffer_m, int buffer_n) {
   assert(blocksize_hierarchical > 0);
   assert(buffer_m > 0);
   assert(buffer_n > 0);

   const int buffer_m_blocked = sc_compute_blockwise_hierarchical ? PIPSgetnOMPthreads() * blocksize_hierarchical
      : buffer_m;

   if (!buffer_blocked_hierarchical)
      buffer_blocked_hierarchical = std::make_unique<DenseMatrix>(buffer_m_blocked, buffer_n);
   else {
      const auto[mbuf, nbuf] = buffer_blocked_hierarchical->n_rows_columns();
      if (mbuf < buffer_m_blocked || nbuf < buffer_n)
         buffer_blocked_hierarchical = std::make_unique<DenseMatrix>(buffer_m_blocked, buffer_n);
   }
   buffer_blocked_hierarchical->putZeros();

   return buffer_m_blocked;
}