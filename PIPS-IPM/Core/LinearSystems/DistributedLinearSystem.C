/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#include "StochOptions.h"
#include "BorderedSymMatrix.h"
#include "DistributedLinearSystem.h"
#include "DistributedTree.h"
#include "DistributedFactory.h"
#include "DistributedQP.hpp"
#include "SparseLinearAlgebraPackage.h"
#include "math.h"

#include "pipsport.h"
#include "omp.h"

DistributedLinearSystem::DistributedLinearSystem(DistributedFactory* factory_, DistributedQP* problem, bool is_hierarchy_root) : LinearSystem(
      factory_, problem), data{problem}, computeBlockwiseSC(pips_options::getBoolParameter("SC_COMPUTE_BLOCKWISE")),
      blocksizemax(pips_options::getIntParameter("SC_BLOCKWISE_BLOCKSIZE_MAX")), is_hierarchy_root(is_hierarchy_root),
      blocksize_hierarchical(pips_options::getIntParameter("SC_BLOCKSIZE_HIERARCHICAL")),
      sc_compute_blockwise_hierarchical{pips_options::getBoolParameter("SC_HIERARCHICAL_COMPUTE_BLOCKWISE")}, stochNode{factory_->tree} {
   if (sc_compute_blockwise_hierarchical && PIPS_MPIgetRank() == 0)
      std::cout << "Computing hierarchical Schur complements blockwise with buffersize " << blocksize_hierarchical
                << " (times # of available OMP threads)\n";

   if (pips_options::getBoolParameter("HIERARCHICAL"))
      assert(is_hierarchy_root);

   problem->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);

   //get the communicator from one of the vectors
   const auto& dds = dynamic_cast<const DistributedVector<double>&>(*primal_diagonal);
   this->mpiComm = dds.mpiComm;
   this->iAmDistrib = dds.iAmDistrib;
}

DistributedLinearSystem::DistributedLinearSystem(DistributedFactory* factory_, DistributedQP* problem, Vector<double>* dd_, Vector<double>* dq_,
      Vector<double>* nomegaInv_, Vector<double>* primal_reg_, Vector<double>* dual_y_reg_, Vector<double>* dual_z_reg_, Vector<double>* rhs_,
      bool create_iter_ref_vecs) : LinearSystem(factory_, problem, dd_, dq_, nomegaInv_, primal_reg_, dual_y_reg_, dual_z_reg_, rhs_,
      create_iter_ref_vecs), data{problem}, computeBlockwiseSC(pips_options::getBoolParameter("SC_COMPUTE_BLOCKWISE")),
      blocksizemax(pips_options::getIntParameter("SC_BLOCKWISE_BLOCKSIZE_MAX")),
      blocksize_hierarchical(pips_options::getIntParameter("SC_BLOCKSIZE_HIERARCHICAL")),
      sc_compute_blockwise_hierarchical{pips_options::getBoolParameter("SC_HIERARCHICAL_COMPUTE_BLOCKWISE")}, stochNode{factory_->tree} {
   problem->getLocalSizes(locnx, locmy, locmz, locmyl, locmzl);

   if (primal_diagonal) {
      const auto& primal_diagonal_stoch = dynamic_cast<const DistributedVector<double>&>(*primal_diagonal);
      mpiComm = primal_diagonal_stoch.mpiComm;
      iAmDistrib = primal_diagonal_stoch.iAmDistrib;
   }
   else {
      mpiComm = MPI_COMM_NULL;
      iAmDistrib = false;
   }

   useRefs = 1;
}


void DistributedLinearSystem::joinRHS(Vector<double>& rhs_in, const Vector<double>& rhs1_in, const Vector<double>& rhs2_in,
      const Vector<double>& rhs3_in) const {
   DistributedVector<double>& rhs = dynamic_cast<DistributedVector<double>&>(rhs_in);
   const DistributedVector<double>& rhs1 = dynamic_cast<const DistributedVector<double>&>(rhs1_in);
   const DistributedVector<double>& rhs2 = dynamic_cast<const DistributedVector<double>&>(rhs2_in);
   const DistributedVector<double>& rhs3 = dynamic_cast<const DistributedVector<double>&>(rhs3_in);

   rhs.jointCopyFrom(rhs1, rhs2, rhs3);
}

void DistributedLinearSystem::separateVars(Vector<double>& x_in, Vector<double>& y_in, Vector<double>& z_in, const Vector<double>& vars_in) const {
   DistributedVector<double>& x = dynamic_cast<DistributedVector<double>&>(x_in);
   DistributedVector<double>& y = dynamic_cast<DistributedVector<double>&>(y_in);
   DistributedVector<double>& z = dynamic_cast<DistributedVector<double>&>(z_in);
   const DistributedVector<double>& vars = dynamic_cast<const DistributedVector<double>&>(vars_in);

   vars.jointCopyTo(x, y, z);
}

void DistributedLinearSystem::factorize(Problem* problem_, Variables* vars) {
#ifdef TIMING
   double tTot = MPI_Wtime();
#endif
   // the call to the the parent's method takes care of all necessary updates
   // to the KKT system (updating diagonals mainly). This is done recursively,
   // we don't have to worry about it anymore.
   LinearSystem::factorize(problem_, vars);

   // now DO THE LINEAR ALGEBRA!

   DistributedQP* problem = dynamic_cast<DistributedQP*>(problem_);
   // in order to avoid a call to QpGenLinsys::factor, call factor2 method.
   factor2(problem, vars);
//  assembleKKT(problem, vars);
//  allreduceAndFactorKKT(problem, vars);

#ifdef TIMING
   tTot = MPI_Wtime() - tTot;
   MPI_Barrier(MPI_COMM_WORLD);
   const int rank == PIPS_MPIgetRANK(MPI_COMM_WORLD);
   //if( 128 * ( myRank / 128 ) == 0 )
   if( 0 == myRank )
       std::cout << "Outer fact. total time " << tTot << std::endl;
#endif
}


/**
 * Computes U = Li\Gi^T.
 *        [ 0 0 0 ]
 * Gi^T = [ A 0 0 ]
 *        [ C 0 0]
 *
 * We have special structure here:
 *             [ 0 ]
 *   U   = Li\ [ A ] ,   U is (nx+my+mz)-by-(np)
 *             [ C ]
 *
 *   V = Di\U
 */
void DistributedLinearSystem::computeU_V(DistributedQP* problem, DenseGenMatrix* U, DenseGenMatrix* V) {
   U->scalarMult(0.0);
   V->scalarMult(0.0);
   assert(false); //need code to deal with cross Hessian term
   SparseGenMatrix& A = problem->getLocalA();
   SparseGenMatrix& C = problem->getLocalC();

   int N, nxP;
   A.getSize(N, nxP);
   assert(N == locmy);

   N = locnx + locmy + locmz;
   SimpleVector<double> uCol(N);

   for (int it = 0; it < nxP; it++) {

      double* p = &uCol[0];
      for (int it1 = 0; it1 < locnx; it1++)
         p[it1] = 0.0;

      A.fromGetDense(0, it, &uCol[locnx], 1, locmy, 1);
      C.fromGetDense(0, it, &uCol[locnx + locmy], 1, locmz, 1);

      solver->Lsolve(uCol);
      U->atPutDense(0, it, &uCol[0], 1, N, 1);

      solver->Dsolve(uCol);
      V->atPutDense(0, it, &uCol[0], 1, N, 1);

   }
}

void DistributedLinearSystem::factorize_with_correct_inertia() {
   assert(false);
   regularization_strategy->notify_new_step();

   /* factor once without applying regularization */
   solver->matrixChanged();
   if (!solver->reports_inertia()) {
      return;
   }

   double last_primal_regularization{0.0};
   double last_dual_equality_regularization{0.0};
   double last_dual_inequality_regularization{0.0};

   while (!regularization_strategy->is_inertia_correct(solver->get_inertia())) {
      auto[primal_regularization_value, dual_equality_regularization_value, dual_inequality_regularization_value] =
      this->regularization_strategy->get_regularization_parameters(solver->get_inertia(), barrier_parameter_current_iterate);

      assert(primal_regularization_value >= last_primal_regularization);
      assert(dual_equality_regularization_value >= last_dual_equality_regularization);
      assert(dual_inequality_regularization_value >= last_dual_inequality_regularization);

      this->add_regularization_local_kkt(primal_regularization_value - last_primal_regularization,
            dual_equality_regularization_value - last_dual_equality_regularization,
            dual_inequality_regularization_value - last_dual_inequality_regularization);
      solver->matrixChanged();
   }
}

void DistributedLinearSystem::allocU(DenseGenMatrix** U, int n0) {
   int lines, cols;
   if (*U == nullptr) {
      *U = new DenseGenMatrix(locnx + locmy + locmz, n0);
   }
   else {
      (*U)->getSize(lines, cols);

      if (lines != locnx + locmy + locmz || cols != n0) {

         delete (*U);
         *U = new DenseGenMatrix(locnx + locmy + locmz, n0);
      }
   }
}

void DistributedLinearSystem::allocV(DenseGenMatrix** V, int n0) {
   int lines, cols;
   if (*V == nullptr)
      *V = new DenseGenMatrix(locnx + locmy + locmz, n0);
   else {
      (*V)->getSize(lines, cols);

      if (lines != locnx + locmy + locmz || cols != n0) {

         delete (*V);
         *V = new DenseGenMatrix(locnx + locmy + locmz, n0);
      }
   }
}

/**
 *       [ R^i^T Ai^T Ci^T ]          [    ]
 * z0 -= [ 0      0   0    ] * Li\Di\ [ zi ]
 *       [ 0      0   0    ]          [    ]
 *
 * 
 */
void DistributedLinearSystem::addLnizi(DistributedQP* problem, Vector<double>& z0_, Vector<double>& zi_) {
   SimpleVector<double>& z0 = dynamic_cast<SimpleVector<double>&>(z0_);
   SimpleVector<double>& zi = dynamic_cast<SimpleVector<double>&>(zi_);

   solver->Dsolve(zi);
   solver->Ltsolve(zi);

   SparseGenMatrix& A = problem->getLocalA();
   SparseGenMatrix& C = problem->getLocalC();
   SparseGenMatrix& R = problem->getLocalCrossHessian();

   //get n0= nx(parent)= #cols of A or C
   int dummy, n0;
   A.getSize(dummy, n0);

   // zi2 and zi3 are just references to fragments of zi
   SimpleVector<double> zi1(&zi[0], locnx);
   SimpleVector<double> zi2(&zi[locnx], locmy);
   SimpleVector<double> zi3(&zi[locnx + locmy], locmz);
   // same for z01 (only the first n0 entries in the output z0 are computed)
   SimpleVector<double> z01(&z0[0], n0);

   R.transMult(1.0, z01, -1.0, zi1);
   A.transMult(1.0, z01, -1.0, zi2);
   C.transMult(1.0, z01, -1.0, zi3);
}

void
DistributedLinearSystem::finalizeDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseGenMatrix& result, int begin_rows, int end_rows) {
   /* compute BiT_buffer += X_j^T Bmodj for all j */
   for (auto& border_mod_block : border_mod) {
      if (border_mod_block.border.isEmpty())
         continue;
      finalizeDenseBorderBlocked(border_mod_block.border, border_mod_block.multiplier, result, begin_rows, end_rows);
   }
}

void
DistributedLinearSystem::multRightDenseBorderModBlocked(std::vector<BorderMod>& border_mod, DenseGenMatrix& result, int begin_cols, int end_cols) {
   /* compute BiT_buffer += X_j^T Bmodj for all j */
   for (auto& border_mod_block : border_mod) {
      std::unique_ptr<BorderBiBlock> BiT_mod{};

      BorderLinsys& border = border_mod_block.border;

      if (border.isEmpty())
         continue;

      if (border.use_local_RAC)
         BiT_mod.reset(
               new BorderBiBlock(data->getLocalCrossHessian().getTranspose(), data->getLocalA().getTranspose(), data->getLocalC().getTranspose(),
                     border.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
      else if (border.has_RAC)
         BiT_mod.reset(new BorderBiBlock(dynamic_cast<SparseGenMatrix&>(*border.R.mat).getTranspose(),
               dynamic_cast<SparseGenMatrix&>(*border.A.mat).getTranspose(), dynamic_cast<SparseGenMatrix&>(*border.C.mat).getTranspose(),
               border.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat)));
      else
         BiT_mod.reset(
               new BorderBiBlock(border.n_empty_rows, dynamic_cast<SparseGenMatrix&>(*border.F.mat), dynamic_cast<SparseGenMatrix&>(*border.G.mat),
                     false));

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
DistributedLinearSystem::finalizeDenseBorderBlocked(BorderLinsys& B, const DenseGenMatrix& X, DenseGenMatrix& result, int begin_rows, int end_rows) {
   const bool has_RAC = B.has_RAC;

   if (!B.has_RAC && !B.use_local_RAC)
      return;

   SparseGenMatrix* F0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.F.mat) : nullptr;
   SparseGenMatrix* G0cons_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.G.mat) : nullptr;

   SparseGenMatrix* A0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.A.mat) : nullptr;
   SparseGenMatrix* C0_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.C.mat) : nullptr;
   SparseGenMatrix* F0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.A.mat_link) : &data->getLocalF();
   SparseGenMatrix* G0vec_border = has_RAC ? dynamic_cast<SparseGenMatrix*>(B.C.mat_link) : &data->getLocalG();

   assert(F0vec_border);
   assert(G0vec_border);
   const int n_rows = end_rows - begin_rows;

   int mA0{0};
   int nA0{0};
   if (A0_border) {
      A0_border->getSize(mA0, nA0);
   }

   int mC0{0};
   int nC0{0};
   if (C0_border) {
      C0_border->getSize(mC0, nC0);
   }

   int mF0C{0};
   int nF0C{0};
   if (F0cons_border) {
      F0cons_border->getSize(mF0C, nF0C);
   }

   int mF0V{0};
   int nF0V{0};
   F0vec_border->getSize(mF0V, nF0V);

   int mG0V{0};
   int nG0V{0};
   G0vec_border->getSize(mG0V, nG0V);

   if (!has_RAC && nF0V == 0 && nG0V == 0) {
      return;
   }

#ifndef NDEBUG
   int mX0, nX0;
   X.getSize(mX0, nX0);
   int mRes, nRes;
   result.getSize(mRes, nRes);
   assert(n_rows <= mX0 && n_rows <= mRes);

   int mG0C{0};
   int nG0C{0};
   if (G0cons_border) {
      G0cons_border->getSize(mG0C, nG0C);
   }

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
void DistributedLinearSystem::multRightDenseBorderBlocked(BorderBiBlock& BT, const DenseGenMatrix& X, DenseGenMatrix& result, int begin_rows,
      int end_rows) {
   /*
    *        [  RiT   AiT   CiT ]
    * Bi^T = [   0     0     0  ]
    *        [   Fi    0     0  ]
    *        [   Gi    0     0  ]
    */
   const bool with_RAC = BT.has_RAC;
   int mX, nX;
   X.getSize(mX, nX);
   int mR, nR;
   BT.R.getSize(mR, nR);
   int mA, nA;
   BT.A.getSize(mA, nA);
   int mF, nF;
   BT.F.getSize(mF, nF);
   int mG, nG;
   BT.G.getSize(mG, nG);


#ifndef NDEBUG
   const int n_empty_rows = BT.n_empty_rows;
   int mRes, nRes;
   result.getSize(mRes, nRes);

   assert(nF == nG);
   if (with_RAC) {
      int mC, nC;
      BT.C.getSize(mC, nC);
      assert(mR == mA);
      assert(mR == mC);
      assert(nR == nF);

      assert(nR + nA + nC == nRes);
      assert(mR + n_empty_rows + mF + mG == nX);
   }
   else {
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

void DistributedLinearSystem::putBiTBorder(DenseGenMatrix& res, const BorderBiBlock& BiT, int begin_rows, int end_rows) const {
   /* add (Bri)^T to res
    *
    *                [ RiT AiT CiT ]
    *                [  0   0   0  ]
    * Bri^T        = [  Fi  0   0  ]
    *                [  Gi  0   0  ]
    */

   int mRt, nRt;
   BiT.R.getSize(mRt, nRt);
   int mAt, nAt;
   BiT.A.getSize(mAt, nAt);
   int mF, nF;
   BiT.F.getSize(mF, nF);
   int mG, nG;
   BiT.G.getSize(mG, nG);

   const int n_empty_rows = BiT.n_empty_rows;

#ifndef NDEBUG
   const int m_border = BiT.has_RAC ? mRt + n_empty_rows + mF + mG : n_empty_rows + mF + mG;

   int mres, nres;
   res.getSize(mres, nres);
   int mCt, nCt;
   BiT.C.getSize(mCt, nCt);
   assert(mres >= end_rows - begin_rows);
   assert(nF == nG);
   if (BiT.has_RAC) {
      assert(nF == nRt);
      assert(nRt + nAt + nCt == nres);
   }
   else
      assert(nres >= nF);

   assert(0 <= begin_rows && begin_rows <= end_rows && end_rows <= m_border);
#endif
   if (BiT.isEmpty())
      return;

   const int length_col = dynamic_cast<SparseSymMatrix&>(*kkt).size();

   const int end_RAC = std::max(0, mRt);
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
   DistributedVector<double>& rhs = dynamic_cast<DistributedVector<double>&>(rhs_);
#ifdef TIMING
   //double tTot=MPI_Wtime();
#endif
   Lsolve(data, rhs);
   Dsolve(data, rhs);
   Ltsolve(data, rhs);
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
void DistributedLinearSystem::LniTransMult(DistributedQP* problem, SimpleVector<double>& y, double alpha, SimpleVector<double>& x) {
   SparseGenMatrix& A = problem->getLocalA();
   int N{0}, nx0{0};

   //get nx(parent) from the number of cols of A (or C). Use N as dummy
   if (data->hasRAC())
      A.getSize(N, nx0);
   // a mild assert
   assert(nx0 <= x.length());

   N = locnx + locmy + locmz;
   assert(y.length() == N);

   //!memopt
   SimpleVector<double> LniTx(N);

   // shortcuts
   SimpleVector<double> x1(&x[0], nx0);
   SimpleVector<double> LniTx1(&LniTx[0], locnx);

   LniTx1.setToZero();
   if (data->hasRAC()) {
      SimpleVector<double> LniTx2(&LniTx[locnx], locmy);
      SimpleVector<double> LniTx3(&LniTx[locnx + locmy], locmz);

      SparseGenMatrix& C = problem->getLocalC();
      SparseGenMatrix& R = problem->getLocalCrossHessian();
      R.mult(0.0, LniTx1, 1.0, x1);
      A.mult(0.0, LniTx2, 1.0, x1);
      C.mult(0.0, LniTx3, 1.0, x1);

   }

   if (locmyl > 0) {
      int nxMyMzP = x.length() - locmyl - locmzl;

      SparseGenMatrix& F = problem->getLocalF();
      SimpleVector<double> xlink(&x[nxMyMzP], locmyl);

      F.transMult(1.0, LniTx1, 1.0, xlink);
   }

   if (locmzl > 0) {
      int nxMyMzMylP = x.length() - locmzl;

      SparseGenMatrix& G = problem->getLocalG();
      SimpleVector<double> xlink(&x[nxMyMzMylP], locmzl);

      G.transMult(1.0, LniTx1, 1.0, xlink);
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

void DistributedLinearSystem::addTermToSchurResidual(DistributedQP* problem, SimpleVector<double>& res, SimpleVector<double>& x) {
   SparseGenMatrix& A = problem->getLocalA();
   SparseGenMatrix& C = problem->getLocalC();
   SparseGenMatrix& F = problem->getLocalF();
   SparseGenMatrix& G = problem->getLocalG();
   SparseGenMatrix& R = problem->getLocalCrossHessian();

   int nxP, aux;
   A.getSize(aux, nxP);
   assert(aux == locmy);
   C.getSize(aux, nxP);
   assert(aux == locmz);
   F.getSize(aux, nxP);
   assert(aux == locmyl);
   G.getSize(aux, nxP);
   assert(aux == locmzl);
   R.getSize(aux, nxP);
   assert(aux == locnx);

   // res contains mz buffer part
   assert(res.length() >= x.length());
   assert(x.length() >= nxP);

   int N = locnx + locmy + locmz;
   SimpleVector<double> y(N);

   R.mult(0.0, &y[0], 1, 1.0, &x[0], 1);
   A.mult(0.0, &y[locnx], 1, 1.0, &x[0], 1);
   C.mult(0.0, &y[locnx + locmy], 1, 1.0, &x[0], 1);

   if (locmyl > 0) {
      assert(res.length() == x.length());
      F.transMult(1.0, &y[0], 1, 1.0, &x[x.length() - locmyl - locmzl], 1);
   }

   if (locmzl > 0) {
      assert(res.length() == x.length());
      G.transMult(1.0, &y[0], 1, 1.0, &x[x.length() - locmzl], 1);
   }

   //cout << "4 - y norm:" << y.twonorm() << endl;
   //printf("%g  %g  %g  %g\n", y[locnx+locmy+0], y[locnx+locmy+1], y[locnx+locmy+2], y[locnx+locmy+3]);
   solver->solve(y);

   R.transMult(1.0, &res[0], 1, 1.0, &y[0], 1);
   A.transMult(1.0, &res[0], 1, 1.0, &y[locnx], 1);
   C.transMult(1.0, &res[0], 1, 1.0, &y[locnx + locmy], 1);

   if (locmyl > 0)
      F.mult(1.0, &res[res.length() - locmyl - locmzl], 1, 1.0, &y[0], 1);

   if (locmzl > 0)
      G.mult(1.0, &res[res.length() - locmzl], 1, 1.0, &y[0], 1);
}

void DistributedLinearSystem::addTermToDenseSchurCompl(DistributedQP* problem, DenseSymMatrix& SC) {
   SparseGenMatrix& A = problem->getLocalA();
   SparseGenMatrix& C = problem->getLocalC();
   SparseGenMatrix& F = problem->getLocalF();
   SparseGenMatrix& G = problem->getLocalG();
   SparseGenMatrix& R = problem->getLocalCrossHessian();

   int N, nxP;

   R.getSize(N, nxP);
   const bool withR = (nxP != -1);

   A.getSize(N, nxP);
   const bool withA = (nxP != -1);

   assert(N == locmy);
   assert(locmyl >= 0);
   assert(locmzl >= 0);

   const int NP = SC.size();
   assert(NP >= nxP);

   const int nxMyP = NP - locmyl - locmzl;
   const int nxMyMzP = NP - locmzl;

   if (nxP == -1)
      C.getSize(N, nxP);

   int N2, nxP2;
   C.getSize(N2, nxP2);
   const bool withC = (nxP2 != -1);

   if (nxP == -1)
      nxP = NP;

   N = locnx + locmy + locmz;

   SimpleVector<double> col(N);
   SimpleVector<int> nnzPerColRAC(nxP);

   if (withR)
      R.addNnzPerCol(nnzPerColRAC);

   if (withA)
      A.addNnzPerCol(nnzPerColRAC);

   if (withC)
      C.addNnzPerCol(nnzPerColRAC);

   const int withMyl = (locmyl > 0);
   const int withMzl = (locmzl > 0);

   for (int it = 0; it < nxP; it++) {
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
      R.transMult(1.0, &SC[it][0], 1, -1.0, &col[0], 1);

      // SC+=At*y
      A.transMult(1.0, &SC[it][0], 1, -1.0, &col[locnx], 1);

      // SC+=Ct*z
      C.transMult(1.0, &SC[it][0], 1, -1.0, &col[locnx + locmy], 1);

      // do we have linking equality constraints? If so, set SC+=F*x
      if (withMyl)
         F.mult(1.0, &SC[it][nxMyP], 1, -1.0, &col[0], 1);

      // do we have linking inequality constraints? If so, set SC+=G*x
      if (withMzl)
         G.mult(1.0, &SC[it][nxMyMzP], 1, -1.0, &col[0], 1);
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

         R.transMult(1.0, &SC[it + nxMyP][0], 1, -1.0, &col[0], 1);
         A.transMult(1.0, &SC[it + nxMyP][0], 1, -1.0, &col[locnx], 1);
         C.transMult(1.0, &SC[it + nxMyP][0], 1, -1.0, &col[locnx + locmy], 1);

         // here we have colGi = inv(H_i)* (it + locnx + locmy)-th col of Gi^t
         // now do colSC = Gi * inv(H_i)* (it + locnx + locmy)-th col of Gi^t

         F.mult(1.0, &SC[it + nxMyP][nxMyP], 1, -1.0, &col[0], 1);

         if (withMzl)
            G.mult(1.0, &SC[it + nxMyP][nxMyMzP], 1, -1.0, &col[0], 1);
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

         R.transMult(1.0, &SC[it + nxMyMzP][0], 1, -1.0, &col[0], 1);
         A.transMult(1.0, &SC[it + nxMyMzP][0], 1, -1.0, &col[locnx], 1);
         C.transMult(1.0, &SC[it + nxMyMzP][0], 1, -1.0, &col[locnx + locmy], 1);

         // here we have colGi = inv(H_i)* (it + locnx + locmy + locmyl)-th col of Gi^t
         // now do colSC = Gi * inv(H_i)* (it + locnx + locmy + locmyl)-th col of Gi^t

         if (withMyl)
            F.mult(1.0, &SC[it + nxMyMzP][nxMyP], 1, -1.0, &col[0], 1);

         G.mult(1.0, &SC[it + nxMyMzP][nxMyMzP], 1, -1.0, &col[0], 1);
      }
   }
}

/* res += [ Bl^T Ki^{-1} BT ]^T for cols begin_rows to end_rows in res */
void DistributedLinearSystem::addBiTLeftKiDenseToResBlockedParallelSolvers(bool sparse_res, bool sym_res, const BorderBiBlock& BlT,
      /* const */ DenseGenMatrix& BT, DoubleMatrix& result, int begin_rows_res, int end_rows_res) {
#ifndef NDEBUG
   int m_res, n_res;
   result.getSize(m_res, n_res);
   assert(m_res >= 0 && n_res >= 0);
   if (sym_res)
      assert(m_res == n_res);
#endif

   if (BlT.isEmpty())
      return;

   int mB, nB;
   BT.getSize(mB, nB);
   assert(0 <= begin_rows_res && begin_rows_res <= end_rows_res);
   const int n_cols = end_rows_res - begin_rows_res;
   assert(n_cols <= mB);
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

      addLeftBorderTimesDenseColsToResTransp(BlT, colsBlockDense_loc, colId.data(), nB, actual_blocksize, sparse_res, sym_res, result);
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
void DistributedLinearSystem::addBiTLeftKiBiRightToResBlockedParallelSolvers(bool sparse_res, bool sym_res, const BorderBiBlock& border_left_transp,
      /* const */ BorderBiBlock& border_right, DoubleMatrix& result, int begin_cols, int end_cols, int begin_rows_res, int end_rows_res) {
   if (sparse_res)
      assert(sym_res);

   int mF_right, nF_right;
   border_right.F.getSize(mF_right, nF_right);
   int mG_right, nG_right;
   border_right.G.getSize(mG_right, nG_right);

   int nR_r, mR_r;
   border_right.R.getSize(mR_r, nR_r);
   int nA_r, mA_r;
   border_right.A.getSize(mA_r, nA_r);

   const bool with_RAC = border_right.has_RAC;
   const bool withF = (nF_right > 0);
   const bool withG = (nG_right > 0);
   const int length_col = dynamic_cast<SparseSymMatrix&>(*kkt).size();

#ifndef NDEBUG
   int m_res_tp, n_res_tp;
   result.getSize(n_res_tp, m_res_tp);

   assert(0 <= begin_cols && begin_cols <= end_cols);
   assert(end_cols - begin_cols <= n_res_tp);

   int mF_left, nF_left;
   border_left_transp.F.getSize(mF_left, nF_left);
   int mG_left, nG_left;
   border_left_transp.G.getSize(mG_left, nG_left);
   if (border_left_transp.has_RAC) {
      int mR_left, nR_left;
      border_left_transp.R.getSize(mR_left, nR_left);
      int mA_left, nA_left;
      border_left_transp.A.getSize(mA_left, nA_left);
      int mC_left, nC_left;
      border_left_transp.C.getSize(mC_left, nC_left);
      assert(mR_left == mA_left);
      assert(mR_left == mC_left);
      assert(nR_left == nF_left);
      assert(nR_left == nG_left);

      assert(mR_left + border_left_transp.n_empty_rows + mF_left + mG_left == m_res_tp);
      assert(nR_left + nA_left + nC_left == length_col);
   }
   else {
      assert(nF_left == nG_left);
      assert(nF_left <= length_col);
      assert(border_left_transp.n_empty_rows + mF_left + mG_left == m_res_tp);
   }

   if (with_RAC) {
      int mC_right, nC_right;
      border_right.C.getSize(mC_right, nC_right);
      assert(nR_r == nA_r);
      assert(nR_r == nC_right);
      assert(mR_r == mF_right);
      assert(mR_r == mG_right);

      assert(end_cols <= nR_r + border_right.n_empty_rows + nF_right + nG_right);
      assert(mR_r + mA_r + mC_right == length_col);
   }
   else {
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
         const int end_block_RAC = std::min(end_cols, nR_r);

         const int n_cols = end_block_RAC - begin_block_RAC;
         // TODO : add buffer for nonzeros and do not reallocate all the time
         SimpleVector<int> nnzPerColRAC(n_cols);

         border_right.R.addNnzPerCol(nnzPerColRAC, begin_block_RAC, end_block_RAC);
         border_right.A.addNnzPerCol(nnzPerColRAC, begin_block_RAC, end_block_RAC);
         border_right.C.addNnzPerCol(nnzPerColRAC, begin_block_RAC, end_block_RAC);

         const int chunks_RAC = std::ceil(static_cast<double>(n_cols) / chunk_length);

         for (int i = 0; i < chunks_RAC; i++) {
            assert(i * chunk_length + begin_block_RAC <= end_block_RAC);

            const int actual_blocksize = std::min((i + 1) * chunk_length + begin_block_RAC, end_block_RAC) - (i * chunk_length + begin_block_RAC);
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
            border_right.A.fromGetColsBlock(colId.data(), nrhs, length_col, mR_r, colsBlockDense.data(), colSparsity_ptr);
            border_right.C.fromGetColsBlock(colId.data(), nrhs, length_col, (mR_r + mA_r), colsBlockDense.data(), colSparsity_ptr);

            solver->solve(nrhs, colsBlockDense.data(), colSparsity_ptr);

            /* map indices back to buffer */
            for (int j = 0; j < nrhs; ++j) {
               colId[j] -= begin_block_RAC + begin_rows_res;
               assert(colId[j] < n_res_tp);
            }

            addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense.data(), colId.data(), length_col, nrhs, sparse_res, sym_res,
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
            const int actual_blocksize = std::min((i + 1) * chunk_length + begin_block_F, end_block_F) - (i * chunk_length + begin_block_F);

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

            addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense.data(), colId.data(), length_col, nrhs, sparse_res, sym_res,
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
            const int actual_blocksize = std::min((i + 1) * chunk_length + begin_block_G, end_block_G) - (i * chunk_length + begin_block_G);

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

            addLeftBorderTimesDenseColsToResTransp(border_left_transp, colsBlockDense.data(), colId.data(), length_col, nrhs, sparse_res, sym_res,
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
      SC.writeToStream(myfile); // todo write out in each iteration with global counter and MPI rank!
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
DistributedLinearSystem::addLeftBorderTimesDenseColsToResTranspSparse(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col,
      int n_cols, SparseSymMatrix& res) const {
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * cols = border_left * cols
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *  the size of the zero rows in border_left is determined by res and can be zero
    */
   int mF, nF;
   Bl.F.getSize(mF, nF);
   int mG, nG;
   Bl.G.getSize(mG, nG);
   int mRes, nRes;
   res.getSize(mRes, nRes);

   const bool with_RAC = Bl.has_RAC;
   const bool with_F = mF > 0;
   const bool with_G = mG > 0;

#ifndef NDEBUG
   assert(mRes == nRes);
   if (with_RAC) {
      int mR, nR;
      Bl.R.getSize(mR, nR);
      int mA, nA;
      Bl.A.getSize(mA, nA);
      int mC, nC;
      Bl.C.getSize(mC, nC);
      assert(nF == nG && nF == nR);
      assert(length_col == nR + nA + nC);
      assert(nRes >= mR + mF + mG);
      assert(mR == mA && mA == mC);
   }
   else
      /* >= since there could be 0linkvars */
      assert(nRes >= mF + mG);
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
         int mR, nR;
         Bl.R.getSize(mR, nR);
         int mA, nA;
         Bl.A.getSize(mA, nA);

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
DistributedLinearSystem::addLeftBorderTimesDenseColsToResTranspDense(const BorderBiBlock& Bl, const double* cols, const int* cols_id, int length_col,
      int n_cols, int n_cols_res, double** res) const {
   /*                  [ R A C ]
    * compute res^T += [ 0 0 0 ] * colsBlockDense = border_left * colsBlockDense
    *                  [ F 0 0 ]
    *                  [ G 0 0 ]
    *
    *  cols lie as rows in storage
    *  the size of the zero rows in border_left is determined by res and can be zero
    */
   int mF, nF;
   Bl.F.getSize(mF, nF);
   int mG, nG;
   Bl.G.getSize(mG, nG);

   const bool with_RAC = Bl.has_RAC;
   const bool with_F = mF > 0;
   const bool with_G = mG > 0;

#ifndef NDEBUG
   if (with_RAC) {
      int mR, nR;
      Bl.R.getSize(mR, nR);
      int mA, nA;
      Bl.A.getSize(mA, nA);
      int mC, nC;
      Bl.C.getSize(mC, nC);
      assert(mR == mA && mA == mC);
      assert(nF == nG && nF == nR);
      assert(length_col == nR + nA + nC);
      assert(n_cols_res >= mR + mF + mG);
   }
   else
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
         int mR, nR;
         Bl.R.getSize(mR, nR);
         int mA, nA;
         Bl.A.getSize(mA, nA);
         Bl.R.mult(1.0, &res[row_res][0], 1, -1.0, &col[0], 1);
         Bl.A.mult(1.0, &res[row_res][0], 1, -1.0, &col[nR], 1);
         Bl.C.mult(1.0, &res[row_res][0], 1, -1.0, &col[nR + nA], 1);
      }

      if (with_F)
         Bl.F.mult(1.0, &res[row_res][n_cols_res - mF - mG], 1, -1.0, &col[0], 1);
      if (with_G)
         Bl.G.mult(1.0, &res[row_res][n_cols_res - mG], 1, -1.0, &col[0], 1);
   }
}

void DistributedLinearSystem::addLeftBorderTimesDenseColsToResTransp(const BorderBiBlock& border_left, const double* cols, const int* cols_id,
      int length_col, int blocksize, bool sparse_res, bool sym_res, DoubleMatrix& res) const {
   if (border_left.isEmpty())
      return;

   assert(cols_id && cols);

   if (sparse_res) {
      assert(sym_res);
      SparseSymMatrix& res_sparse = dynamic_cast<SparseSymMatrix&>(res);
      addLeftBorderTimesDenseColsToResTranspSparse(border_left, cols, cols_id, length_col, blocksize, res_sparse);
   }
   else {
      double** res_array;
      int res_ncols;

#ifndef NDEBUG
      int res_mrows;
      if (sym_res)
         res_mrows = dynamic_cast<DenseSymMatrix&>(res).size();
      else
         res_mrows = dynamic_cast<DenseGenMatrix&>(res).mStorage->m;
      for (int i = 0; i < blocksize; ++i)
         assert(cols_id[i] < res_mrows);
#endif

      if (sym_res) {
         DenseSymMatrix& res_dense = dynamic_cast<DenseSymMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.size();
      }
      else {
         DenseGenMatrix& res_dense = dynamic_cast<DenseGenMatrix&>(res);
         res_array = res_dense.mStorage->M;
         res_ncols = res_dense.mStorage->n;
      }

      addLeftBorderTimesDenseColsToResTranspDense(border_left, cols, cols_id, length_col, blocksize, res_ncols, res_array);
   }
}

int DistributedLinearSystem::allocateAndZeroBlockedComputationsBuffer(int buffer_m, int buffer_n) {
   assert(blocksize_hierarchical > 0);
   assert(buffer_m > 0);
   assert(buffer_n > 0);

   const int buffer_m_blocked = sc_compute_blockwise_hierarchical ? PIPSgetnOMPthreads() * blocksize_hierarchical : buffer_m;

   if (!buffer_blocked_hierarchical)
      buffer_blocked_hierarchical.reset(new DenseGenMatrix(buffer_m_blocked, buffer_n));
   else {
      int mbuf, nbuf;
      buffer_blocked_hierarchical->getSize(mbuf, nbuf);
      if (mbuf < buffer_m_blocked || nbuf < buffer_n)
         buffer_blocked_hierarchical.reset(new DenseGenMatrix(buffer_m_blocked, buffer_n));
   }
   buffer_blocked_hierarchical->putZeros();

   return buffer_m_blocked;
}

template<>
bool DistributedLinearSystem::BorderLinsys::isEmpty() const {
   if (use_local_RAC)
      return false;
   else {
      if (F.numberOfNonZeros() == 0 && G.numberOfNonZeros() == 0) {
         if (has_RAC)
            return R.numberOfNonZeros() == 0 && A.numberOfNonZeros() == 0 && C.numberOfNonZeros() == 0;
         else
            return true;
      }
      else
         return false;
   }
};

template<>
bool DistributedLinearSystem::BorderBiBlock::isEmpty() const {
   if (use_local_RAC)
      return false;
   else {
      if (F.numberOfNonZeros() == 0 && G.numberOfNonZeros() == 0) {
         if (has_RAC)
            return R.numberOfNonZeros() == 0 && A.numberOfNonZeros() == 0 && C.numberOfNonZeros() == 0;
         else
            return true;
      }
      else
         return false;
   }
};
