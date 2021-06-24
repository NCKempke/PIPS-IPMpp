/* PIPS
   Authors: Cosmin Petra and Miles Lubin
   See license and copyright information in the documentation */

#include <utility>

#include "DistributedRootLinearSystem.h"
#include "DistributedFactory.hpp"
#include "DistributedProblem.hpp"
#include "DistributedDummyLinearSystem.h"
#include "DistributedLeafLinearSystem.h"
#include "PIPSIPMppOptions.h"

/*********************************************************************/
/************************** ROOT *************************************/
/*********************************************************************/

#ifdef STOCH_TESTING
double g_scenNum;
#endif

DistributedRootLinearSystem::DistributedRootLinearSystem(const DistributedFactory& factory_, DistributedProblem* prob_,
   bool is_hierarchy_root) : DistributedLinearSystem(factory_, prob_, is_hierarchy_root) {
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      assert(is_hierarchy_root);
   init();
}

DistributedRootLinearSystem::DistributedRootLinearSystem(const DistributedFactory& factory_, DistributedProblem* prob_,
   std::shared_ptr<Vector<double>> dd_, std::shared_ptr<Vector<double>> dq_, std::shared_ptr<Vector<double>> nomegaInv_,
   std::shared_ptr<Vector<double>> primal_reg_, std::shared_ptr<Vector<double>> dual_y_reg_,
   std::shared_ptr<Vector<double>> dual_z_reg_, std::shared_ptr<Vector<double>> rhs_) : DistributedLinearSystem(
   factory_,
   prob_, std::move(dd_), std::move(dq_), std::move(nomegaInv_), std::move(primal_reg_), std::move(dual_y_reg_), std::move(dual_z_reg_), std::move(rhs_), true) {
   init();
}

void DistributedRootLinearSystem::init() {
   createChildren();

   precondSC = SCsparsifier(mpiComm);
   usePrecondDist = pipsipmpp_options::get_bool_parameter("PRECONDITION_DISTRIBUTED");

   // use sparse KKT if (enough) 2 links are present
   hasSparseKkt = data->exploitingLinkStructure();
   allreduce_kkt = pipsipmpp_options::get_bool_parameter("ALLREDUCE_SCHUR_COMPLEMENT");
   if (pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      assert(allreduce_kkt);

   usePrecondDist = usePrecondDist && hasSparseKkt && iAmDistrib;
   MatrixEntryTriplet_mpi = MPI_DATATYPE_NULL;

   initProperChildrenRange();
}

DistributedRootLinearSystem::~DistributedRootLinearSystem() {
   delete kktDist;
   delete[] sparseKktBuffer;
}

void DistributedRootLinearSystem::assembleKKT() {
   if (is_hierarchy_root)
      assert(children.size() == 1);

   /* set kkt to zero */
   initializeKKT();

   /* important that int separate loops! else block in Allreduce might occur */
   for (auto& c : children)
      c->assembleKKT();
   for (auto& c : children)
      c->allreduceAndFactorKKT();

   /* build KKT from local children */
   assembleLocalKKT();
}

void DistributedRootLinearSystem::allreduceAndFactorKKT() {
   reduceKKT();

   finalizeKKT();

   factorizeKKT();
}

void DistributedRootLinearSystem::factor2() {
   if (PIPS_MPIgetRank(mpiComm) == 0) {
      if (is_hierarchy_root) {
         assert(children.size() == 1);
         std::cout << "\n Building dense outer Schur Complement...\n\n";
      } else if (data->isHierarchyInnerRoot())
         std::cout << " Building sparse top level Schur Complement...\n\n";
   }

   /* set kkt to zero */
   initializeKKT();

   // First tell children to factorize.
   for (auto& c : children)
      c->factor2();

   /* build KKT from local children */
   assembleLocalKKT();

#ifdef TIMING
   MPI_Barrier(MPI_COMM_WORLD);
   stochNode->resMon.recReduceTmLocal_start();
#endif

   reduceKKT();
#ifdef TIMING
   stochNode->resMon.recReduceTmLocal_stop();
#endif

   finalizeKKT();

   if (PIPS_MPIgetRank(mpiComm) == 0) {
      if (is_hierarchy_root)
         std::cout << " Dense Schur Complement factorization is starting... ";
      else if (data->isHierarchyInnerRoot())
         std::cout << " Sparse top level Schur Complement factorization is starting... ";
   }

   factorizeKKT();

   if (PIPS_MPIgetRank(mpiComm) == 0) {
      if (is_hierarchy_root || data->isHierarchyInnerRoot())
         std::cout << "done\n\n";
   }

#ifdef TIMING
   afterFactor();
#endif
}

#ifdef TIMING
void sLinsysRoot::afterFactor()
{
  int mype; MPI_Comm_rank(mpiComm, &mype);

  if( (mype/256)*256==mype) {
    for (size_t c=0; c<children.size(); c++) {
      if (children[c]->mpiComm == MPI_COMM_NULL) continue;
      
      printf("  rank %d NODE %4zu SPFACT %g BACKSOLVE %g SEC ITER %d\n", mype, c,
        children[c]->stochNode->resMon.eFact.tmLocal,
        children[c]->stochNode->resMon.eFact.tmChildren, (int)g_iterNumber);
    }
  }
  if( (mype/1024)*1024==mype) {
    for (size_t c=0; c<children.size(); c++) {
      if (children[c]->mpiComm == MPI_COMM_NULL) continue;

      double redall = stochNode->resMon.eReduce.tmLocal;
      double redscat = stochNode->resMon.eReduceScatter.tmLocal;
      printf("  rank %d REDUCE %g SEC ITER %d REDSCAT %g DIFF %g\n", mype, redall,
        (int)g_iterNumber, redscat, redall-redscat);
    }
  }
}
#endif

/* compute
 *             locnx locmy locmz locmyl locmzl
 *           [   0    A0T   C0T   F0VT   G0VT ] nx_border
 * buffer += [  F0C    0     0     0      0   ] myl_border
 *           [  G0C    0     0     0      0   ] mzl_border
 *
 * [  0 F0C^T  G0C^T ]^T
 * [ A0   0     0    ]
 * [ C0   0     0    ]
 * [ F0V  0     0    ]   + buffer
 * [ G0V  0     0    ]
 */
// TODO : move to aug..
void
DistributedRootLinearSystem::finalizeZ0Hierarchical(DenseMatrix& buffer, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, int begin_rows, int end_rows) {
   assert(0 <= begin_rows && begin_rows <= end_rows);
   // TODO : parallelize over MPI procs?
   finalizeDenseBorderModBlocked(Br_mod_border, buffer, begin_rows, end_rows);

   if (sc_compute_blockwise_hierarchical)
      assert(end_rows - begin_rows <= buffer.n_rows());
   else
      assert(end_rows <= buffer.n_rows());

   if (!Br.has_RAC && !Br.use_local_RAC)
      return;

   bool has_RAC = Br.has_RAC;

   SparseMatrix* F0cons_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.F.first.get()) : nullptr;
   SparseMatrix* G0cons_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.G.first.get()) : nullptr;

   SparseMatrix* A0_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.A.first.get()) : nullptr;
   SparseMatrix* C0_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.C.first.get()) : nullptr;

   const SparseMatrix* F0vec_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.A.last.get()) : &data->getLocalF();
   const SparseMatrix* G0vec_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.C.last.get()) : &data->getLocalG();

   if (has_RAC)
      assert(F0cons_border && G0cons_border && A0_border && C0_border);
   assert(F0vec_border);
   assert(G0vec_border);

   int mA0{0};
   int nA0{0};
   if (A0_border)
      std::tie(mA0, nA0) = A0_border->n_rows_columns();

   int mC0 = C0_border ? C0_border->n_rows() : 0;

   int mF0C{0};
   int nF0C{0};
   if (F0cons_border)
      std::tie(mF0C, nF0C) = F0cons_border->n_rows_columns();

   int mG0C = G0cons_border ? G0cons_border->n_rows() : 0;
   const auto[mF0V, nF0V] = F0vec_border->n_rows_columns();
   int mG0V = G0vec_border->n_rows();

#ifndef NDEBUG
   long long nC0 = C0_border ? C0_border->n_columns() : 0;
   long long nG0C = G0cons_border ? G0cons_border->n_columns() : 0;
   long long nG0V = G0vec_border->n_columns();

   assert(nA0 == nC0);
   assert(nF0V == nG0V);

   if (has_RAC)
      assert(nA0 == nF0V);

   assert(nF0C == nG0C);
   if (mA0 != 0)
      assert(nF0C + mA0 + mC0 + mF0V + mG0V == buffer.n_columns());

   if (!sc_compute_blockwise_hierarchical) {
      if (has_RAC)
         assert(buffer.n_rows() >= nF0V + mF0C + mG0C);
      else
         assert(buffer.n_rows() >= nF0V);
   }
#endif

   /* add A0^T, C0^T, F0V^T, G0V^T */
   if (begin_rows < nF0V) {
      const int end_f0vblock = std::min((int) nF0V, end_rows);

      /* A0^T */
      if (mA0 > 0)
         buffer.addMatAt(A0_border->getTranspose(), begin_rows, end_f0vblock, 0, nF0C);

      /* C0^T */
      if (mC0 > 0)
         buffer.addMatAt(C0_border->getTranspose(), begin_rows, end_f0vblock, 0, nF0C + mA0);

      /* F0V^T */
      if (mF0V > 0)
         buffer.addMatAt(F0vec_border->getTranspose(), begin_rows, end_f0vblock, 0, nF0C + mA0 + mC0);

      /* G0V^T */
      if (mG0V > 0)
         buffer.addMatAt(G0vec_border->getTranspose(), begin_rows, end_f0vblock, 0, nF0C + mA0 + mC0 + mF0V);
   }

   /* F0C */
   {
      const int start_F0C_block = nA0;
      const int end_F0C_block = nA0 + mF0C;

      if (mF0C > 0 && begin_rows < end_F0C_block && start_F0C_block <= end_rows) {
         const int start_F0C_mat = std::max(begin_rows, start_F0C_block) - start_F0C_block;
         const int end_F0C_mat = std::min(end_rows, end_F0C_block) - start_F0C_block;
         assert(0 <= start_F0C_mat && start_F0C_mat <= end_F0C_mat);
         buffer.addMatAt(*F0cons_border, start_F0C_mat, end_F0C_mat, std::max(0, start_F0C_block - begin_rows), 0);
      }
   }

   /* G0C */
   {
      const int start_G0C_block = nA0 + mF0C;
      const int end_G0C_block = nA0 + mF0C + mG0C;

      if (mG0C > 0 && begin_rows < end_G0C_block && start_G0C_block <= end_rows) {
         const int start_G0C_mat = std::max(begin_rows, start_G0C_block) - start_G0C_block;
         const int end_G0C_mat = std::min(end_rows, end_G0C_block) - start_G0C_block;
         assert(0 <= start_G0C_mat && start_G0C_mat <= end_G0C_mat);
         buffer.addMatAt(*G0cons_border, start_G0C_mat, end_G0C_mat, std::max(0, start_G0C_block - begin_rows), 0);
      }
   }
}

/* result -= [ Br0^T X0 ]^T = X0^T Br0
 *
 * compute result -= Br0^T X0
 *         [  0  A0T C0T F0VT G0VT ]
 * Br0^T = [ F0C  0   0   0    0   ]
 *         [ G0C  0   0   0    0   ]
 *
 * result and X0 are stored in transposed form
 *
 * result -= X0 Br0 instead
 *
 * Br0 = [  0  F0CT G0CT ]
 *       [  A   0    0   ]
 *       [  C   0    0   ]
 *       [ F0V  0    0   ]
 *       [ G0V  0    0   ]
 *
 */
void DistributedRootLinearSystem::finalizeInnerSchurComplementContribution(AbstractMatrix& result, DenseMatrix& X0,
   BorderLinsys& Br, bool is_sym, bool is_sparse,
   int begin_rows, int end_rows) {

   if (is_sparse)
      assert(is_sym);

   const bool has_RAC = Br.has_RAC;
   if (!has_RAC && !Br.use_local_RAC)
      return;

   SparseMatrix* F0cons_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.F.first.get()) : nullptr;
   SparseMatrix* G0cons_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.G.first.get()) : nullptr;

   SparseMatrix* A0_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.A.first.get()) : nullptr;
   SparseMatrix* C0_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.C.first.get()) : nullptr;
   const SparseMatrix* F0vec_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.A.last.get()) : &data->getLocalF();
   const SparseMatrix* G0vec_border = has_RAC ? dynamic_cast<SparseMatrix*>(Br.C.last.get()) : &data->getLocalG();

   assert(F0vec_border);
   assert(G0vec_border);

   if (!has_RAC && F0vec_border->n_columns() == 0 && G0vec_border->n_columns() == 0)
      return;

#ifndef NDEBUG
   const auto[m_result, n_result] = result.n_rows_columns();
   const auto[mX0, nX0] = X0.n_rows_columns();
   const auto[mF0V, nF0V] = F0vec_border->n_rows_columns();
   const auto[mG0V, nG0V] = G0vec_border->n_rows_columns();

   int mF0C{0};
   int nF0C{0};
   if (F0cons_border) {
      std::tie(mF0C, nF0C) = F0cons_border->n_rows_columns();
   }

   int mG0C{0};
   int nG0C{0};
   if (G0cons_border) {
      std::tie(mG0C, nG0C) = G0cons_border->n_rows_columns();
   }

   int mA0{0};
   int nA0{0};
   if (A0_border) {
      std::tie(mA0, nA0) = A0_border->n_rows_columns();
   }

   int mC0{0};
   int nC0{0};
   if (C0_border) {
      std::tie(mC0, nC0) = C0_border->n_rows_columns();
   }

   const int n_rows = end_rows - begin_rows;

   assert(n_rows <= mX0);
   assert(0 <= begin_rows && begin_rows <= end_rows);
   if (sc_compute_blockwise_hierarchical) { ;
   } else
      assert(end_rows <= m_result && begin_rows == 0);

   assert(nA0 == nC0);
   assert(nF0V == nG0V);

   if (has_RAC)
      assert(nA0 == nF0V);

   assert(nF0C + mA0 + mC0 + mF0V + mG0V == nX0);
   assert(nF0C == nG0C);

   if (has_RAC)
      assert(n_result == nF0V + mF0C + mG0C);
   else
      assert(n_result >= nF0V);
#endif

   if (is_sparse)
      finalizeInnerSchurComplementContributionSparse(result, X0, A0_border, C0_border, F0vec_border, G0vec_border,
         F0cons_border, G0cons_border,
         begin_rows, end_rows);
   else
      finalizeInnerSchurComplementContributionDense(result, X0, A0_border, C0_border, F0vec_border, G0vec_border,
         F0cons_border, G0cons_border,
         is_sym, begin_rows, end_rows);
}

/* SC and X0 stored in transposed form */
void
DistributedRootLinearSystem::finalizeInnerSchurComplementContributionSparse(AbstractMatrix& SC_, const DenseMatrix& X0,
   const SparseMatrix* A0_border,
   const SparseMatrix* C0_border, const SparseMatrix* F0vec_border, const SparseMatrix* G0vec_border,
   const SparseMatrix* F0cons_border,
   const SparseMatrix* G0cons_border, int begin_rows, int end_rows) {
   assert(F0vec_border);
   assert(G0vec_border);

   auto& SC = dynamic_cast<SparseSymmetricMatrix&>(SC_);

   assert(0 <= begin_rows && begin_rows <= end_rows && end_rows - begin_rows <= X0.n_rows());

   int mA0{0};
   int nA0{0};
   if (A0_border)
      std::tie(mA0, nA0) = A0_border->n_rows_columns();

   int mC0 = C0_border ? C0_border->n_rows() : 0;

   int nF0C{0};
   int mF0C{0};
   if (F0cons_border)
      std::tie(mF0C, nF0C) = F0cons_border->n_rows_columns();

   const int mF0V = F0vec_border->n_rows();
   const int mG0V = G0vec_border->n_rows();

   // multiply each column with B_{outer]}^T and add it to res
   // todo: #pragma omp parallel for schedule(dynamic, 10)
   for (int i = begin_rows; i < end_rows; i++) {
      const double* const col = X0[i];
      if (A0_border)
         A0_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C], i, 0);

      if (C0_border)
         C0_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C + mA0], i, 0);

      if (mF0V > 0)
         F0vec_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C + mA0 + mC0], i, 0);

      if (mG0V > 0)
         G0vec_border->transmultMatSymUpper(1.0, SC, -1.0, &col[nF0C + mA0 + mC0 + mF0V], i, 0);

      if (F0cons_border)
         F0cons_border->multMatSymUpper(1.0, SC, -1.0, &col[0], i, nA0);

      if (G0cons_border)
         G0cons_border->multMatSymUpper(1.0, SC, -1.0, &col[0], i, nA0 + mF0C);
   }
}

/* SC and X0 stored in transposed form */
void
DistributedRootLinearSystem::finalizeInnerSchurComplementContributionDense(AbstractMatrix& SC_, const DenseMatrix& X0,
   const SparseMatrix* A0_border,
   const SparseMatrix* C0_border, const SparseMatrix* F0vec_border, const SparseMatrix* G0vec_border,
   const SparseMatrix* F0cons_border,
   const SparseMatrix* G0cons_border, bool is_sym, int begin_rows, int end_rows) {
   assert(F0vec_border);
   assert(G0vec_border);

   double** SC = is_sym ? dynamic_cast<DenseSymmetricMatrix&>(SC_).Mat() : dynamic_cast<DenseMatrix&>(SC_).Mat();

   const int n_rows = end_rows - begin_rows;
   assert(0 <= begin_rows && begin_rows <= end_rows && n_rows <= X0.n_rows());
   assert(n_rows <= SC_.n_rows());

   int mA0{0};
   int nA0{0};
   if (A0_border)
      std::tie(mA0, nA0) = A0_border->n_rows_columns();

   int mC0 = C0_border ? C0_border->n_rows() : 0;

   int nF0C{0};
   int mF0C{0};
   if (F0cons_border)
      std::tie(mF0C, nF0C) = F0cons_border->n_rows_columns();

   int mF0V = F0vec_border->n_rows();

   // multiply each column with B_{outer]}^T and add it to res
   // todo: #pragma omp parallel for schedule(dynamic, 10)
   for (int i = 0; i < n_rows; ++i) {
      const int row = is_sym ? i + begin_rows : i;
      assert(i < X0.n_rows() && row < SC_.n_rows());

      const double* const col = X0[i];

      if (A0_border)
         A0_border->getStorage().transMult(1.0, &SC[row][0], -1.0, &col[nF0C]);

      if (C0_border)
         C0_border->getStorage().transMult(1.0, &SC[row][0], -1.0, &col[nF0C + mA0]);

      F0vec_border->getStorage().transMult(1.0, &SC[row][0], -1.0, &col[nF0C + mA0 + mC0]);

      G0vec_border->getStorage().transMult(1.0, &SC[row][0], -1.0, &col[nF0C + mA0 + mC0 + mF0V]);

      if (F0cons_border)
         F0cons_border->getStorage().mult(1.0, &SC[row][nA0], -1.0, &col[0]);

      if (G0cons_border)
         G0cons_border->getStorage().mult(1.0, &SC[row][nA0 + mF0C], -1.0, &col[0]);
   }
}

/* compute result += [ SUM_i Bi_{inner}^T Ki^{-1} (Bri - SUM_j Bmodij Xij) ]^T += SUM_i (Bri - SUM_j Xij^T Bmodij^T) Ki^{-1} Bi_{inner}^T */
void DistributedRootLinearSystem::LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool use_local_RAC,
   bool two_link_border, int begin_cols, int end_cols) {
   assert(children.size() == Br.F.children.size());

   /* get contribution to schur_complement from each child */
   /* for a pure 2-link border this has only be done for the last and first process in a communicator - the rest is zero */
   for (size_t it = 0; it < children.size(); it++) {
      BorderLinsys border_child = getChild(Br, it);

      std::vector<BorderMod> Br_mod_border_child;

      for (auto& br_mod : Br_mod_border) {
         auto child = getChild(br_mod, it);
         if (!child.border.isEmpty())
            Br_mod_border_child.push_back(child);
      }

      if (border_child.isEmpty() && Br_mod_border.empty())
         continue;

      const bool sparse_res = false;
      const bool sym_res = false;
      children[it]->addInnerBorderKiInvBrToRes(result, border_child, Br_mod_border_child, use_local_RAC, sparse_res,
         sym_res, begin_cols, end_cols,
         locmy + locmz);
   }

   /* allreduce the result */
   // TODO : optimize -> do not reduce A_0 part ( all zeros... )
   /* in a two-link border we have at most 2 contributions and these are guaranteed disjunct */
   if (iAmDistrib)
      allreduceMatrix(result, false, false, mpiComm);
}

/* compute SUM_i Bli^T X_i = SUM_i Bli^T Ki^-1 (( Bri - sum_j Bmodij Xij ) - Bi_{inner} X0) */
void DistributedRootLinearSystem::LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl,
   BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool sym_res, bool sparse_res, bool use_local_RAC, int begin_cols,
   int end_cols) {
   assert(!is_hierarchy_root);

   /* X0 is still in transposed form */
   assert(children.size() == Br.F.children.size());
   assert(children.size() == Bl.F.children.size());
   assert(!Br.F.is_a(kStringGenDummyMatrix));
   assert(!Bl.F.is_a(kStringGenDummyMatrix));

   /* for every child - add Bi_{outer}^T Ki^-1 (Bi_{outer} - Bi_{inner} X0) */
   for (size_t it = 0; it < children.size(); it++) {
      if (getChild(Bl, it).isEmpty())
         continue;

      BorderLinsys bl_child = getChild(Bl, it);
      BorderLinsys br_child = getChild(Br, it);

      std::vector<BorderMod> border_mod_child;
      for (auto& bm : Br_mod_border) {
         const BorderMod child = getChild(bm, it);

         if (!child.border.isEmpty())
            border_mod_child.push_back(getChild(bm, it));
      }

      children[it]->LniTransMultHierarchyBorder(res, X0, bl_child, br_child, border_mod_child, sparse_res, sym_res,
         use_local_RAC, begin_cols,
         end_cols, locmy + locmz);
   }
}

void DistributedRootLinearSystem::addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0,
   BorderLinsys& border) {
   assert(rhs.children.size() == children.size());
   assert(border.A.children.size() == children.size());

   for (size_t i = 0; i < children.size(); ++i) {
      BorderLinsys child_border = getChild(border, i);
      if (child_border.isEmpty())
         continue;

      children[i]->addBorderX0ToRhs(*rhs.children[i], x0, child_border);
   }

   /* add schur complement part */
   assert(border.A.first);
   assert(border.C.first);

   auto& A0_border = dynamic_cast<SparseMatrix&>(*border.A.first);
   const auto[mA0, nA0] = A0_border.n_rows_columns();

   auto& C0_border = dynamic_cast<SparseMatrix&>(*border.C.first);
   const auto mC0 = C0_border.n_rows();

   assert(border.F.first);
   assert(border.A.last);

   auto& F0vec_border = dynamic_cast<SparseMatrix&>(*border.A.last);
   auto& G0vec_border = dynamic_cast<SparseMatrix&>(*border.C.last);
   const auto mF0V = F0vec_border.n_rows();

   auto& F0cons_border = dynamic_cast<SparseMatrix&>(*border.F.first);
   auto& G0cons_border = dynamic_cast<SparseMatrix&>(*border.G.first);
   const auto[mF0C, nF0C] = F0cons_border.n_rows_columns();

#ifndef NDEBUG
   assert(border.C.last);
   assert(border.G.first);


   assert(rhs.first);
   assert(rhs.first->length() == nF0C + mA0 + mC0 + mF0V + G0vec_border.n_rows());
   assert(x0.length() == nA0 + mF0C + G0cons_border.n_rows());
#endif

   auto& rhs0 = dynamic_cast<SimpleVector<double>&>(*rhs.first);

   double* rhs01 = &rhs0[0];
   double* rhs02 = &rhs0[nF0C];
   double* rhs03 = &rhs0[nF0C + mA0];
   double* rhs04 = &rhs0[nF0C + mA0 + mC0];
   double* rhs05 = &rhs0[nF0C + mA0 + mC0 + mF0V];

   const double* x01 = &x0[0];
   const double* x02 = &x0[nA0];
   const double* x03 = &x0[nA0 + mF0C];

   A0_border.getStorage().mult(1.0, rhs02, -1.0, x01);
   C0_border.getStorage().mult(1.0, rhs03, -1.0, x01);
   F0vec_border.getStorage().mult(1.0, rhs04, -1.0, x01);
   G0vec_border.getStorage().mult(1.0, rhs05, -1.0, x01);

   F0cons_border.getStorage().transMult(1.0, rhs01, -1.0, x02);
   G0cons_border.getStorage().transMult(1.0, rhs01, -1.0, x03);
}

void DistributedRootLinearSystem::addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0,
   BorderLinsys& border) {
   assert(rhs.children.size() == children.size());
   assert(border.A.children.size() == children.size());

   for (size_t i = 0; i < children.size(); ++i) {
      BorderLinsys child_border = getChild(border, i);
      if (child_border.isEmpty())
         continue;

      children[i]->addBorderTimesRhsToB0(*rhs.children[i], b0, child_border);
   }

   /* add schur complement part */
   if (PIPS_MPIgetSize(mpiComm) == 0 || PIPS_MPIgetRank(mpiComm) == 0) {
      assert(border.A.first);
      assert(border.C.first);

      auto& A0_border = dynamic_cast<SparseMatrix&>(*border.A.first);
      const auto[mA0, nA0] = A0_border.n_rows_columns();

      auto& C0_border = dynamic_cast<SparseMatrix&>(*border.C.first);
      const auto mC0 = C0_border.n_rows();

      assert(border.F.first);
      assert(border.A.last);
      auto& F0vec_border = dynamic_cast<SparseMatrix&>(*border.A.last);
      const auto mF0V = F0vec_border.n_rows();
      auto& F0cons_border = dynamic_cast<SparseMatrix&>(*border.F.first);
      const auto[mF0C, nF0C]  = F0cons_border.n_rows_columns();

      assert(border.C.last);
      assert(border.G.first);
      auto& G0vec_border = dynamic_cast<SparseMatrix&>(*border.C.last);
      const auto mG0V = G0vec_border.n_rows();
      auto& G0cons_border = dynamic_cast<SparseMatrix&>(*border.G.first);
      const auto mG0C = G0cons_border.n_rows();

      assert(rhs.first);
      assert(rhs.first->length() == nF0C + mA0 + mC0 + mF0V + mG0V);
      assert(b0.length() == nA0 + mF0C + mG0C);

      auto& zi = dynamic_cast<SimpleVector<double>&>(*rhs.first);

      SimpleVector<double> zi1(&zi[0], nF0C);
      SimpleVector<double> zi2(&zi[nF0C], mA0);
      SimpleVector<double> zi3(&zi[nF0C + mA0], mC0);
      SimpleVector<double> zi4(&zi[nF0C + mA0 + mC0], mF0V);
      SimpleVector<double> zi5(&zi[nF0C + mA0 + mC0 + mF0V], mG0V);

      SimpleVector<double> b1(&b0[0], nA0);
      SimpleVector<double> b2(&b0[nA0], mF0C);
      SimpleVector<double> b3(&b0[nA0 + mF0C], mG0C);

      A0_border.transMult(1.0, b1, -1.0, zi2);
      C0_border.transMult(1.0, b1, -1.0, zi3);
      F0vec_border.transMult(1.0, b1, -1.0, zi4);
      G0vec_border.transMult(1.0, b1, -1.0, zi5);

      F0cons_border.mult(1.0, b2, -1.0, zi1);
      G0cons_border.mult(1.0, b3, -1.0, zi1);
   }
}

void DistributedRootLinearSystem::Ltsolve2(DistributedVector<double>& x, SimpleVector<double>& x0, bool) {
   assert(false && "not in use");
   assert(pipsipmpp_options::get_bool_parameter("HIERARCHICAL"));
   assert(children.size() == x.children.size());

   auto& b = dynamic_cast<DistributedVector<double>&>(x);

   for (size_t i = 0; i < children.size(); ++i) {
      children[i]->computeInnerSystemRightHandSide(*b.children[i], x0, true);
      children[i]->solveCompressed(*x.children[i]);
   }
}

void DistributedRootLinearSystem::createChildren() {
   std::unique_ptr<DistributedLinearSystem> child{};
   assert(primal_diagonal && dq && nomegaInv && rhs);

   auto& primal_diagonalst = dynamic_cast<DistributedVector<double>&>(*primal_diagonal);
   auto& dqst = dynamic_cast<DistributedVector<double>&>(*dq);
   auto& nomegaInvst = dynamic_cast<DistributedVector<double>&>(*nomegaInv);
   auto& regPst = dynamic_cast<DistributedVector<double>&>(*primal_regularization_diagonal);
   auto& regDyst = dynamic_cast<DistributedVector<double>&>(*dual_equality_regularization_diagonal);
   auto& regDzst = dynamic_cast<DistributedVector<double>&>(*dual_inequality_regularization_diagonal);
   auto& rhsst = dynamic_cast<DistributedVector<double>&>(*rhs);

   const auto& distributed_factory = dynamic_cast<const DistributedFactory&>(factory);

   for (size_t it = 0; it < data->children.size(); it++) {
      assert(primal_diagonalst.children[it]);

      if (MPI_COMM_NULL == primal_diagonalst.children[it]->mpiComm) {
         child = std::make_unique<DistributedDummyLinearSystem>(distributed_factory, data->children[it]);
      } else {
         assert(data->children[it]);
         if (is_hierarchy_root) {
            assert(data->isHierarchyRoot());
            assert(data->children.size() == 1);
            assert(data->children[0]);
            assert(primal_diagonalst.children.size() == 1 && dqst.children.size() == 1 &&
               nomegaInvst.children.size() == 1 &&
               rhsst.children.size() == 1);
            assert(MPI_COMM_NULL != primal_diagonalst.children[0]->mpiComm);
         }

         if (data->children[it]->children.empty()) {
            child = distributed_factory.make_linear_system_leaf(data->children[it], primal_diagonalst.children[it],
               dqst.children[it],
               nomegaInvst.children[it], regPst.children[it], regDyst.children[it], regDzst.children[it],
               rhsst.children[it]);
         } else {
            assert(data->children[it]);
            child = distributed_factory.make_linear_system_root(data->children[it], primal_diagonalst.children[it],
               dqst.children[it],
               nomegaInvst.children[it], regPst.children[it], regDyst.children[it], regDzst.children[it],
               rhsst.children[it]);
         }
      }
      assert(child);
      AddChild(std::move(child));
   }
}

void DistributedRootLinearSystem::deleteChildren() {
   children.clear();
}

void DistributedRootLinearSystem::initProperChildrenRange() {
   assert(!children.empty());

   int childStart = -1;
   int childEnd = -1;
   for (size_t it = 0; it < children.size(); it++) {
      if (childEnd != -1)
         assert(children[it]->isDummy());

      if (children[it]->isDummy()) {
         // end of range?
         if (childStart != -1 && childEnd == -1)
            childEnd = int(it);

         continue;
      }

      // start of range?
      if (childStart == -1)
         childStart = int(it);
   }

   assert(childStart >= 0);

   if (childEnd == -1) {
      assert(!children[children.size() - 1]->isDummy());
      childEnd = int(children.size());
   }

   assert(childStart < childEnd && childEnd <= int(children.size()));

   childrenProperStart = childStart;
   childrenProperEnd = childEnd;
}

void DistributedRootLinearSystem::put_primal_diagonal() {
   assert(primal_diagonal);
   const auto& primal_diagonal_stoch = dynamic_cast<const DistributedVector<double>&>(*primal_diagonal);
   assert(children.size() == primal_diagonal_stoch.children.size());

   xDiag = primal_diagonal_stoch.first.get();

   for (auto& it : children)
      it->put_primal_diagonal();
}

void DistributedRootLinearSystem::clear_dual_equality_diagonal() {
   /* we don't have to do anything for our own kkt since the kkt will be reset completely for each iteration */
   for (auto& child : children) {
      child->clear_dual_equality_diagonal();
   }
}

void DistributedRootLinearSystem::put_dual_inequalites_diagonal() {
   assert(nomegaInv);
   const auto& nomegaInv_stoch = dynamic_cast<const DistributedVector<double>&>(*nomegaInv);
   assert(children.size() == nomegaInv_stoch.children.size());

   //kkt->atPutDiagonal( locnx+locmy, *zdiag.first );
   zDiag = nomegaInv_stoch.first.get();
   zDiagLinkCons = nomegaInv_stoch.last.get();

   for (auto& it : children)
      it->put_dual_inequalites_diagonal();
}

void DistributedRootLinearSystem::put_barrier_parameter(double barrier) {
   this->barrier_parameter_current_iterate = barrier;

   for (auto& child : children)
      child->put_barrier_parameter(barrier);
}

void DistributedRootLinearSystem::AddChild(std::unique_ptr<DistributedLinearSystem> child) {
   children.push_back(std::move(child));
}

///////////////////////////////////////////////////////////
// ATOMS of FACTOR 2
//////////////////////////////////////////////////////////
/* Atoms methods of FACTOR2 for a non-leaf linear system */
void DistributedRootLinearSystem::initializeKKT() {
   if (hasSparseKkt)
      dynamic_cast<SparseSymmetricMatrix*>(kkt.get())->symPutZeroes();
   else {
      auto* kktd = dynamic_cast<DenseSymmetricMatrix*>(kkt.get());
      myAtPutZeros(kktd);
   }
}

void DistributedRootLinearSystem::reduceKKT() {
   if (usePrecondDist)
      reduceKKTdist();
   else if (hasSparseKkt)
      reduceKKTsparse();
   else
      reduceKKTdense();
}


/* collects (reduces) lower left part of dense global symmetric Schur complement */
void DistributedRootLinearSystem::reduceKKTdense() {
   auto* const kktd = dynamic_cast<DenseSymmetricMatrix*>(kkt.get());

   // parallel communication
   if (iAmDistrib) {
      if (locnx > 0)
         submatrixAllReduceDiagLower(kktd, 0, locnx, mpiComm);

      if (locmyl > 0 || locmzl > 0) {
         const int locNxMy = locnx + locmy;
         assert(kktd->size() == locnx + locmy + locmyl + locmzl);

         // reduce lower left part
         if (locnx > 0)
            submatrixAllReduceFull(kktd, locNxMy, 0, locmyl + locmzl, locnx, mpiComm);

         // reduce lower diagonal linking part
         submatrixAllReduceDiagLower(kktd, locNxMy, locmyl + locmzl, mpiComm);
      }
   }
}


// collects sparse global Schur complement
void DistributedRootLinearSystem::reduceKKTsparse() {
   if (!iAmDistrib)
      return;
   assert(kkt);

   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;
   const int nnzKkt = krowKkt[sizeKkt];

   assert(kkts.size() == sizeKkt);
   assert(!kkts.is_lower());

   if (allreduce_kkt)
      reduceToAllProcs(nnzKkt, MKkt);
   else
      reduceToProc0(nnzKkt, MKkt);
}

#define CHUNK_SIZE (1024*1024*64) //doubles = 128 MBytes (maximum)

void DistributedRootLinearSystem::reduceToAllProcs(int size, double* values) {
   assert(values && values != sparseKktBuffer);
   assert(size > 0);

   if (sparseKktBuffer == nullptr)
      sparseKktBuffer = new double[CHUNK_SIZE];

   const int reps = size / CHUNK_SIZE;
   const int res = size - CHUNK_SIZE * reps;
   assert(res >= 0 && res < CHUNK_SIZE);

   for (int i = 0; i < reps; i++) {
      double* const start = &values[i * CHUNK_SIZE];
      MPI_Allreduce(start, sparseKktBuffer, CHUNK_SIZE, MPI_DOUBLE, MPI_SUM, mpiComm);

      memcpy(start, sparseKktBuffer, size_t(CHUNK_SIZE) * sizeof(double));
   }

   if (res > 0) {
      double* const start = &values[reps * CHUNK_SIZE];
      MPI_Allreduce(start, sparseKktBuffer, res, MPI_DOUBLE, MPI_SUM, mpiComm);

      memcpy(start, sparseKktBuffer, size_t(res) * sizeof(double));
   }
}

#define CHUNK_SIZE (1024*1024*64) //doubles = 128 MBytes (maximum)

void DistributedRootLinearSystem::reduceToProc0(int size, double* values) {
   assert(values && values != sparseKktBuffer);
   assert(size > 0);

   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);

   if (myRank == 0 && sparseKktBuffer == nullptr)
      sparseKktBuffer = new double[CHUNK_SIZE];

   const int reps = size / CHUNK_SIZE;
   const int res = size - CHUNK_SIZE * reps;
   assert(res >= 0 && res < CHUNK_SIZE);

   for (int i = 0; i < reps; i++) {
      double* const start = &values[i * CHUNK_SIZE];
      MPI_Reduce(start, sparseKktBuffer, CHUNK_SIZE, MPI_DOUBLE, MPI_SUM, 0, mpiComm);

      if (myRank == 0)
         memcpy(start, sparseKktBuffer, size_t(CHUNK_SIZE) * sizeof(double));
   }

   if (res > 0) {
      double* const start = &values[reps * CHUNK_SIZE];
      MPI_Reduce(start, sparseKktBuffer, res, MPI_DOUBLE, MPI_SUM, 0, mpiComm);

      if (myRank == 0)
         memcpy(start, sparseKktBuffer, size_t(res) * sizeof(double));
   }
}


void DistributedRootLinearSystem::registerMatrixEntryTripletMPI() {
   assert(MatrixEntryTriplet_mpi == MPI_DATATYPE_NULL);

   const int nitems = 3;
   int blocklengths[3] = {1, 1, 1};
   MPI_Datatype Types[3] = {MPI_DOUBLE, MPI_INT, MPI_INT};
   MPI_Aint offsets[3];

   offsets[0] = offsetof(MatrixEntryTriplet, val);
   offsets[1] = offsetof(MatrixEntryTriplet, row);
   offsets[2] = offsetof(MatrixEntryTriplet, col);

   MPI_Type_create_struct(nitems, blocklengths, offsets, Types, &MatrixEntryTriplet_mpi);
   MPI_Type_commit(&MatrixEntryTriplet_mpi);
}

void DistributedRootLinearSystem::syncKKTdistLocalEntries() {
   if (!iAmDistrib)
      return;

   assert(kkt && hasSparseKkt);

   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();

   const int childStart = childrenProperStart;
   const int childEnd = childrenProperEnd;
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   int size;
   MPI_Comm_size(mpiComm, &size);

   assert(size > 1);

   // MPI matrix entries triplet not registered yet?
   if (MatrixEntryTriplet_mpi == MPI_DATATYPE_NULL)
      registerMatrixEntryTripletMPI();

   // pack the entries that will be send below
   std::vector<MatrixEntryTriplet> prevEntries = this->packKKTdistOutOfRangeEntries(childStart, childEnd);
   std::vector<MatrixEntryTriplet> myEntries(0);

   assert(!prevEntries.empty() && prevEntries[0].row == -1 && prevEntries[0].col == -1);

   // odd processes send first (one process back)
   if (myRank % 2 != 0) {
      this->sendKKTdistLocalEntries(prevEntries);
   }

   // even processes (except last) receive first
   if (myRank % 2 == 0 && myRank != size - 1) {
      assert(myEntries.empty());
      myEntries = this->receiveKKTdistLocalEntries();
   }

   // even processes (except first) send
   if (myRank % 2 == 0 && myRank > 0) {
      this->sendKKTdistLocalEntries(prevEntries);
   }

   // odd processes (except last) receive
   if (myRank % 2 != 0 && myRank != size - 1) {
      assert(myEntries.empty());
      myEntries = this->receiveKKTdistLocalEntries();
   }

   assert(!myEntries.empty() || myRank == size - 1);

   int lastRow = 0;
   int lastC = -1;

#ifndef NDEBUG
   const std::vector<bool>& rowIsLocal = data->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = data->getSCrowMarkerMyLocal();
#endif

   // finally, put received data into Schur complement matrix
   for (size_t i = 1; i < myEntries.size(); i++) {
      const double val = myEntries[i].val;
      const int row = myEntries[i].row;
      const int col = myEntries[i].col;

      assert(myRank != size - 1);
      assert(row >= 0 && row < locnx + locmy + locmyl + locmzl);
      assert(col >= row && col < locnx + locmy + locmyl + locmzl);
      assert(val == val); // catch NaNs

      assert(rowIsMyLocal[row] || (rowIsMyLocal[col] && !rowIsLocal[row]));

      int c;

      // continue from last position?
      if (row == lastRow) {
         for (c = lastC + 1; c < krowKkt[row + 1]; c++) {
            const int colKkt = jColKkt[c];

            if (colKkt == col) {
               MKkt[c] += val;
               break;
            }
         }

         // found the correct entry in last row?
         if (c != krowKkt[row + 1]) {
            assert(c < krowKkt[row + 1]);
            lastRow = row;
            lastC = c;
            continue;
         }
      }

      c = krowKkt[row];
      assert(col >= jColKkt[c]);

      for (; c < krowKkt[row + 1]; c++) {
         const int colKkt = jColKkt[c];

         if (colKkt == col) {
            MKkt[c] += val;
            break;
         }
      }

      assert(c != krowKkt[row + 1]);

      lastRow = row;
      lastC = c;
   }
}


std::vector<DistributedRootLinearSystem::MatrixEntryTriplet>
DistributedRootLinearSystem::receiveKKTdistLocalEntries() const {
   assert(kkt && hasSparseKkt);
   assert(MatrixEntryTriplet_mpi != MPI_DATATYPE_NULL);

   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   int size;
   MPI_Comm_size(mpiComm, &size);
   const int nextRank = myRank + 1;
   assert(nextRank < size);

   // receive data from next process

   MPI_Status status;
   MPI_Probe(nextRank, 0, mpiComm, &status);

   int nEntries;
   MPI_Get_count(&status, MatrixEntryTriplet_mpi, &nEntries);

   assert(nEntries >= 1);

   std::vector<MatrixEntryTriplet> entries(nEntries);

   MPI_Recv((void*) &entries[0], nEntries, MatrixEntryTriplet_mpi, nextRank, 0, mpiComm, MPI_STATUS_IGNORE);

   assert(entries[0].row == -1 && entries[0].col == -1); // dummy check

   PIPSdebugMessage("myRank=%d received %d \n", myRank, nEntries);

   return entries;
}


void DistributedRootLinearSystem::sendKKTdistLocalEntries(const std::vector<MatrixEntryTriplet>& prevEntries) const {
   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   const int prevRank = myRank - 1;
   const int nEntries = int(prevEntries.size());

   assert(myRank >= 0);
   assert(nEntries > 0);
   assert(MatrixEntryTriplet_mpi != MPI_DATATYPE_NULL);

   PIPSdebugMessage("myRank=%d sends %d \n", myRank, nEntries);
   MPI_Send(&prevEntries[0], nEntries, MatrixEntryTriplet_mpi, prevRank, 0, mpiComm);
}

std::vector<DistributedRootLinearSystem::MatrixEntryTriplet>
DistributedRootLinearSystem::packKKTdistOutOfRangeEntries(int childStart, int) const {
   assert(kkt && hasSparseKkt);

   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);

   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);
   const std::vector<bool>& rowIsLocal = data->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = data->getSCrowMarkerMyLocal();
   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;

   std::vector<MatrixEntryTriplet> packedEntries(0);

   // add dummy value
   const MatrixEntryTriplet entry_zero = {-1.0, -1, -1};
   packedEntries.push_back(entry_zero);

   if (childStart > 0) {
      assert(myRank > 0);

      // pack data
      for (int r = 0; r < sizeKkt; r++) {
         const bool rIsLocal = rowIsLocal[r];

         if (rIsLocal && !rowIsMyLocal[r]) {
            for (int c = krowKkt[r]; c < krowKkt[r + 1]; c++) {
               const int col = jColKkt[c];

               if (!rowIsMyLocal[col]) {
                  const double val = MKkt[c];

                  if (PIPSisZero(val))
                     continue;

                  const MatrixEntryTriplet entry = {val, r, col};
                  packedEntries.push_back(entry);
               }
            }
         }

         if (rIsLocal)
            continue;

         for (int c = krowKkt[r]; c < krowKkt[r + 1]; c++) {
            const int col = jColKkt[c];

            if (rowIsLocal[col] && !rowIsMyLocal[col]) {
               const double val = MKkt[c];

               if (PIPSisZero(val))
                  continue;

               const MatrixEntryTriplet entry = {val, r, col};
               packedEntries.push_back(entry);
            }
         }
      }
   }

   return packedEntries;
}


void DistributedRootLinearSystem::reduceKKTdist() {
   assert(iAmDistrib);
   assert(kkt);

   const std::vector<bool>& rowIsLocal = data->getSCrowMarkerLocal();
   const std::vector<bool>& rowIsMyLocal = data->getSCrowMarkerMyLocal();

   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);

   int* const krowKkt = kkts.krowM();
   int* const jColKkt = kkts.jcolM();
   double* const MKkt = kkts.M();
   const int sizeKkt = locnx + locmy + locmyl + locmzl;
   int nnzDistMyLocal = 0;
   int nnzDistShared = 0;
   int nnzDistLocal;

   std::vector<int> rowSizeMyLocal(sizeKkt, 0);
   std::vector<int> rowSizeShared(sizeKkt, 0);
   std::vector<int> rowSizeLocal(sizeKkt, 0);
   std::vector<int> rowIndexMyLocal(0);
   std::vector<int> colIndexMyLocal(0);

   assert(int(rowIsLocal.size()) == sizeKkt);

   // add up locally owned entries
   syncKKTdistLocalEntries();

   // add B_0, F_0, G_0 and diagonals (all scattered)
   finalizeKKTdist();

   precondSC.updateDiagDomBound();
   precondSC.unmarkDominatedSCdistLocals(*data, kkts);

   // compute row lengths
   for (int r = 0; r < sizeKkt; r++) {
      if (rowIsMyLocal[r]) {
         for (int c = krowKkt[r]; c < krowKkt[r + 1]; c++) {
            const int col = jColKkt[c];

            if (col < 0)
               continue;

            nnzDistMyLocal++;
            rowSizeMyLocal[r]++;
            rowIndexMyLocal.push_back(r);
            colIndexMyLocal.push_back(col);
         }

         continue;
      }

      const bool rIsLocal = rowIsLocal[r];

      for (int c = krowKkt[r]; c < krowKkt[r + 1]; c++) {
         const int col = jColKkt[c];

         if (col < 0) {
            assert(-col - 1 >= 0 && -col - 1 < sizeKkt);
            assert((!(!rIsLocal && !rowIsLocal[-col - 1])));
            continue;
         }

         // is (r, col) a shared entry?
         if (!rIsLocal && !rowIsLocal[col]) {
            nnzDistShared++;
            rowSizeShared[r]++;
            assert(!rowIsMyLocal[col]);
         }

         if (rowIsMyLocal[col]) {
            nnzDistMyLocal++;
            rowSizeMyLocal[r]++;
            rowIndexMyLocal.push_back(r);
            colIndexMyLocal.push_back(col);
         }
      }
   }

   assert(int(rowIndexMyLocal.size()) == nnzDistMyLocal);

   // sum up local sizes
   MPI_Allreduce(&nnzDistMyLocal, &nnzDistLocal, 1, MPI_INT, MPI_SUM, mpiComm);
   MPI_Allreduce(&rowSizeMyLocal[0], &rowSizeLocal[0], sizeKkt, MPI_INT, MPI_SUM, mpiComm);

#ifndef NDEBUG
   {
      int nnzDistSharedMax;
      std::vector<int> rowSizeSharedMax(sizeKkt, 0);

      MPI_Allreduce(&nnzDistShared, &nnzDistSharedMax, 1, MPI_INT, MPI_MAX, mpiComm);
      MPI_Allreduce(&rowSizeShared[0], &rowSizeSharedMax[0], sizeKkt, MPI_INT, MPI_MAX, mpiComm);

      assert(nnzDistSharedMax == nnzDistShared);
      for (int i = 0; i < sizeKkt; i++)
         assert(rowSizeShared[i] == rowSizeSharedMax[i]);
   }
#endif

   int localGatheredMyStart;
   int localGatheredMyEnd;

   std::vector<int> rowIndexGathered = PIPSallgathervInt(rowIndexMyLocal, mpiComm);
   std::vector<int> colIndexGathered = PIPSallgathervInt(colIndexMyLocal, mpiComm, localGatheredMyStart,
      localGatheredMyEnd);

#ifndef NDEBUG
   assert(int(rowIndexGathered.size()) == nnzDistLocal);
   assert(int(colIndexGathered.size()) == nnzDistLocal);
   assert(localGatheredMyEnd - localGatheredMyStart == nnzDistMyLocal);

   for (int i = 0; i < int(rowIndexMyLocal.size()); i++) {
      assert(rowIndexMyLocal[i] == rowIndexGathered[i + localGatheredMyStart]);
      assert(colIndexMyLocal[i] == colIndexGathered[i + localGatheredMyStart]);
   }
#endif

   const int nnzDist = nnzDistLocal + nnzDistShared;

   assert(!kktDist || !kktDist->is_lower());

   delete kktDist;
   kktDist = new SparseSymmetricMatrix(sizeKkt, nnzDist, false);

   int* const krowDist = kktDist->krowM();
   int* const jColDist = kktDist->jcolM();
   double* const MDist = kktDist->M();

   assert(krowDist[0] == 0);
   assert(sizeKkt > 0 && krowDist[1] == 0);

   memset(MDist, 0, nnzDist * sizeof(double));

   for (int r = 1; r < sizeKkt; r++)
      krowDist[r + 1] = krowDist[r] + rowSizeLocal[r - 1] + rowSizeShared[r - 1];

   // fill in global and locally owned positions and values
   for (int r = 0; r < sizeKkt; r++) {
      if (rowIsMyLocal[r]) {
         for (int c = krowKkt[r]; c < krowKkt[r + 1]; c++) {
            assert(krowDist[r + 1] < nnzDist);

            const int col = jColKkt[c];

            if (col < 0)
               continue;

            const double val = MKkt[c];

            MDist[krowDist[r + 1]] = val;
            jColDist[krowDist[r + 1]++] = col;
         }

         continue;
      }

      for (int c = krowKkt[r]; c < krowKkt[r + 1]; c++) {
         const int col = jColKkt[c];

         if (col < 0) {
            assert(-col - 1 >= 0 && -col - 1 < sizeKkt);
            assert(!(!rowIsLocal[r] && !rowIsLocal[-col - 1]));
            continue;
         }

         // is (r, col) a shared entry or locally owned?
         if ((!rowIsLocal[r] && !rowIsLocal[col]) || rowIsMyLocal[col]) {
            assert(krowDist[r + 1] < nnzDist);

            const double val = MKkt[c];

            MDist[krowDist[r + 1]] = val;
            jColDist[krowDist[r + 1]++] = col;
         }
      }
   }

   precondSC.resetSCdistEntries(kkts);

   // fill in gathered local pairs not inserted yet
   for (int i = 0; i < nnzDistLocal; i++) {
      const int row = rowIndexGathered[i];
      const int col = colIndexGathered[i];

      assert(row >= 0 && row < sizeKkt);
      assert(col >= row && col < sizeKkt);

      // pair already added?
      if (i >= localGatheredMyStart && i < localGatheredMyEnd)
         continue;

      assert(krowDist[row + 1] < nnzDist);
      assert(MDist[krowDist[row + 1]] == 0.0);
      jColDist[krowDist[row + 1]++] = col;
   }

#ifndef NDEBUG
   assert(krowDist[0] == 0);
   assert(krowDist[sizeKkt] == nnzDist);

   for (int r = 0; r < sizeKkt; r++) {
      assert(krowDist[r + 1] == krowDist[r] + rowSizeLocal[r] + rowSizeShared[r]);
      assert(krowDist[r + 1] >= krowDist[r]);
   }
#endif

   kktDist->getStorage().sortCols();

   assert(kktDist->getStorage().isValid());

   if (allreduce_kkt)
      reduceToAllProcs(nnzDist, MDist);
   else
      reduceToProc0(nnzDist, MDist);

   assert(kktDist->getStorage().isValid());
   assert(kktDist->getStorage().isSorted());
}

void DistributedRootLinearSystem::factorizeKKT() {
   //stochNode->resMon.recFactTmLocal_start();
#ifdef TIMING
   MPI_Barrier(mpiComm);
   extern double g_iterNumber;
   double st=MPI_Wtime();
#endif
   if (is_hierarchy_root)
      assert(!usePrecondDist);

   if (usePrecondDist) {
      const int myRank = PIPS_MPIgetRank(mpiComm);

      assert(kktDist);

      if (allreduce_kkt || myRank == 0)
         precondSC.getSparsifiedSC_fortran(*data, *kktDist);

      // todo do that properly
      precondSC.updateStats();

#if 0
      {
         ofstream myfile;
         int mype; MPI_Comm_rank(mpiComm, &mype);

         if( mype == 0 )
         {
            printf("\n\n ...WRITE OUT kktDist! \n\n");
            myfile.open("../ADist.txt");
            int* ia = kktDist->krowM(); int* ja = kktDist->jcolM(); double* a = kktDist->M();

            for( int i = 0; i < kktDist->size(); i++ )
               for( int k = ia[i]; k < ia[i + 1]; k++ )
                  myfile << i << '\t' << ja[k - 1] << '\t' << a[k - 1] << endl;

            myfile.close();
         }

         MPI_Barrier(mpiComm);
         printf("...exiting (root) \n");
         exit(1);
      }
#endif
      if (apply_regularization) {
         assert(false && "TODO: implement");
         solver->matrixRebuild(*kktDist);
      } else {
         solver->matrixRebuild(*kktDist);
      }
   } else {
      if (apply_regularization) {
         factorize_with_correct_inertia();
      } else {
         solver->matrixChanged();
      }
   }

   //stochNode->resMon.recFactTmLocal_stop();
#ifdef TIMING
   st = MPI_Wtime()-st;
   MPI_Barrier(mpiComm);
   int mype; MPI_Comm_rank(mpiComm, &mype);
   // note, this will include noop scalapack processors
   if( (mype/512)*512==mype )
     printf("  rank %d 1stSTAGE FACT %g SEC ITER %d\n", mype, st, (int)g_iterNumber);
#endif
}

//faster than DenseSymmetricMatrix::atPutZeros
void
DistributedRootLinearSystem::myAtPutZeros(DenseSymmetricMatrix* mat, int row, int col, int rowExtent, int colExtent) {
   assert(row >= 0 && row + rowExtent <= mat->size());
   assert(col >= 0 && col + colExtent <= mat->size());

   double** M = mat->getStorage().M;

   for (int j = col; j < col + colExtent; j++) {
      M[row][j] = 0.0;
   }

   unsigned long nToCopy = colExtent * sizeof(double);

   for (int i = row + 1; i < row + rowExtent; i++) {
      memcpy(M[i] + col, M[row] + col, nToCopy);
   }
}

void DistributedRootLinearSystem::myAtPutZeros(DenseSymmetricMatrix* mat) {
   int n = static_cast<int>(mat->size());
   myAtPutZeros(mat, 0, 0, n, n);
}

void DistributedRootLinearSystem::addTermToSchurCompl(size_t childindex, bool use_local_RAC) {
   assert(childindex < data->children.size());

   if (computeBlockwiseSC) {
      const int n_empty_rows_border = use_local_RAC ? locmy : locnx + locmy;
      children[childindex]->addTermToSchurComplBlocked(hasSparseKkt, *kkt, use_local_RAC, n_empty_rows_border);
   } else {
      if (hasSparseKkt) {
         auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);

         children[childindex]->addTermToSparseSchurCompl(kkts);
      } else {
         auto& kktd = dynamic_cast<DenseSymmetricMatrix&>(*kkt);
         children[childindex]->addTermToDenseSchurCompl(kktd);
      }
   }
}

void DistributedRootLinearSystem::submatrixAllReduce(DenseSymmetricMatrix* A, int startRow, int startCol, int nRows,
   int nCols, MPI_Comm comm) {
   double** M = A->mStorage->M;
   int n = A->mStorage->n;

   assert(nRows > 0);
   assert(nCols > 0);
   assert(startRow >= 0);
   assert(startCol >= 0);

   int endRow = startRow + nRows;
   int endCol = startCol + nCols;

   assert(n >= endRow);
   assert(n >= endCol);

   int chunk_size = (CHUNK_SIZE / n) * n;
   chunk_size = std::min(chunk_size, n * nRows);

   auto* chunk = new double[chunk_size];

   int rows_in_chunk = chunk_size / n;

   int iRow = startRow;

   // main loop
   do {

      if (iRow + rows_in_chunk > endRow)
         rows_in_chunk = endRow - iRow;

      assert(rows_in_chunk > 0);
#ifndef NDEBUG
      const int iErr = MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk * n, MPI_DOUBLE, MPI_SUM, comm);
      assert(iErr == MPI_SUCCESS);
#else
      MPI_Allreduce(&M[iRow][0], chunk, rows_in_chunk*n, MPI_DOUBLE, MPI_SUM, comm);
#endif

      int shift = 0;

      // copy into M
      for (int i = iRow; i < iRow + rows_in_chunk; i++) {
         for (int j = startCol; j < endCol; j++)
            M[i][j] = chunk[shift + j];

         // shift one row forward
         shift += n;
      }
      iRow += rows_in_chunk;

   } while (iRow < endRow);

   delete[] chunk;
}

// TODO: move all this to the respective matrix and storages.......
void DistributedRootLinearSystem::allreduceMatrix(AbstractMatrix& mat, bool is_sparse, bool is_sym, MPI_Comm comm) {
   const auto[m, n] = mat.n_rows_columns();

   if (is_sparse) {
      if (is_sym) {
         auto& matsp = dynamic_cast<SparseSymmetricMatrix&>(mat);

         int* const krowKkt = matsp.krowM();
         double* const MKkt = matsp.M();
         const int nnzKkt = krowKkt[m];

         assert(!matsp.is_lower());

         reduceToAllProcs(nnzKkt, MKkt);
      } else {
         auto& matsp = dynamic_cast<SparseMatrix&>(mat);

         int* const krowKkt = matsp.krowM();
         double* const MKkt = matsp.M();
         const int nnzKkt = krowKkt[m];

         reduceToAllProcs(nnzKkt, MKkt);
      }
   } else {
      // TODO : these seem to be not proper handling of the symmetric dense schur complement - need ot check maths first
      if (is_sym)
         submatrixAllReduceFull(&dynamic_cast<DenseSymmetricMatrix&>(mat), 0, 0, m, n, comm);
      else
         submatrixAllReduceFull(&dynamic_cast<DenseMatrix&>(mat), 0, 0, m, n, comm);
   }
}

void DistributedRootLinearSystem::submatrixAllReduceFull(DenseSymmetricMatrix* A, int startRow, int startCol, int nRows,
   int nCols, MPI_Comm comm) {
   assert(A->n_rows() >= startRow + nRows);
   assert(A->n_columns() >= startCol + nCols);

   submatrixAllReduceFull(A->mStorage->M, startRow, startCol, nRows, nCols, comm);
}

void
DistributedRootLinearSystem::submatrixAllReduceFull(DenseMatrix* A, int startRow, int startCol, int nRows, int nCols,
   MPI_Comm comm) {
   assert(A->n_rows() >= startRow + nRows);
   assert(A->n_columns() >= startCol + nCols);

   submatrixAllReduceFull(A->mStorage->M, startRow, startCol, nRows, nCols, comm);
}

void DistributedRootLinearSystem::submatrixAllReduceFull(double** A, int startRow, int startCol, int nRows, int nCols,
   MPI_Comm comm) {
   assert(nRows >= 0);
   assert(nCols >= 0);
   if (nRows == 0 || nCols == 0)
      return;
   assert(startRow >= 0);
   assert(startCol >= 0);

   const int endRow = startRow + nRows;
   const int buffersize = nRows * nCols;

   auto* const bufferSend = new double[buffersize];
   auto* const bufferRecv = new double[buffersize];

   // copy into send buffer
   int counter = 0;
   const size_t nColBytes = nCols * sizeof(double);

   for (int r = startRow; r < endRow; r++) {
      memcpy(&bufferSend[counter], &A[r][startCol], nColBytes);
      counter += nCols;
   }

   assert(counter == buffersize);

#ifndef NDEBUG
   const int iErr = MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
   assert(iErr == MPI_SUCCESS);
#else
   MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
#endif

   // copy back
   counter = 0;
   for (int r = startRow; r < endRow; r++) {
      memcpy(&A[r][startCol], &bufferRecv[counter], nColBytes);
      counter += nCols;
   }

   assert(counter == buffersize);

   delete[] bufferRecv;
   delete[] bufferSend;
}


void DistributedRootLinearSystem::submatrixAllReduceDiagLower(DenseSymmetricMatrix* A, int substart, int subsize,
   MPI_Comm comm) {
   double** const M = A->mStorage->M;

   assert(subsize >= 0);
   assert(substart >= 0);

   if (subsize == 0)
      return;

   const int subend = substart + subsize;
   assert(A->mStorage->n >= subend);

   // number of elements in lower matrix triangle (including diagonal)
   const int buffersize = (subsize * subsize + subsize) / 2;
   assert(buffersize > 0);

   auto* const bufferSend = new double[buffersize];
   auto* const bufferRecv = new double[buffersize];

   int counter = 0;

   for (int i = substart; i < subend; i++)
      for (int j = substart; j <= i; j++) {
         assert(counter < buffersize);
         bufferSend[counter++] = M[i][j];
      }

   assert(counter == buffersize);

#ifndef NDEBUG
   const int iErr = MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
   assert(iErr == MPI_SUCCESS);
#else
   MPI_Allreduce(bufferSend, bufferRecv, buffersize, MPI_DOUBLE, MPI_SUM, comm);
#endif

   counter = 0;
   for (int i = substart; i < subend; i++)
      for (int j = substart; j <= i; j++) {
         assert(counter < buffersize);
         M[i][j] = bufferRecv[counter++];
      }

   delete[] bufferSend;
   delete[] bufferRecv;
}


#ifdef STOCH_TESTING
void sLinsysRoot::dumpMatrix(int scen, int proc, const char* nameToken, DenseSymmetricMatrix& M)
{
  int n = M.size();
  char szNumber[30];
  string strBuffer="";

  //assert(false);

  int iter = g_iterNumber;

  if(iter!=1 && iter!=5 && iter!=15 && iter!=25 && iter!=35 && iter!=45) return;


  char szFilename[256];
  if(scen==-1)
    sprintf(szFilename, "%s_%d__%d.mat", nameToken, n, iter);
  else 
    sprintf(szFilename, "%s_%03d_%d__%d.mat", nameToken, scen+1, n, iter);
  FILE* file = fopen(szFilename, "w");
  assert(file);
  

  for(int j=0; j<n; j++) {
    for(int i=0; i<n; i++) {
      sprintf(szNumber, "%22.16f ", M[i][j]);
      strBuffer += szNumber;
    }
    strBuffer += "\n";
    
    if(strBuffer.length()>1250000) {
      fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
      strBuffer = "";
    }
  }

  if(strBuffer.length()>0) {
    fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
  }
  
  fclose(file);
}

void sLinsysRoot::dumpRhs(int proc, const char* nameToken,  SimpleVector<double>& rhs)
{
  int n = rhs.length();
  char szNumber[30];
  string strBuffer="";


  int iter = g_iterNumber;
  if(iter!=0 && iter!=2 && iter!=20 && iter!=25 && iter!=55) return;

  char ipmPhase[4];
  if(g_iterNumber-iter>0) strcpy(ipmPhase, "co");
  else strcpy(ipmPhase, "pr");

  char szFilename[256];
  sprintf(szFilename, "%s_%s_%d__%d.mat", nameToken,ipmPhase, n, iter);


  for(int i=0; i<n; i++) {
    sprintf(szNumber, "%22.16f ", rhs[i]);
    strBuffer += szNumber;
  }

  FILE* file = fopen(szFilename, "w");
  assert(file);

  fwrite(strBuffer.c_str(), 1, strBuffer.length(), file);
  fclose(file);
}

#endif
