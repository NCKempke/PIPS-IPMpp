/* PIPS
   Authors: Cosmin Petra
   See license and copyright information in the documentation */

#include "sLinsysRootAug.h"
#include "DistributedQP.hpp"
#include "BorderedSymmetricMatrix.h"
#include "PIPSIPMppOptions.h"
#include <memory>
#include <unistd.h>
#include <cmath>

//#define DUMPKKT
#ifdef DUMPKKT
#include <iostream>
#include <fstream>
#endif
#ifdef STOCH_TESTING
extern double g_iterNumber;
#endif
extern int gInnerBiCGIter;
extern int gInnerBiCGFails;


static void biCGStabCommunicateStatus(int flag, int it) {
   gInnerBiCGIter = it;

   if (flag != 0)
      gInnerBiCGFails++;
}

sLinsysRootAug::sLinsysRootAug(DistributedFactory* factory_, DistributedQP* prob_, bool is_hierarchy_root) : DistributedRootLinearSystem(
   factory_, prob_, is_hierarchy_root) {

   if (!is_hierarchy_root) {
      static bool printed = false;
      if (!printed && PIPS_MPIgetRank() == 0) {
         print_solver_regularization_and_sc_info("sLinsysRootAug");
         printed = true;
      }
   }

   redRhs = std::make_unique<SimpleVector<double>>(locnx + locmy + locmz + locmyl + locmzl);
}

sLinsysRootAug::sLinsysRootAug(DistributedFactory* factory_, DistributedQP* prob_, std::shared_ptr<Vector<double>> dd_,
   std::shared_ptr<Vector<double>> dq_,
   std::shared_ptr<Vector<double>> nomegaInv_, std::shared_ptr<Vector<double>> regP,
   std::shared_ptr<Vector<double>> regDy, std::shared_ptr<Vector<double>> regDz, std::shared_ptr<Vector<double>> rhs_,
   bool create_sub_root_solver)
   : DistributedRootLinearSystem(factory_, prob_, std::move(dd_), std::move(dq_), std::move(nomegaInv_),
   std::move(regP), std::move(regDy), std::move(regDz), std::move(rhs_), create_sub_root_solver) {
   assert(pipsipmpp_options::get_bool_parameter("HIERARCHICAL"));
   assert(computeBlockwiseSC);

   std::cout << "slinsysrootaug : " << locnx << " " << locmy << " " << locmz << " " << locmyl << " "
             << locmzl << std::endl;

   redRhs = std::make_unique<SimpleVector<double>>(locnx + locmy + locmz + locmyl + locmzl);
}

/* Add corner block
 * [ Q0 A0T C0T F0T G0T  ]
 * [ A0  0   0   0   0   ]
 * [ C0  0  Om0  0   0   ]
 * [ F0  0   0   0   0   ]
 * [ G0  0   0   0 OmN+1 ]
 */
void sLinsysRootAug::finalizeKKT() {

   resource_monitor->recFactTmLocal_start();
   resource_monitor->recSchurMultLocal_start();

   if (usePrecondDist) {
      // don't do anything, already done previously
   } else {
      if (hasSparseKkt)
         finalizeKKTsparse();
      else
         finalizeKKTdense();
   }

   resource_monitor->recSchurMultLocal_stop();
   resource_monitor->recFactTmLocal_stop();
}


void sLinsysRootAug::finalizeKKTdist() {
   assert(kkt && hasSparseKkt);

   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);

   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   int mpiCommSize;
   MPI_Comm_size(mpiComm, &mpiCommSize);
   const bool iAmLastRank = (myRank == mpiCommSize - 1);
   const int childStart = childrenProperStart;
   const int childEnd = childrenProperEnd;

#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();

   assert(childStart >= 0 && childStart < childEnd);
   assert(kkts.size() == locnx + locmy + locmyl + locmzl);
   assert(!kkts.is_lower());
   assert(locmyl >= 0 && locmzl >= 0);
   assert(data->getLocalQ().krowM()[locnx] == 0 && "Q currently not supported for dist. sparse kkt");

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////
   if (xDiag && iAmLastRank) {
      // TODO : use addDiag method
      const auto& sxDiag = dynamic_cast<const SimpleVector<double>&>(*xDiag);

      for (int i = 0; i < locnx; i++) {
         const int diagIdx = krowKkt[i];
         assert(jcolKkt[diagIdx] == i);

         MKkt[diagIdx] += sxDiag[i];
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if (locmz > 0 && iAmLastRank) {
      assert(zDiag);
      SymmetricMatrix* CtDCptr = CtDC ? CtDC.get() : nullptr;

      compute_CtDC_and_add_to_Schur_complement(CtDCptr, *zDiag);
      if (!CtDC) {
         CtDC.reset(CtDCptr);
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with At (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if (locmy > 0 && iAmLastRank) {
      const SparseMatrix& At = data->getLocalB().getTranspose(); // yes, B
      const double* MAt = At.M();
      const int* krowAt = At.krowM();

      for (int i = 0; i < locnx; ++i) {
         const int pstart = krowAt[i];
         const int pend = krowAt[i + 1];

         // get start position of sparse kkt block
         const int blockStart = krowKkt[i] + locnx - i;

         assert(blockStart <= krowKkt[i + 1]);

         for (int p = pstart; p < pend; ++p) {
            assert(At.jcolM()[p] < locmy);
            assert(blockStart + (p - pstart) <= krowKkt[i + 1]);
            assert(jcolKkt[blockStart + (p - pstart)] == (locnx + At.jcolM()[p]));

            MKkt[blockStart + (p - pstart)] += MAt[p];
         }
      }
   }

   int local2linksStartEq;
   int local2linksEndEq;
   int local2linksStartIneq;
   int local2linksEndIneq;

   data->getSCrangeMarkersMy(childStart, childEnd, local2linksStartEq, local2linksEndEq, local2linksStartIneq,
      local2linksEndIneq);

   const int n2linksRowsLocalEq = local2linksEndEq - local2linksStartEq;

   PIPSdebugMessage("rank %d FT local columns: %d-%d \n", myRank, local2linksStartEq, local2linksEndEq);
   PIPSdebugMessage("rank %d GT local columns: %d-%d \n", myRank, local2linksStartIneq, local2linksEndIneq);
   PIPSdebugMessage("rank %d FT local columns: %d-%d \n", myRank, local2linksStartEq, local2linksEndEq);

   /////////////////////////////////////////////////////////////
   // update the KKT with Ft
   /////////////////////////////////////////////////////////////
   if (locmyl > 0) {
      const SparseMatrix& Ft = data->getLocalF().getTranspose();

      // add locally owned sparse part of Ft
      addLinkConsBlock0Matrix(Ft, locnx + locmy, 0, local2linksStartEq, local2linksEndEq);

      const int n2linksRowsEq = data->n2linkRowsEq();
      const int bordersizeEq = locmyl - n2linksRowsEq;
      const int borderstartEq = locnx + locmy + n2linksRowsEq;
      if (myRank == 0) {
         PIPSdebugMessage("rank %d FT border columns: %d-%d \n", myRank, borderstartEq, borderstartEq + bordersizeEq);

         // add (shared) border part of Ft
         addLinkConsBlock0Matrix(Ft, locnx + locmy, n2linksRowsLocalEq, borderstartEq, borderstartEq + bordersizeEq);
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Gt and add z diagonal
   /////////////////////////////////////////////////////////////
   if (locmzl > 0) {
      const SparseMatrix& Gt = data->getLocalG().getTranspose();
      const int n2linksRowsIneq = data->n2linkRowsIneq();
      const int bordersizeIneq = locmzl - n2linksRowsIneq;
      const int borderstartIneq = locnx + locmy + locmyl + n2linksRowsIneq;

      // add locally owned sparse part of Gt
      addLinkConsBlock0Matrix(Gt, locnx + locmy + locmyl, n2linksRowsLocalEq, local2linksStartIneq, local2linksEndIneq);

      if (myRank == 0) {
         const int n2linksRowsLocalIneq = local2linksEndIneq - local2linksStartIneq;
         PIPSdebugMessage("rank %d GT border columns: %d-%d\n", myRank, borderstartIneq,
               borderstartIneq + bordersizeIneq);

         // add (shared) border part of Gt
         addLinkConsBlock0Matrix(Gt, locnx + locmy + locmyl, n2linksRowsLocalEq + n2linksRowsLocalIneq, borderstartIneq,
            borderstartIneq + bordersizeIneq);
      }

      assert(zDiagLinkCons);
      const auto& szDiagLinkCons = dynamic_cast<const SimpleVector<double>&>(*zDiagLinkCons);

      assert(local2linksStartIneq >= locnx + locmy + locmyl);
      assert(local2linksEndIneq <= locnx + locmy + locmyl + locmzl);

      const int szDiagLocalStart = local2linksStartIneq - (locnx + locmy + locmyl);
      assert(szDiagLocalStart >= 0);
      assert(szDiagLocalStart < locmzl || (szDiagLocalStart == locmzl && local2linksStartIneq == local2linksEndIneq));

      // add locally owned part of z diagonal
      for (int i = szDiagLocalStart, iKkt = local2linksStartIneq; iKkt < local2linksEndIneq; ++i, ++iKkt) {
         const int idx = krowKkt[iKkt];
         assert(jcolKkt[idx] == iKkt);
         assert(i < locmzl);

         MKkt[idx] += szDiagLinkCons[i];
      }

      if (myRank == 0) {
         const int szDiagBorderStart = borderstartIneq - (locnx + locmy + locmyl);

         assert(szDiagBorderStart >= 0 && szDiagBorderStart <= locmzl);
         assert(szDiagBorderStart + bordersizeIneq == locmzl);

         // add border part of diagonal
         for (int i = szDiagBorderStart, iKkt = borderstartIneq; iKkt < borderstartIneq + bordersizeIneq; ++i, ++iKkt) {
            const int idx = krowKkt[iKkt];
            assert(jcolKkt[idx] == iKkt);
            assert(i < locmzl);

            MKkt[idx] += szDiagLinkCons[i];
         }
      }
   }
}

void sLinsysRootAug::assembleLocalKKT() {
   const bool is_layer_only_twolinks = data->isHierarchySparseTopLayerOnlyTwolinks();
   if (!pipsipmpp_options::get_bool_parameter("HIERARCHICAL"))
      assert(!is_layer_only_twolinks);

   for (size_t c = 0; c < children.size(); ++c) {
      if (children[c]->mpiComm == MPI_COMM_NULL)
         continue;

      children[c]->resource_monitor->recFactTmChildren_start();

      //---------------------------------------------
      addTermToSchurCompl(c, !is_layer_only_twolinks);

      //---------------------------------------------
      children[c]->resource_monitor->recFactTmChildren_stop();
   }
}

/* compute Schur rhs b0 - sum Bi^T Ki^-1 bi for all children */
void sLinsysRootAug::Lsolve(Vector<double>& x) {
   assert(!is_hierarchy_root);

   auto& b = dynamic_cast<DistributedVector<double>&>(x);
   assert(children.size() == b.children.size());

   auto& b0 = dynamic_cast<SimpleVector<double>&>(*b.first);
   assert(!b.last);

   if (iAmDistrib && PIPS_MPIgetRank(mpiComm) > 0)
      b0.setToZero();

   // compute Bi^T Ki^-1 rhs_i and sum it up
   for (size_t it = 0; it < children.size(); it++) {
      children[it]->addLniziLinkCons(b0, *b.children[it], true);

   }
   if (iAmDistrib)
      PIPS_MPIsumArrayInPlace(b0.elements(), b0.length(), mpiComm);

   //dumpRhs(0, "rhs",  b0);
}

/* does Schur Complement solve */
void sLinsysRootAug::Dsolve(Vector<double>& x) {
   /* Ki^-1 bi has already been computed in Lsolve */

   /* children have already computed Li^T\Di\Li\bi in Lsolve() */
   auto& b = dynamic_cast<DistributedVector<double>&>(x);
   auto& b0 = dynamic_cast<SimpleVector<double>&>(*b.first);
   solveReducedLinkCons(b0);
}

void sLinsysRootAug::Ltsolve(Vector<double>& x) {
   auto& b = dynamic_cast<DistributedVector<double>&>(x);
   auto& b0 = dynamic_cast<SimpleVector<double>&>(*b.first);

   //dumpRhs(0, "sol",  b0);
   SimpleVector<double>& z0 = b0; //just another name, for clarity

   for (size_t it = 0; it < children.size(); it++)
      children[it]->Ltsolve2(*b.children[it], z0, true);
}

/* gets called for computing the dense schur complement*/
void sLinsysRootAug::LsolveHierarchyBorder(DenseMatrix& result, BorderLinsys& Br, std::vector<BorderMod>& Br_mod_border,
   bool two_link_border,
   int begin_cols, int end_cols) {
   LsolveHierarchyBorder(result, Br, Br_mod_border, true, two_link_border, begin_cols, end_cols);
}

void
sLinsysRootAug::LtsolveHierarchyBorder(AbstractMatrix& res, const DenseMatrix& X0, BorderLinsys& Bl, BorderLinsys& Br,
   std::vector<BorderMod>& br_mod_border, bool sym_res, bool sparse_res, int begin_cols, int end_cols) {
   if (Bl.isEmpty() || (Br.isEmpty() && br_mod_border.empty()))
      return;

   LtsolveHierarchyBorder(res, X0, Bl, Br, br_mod_border, sym_res, sparse_res, true, begin_cols, end_cols);
}

extern int gLackOfAccuracy;

void sLinsysRootAug::solveReducedLinkCons(SimpleVector<double>& b_vec) {
   assert(locmyl >= 0 && locmzl >= 0);
   assert(locnx + locmy + locmz + locmyl + locmzl == b_vec.length());

   ///////////////////////////////////////////////////////////////////////
   // LOCAL SOLVE WITH SCHUR COMPLEMENT and i-th rhs from buffer
   ///////////////////////////////////////////////////////////////////////
   const SparseMatrix& C = data->getLocalD();

   double* rhs_reduced = redRhs->elements();
   assert(redRhs->length() >= b_vec.length());

   ///////////////////////////////////////////////////////////////////////
   // b = [b1; b2; b3; b4; b5] is a locnx + locmy + locmz + locmyl + locmz vector
   // the new rhs should be
   //           r = [b1-C^T*(zDiagReg)^{-1}*b3; b2; b4; b5]
   ///////////////////////////////////////////////////////////////////////
   double* b = b_vec.elements();

   //copy all elements from b into r except for the the residual values corresponding to z0 = b3
   // copy b1, b2
   std::copy(b, b + locnx + locmy, rhs_reduced);
   // copy b4, b5
   std::copy(b + locnx + locmy + locmz, b + locnx + locmy + locmz + locmyl + locmzl, rhs_reduced + locnx + locmy);
   // copy b3 to the end - used as buffer for reduction computations
   std::copy(b + locnx + locmy, b + locnx + locmy + locmz, rhs_reduced + locnx + locmy + locmyl + locmzl);
   // rhs_reduced now : [ b1; b2; b4; b5; b3]

   // alias to r1 part (no mem allocations)
   SimpleVector<double> rhs1(rhs_reduced, locnx);
   SimpleVector<double> rhs_reduced_b3(rhs_reduced + locnx + locmy + locmyl + locmzl, locmz);

   ///////////////////////////////////////////////////////////////////////
   // compute r1 = b1 - C^T * (zDiag)^{-1} * rhs_reduced_b3
   ///////////////////////////////////////////////////////////////////////
   // if we have C part
   if (locmz > 0) {
      assert(zDiag);
      assert(rhs_reduced_b3.length() == zDiag->length());

      rhs_reduced_b3.componentDiv(*zDiag);
      C.transMult(1.0, rhs1, -1.0, rhs_reduced_b3);
   }

   ///////////////////////////////////////////////////////////////////////
   // rhs_reduced now contains all components -> solve for it
   ///////////////////////////////////////////////////////////////////////

   // we do not need the last locmz elements of r since they were buffer only
   SimpleVector<double> rhs_short(rhs_reduced, locnx + locmy + locmyl + locmzl);

   if (innerSCSolve == 0) {
      // Option 1. - solve with the factors
      solver->solveSynchronized(rhs_short);
   } else if (innerSCSolve == 1) {
      // Option 2 - solve with the factors and perform iter. ref.
      solveWithIterRef(rhs_short);
   } else {
      assert(innerSCSolve == 2);
      // Option 3 - use the factors as preconditioner and apply BiCGStab
      solveWithBiCGStab(rhs_short);
   }

   ///////////////////////////////////////////////////////////////////////
   // rhs_small is now the solution to the reduced system
   // the solution to the augmented system can now be computed as
   //      x = [rhs1; rhs2; zDiag^{-1} * (b3 - C * r1); rhs3; rhs4]
   ///////////////////////////////////////////////////////////////////////

   // copy the solution components and calculate r3
   // copy rhs1 and rhs2
   std::copy(rhs_reduced, rhs_reduced + locnx + locmy, b);
   // compute x3
   if (locmz > 0) {
      SimpleVector<double> b3(b + locnx + locmy, locmz);
      C.mult(1.0, b3, -1.0, rhs1);
      b3.componentDiv(*zDiag);
   }

   // copy rhs3 and rhs4
   std::copy(rhs_reduced + locnx + locmy, rhs_reduced + locnx + locmy + locmyl + locmzl, b + locnx + locmy + locmz);
}

void sLinsysRootAug::solveReducedLinkConsBlocked(DenseMatrix& rhs_mat_transp, int rhs_start, int n_rhs) {
#ifdef TIMING
   t_start = MPI_Wtime();
   troot_total = tchild_total = tcomm_total = 0.0;
#endif

   const int length_rhs = rhs_mat_transp.n_columns();

   assert(locmyl >= 0 && locmzl >= 0);
   assert(locnx + locmy + locmz + locmyl + locmzl == length_rhs);

   const int length_reduced = locnx + locmy + locmyl + locmzl;
   if (reduced_rhss_blocked.size() <= static_cast<unsigned int>(n_rhs * length_reduced))
      reduced_rhss_blocked.resize(n_rhs * length_reduced);

   ///////////////////////////////////////////////////////////////////////
   // LOCAL SOLVE WITH SCHUR COMPLEMENT and set of buffer rhs
   ///////////////////////////////////////////////////////////////////////
   const SparseMatrix& C = data->getLocalD();

#pragma omp parallel for schedule(dynamic, 1)
   for (int rhs_i = rhs_start; rhs_i < rhs_start + n_rhs; ++rhs_i) {
      assert(rhs_i < rhs_mat_transp.n_rows());

      double* rhs_reduced = reduced_rhss_blocked.data() + (rhs_i - rhs_start) * length_reduced;

      SimpleVector<double> b_vec(rhs_mat_transp[rhs_i], length_rhs);

      double* b = b_vec.elements();

      ///////////////////////////////////////////////////////////////////////
      // b = [b1; b2; b3; b4; b5] is a locnx + locmy + locmz + locmyl + locmz vector
      // the new rhs should be
      //           r = [b1-C^T*(zDiag)^{-1}*b3; b2; b4; b5]
      ///////////////////////////////////////////////////////////////////////

      //copy all elements from b into r except for the the residual values corresponding to z0 = b3
      // copy b1, b2
      std::copy(b, b + locnx + locmy, rhs_reduced);

      // copy b4, b5
      std::copy(b + locnx + locmy + locmz, b + locnx + locmy + locmz + locmyl + locmzl, rhs_reduced + locnx + locmy);

      // rhs_reduced now : [ b1; b2; b4; b5 ]

      // alias to r1 part (no mem allocations)
      SimpleVector<double> rhs1(rhs_reduced, locnx);
      SimpleVector<double> b3(b + locnx + locmy, locmz);

      ///////////////////////////////////////////////////////////////////////
      // compute r1 = b1 - C^T * (zDiag)^{-1} * b3
      ///////////////////////////////////////////////////////////////////////
      // if we have C part
      if (locmz > 0) {
         assert(zDiag);
         assert(b3.length() == zDiag->length());
         C.transMultD(1.0, rhs1, -1.0, b3, *zDiag);
      }
   }


   ///////////////////////////////////////////////////////////////////////
   // reduced_rhss_blocked now contains all reduced rhs -> solve for it
   ///////////////////////////////////////////////////////////////////////
   if (innerSCSolve == 0)
      solver->solve(n_rhs, reduced_rhss_blocked.data(), nullptr);
   else
      assert(false && "bicg and iterref not available for blocked solution");

   ///////////////////////////////////////////////////////////////////////
   // rhs_small is now the solution to the reduced system
   // the solution to the augmented system can now be computed as
   //      x = [rhs1; rhs2; zDiag^{-1} * (b3 - C * r1); rhs3; rhs4]
   ///////////////////////////////////////////////////////////////////////


   // copy the solution components and calculate r3
   // copy rhs1 and rhs2
#pragma omp parallel for schedule(dynamic, 1)
   for (int rhs_i = rhs_start; rhs_i < rhs_start + n_rhs; ++rhs_i) {
      double* rhs_reduced = reduced_rhss_blocked.data() + (rhs_i - rhs_start) * length_reduced;

      SimpleVector<double> b_vec(rhs_mat_transp[rhs_i], length_rhs);
      double* b = b_vec.elements();

      std::copy(rhs_reduced, rhs_reduced + locnx + locmy, b);
      // compute x3
      if (locmz > 0) {
         SimpleVector<double> rhs1(rhs_reduced, locnx);
         SimpleVector<double> b3(b + locnx + locmy, locmz);
         C.mult(1.0, b3, -1.0, rhs1);
         b3.componentDiv(*zDiag);
      }

      // copy rhs3 and rhs4
      std::copy(rhs_reduced + locnx + locmy, rhs_reduced + locnx + locmy + locmyl + locmzl, b + locnx + locmy + locmz);
   }
#ifdef TIMING
   // TODO
  if( myRank == 0 && innerSCSolve >= 1 )
    std::cout << "Root - Refin times: child=" << tchild_total << " root=" << troot_total
       << " comm=" << tcomm_total << " total=" << MPI_Wtime()-t_start << "\n";
#endif
}


/** Ht should be either Ft or Gt */
void
sLinsysRootAug::addLinkConsBlock0Matrix(const SparseMatrix& Ht, int nHtOffsetCols, int nKktOffsetCols, int startCol,
   int endCol) {
   assert(startCol >= 0 && startCol <= endCol && nKktOffsetCols >= 0 && nKktOffsetCols <= startCol);

   if (startCol == endCol)
      return;

   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);

   int* const jcolKkt = kkts.jcolM();
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const double* const MHt = Ht.M();
   const int* const krowHt = Ht.krowM();
   const int* const jcolHt = Ht.jcolM();
   const int n0Links = data->getN0LinkVars();

   /* main loop going over all rows of Ht */
   for (int i = 0; i < locnx; ++i) {
      const bool sparseRow = (i >= locnx - n0Links);

      // note: upper left block ignores 0-link variables pattern, since CtC pattern is not implemented
      int pKkt = krowKkt[i] + locnx - i;

      if (!sparseRow)
         pKkt += nKktOffsetCols;

      assert(pKkt <= krowKkt[i + 1]);
      assert(sparseRow || pKkt == krowKkt[i + 1] || jcolKkt[pKkt] <= startCol);

      if (jcolKkt[pKkt] >= endCol) {
#ifndef NDEBUG
         // make sure that there is no entry of Ht in the given range
         int pHt;

         for (pHt = krowHt[i]; pHt < krowHt[i + 1]; pHt++) {
            const int colHt = jcolHt[pHt] + nHtOffsetCols;
            if (colHt >= startCol && colHt < endCol)
               break;
         }

         assert(pHt == krowHt[i + 1]);
#endif
         return;
      }

      bool hit = false;

      // get first in-range entry of Kkt
      for (; pKkt < krowKkt[i + 1]; pKkt++) {
         const int colKkt = jcolKkt[pKkt];
         if (colKkt >= startCol && colKkt < endCol) {
            hit = true;
            break;
         }

         if (colKkt >= endCol)
            break;
      }

      // no entry of Kkt in range?
      if (!hit) {
         assert(startCol == endCol || sparseRow);
         continue;
      }

      assert(pKkt < krowKkt[i + 1]);

      int pHt;
      int colHt = -1;
      hit = false;

      // get first in-range entry of Ht
      for (pHt = krowHt[i]; pHt < krowHt[i + 1]; pHt++) {
         colHt = jcolHt[pHt] + nHtOffsetCols;
         if (colHt >= startCol && colHt < endCol) {
            hit = true;
            break;
         }

         if (colHt >= endCol)
            break;
      }

      // no entry of Ht in range?
      if (!hit)
         continue;

      assert(colHt >= startCol && colHt < endCol);

      // add in-range entries of Ht to Kkt
      for (; pKkt < krowKkt[i + 1]; pKkt++) {
         const int colKkt = jcolKkt[pKkt];

         if (colKkt >= endCol)
            break;

         if (colKkt == colHt) {
            assert(pHt < krowHt[i + 1]);

            MKkt[pKkt] += MHt[pHt++];

            // end of Ht row reached?
            if (pHt == krowHt[i + 1])
               break;

            colHt = jcolHt[pHt] + nHtOffsetCols;
         }
      }

      assert(
         pHt == krowHt[i + 1] || jcolHt[pHt] + nHtOffsetCols >= endCol); // asserts that no entry of Ht has been missed
   }
}

// TODO : should be const
/** rxy = beta*rxy + alpha * SC * x */
void sLinsysRootAug::SCmult(double beta, SimpleVector<double>& rxy, double alpha, SimpleVector<double>& x) {
   assert(false && "TODO: implement regularization");
   //if (iAmDistrib) {
   //only one process subtracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy + F'*xxl + G'*xyl ] from r
   //                           [  A*xx  + xReg*xy                              ]
   //                           [  F*xx                 yRegLink*xxl            ]
   //                           [  G*xx                           + Omega * xyl ]

   const int myRank = PIPS_MPIgetRank(mpiComm);
   const int mpiCommSize = PIPS_MPIgetSize(mpiComm);
   const bool iAmLastRank = (myRank == mpiCommSize - 1);
   assert(mpiCommSize >= 1);

   if (iAmLastRank) { // needs to be the last rank because only this rank is guaranteed to have CtDC
      assert(rxy.length() == locnx + locmy + locmyl + locmzl);

      //only this proc subtracts from rxy
      rxy.scalarMult(beta);
      const SparseSymmetricMatrix& Q = data->getLocalQ();
      Q.getStorage().mult(1.0, &rxy[0], alpha, &x[0]);

      if (locmz > 0) {
         auto* CtDC_sp = dynamic_cast<SparseSymmetricMatrix*>(CtDC.get());
         assert(CtDC_sp);

         CtDC_sp->getStorage().multSym(1.0, &rxy[0], -alpha, &x[0]);
      }

      auto& xDiagv = dynamic_cast<SimpleVector<double>&>(*xDiag);
      assert(xDiagv.length() == locnx);
      for (int i = 0; i < xDiagv.length(); i++)
         rxy[i] += alpha * xDiagv[i] * x[i];

      const SparseMatrix& A = data->getLocalB();
      A.getStorage().transMult(1.0, &rxy[0], alpha, &x[locnx]);
      A.getStorage().mult(1.0, &rxy[locnx], alpha, &x[0]);

      assert(locmyl >= 0 && locmzl >= 0);

      if (locmyl > 0) {
         const SparseMatrix& F = data->getLocalF();
         F.getStorage().transMult(1.0, &rxy[0], alpha, &x[locnx + locmy]);
         F.getStorage().mult(1.0, &rxy[locnx + locmy], alpha, &x[0]);
      }

      if (locmzl > 0) {
         const SparseMatrix& G = data->getLocalG();
         G.getStorage().transMult(1.0, &rxy[0], alpha, &x[locnx + locmy + locmyl]);
         G.getStorage().mult(1.0, &rxy[locnx + locmy + locmyl], alpha, &x[0]);

         auto& zDiagLinkConsv = dynamic_cast<SimpleVector<double>&>(*zDiagLinkCons);
         assert(zDiagLinkConsv.length() == locmzl);
         const int shift = locnx + locmy + locmyl;
         for (int i = 0; i < zDiagLinkConsv.length(); i++)
            rxy[i + shift] += alpha * zDiagLinkConsv[i] * x[i + shift];
      }
   } else {
      //other processes set r to zero since they will get this portion from process 0
      rxy.setToZero();
   }

#ifdef TIMING
   taux=MPI_Wtime();
#endif
   // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
   SimpleVector<double> xx((locmyl || locmzl) ? (locnx + locmy + locmyl + locmzl) : locnx);
   xx.copyFromArray(x.elements());
   xx.scalarMult(-alpha);

   for (size_t it = 0; it < children.size(); it++) {
      children[it]->addTermToSchurResidual(rxy, xx);
   }

#ifdef TIMING
   tchild_total +=  (MPI_Wtime()-taux);
#endif
   //~done computing residual

#ifdef TIMING
   taux=MPI_Wtime();
#endif
   //all-reduce residual
   if (iAmDistrib) {
      PIPS_MPIsumArrayInPlace(rxy.elements(), locnx + locmy + locmyl + locmzl, mpiComm);
   }
#ifdef TIMING
   tcomm_total += (MPI_Wtime()-taux);
#endif

}


void sLinsysRootAug::solveWithIterRef(SimpleVector<double>& r) {
   assert(false && " TODO : not sure if working correctly...");
   SimpleVector<double> r2(&r[locnx], locmy);
   SimpleVector<double> r1(&r[0], locnx);

   //SimpleVector<double> realRhs(&r[0], locnx+locmy);
#ifdef TIMING
   taux=MPI_Wtime();
#endif

   double rhsNorm = r.two_norm(); //r== the initial rhs of the reduced system here

   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);
   SimpleVector<double> rxy(locnx + locmy);
   rxy.copyFrom(r);
   SimpleVector<double> x(locnx + locmy);
   x.setToZero(); //solution
   SimpleVector<double> dx(locnx + locmy);                //update from iter refinement
   SimpleVector<double> x_prev(locnx + locmy);
   int refinSteps = 0;
   std::vector<double> histResid;
   int maxRefinSteps = (gLackOfAccuracy > 0 ? 9 : 8);
   do { //iterative refinement
#ifdef TIMING
      taux=MPI_Wtime();
#endif

      x_prev.copyFrom(x);
      //dx = Ainv * r
      dx.copyFrom(rxy);
      solver->Dsolve(dx);
      //update x
      x.axpy(1.0, dx);

#ifdef TIMING
      troot_total += (MPI_Wtime()-taux);
#endif

      if (gLackOfAccuracy < 0)
         break;
      if (refinSteps == maxRefinSteps)
         break;

      //////////////////////////////////////////////////////////////////////
      //iterative refinement
      //////////////////////////////////////////////////////////////////////
      //compute residual

      //if (iAmDistrib) {
      //only one process substracts [ (Q+Dx0+C'*Dz0*C)*xx + A'*xy ] from r
      //                            [  A*xx                       ]
      if (myRank == 0) {
         rxy.copyFrom(r);
         if (locmz > 0) {
            auto* CtDC_sp = dynamic_cast<SparseSymmetricMatrix*>(CtDC.get());
            CtDC_sp->getStorage().multSym(1.0, &rxy[0], 1.0, &x[0]);
         }
         const SparseSymmetricMatrix& Q = data->getLocalQ();
         Q.getStorage().mult(1.0, &rxy[0], -1.0, &x[0]);

         auto& xDiagv = dynamic_cast<SimpleVector<double>&>(*xDiag);
         assert(xDiagv.length() == locnx);
         for (int i = 0; i < xDiagv.length(); i++)
            rxy[i] -= xDiagv[i] * x[i];

         const SparseMatrix& A = data->getLocalB();
         A.getStorage().transMult(1.0, &rxy[0], -1.0, &x[locnx]);
         A.getStorage().mult(1.0, &rxy[locnx], -1.0, &x[0]);
      } else {
         //other processes set r to zero since they will get this portion from process 0
         rxy.setToZero();
      }

#ifdef TIMING
      taux=MPI_Wtime();
#endif
      // now children add [0 A^T C^T ]*inv(KKT)*[0;A;C] x
      SimpleVector<double> xx(&x[0], locnx);
      for (auto & it : children) {
         it->addTermToSchurResidual(rxy, xx);
      }
#ifdef TIMING
      tchild_total +=  (MPI_Wtime()-taux);
#endif
      //~done computing residual

#ifdef TIMING
      taux=MPI_Wtime();
#endif
      //all-reduce residual
      if (iAmDistrib) {
         dx.setToZero(); //we use dx as the recv buffer
         MPI_Allreduce(rxy.elements(), dx.elements(), locnx + locmy, MPI_DOUBLE, MPI_SUM, mpiComm);
         rxy.copyFrom(dx);
      }
#ifdef TIMING
      tcomm_total += (MPI_Wtime()-taux);
#endif

      double relResNorm = rxy.two_norm() / rhsNorm;

      if (relResNorm < 1.0e-10) {
         break;
      } else {
         double prevRelResNorm = 1.0e10;
         if (!histResid.empty())
            prevRelResNorm = histResid[histResid.size() - 1];

         //check for stop, divergence or slow convergence conditions
         if (relResNorm > prevRelResNorm) {
            // diverging; restore iteration
            if (myRank == 0) {
               std::cout << "1st stg - iter refinement diverges relResNorm=" << relResNorm << "  before was "
                         << prevRelResNorm << "\n";
               std::cout << "Restoring iterate.\n";
            }
            x.copyFrom(x_prev);
            break;
         } else {
            //check slow convergence for the last xxx iterates.
            // xxx is 1 for now
            //if(relResNorm>0.*prevRelResNorm) {

            //  if(myRank==0) {
            //    cout << "1st stg - iter refinement stuck relResNorm=" << relResNorm
            //	 << "  before was " << prevRelResNorm << endl;
            //    cout << "exiting refinement." << endl;
            //  }
            //  break;
            //
            //} else {
            //  //really nothing, continue
            //}
         }
         histResid.push_back(relResNorm);
         if (myRank == 0)
            std::cout << "1st stg - sol does NOT  have enough accuracy (" << relResNorm << ") after " << refinSteps
                      << " refinement steps\n";
      }
      refinSteps++;
   } while (refinSteps <= maxRefinSteps);

#ifdef TIMING
   taux = MPI_Wtime();
#endif

   r1.copyFrom(x);
   r2.copyFromArray(&x[locnx]);

#ifdef TIMING
   troot_total += (MPI_Wtime()-taux);
#endif
}

void sLinsysRootAug::solveWithBiCGStab(SimpleVector<double>& b) {
   assert(false && "TODO: not sure if working correctly");
   int n = b.length();

   const int maxit = 75; //500
   const double tol = 1e-10, EPS = 1e-15; // EPS=2e-16

   int myRank;
   MPI_Comm_rank(mpiComm, &myRank);

   SimpleVector<double> r(n);           //residual
   SimpleVector<double> s(n);           //residual associated with half iterate
   SimpleVector<double> rt(n);          //shadow residual
   SimpleVector<double> xmin(n);        //minimal residual iterate
   SimpleVector<double> x(n);           //iterate
   SimpleVector<double> xhalf(n);       // half iterate of BiCG
   SimpleVector<double> p(n), paux(n);
   SimpleVector<double> v(n), t(n);
   int flag;
   double n2b;                  //norm of b
   double normr, normrmin;      //norm of the residual and norm of residual at min-resid iterate
   double normr_act;
   double tolb;                 //relative tolerance
   double rho, omega;
   double alpha = -1;
   int stag, maxmsteps, maxstagsteps, moresteps;
   //double imin;
   //maxit = n/2+1;

   //////////////////////////////////////////////////////////////////
   //  Problem Setup and initialization
   //////////////////////////////////////////////////////////////////

   n2b = b.two_norm();
   tolb = n2b * tol;

   tolb = std::max(tolb, EPS);

#ifdef TIMING
   double relres;
   double iter=0.0;
   if( myRank == 0 )
      std::cout << "initial norm of b " << n2b << "\n";
   taux = MPI_Wtime();
#endif
   //initial guess
   x.copyFrom(b);

   solver->Dsolve(x);
   //initial residual
   r.copyFrom(b);

#ifdef TIMING
   troot_total += (MPI_Wtime()-taux);
   taux = MPI_Wtime();
#endif

   //applyA(1.0, r, -1.0, x);
   SCmult(1.0, r, -1.0, x);

#ifdef TIMING
   tchild_total +=  (MPI_Wtime()-taux);
#endif

   normr = r.two_norm();
   normr_act = normr;

   if (normr <= tolb) {
      //initial guess is good enough
      b.copyFrom(x);
      flag = 0;
      return;
   }

   if (myRank == 0)
      std::cout << "innerBICG starts: " << normr << " > " << tolb << "\n";

   rt.copyFrom(r); //Shadow residual
   auto* resvec = new double[2 * maxit + 1];
   resvec[0] = normr;
   normrmin = normr;
   rho = 1.0;
   omega = 1.0;
   stag = 0;
   maxmsteps = std::min(std::min(n / 50, 5), n - maxit);
   maxstagsteps = 3;
   moresteps = 0;

   //////////////////////////////////////////////////////////////////
   // loop over maxit iterations
   //////////////////////////////////////////////////////////////////
   int ii = 0;
   while (ii < maxit) {
      //cout << ii << " ";
      flag = -1;
      ///////////////////////////////
      // First half of the iterate
      ///////////////////////////////
      double rho1 = rho;
      double beta;
      rho = rt.dotProductWith(r);
      //printf("rho=%g\n", rho);
      if (0.0 == rho) {
         flag = 4;
         break;
      }

      if (ii == 0)
         p.copyFrom(r);
      else {
         beta = (rho / rho1) * (alpha / omega);
         if (beta == 0.0) {
            flag = 4;
            break;
         }

         //-------- p = r + beta*(p - omega*v) --------
         p.axpy(-omega, v);
         p.scale(beta);
         p.axpy(1.0, r);
      }

#ifdef TIMING
      taux = MPI_Wtime();
#endif
      //------ v = A*(M2inv*(M1inv*p)) and ph=M2inv*(M1inv*p)
      //first use v as temp storage
      //applyM1(0.0, v,    1.0, p);
      //applyM2(0.0, paux, 1.0, v);
      //applyA (0.0, v,    1.0, paux);
      paux.copyFrom(p);
      solver->solve(paux);

#ifdef TIMING
      troot_total += (MPI_Wtime()-taux);
#endif

      SCmult(0.0, v, 1.0, paux);

      SimpleVector<double>& ph = paux;

      double rtv = rt.dotProductWith(v);
      if (rtv == 0.0) {
         flag = 4;
         break;
      }

      alpha = rho / rtv;
      if (fabs(alpha) * ph.two_norm() < EPS * x.two_norm())
         stag++;
      else
         stag = 0;

      // xhalf = x + alpha*ph and the associated residual
      xhalf.copyFrom(x);
      xhalf.axpy(alpha, ph);
      s.copyFrom(r);
      s.axpy(-alpha, v);
      normr = s.two_norm();
      normr_act = normr;
      resvec[2 * ii] = normr;

      //printf("iter %g normr=%g\n", ii+0.5, normr);
      //-------- check for convergence in the middle of the iterate.  --------
      if (normr <= tolb || stag >= maxstagsteps || moresteps) {
         s.copyFrom(b);
         //applyA(1.0, s, -1.0, xhalf); // s=b-Ax
         SCmult(1.0, s, -1.0, xhalf);
         normr_act = s.two_norm();

         if (normr <= tolb) {
            //converged
            x.copyFrom(xhalf);
            flag = 0;
#ifdef TIMING
            iter = 0.5+ii;
#endif
            break;
         } else {
            if (stag >= maxstagsteps && moresteps == 0) {
               stag = 0;
            }
            moresteps++;
            if (moresteps >= maxmsteps) {
               //method stagnated
               flag = 3;
               x.copyFrom(xhalf);
               break;
            }
         }
      }
      if (stag >= maxstagsteps) {
         flag = 3;
         break;
      } //stagnation

      //update quantities related to minimal norm iterate
      if (normr_act < normrmin) {
         xmin.copyFrom(xhalf);
         normrmin = normr_act;
         //imin=0.5+ii;
      }

#ifdef TIMING
      taux = MPI_Wtime();
#endif
      ///////////////////////////////
      // Second half of the iterate
      //////////////////////////////
      //applyM1(0.0, t,    1.0, s); //applyM1(s,     stemp);
      //applyM2(0.0, paux, 1.0, t); //applyM2(stemp, sh);
      //applyA (0.0, t,    1.0, paux); //applyA (sh, t);
      //kkt->mult(0.0,paux, 1.0,s);
      paux.copyFrom(s);
      solver->solve(paux);
#ifdef TIMING
      troot_total += (MPI_Wtime()-taux);
#endif

      SCmult(0.0, t, 1.0, paux);

      SimpleVector<double>& sh = paux;
      double tt = t.dotProductWith(t);
      if (tt == 0.0) {
         flag = 4;
         break;
      }

      omega = t.dotProductWith(s);
      omega /= tt;

      if (fabs(omega) * sh.two_norm() < EPS * xhalf.two_norm())
         stag++;
      else
         stag = 0;

      x.copyFrom(xhalf);
      x.axpy(omega, sh); // x=xhalf+omega*sh
      r.copyFrom(s);
      r.axpy(-omega, t); // r=s-omega*t

      normr = r.two_norm();
      normr_act = normr;
      resvec[2 * ii + 1] = normr;

      //printf("stag=%d  maxstagsteps=%d moresteps=%d  normr=%g\n",
      //	   stag, maxstagsteps, moresteps, normr);

      //-------- check for convergence at the end of the iterate.  --------
      if (normr <= tolb || stag >= maxstagsteps || moresteps) {
         r.copyFrom(b);
         //applyA(1.0, r, -1.0, x); //r=b-Ax
         SCmult(1.0, r, -1.0, x);
         normr_act = r.two_norm();

         if (normr <= tolb) {
            flag = 0;
#ifdef TIMING
            iter = 1.0+ii;
#endif
            break;
         } else {
            if (stag >= maxstagsteps && moresteps == 0) {
               stag = 0;
            }
            moresteps++;
            if (moresteps >= maxmsteps) {
               //method stagnated
               flag = 3;
               break;
            }
         }
      } // end convergence check
      if (stag >= maxstagsteps) {
         flag = 3;
         break;
      } //stagnation

      //update quantities related to minimal norm iterate
      if (normr_act < normrmin) {
         xmin.copyFrom(x);
         normrmin = normr_act;
         //imin=1.5+ii;
      }
      //printf("iter %g normr=%g\n", ii+1.0, normr);
      ///////////////////////////////
      // Next iterate
      ///////////////////////////////
      ii++;

   }//end while

   if (ii >= maxit) {
#ifdef TIMING
      iter=ii;
#endif
      flag = 10;
   }

   if (flag == 0 || flag == -1) {
#ifdef TIMING
      relres = normr_act/n2b;
      if(myRank==0) {
        printf("INNER BiCGStab converged: normResid=%g relResid=%g iter=%g\n",
             normr_act, relres, iter);
      }
#endif
   } else {
      if (ii == maxit)
         flag = 10;//aaa
      //FAILURE -> return minimum resid-norm iterate
      r.copyFrom(b);
      //applyA(1.0, r, -1.0, xmin);
      SCmult(1.0, r, -1.0, xmin);

      normr = r.two_norm();
      if (normr >= normr_act) {
         x.copyFrom(xmin);
         //iter=imin;
#ifdef TIMING
         relres=normr/n2b;
#endif
      } else {
#ifdef TIMING
         iter=1.0+ii;
         relres = normr/n2b;
#endif
      }

#ifdef TIMING
      if(myRank==0) {
        printf("INNERBiCGStab did not NOT converged after %g[%d] iterations.\n", iter,ii);
        printf("\t - Error code %d\n\t - Act res=%g\n\t - Rel res=%g %g\n\n",
          flag, normr, relres, normrmin);
      }
#endif
   }

   if (myRank == 0)
      std::cout << "innerBICG: " << "ii=" << ii << " flag=" << flag << " normr=" << normr << " normr_act=" << normr_act
                << " tolb=" << tolb << "\n";

   biCGStabCommunicateStatus(flag, ii);

   b.copyFrom(x);
   delete[] resvec;
}

void sLinsysRootAug::add_CtDC_to_dense_schur_complement(const SymmetricMatrix& CtDC_loc) {
   auto* const kktd = dynamic_cast<DenseSymmetricMatrix*>(kkt.get());
   double** const dKkt = kktd->Mat();

   const auto& CtDC_sparse = dynamic_cast<const SparseSymmetricMatrix&>(CtDC_loc);

   const int* krow_CtDC = CtDC_sparse.krowM();
   const int* jcol_CtDC = CtDC_sparse.jcolM();
   const double* M_CtDC = CtDC_sparse.M();

   for (int i = 0; i < locnx; i++) {
      for (int p = krow_CtDC[i]; p < krow_CtDC[i + 1]; p++) {
         const int j = jcol_CtDC[p];

         if (j <= i)
            dKkt[i][j] -= M_CtDC[p];
      }
   }
}

void sLinsysRootAug::add_CtDC_to_sparse_schur_complement(const SymmetricMatrix& CtDC_loc) {
   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);
#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();

   const auto& CtDC_sparse = dynamic_cast<const SparseSymmetricMatrix&>(CtDC_loc);
   const int* krow_CtDC = CtDC_sparse.krowM();
   const int* jcol_CtDC = CtDC_sparse.jcolM();
   const double* M_CtDC = CtDC_sparse.M();

   for (int i = 0; i < locnx; i++) {
      for (int p = krow_CtDC[i]; p < krow_CtDC[i + 1]; p++) {
         const int col = jcol_CtDC[p];

         if (col >= i) {
            // get start position of dense kkt block
            const int blockStart = krowKkt[i];
            assert(col < locnx && jcolKkt[blockStart + col - i] == col);

            MKkt[blockStart + col - i] -= M_CtDC[p];
         }
      }
   }
}

void
sLinsysRootAug::compute_CtDC_and_add_to_Schur_complement(SymmetricMatrix*& CtDC_loc, const Vector<double>& diagonal) {
   assert(diagonal.all_of([](auto& v) { return v <= 0.0; }));

   const SparseMatrix& C = data->getLocalD();
   C.matTransDinvMultMat(diagonal, &CtDC_loc);
   assert(CtDC_loc->size() == locnx);
   assert(CtDC_loc);

   if (this->hasSparseKkt) {
      add_CtDC_to_sparse_schur_complement(*CtDC_loc);
   } else {
      add_CtDC_to_dense_schur_complement(*CtDC_loc);
   }
}

void sLinsysRootAug::clear_CtDC_from_dense_schur_complement(const SymmetricMatrix& CtDC_loc) {
   auto* const kktd = dynamic_cast<DenseSymmetricMatrix*>(kkt.get());
   double** const dKkt = kktd->Mat();

   const auto& CtDC_sparse = dynamic_cast<const SparseSymmetricMatrix&>(CtDC_loc);

   const int* krow_CtDC = CtDC_sparse.krowM();
   const int* jcol_CtDC = CtDC_sparse.jcolM();
   const double* M_CtDC = CtDC_sparse.M();

   for (int i = 0; i < locnx; i++) {
      for (int p = krow_CtDC[i]; p < krow_CtDC[i + 1]; p++) {
         const int j = jcol_CtDC[p];

         if (j <= i)
            dKkt[i][j] += M_CtDC[p];
      }
   }
}

void sLinsysRootAug::clear_CtDC_from_sparse_schur_complement(const SymmetricMatrix& CtDC_loc) {
   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);
#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();

   const auto& CtDC_sparse = dynamic_cast<const SparseSymmetricMatrix&>(CtDC_loc);
   const int* krow_CtDC = CtDC_sparse.krowM();
   const int* jcol_CtDC = CtDC_sparse.jcolM();
   const double* M_CtDC = CtDC_sparse.M();

   for (int i = 0; i < locnx; i++) {
      for (int p = krow_CtDC[i]; p < krow_CtDC[i + 1]; p++) {
         const int col = jcol_CtDC[p];

         if (col >= i) {
            // get start position of dense kkt block
            const int blockStart = krowKkt[i];
            assert(col < locnx && jcolKkt[blockStart + col - i] == col);

            MKkt[blockStart + col - i] += M_CtDC[p];
         }
      }
   }
}

void sLinsysRootAug::clear_CtDC_from_schur_complement(const SymmetricMatrix& CtDC_loc) {
   assert(CtDC_loc.size() == locnx);

   if (this->hasSparseKkt) {
      clear_CtDC_from_sparse_schur_complement(CtDC_loc);
   } else {
      clear_CtDC_from_dense_schur_complement(CtDC_loc);
   }
}

//   * Br0 = [  0  A00T F0CT G0CT ]
//   *       [[A0]   0   0    0   ]
//   *       [  C    0   0    0   ]
//   *       [ F0V   0   0    0   ]
//   *       [ G0V   0   0    0   ]
void sLinsysRootAug::addBorderX0ToRhs(DistributedVector<double>& rhs, const SimpleVector<double>& x0,
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

   const auto* A0_border = distributed_tree->was_A0_moved_to_border() ? nullptr : dynamic_cast<const SparseMatrix*>(border.A.first.get());
   long long mA0 = A0_border ? A0_border->n_rows() : 0;

   const auto* A00_border = distributed_tree->was_A0_moved_to_border() ? &data->getLocalB() : nullptr;
   long long mA00 = A00_border ? A00_border->n_rows() : 0;

   auto& C0_border = dynamic_cast<SparseMatrix&>(*border.C.first);
   const auto mC0 = C0_border.n_rows();

   assert(border.F.first);
   assert(border.A.last);

   auto& F0vec_border = dynamic_cast<SparseMatrix&>(*border.A.last);
   auto& G0vec_border = dynamic_cast<SparseMatrix&>(*border.C.last);
   const auto [mF0V, nF0V] = F0vec_border.n_rows_columns();

   auto& F0cons_border = dynamic_cast<SparseMatrix&>(*border.F.first);
   auto& G0cons_border = dynamic_cast<SparseMatrix&>(*border.G.first);
   const auto[mF0C, nF0C] = F0cons_border.n_rows_columns();

#ifndef NDEBUG
   assert(border.C.last);
   assert(border.G.first);

   assert(rhs.first);
   assert(rhs.first->length() == nF0C + mA0 + mC0 + mF0V + G0vec_border.n_rows());
   assert(x0.length() == nF0V + mA00 + mF0C + G0cons_border.n_rows());
#endif

   auto& rhs0 = dynamic_cast<SimpleVector<double>&>(*rhs.first);

   double* rhs01 = &rhs0[0];
   double* rhs02 = &rhs0[nF0C];
   double* rhs03 = &rhs0[nF0C + mA0];
   double* rhs04 = &rhs0[nF0C + mA0 + mC0];
   double* rhs05 = &rhs0[nF0C + mA0 + mC0 + mF0V];

   const double* x01 = &x0[0];
   const double* x02 = &x0[nF0V];
   const double* x03 = &x0[nF0V + mA00];
   const double* x04 = &x0[nF0V + mA00 + mF0C];

   if (A0_border)
      A0_border->getStorage().mult(1.0, rhs02, -1.0, x01);
   C0_border.getStorage().mult(1.0, rhs03, -1.0, x01);
   F0vec_border.getStorage().mult(1.0, rhs04, -1.0, x01);
   G0vec_border.getStorage().mult(1.0, rhs05, -1.0, x01);

   if (A00_border)
      A00_border->getStorage().transMult(1.0, rhs01, -1.0, x02);
   F0cons_border.getStorage().transMult(1.0, rhs01, -1.0, x03);
   G0cons_border.getStorage().transMult(1.0, rhs01, -1.0, x04);
}

void sLinsysRootAug::addBorderTimesRhsToB0(DistributedVector<double>& rhs, SimpleVector<double>& b0,
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

      const auto* A0_border = distributed_tree->was_A0_moved_to_border() ? nullptr : dynamic_cast<const SparseMatrix*>(border.A.first.get());
      long long mA0 = A0_border ? A0_border->n_rows() : 0;

      const auto* A00_border = distributed_tree->was_A0_moved_to_border() ? &data->getLocalB() : nullptr;
      long long mA00 = A00_border ? A00_border->n_rows() : 0;

      auto& C0_border = dynamic_cast<SparseMatrix&>(*border.C.first);
      const auto mC0 = C0_border.n_rows();

      assert(border.F.first);
      assert(border.A.last);
      auto& F0vec_border = dynamic_cast<SparseMatrix&>(*border.A.last);
      const auto [mF0V, nF0V] = F0vec_border.n_rows_columns();
      auto& F0cons_border = dynamic_cast<SparseMatrix&>(*border.F.first);
      const auto [mF0C, nF0C]  = F0cons_border.n_rows_columns();

      assert(border.C.last);
      assert(border.G.first);
      auto& G0vec_border = dynamic_cast<SparseMatrix&>(*border.C.last);
      const auto mG0V = G0vec_border.n_rows();
      auto& G0cons_border = dynamic_cast<SparseMatrix&>(*border.G.first);
      const auto mG0C = G0cons_border.n_rows();

      assert(rhs.first);
      assert(rhs.first->length() == nF0C + mA0 + mC0 + mF0V + mG0V);
      assert(b0.length() == nF0V + mA00 + mF0C + mG0C);

      auto& zi = dynamic_cast<SimpleVector<double>&>(*rhs.first);

      SimpleVector<double> zi1(&zi[0], nF0C);
      SimpleVector<double> zi2(&zi[nF0C], mA0);
      SimpleVector<double> zi3(&zi[nF0C + mA0], mC0);
      SimpleVector<double> zi4(&zi[nF0C + mA0 + mC0], mF0V);
      SimpleVector<double> zi5(&zi[nF0C + mA0 + mC0 + mF0V], mG0V);

      SimpleVector<double> b1(&b0[0], nF0V);
      SimpleVector<double> b2(&b0[nF0V], mA00);
      SimpleVector<double> b3(&b0[nF0V + mA00], mF0C);
      SimpleVector<double> b4(&b0[nF0V + mA00 + mF0C], mG0C);

      if (A0_border)
         A0_border->transMult(1.0, b1, -1.0, zi2);
      C0_border.transMult(1.0, b1, -1.0, zi3);
      F0vec_border.transMult(1.0, b1, -1.0, zi4);
      G0vec_border.transMult(1.0, b1, -1.0, zi5);

      if (A00_border)
         A00_border->mult(1.0, b2, -1.0, zi1);
      F0cons_border.mult(1.0, b3, -1.0, zi1);
      G0cons_border.mult(1.0, b4, -1.0, zi1);
   }
}

void sLinsysRootAug::add_regularization_local_kkt(double primal_regularization, double dual_equality_regularization,
   double dual_inequality_regularization) {
   assert(apply_regularization);
   assert(primal_regularization_diagonal);

   assert(dynamic_cast<const DistributedVector<double>*>(this->primal_regularization_diagonal.get()));
   assert(dynamic_cast<const DistributedVector<double>*>(this->dual_equality_regularization_diagonal.get()));
   assert(dynamic_cast<const DistributedVector<double>*>(this->dual_inequality_regularization_diagonal.get()));


   assert(primal_regularization >= 0);
   assert(dual_inequality_regularization >= 0);
   assert(dual_equality_regularization >= 0);

   if (PIPS_MPIgetRank() == 0) {
      std::cout << "regularizing with root " << primal_regularization << " " << dual_equality_regularization << " "
                << dual_inequality_regularization << std::endl;
   }

   /* primal diagonal */
   if (locnx > 0) {
      const auto& primal_regularization_vec = dynamic_cast<DistributedVector<double>&>(*this->primal_regularization_diagonal).first;
      assert(primal_regularization_vec);

      primal_regularization_vec->addConstant(primal_regularization);
      kkt->diagonal_add_constant_from(0, locnx, primal_regularization);
   }

   /* C^T reg^-1 C block */
   if (locmz > 0) {
      const auto& dual_inequality_regularization_vec = dynamic_cast<DistributedVector<double>&>(*this->dual_inequality_regularization_diagonal).first;
      dual_inequality_regularization_vec->addConstant(-dual_inequality_regularization);
      assert(dual_inequality_regularization_vec);

      assert(zDiag);
      assert(CtDC);

      if (!dual_inequality_non_link_diagonal_regularized)
         dual_inequality_non_link_diagonal_regularized.reset(dynamic_cast<SimpleVector<double>*>(zDiag->clone()));

      dual_inequality_non_link_diagonal_regularized->copyFrom(*zDiag);
      dual_inequality_non_link_diagonal_regularized->axpy(1.0, *dual_inequality_regularization_vec);

      SymmetricMatrix* CTDC_ptr = CtDC.get();
      clear_CtDC_from_schur_complement(*CtDC);
      compute_CtDC_and_add_to_Schur_complement(CTDC_ptr, *dual_inequality_non_link_diagonal_regularized);
   }


   /* A0 dual equalities */
   if (locmy > 0) {
      const auto& dual_equality_regularization_vec = dynamic_cast<DistributedVector<double>&>(*this->dual_equality_regularization_diagonal).first;
      assert(dual_equality_regularization_vec);

      dual_equality_regularization_vec->addConstant(-dual_equality_regularization);
      kkt->diagonal_add_constant_from(locnx, locmy, -dual_equality_regularization);
   }

   /* dual linking equalities */
   if (locmyl > 0) {
      const auto& dual_equality_regularization_link_cons = dynamic_cast<DistributedVector<double>&>(*this->dual_equality_regularization_diagonal).last;
      assert(dual_equality_regularization_link_cons);

      dual_equality_regularization_link_cons->addConstant(-dual_equality_regularization);
      kkt->diagonal_add_constant_from(locnx + locmy, locmyl, -dual_equality_regularization);
   }

   if (locmzl > 0) {
      const auto& dual_inequality_regularization_link_cons = dynamic_cast<DistributedVector<double>&>(*this->dual_inequality_regularization_diagonal).last;
      assert(dual_inequality_regularization_link_cons);

      dual_inequality_regularization_link_cons->addConstant(-dual_inequality_regularization);
      kkt->diagonal_add_constant_from(locnx + locmy + locmyl, locmzl, -dual_inequality_regularization);
   }
}

void sLinsysRootAug::finalizeKKTsparse() {
   auto& kkts = dynamic_cast<SparseSymmetricMatrix&>(*kkt);

#ifndef NDEBUG
   int* const jcolKkt = kkts.jcolM();
#endif
   int* const krowKkt = kkts.krowM();
   double* const MKkt = kkts.M();
   const int n0Links = data->getN0LinkVars();

   assert(kkts.size() == locnx + locmy + locmyl + locmzl);
   assert(!kkts.is_lower());
   assert(locmyl >= 0 && locmzl >= 0);

   //////////////////////////////////////////////////////
   // compute Q+diag(xdiag) - C' * diag(zDiag) * C
   // and update the KKT
   //////////////////////////////////////////////////////

   /////////////////////////////////////////////////////////////
   // update the KKT with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   assert(data->getLocalQ().krowM()[locnx] == 0 && "Q currently not supported for sparse kkt");

   /////////////////////////////////////////////////////////////
   // update the KKT with the diagonals
   // xDiag is in fact diag(Q)+X^{-1} S
   /////////////////////////////////////////////////////////////
   if (locnx > 0) {
      assert(xDiag);
      assert(xDiag->all_of([](const auto& d) { return d >= 0.0; }));
      kkt->atAddDiagonal(0, *xDiag);
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if (locmz > 0) {
      assert(zDiag);

      SymmetricMatrix* CtDCptr = CtDC ? CtDC.get() : nullptr;
      compute_CtDC_and_add_to_Schur_complement(CtDCptr, *zDiag);
      if (!CtDC) {
         CtDC.reset(CtDCptr);
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with At (symmetric update forced)
   /////////////////////////////////////////////////////////////
   if (locmy > 0) {
      const SparseMatrix& At = data->getLocalB().getTranspose(); // yes, B
      const double* MAt = At.M();
      const int* krowAt = At.krowM();

      for (int i = 0; i < locnx; ++i) {
         const int pstart = krowAt[i];
         const int pend = krowAt[i + 1];

         // get start position of sparse kkt block
         const int blockStart = krowKkt[i] + locnx - i;

         assert(blockStart <= krowKkt[i + 1]);

         for (int p = pstart; p < pend; ++p) {
            assert(At.jcolM()[p] < locmy);
            assert(blockStart + (p - pstart) <= krowKkt[i + 1]);
            assert(jcolKkt[blockStart + (p - pstart)] == (locnx + At.jcolM()[p]));

            MKkt[blockStart + (p - pstart)] += MAt[p];
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Ft
   /////////////////////////////////////////////////////////////
   if (locmyl > 0) {
      const SparseMatrix& Ft = data->getLocalF().getTranspose();
      const double* MFt = Ft.M();
      const int* krowFt = Ft.krowM();
      const int* jcolFt = Ft.jcolM();

      const int* krowGt = nullptr;

      if (locmzl > 0) {
         const SparseMatrix& Gt = data->getLocalG().getTranspose();
         krowGt = Gt.krowM();
      }

      for (int i = 0; i < locnx; ++i) {
         const bool sparseRow = (i >= locnx - n0Links);
         const int pend = krowFt[i + 1];

         if (sparseRow) {
            int blockStart = krowKkt[i + 1] - (krowFt[i + 1] - krowFt[i]);

            if (locmzl > 0)
               blockStart -= (krowGt[i + 1] - krowGt[i]);

            assert(blockStart >= krowKkt[i]);

            for (int p = krowFt[i], shift = 0; p < pend; ++p, ++shift) {
               assert(jcolFt[p] < locmyl && jcolKkt[blockStart + shift] == (locnx + locmy + jcolFt[p]));

               MKkt[blockStart + shift] += MFt[p];
            }
         } else {
            const int blockStart = krowKkt[i + 1] - locmyl - locmzl;
            assert(blockStart >= krowKkt[i]);

            for (int p = krowFt[i]; p < pend; ++p) {
               const int col = jcolFt[p];
               assert(col < locmyl && jcolKkt[blockStart + col] == (locnx + locmy + col));

               MKkt[blockStart + col] += MFt[p];
            }
         }
      }
   }

   /////////////////////////////////////////////////////////////
   // update the KKT with Gt and add z diagonal
   /////////////////////////////////////////////////////////////
   if (locmzl > 0) {
      const SparseMatrix& Gt = data->getLocalG().getTranspose();
      const double* MGt = Gt.M();
      const int* krowGt = Gt.krowM();
      const int* jcolGt = Gt.jcolM();

      for (int i = 0; i < locnx; ++i) {
         const int pend = krowGt[i + 1];
         const bool sparseRow = (i >= locnx - n0Links);

         if (sparseRow) {
            const int blockStart = krowKkt[i + 1] - (krowGt[i + 1] - krowGt[i]);

            assert(blockStart >= krowKkt[i]);

            for (int p = krowGt[i], shift = 0; p < pend; ++p, ++shift) {
               assert(jcolGt[p] < locmzl && jcolKkt[blockStart + shift] == (locnx + locmy + locmyl + jcolGt[p]));

               MKkt[blockStart + shift] += MGt[p];
            }
         } else {
            const int blockStart = krowKkt[i + 1] - locmzl;

            assert(blockStart >= krowKkt[i]);

            for (int p = krowGt[i]; p < pend; ++p) {
               const int col = jcolGt[p];
               assert(col < locmzl && jcolKkt[blockStart + col] == (locnx + locmy + locmyl + col));

               MKkt[blockStart + col] += MGt[p];
            }
         }
      }

      const auto& szDiagLinkCons = dynamic_cast<const SimpleVector<double>&>(*zDiagLinkCons);
      kkt->atAddDiagonal(locnx + locmy + locmyl, szDiagLinkCons);
   }

#ifdef DUMPKKT
   ofstream myfile;
   myfile.open("../sparsekkt");

   int zerocount = 0;
   const int sizeKkt = locnx + locmy + locmyl + locmzl;

   for( int r = 0; r < sizeKkt; r++ )
   {
      for( int i = krowKkt[r]; i < krowKkt[r + 1]; i++ )
      {
         const double val = MKkt[i];
         const double col = jcolKkt[i];
         if( val != 0.0 )
            myfile << r << " " << col << " " << val << "\n";
         else
            zerocount++;
      }
   }

   std::cout << "zero-count " << zerocount << " of " << krowKkt[sizeKkt] << "\n";

   myfile.close();

   assert(0);

#endif
}

void sLinsysRootAug::finalizeKKTdense() {
   auto& schur_complement = dynamic_cast<DenseSymmetricMatrix&>(*kkt);

   /////////////////////////////////////////////////////////////
   // update the schur complement with Q (DO NOT PUT DIAG)
   /////////////////////////////////////////////////////////////
   // update the schur complement with the diagonals
   // xDiag is in fact diag(Q)+X^{-1}S
   /////////////////////////////////////////////////////////////
   if (locnx > 0) {
      const auto& Q = data->getLocalQ();
      assert(Q.krowM()[Q.size()] == 0 && "Currently only empty Q supported!!");
      schur_complement.add_matrix_at_without_diag(Q, 0, 0);

      schur_complement.atAddDiagonal(0, *xDiag);
   }

   /////////////////////////////////////////////////////////////
   // eliminate C from the schur complement
   // update with   - C' * diag(zDiag) *C
   /////////////////////////////////////////////////////////////
   if (locmz > 0) {
      assert(zDiag);
      SymmetricMatrix* CtDCptr = CtDC ? CtDC.get() : nullptr;

      compute_CtDC_and_add_to_Schur_complement(CtDCptr, *zDiag);

      if (!CtDC) {
         CtDC.reset(CtDCptr);
      }
   }

   /////////////////////////////////////////////////////////////
   // update the schur complement with A
   /////////////////////////////////////////////////////////////
   if (locmy > 0) {
      schur_complement.add_matrix_at(data->getLocalB(), locnx, 0); // yes, B
   }

   /////////////////////////////////////////////////////////////
   // update the schur complement with F
   /////////////////////////////////////////////////////////////
   if (locmyl > 0) {
      schur_complement.add_matrix_at(data->getLocalF(), locnx + locmy, 0);
   }

   /////////////////////////////////////////////////////////////
   // update the schur complement with G and put z diagonal
   /////////////////////////////////////////////////////////////
   if (locmzl > 0) {
      assert(zDiagLinkCons);
      schur_complement.add_matrix_at(data->getLocalG(), locnx + locmy + locmyl, 0);
      schur_complement.atAddDiagonal(locnx + locmy + locmyl, *zDiagLinkCons);
   }

#ifdef DUMPKKT
   const int msize = locnx + locmy + locmyl + locmzl;

   ofstream myfile;
   myfile.open("../densekkt");

   for( int col = 0; col < msize; col++ )
      for( int row = col; row < msize; row++ )
         if( dKkt[row][col] != 0.0 )
            myfile << col << " " << row << " " << dKkt[row][col] << "\n";

   myfile.close();

   assert(0);
#endif
}

void sLinsysRootAug::DsolveHierarchyBorder(DenseMatrix& rhs_mat_transp, int n_cols) {
   /* b holds all rhs in transposed form - C part from schur complement is already missing in b */
   const int my_rank = PIPS_MPIgetRank(mpiComm);
#ifdef TIMING
   // TODO
#endif

   const auto[m, n] = rhs_mat_transp.n_rows_columns();
#ifndef NDEBUG
   assert(locmyl >= 0 && locmzl >= 0);

   assert(n_cols <= m);
   assert(locnx + locmy + locmz + locmyl + locmzl == n);
#endif

   /*
    * for every right hand side one of the processes now does the SC solve operation and puts it at the corresponding
    * position in b
    * Every process has to do n_rhs / n_procs right hand sides while the first few might have to solve with one additional one
    */
   const int size = PIPS_MPIgetSize(mpiComm);
   const int n_blockrhs = static_cast<int>( n_cols / size );
   const int leftover = n_cols % size;

   const int n_rhs = (my_rank < leftover) ? n_blockrhs + 1 : n_blockrhs;
   const int rhs_start =
      my_rank < leftover ? (n_blockrhs + 1) * my_rank : (n_blockrhs + 1) * leftover + (my_rank - leftover) * n_blockrhs;

   assert(rhs_start <= n_cols);
   assert(rhs_start + n_rhs <= n_cols);

   // set rhs contributed by other procs to zero
#pragma omp parallel for schedule(dynamic, 1)
   for (int rhs_i = 0; rhs_i < n_cols; ++rhs_i) {
      if (rhs_start <= rhs_i && rhs_i < rhs_start + n_rhs)
         continue;
      SimpleVector<double> b(rhs_mat_transp[rhs_i], n);
      b.setToZero();
   }

   solveReducedLinkConsBlocked(rhs_mat_transp, rhs_start, n_rhs);

   if (iAmDistrib) {
      // TODO only allreduce relevant part
      // TODO is allreduce even worth it here? Every proc could also compute all its rhs - add if n_rhs big
      submatrixAllReduceFull(rhs_mat_transp, 0, 0, m, n, mpiComm);
   }
}

void sLinsysRootAug::addBlTKiInvBrToRes(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border, bool sym_res,
   bool sparse_res) {
   assert(!is_hierarchy_root);

   if (Bl.isEmpty() || (Br.isEmpty() && Br_mod_border.empty()))
      return;

   const int m_result = result.n_rows();

   assert(m_result > 0);
   assert(blocksize_hierarchical > 0);

   // buffer b0 for blockwise computation of Br0 - SUM_i  Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ), stored in transposed form (for quick access of cols in solve)
   // dense since we have no clue about any structure in the system and Xij are dense
   const int n_buffer = locnx + locmy + locmz + locmyl + locmzl;
   const int m_buffer = allocateAndZeroBlockedComputationsBuffer(m_result, n_buffer);

   const int n_chunks = std::ceil(static_cast<double>(m_result) / m_buffer);

   assert(n_chunks > 0);
   if (!sc_compute_blockwise_hierarchical) {
      assert(n_chunks == 1);
      assert(m_buffer == m_result);
   }

   for (int i = 0; i < n_chunks; ++i) {
      const int begin_chunk = i * m_buffer;
      const int end_chunk = std::min(m_result, (i + 1) * m_buffer);

      assert(end_chunk - begin_chunk <= m_buffer);
      addBlTKiInvBrToResBlockwise(result, Bl, Br, Br_mod_border, sym_res, sparse_res, *buffer_blocked_hierarchical,
         begin_chunk, end_chunk);
   }
}

//TODO: determine start and end of the two-links (exclusively on first and last process);
//TODO: probably separate method for two-links?

/* if Bl^T is a two link border this needs to be done for at most the first and the last child of this communicator (all the other children are guaranteed 0) */
/* compute res += [ Bl^T Ki^-1 (Br - sum_j Bmodj Xj) ]^T = (Br^T - SUM_j Xj^T Bmodj^T) Ki^-1 Bl */
void sLinsysRootAug::addBlTKiInvBrToResBlockwise(AbstractMatrix& result, BorderLinsys& Bl, BorderLinsys& Br,
   std::vector<BorderMod>& Br_mod_border,
   bool sym_res, bool sparse_res, DenseMatrix& buffer_b0, int begin_cols, int end_cols) {
   /* only called on sLinsysRootBordered and sLinsysRootAugHierInner */
   buffer_b0.putZeros();

   /* buffer_b0 is in transposed for so that we can access its cols (or rows in the storage) quickly in Dsolve */
   assert(!is_hierarchy_root);
   assert(0 <= begin_cols && begin_cols <= end_cols);
   assert(buffer_b0.n_columns() == kkt->size() || buffer_b0.n_columns() == kkt->size() + locmz);
   assert(buffer_b0.n_rows() >= end_cols - begin_cols);

   const bool two_link_border_left = !(Bl.has_RAC || Bl.use_local_RAC);
   const bool two_link_border_right = !(Br.has_RAC || Br.use_local_RAC);

   if (two_link_border_right);

   /* compute Schur Complement right hand sides SUM_i Bi_{inner} Ki^-1 ( Bri - sum_j Bmodij Xij )
    * (keep in mind that in Bi_{this} and the SC we projected C0 Omega0 out) */
   // buffer_b0 = - [ SUM_i Bi_{inner}^T Ki^{-1} (Bri - SUM_j Bmodij Xij) ]^T = SUM_i (Bri^T - SUM_j Xij^T Bmodij^T) Ki^{-1} Bi_{inner}
   LsolveHierarchyBorder(buffer_b0, Br, Br_mod_border, two_link_border_left, begin_cols, end_cols);

   // buffer_b0 = (Br0 - sum_j Bmod0J X0j ) - buffer_b0 = Br0 - sum_j Bmod0J X0j - SUM_i Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij )}
   finalizeZ0Hierarchical(buffer_b0, Br, Br_mod_border, begin_cols, end_cols);

   // solve with Schur Complement for B0_{outer} - SUM_i Bi_{inner}^T Ki^{-1} ( Bri - sum_j Bmodij Xij ) (stored in transposed form!)
   // buffer_b0 = SC_{inner}^-1 buffer_b0 = X0
   DsolveHierarchyBorder(buffer_b0, end_cols - begin_cols);

   // compute result += -SUM_i Bli^T Ki^{-1} ( ( Bri - sum_j Bmodij Xij )  - Bi_{inner} X0 ) += -SUM_i Bli^T Xi
   LtsolveHierarchyBorder(result, buffer_b0, Bl, Br, Br_mod_border, sym_res, sparse_res, begin_cols, end_cols);

   // compute result += Bl0^T X0
   if (PIPS_MPIgetRank(mpiComm) == 0)
      finalizeInnerSchurComplementContribution(result, buffer_b0, Bl, sym_res, sparse_res, begin_cols, end_cols);
}
