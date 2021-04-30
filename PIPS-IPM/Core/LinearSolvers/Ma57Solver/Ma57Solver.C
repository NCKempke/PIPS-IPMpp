/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP
 * Modified by Cosmin Petra to perform solves with the factors.
 */

#include "Ma57Solver.h"

#include <algorithm>
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "DenseGenMatrix.h"
#include "pipsport.h"
#include <mpi.h>

extern int gOoqpPrintLevel;


void dumpdata(int* irow, int* jcol, double* M, int, int nnz) {
   printf("======================================================\n");
   for (int i = 0; i < nnz; i++)
      printf("%6d %6d %10.2f\n", irow[i], jcol[i], M[i]);
   printf("\n");
   printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
}

Ma57Solver::Ma57Solver(const SparseSymMatrix* sgm, const std::string& name_) : mat_storage{sgm->getStorageHandle()}, n{mat_storage->n},
      nnz{mat_storage->numberOfNonZeros()}, lkeep{7 * n + nnz + 2 * std::max(n, nnz) + 42}, n_threads{PIPSgetnOMPthreads()}, name(name_) {
   assert(n_threads >= 1);
   x.resize(n * n_threads);
   resid.resize(n * n_threads);
   if (n_threads > 1) {
      info.resize(40 * n_threads);
      rinfo.resize(20 * n_threads);
   }

   init();
}

void Ma57Solver::init() {
   assert(mat_storage->n == mat_storage->m);
   assert(n > 0);

   FNAME(ma57id)(cntl.data(), icntl.data());

   icntl[1] = -1; // don't print warning messages
   // 1 -> use pivot order in KEEP, 0 -> AMD ordering using MC47 (no dense rows), 2 -> AMD using MC47,
   // 3 -> min degree MA27, 4 -> METIS ordering, 5 -> auto choice (default)
   icntl[5] = 4;
   icntl[8] = n_iterative_refinement;
   // default 1 -> scale using symmetrized MC64 version, if != 1 no sclaing
   //icntl[14] = 1;
   // if 0 nothing (default) if 1 remove small entries cntl[1] and place corresponding pivots at the end of factorization -> for highly rank deficient matrices ..
   icntl[15] = 1;

   /* in the solve stage each solver will need its own copy of icntl to switch between solve and iterative refinement */
   icntl.resize(n_threads * icntl.size());
   for (int i = 1; i < n_threads; ++i)
      std::copy(icntl.begin(), icntl.begin() + 20, icntl.begin() + 20 * i);

   // only use if icntl[6] == 1 (default)
   setThresholdPivoting(default_pivoting);
//   setSmallPivot(default_pivtol);
   //   cntl[2] convergence for Ma57DD iterative refinement (norm improvement by factor..) default 0.5
   //   cntl[3] Ma57BD (reuse of old factorization)
   //   cntl[4] static pivoting
}

void Ma57Solver::firstCall() {
   irowM.resize(nnz);
   jcolM.resize(nnz);
   getIndices(irowM, jcolM);

   if (keep.size() < static_cast<size_t>(lkeep))
      keep.resize(lkeep);

   iworkn.resize(std::max(5 * n, n_threads * n * 4));

   FNAME(ma57ad)(&n, &nnz, irowM.data(), jcolM.data(), &lkeep, keep.data(), iworkn.data(), icntl.data(), info.data(), rinfo.data());
   assert(info[0] >= 0);

   lfact = 2 * static_cast<int>(rpessimism * info[8]);
   fact.resize(lfact);

   lifact = static_cast<int>(ipessimism * info[9]);
   ifact.resize(lifact);

   dworkn.resize(n_threads * n * 4);

}

void Ma57Solver::diagonalChanged(int /* idiag */, int /* extent */) {
   matrixChanged();
}

void Ma57Solver::matrixChanged() {
   if (keep.empty())
      firstCall();

   bool done = false;
   int tries = 0;

   do {
      FNAME(ma57bd)(&n, &nnz, mat_storage->M, fact.data(), &lfact, ifact.data(), &lifact, &lkeep, keep.data(), iworkn.data(), icntl.data(),
            cntl.data(), info.data(), rinfo.data());
      done = checkErrorsAndReact();
      ++tries;
   } while (!done && tries < max_tries);

   if (!done && tries > max_tries) {
      std::cerr << "ERROR MA57: " << name << ":could not get factorization of matrix after max " << max_tries << " tries\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
   }
   freshFactor = true;
}

void Ma57Solver::solve(Vector<double>& rhs_in) {
   const int my_id = omp_get_thread_num();
   assert(my_id < n_threads);

   assert(dworkn.size() >= static_cast<unsigned int>((my_id + 1) * n * 4 ));
   assert(iworkn.size() >= static_cast<unsigned int>((my_id + 1) * n ));
   assert(icntl.size() >= static_cast<unsigned int>((my_id + 1) * 20 ));
   assert(info.size() >= static_cast<unsigned int>((my_id + 1) * 40 ));
   assert(rinfo.size() >= static_cast<unsigned int>((my_id + 1) * 20 ));
   assert(x.size() >= static_cast<unsigned int>((my_id + 1) * n ));
   assert(resid.size() >= static_cast<unsigned int>((my_id + 1) * n ));

   double* dworkn_loc = dworkn.data() + my_id * n * 4;
   int* iworkn_loc = iworkn.data() + my_id * n;
   int* icntl_loc = icntl.data() + my_id * 20;
   int* info_loc = info.data() + my_id * 40;
   double* rinfo_loc = rinfo.data() + my_id * 20;

   SimpleVector<double> x_loc(x.data() + my_id * n, n);
   SimpleVector<double> resid_loc(resid.data() + my_id * n, n);

   SimpleVector<double>& rhs = dynamic_cast<SimpleVector<double>&>(rhs_in);

   /* job 1 : perform one step of iterative refinement */
   int job = 0;

   if (freshFactor)
      icntl_loc[8] = 1; // No iterative refinement
   else
      icntl_loc[8] = n_iterative_refinement; // Iterative refinement

   const double rhsnorm = rhs.infnorm();

   bool done = false;

   while (!done) {
      FNAME(ma57dd)(&job, &n, &nnz, mat_storage->M, irowM.data(), jcolM.data(), fact.data(), &lfact, ifact.data(), &lifact, rhs.elements(),
            x_loc.elements(), resid_loc.elements(), dworkn_loc, iworkn_loc, icntl_loc, cntl.data(), info_loc, rinfo_loc);

      done = checkErrorsAndReact();

      if (resid_loc.infnorm() < precision * (1 + rhsnorm))
         done = true;
      else {
         if (thresholdPivoting() < threshold_pivoting_max) {
            // refactor with a higher Threshold Pivoting parameter
            setThresholdPivoting(std::min(thresholdPivoting() * threshold_pivoting_factor, threshold_pivoting_max));

            if (gOoqpPrintLevel >= ooqp_print_level_warnings)
               std::cout << "WARNING MA57 " << name << ": Setting ThresholdPivoting parameter to " << thresholdPivoting()
                         << " for future factorizations\n";
         }

         if (freshFactor) {
            job = 2; // set to iterative refinement
            icntl_loc[8] = 10;

            freshFactor = false;
            done = false;
         }
         else {
            if (gOoqpPrintLevel >= ooqp_print_level_warnings) {
               if (thresholdPivoting() == threshold_pivoting_max)
                  std::cout << "WARNING MA57 " << name << ": Unprecise solution but ThresholdPivoting is already at its max\n";
               else
                  std::cout << "WARNING MA57 " << name << ": Setting ThresholdPivoting to " << thresholdPivoting()
                            << " - refactorization suggested\n";
            }
            done = true;
         }
      }
   }

   rhs.copyFrom(x_loc);
}

void Ma57Solver::solve(int solveType, Vector<double>& rhs_in) {
   if (solveType < 1 || solveType > 4)
      assert("Unknown JOB assigned for use in MA57CD!" && 0);
   else if (solveType == 1) {
      solve(rhs_in);
   }
   else {
      // not thread safe
      int job = solveType;
      int one_rhs = 1;

      SimpleVector<double>& rhs = dynamic_cast<SimpleVector<double>&>(rhs_in);

      double* drhs = rhs.elements();

      FNAME(ma57cd)(&job, &n, fact.data(), &lfact, ifact.data(), &lifact, &one_rhs, drhs, &n, dworkn.data(), &n, iworkn.data(), icntl.data(),
            info.data());
   }
}

// rhs vectors are on the "rows", for continuous memory
void Ma57Solver::solve(GenMatrix& rhs_in) {
   DenseGenMatrix& rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

   int N, NRHS;
   rhs.getSize(NRHS, N);

   assert(n == N);
   solve(NRHS, rhs.elements(), nullptr);
}

void Ma57Solver::solve(int nrhss, double* rhss, int*) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
   /* the multiple rhs (macd) option in MA57 does not allow for iterative refinement */
   for (int i = 0; i < nrhss; i++) {
      SimpleVector<double> v(rhss + i * n, n);

      if (v.isZero())
         continue;
      solve(v);
   }
}

bool Ma57Solver::checkErrorsAndReact() {
   bool error = false;
   const int error_flag = info[0];
   const int error_info = info[1];

   switch (error_flag) {
      case 0 :
         break;
      case -1 : {
         std::cerr << "ERROR MA57 " << name << ": N out of range or < -1: " << error_info << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -2 : {
         std::cerr << "ERROR MA57 " << name << ": NNZ out of range or < -1 : " << error_info << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -3 : {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": insufficient space in fact: " << info[1] << " suggest reset to " << info[16] << "\n";
         rpessimism *= 1.1;

         lfact = std::max(static_cast<int>(rpessimism * lfact), static_cast<int>(rpessimism * info[16]));
         fact.resize(lfact);

         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << " resetting to " << lfact << "\n";

         error = true;
      };
         break;
      case -4 : {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": insufficient factorization space in ifact: " << info[1] << " suggest reset to " << info[17]
                      << "\n";;
         ipessimism *= 1.1;

         lifact = std::max(static_cast<int>(ipessimism * lifact), static_cast<int>(ipessimism * info[17]));
         ifact.resize(lifact);
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << " resetting to " << lifact << "\n";

         error = true;
      }
         break;
      case -5: {
         std::cerr << "WARNING MA57 " << name << ": Small pivot found when pivoting disabled\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -6: {
         std::cerr << "ERROR M527 " << name << ": change of sign of pivots detected at stage even though matrix supposedly definiet" << error_info
                   << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -7: {
         std::cerr << "ERROR MA57 " << name << ": value new int/real array not bigger than old value in MA57ED " << error_info << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -8: {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": iterative refinement failed to converge\n";
         error = true;
      };
         break;
      case -9: {
         std::cerr << "ERROR MA57 " << name << ": error in permutation array pos " << error_info << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -10: {
         std::cerr << "ERROR MA57 " << name << ": icntl[6] " << error_info << " is out of range\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -11: {
         std::cerr << "ERROR MA57 " << name << ": LRHS " << error_info << " smaller than n " << n << " on solve call\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -12: {
         std::cerr << "ERROR MA57 " << name << ": invalid value " << error_info << " for job\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -13: {
         std::cerr << "ERROR MA57 " << name << ": invalid value " << error_info << " for icntl[8]\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -14: {
         std::cerr << "ERROR MA57 " << name << ": MC71AD (icntl > 10) failed inside of MA57DD\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -15: {
         std::cerr << "ERROR MA57 " << name << ": lkeep " << error_info << " less than required" << error_info << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -16: {
         std::cerr << "ERROR MA57 " << name << ": NRHS " << error_info << " less than one" << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -17: {
         std::cerr << "ERROR MA57 " << name << ": lwork too small: " << error_info << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -18: {
         std::cerr << "ERROR MA57 " << name << ": METIS ordering requested but METIS was not linked\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case 1 : {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": detected " << info[2] << " entries out of range in irowM and jcolM; ignored\n";
      };
         break;
      case 2 : {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING M57 " << name << ": detected " << info[3]
                      << " duplicate entries in user supplied matrix detected - summing them up\n";
      }
         break;
      case 3: {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": out of range etries and duplicates detected .. ignoring/summing them up\n";
      };
         break;
      case 4: {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": rank deficient matrix detected; apparent rank is " << info[24] << "\n";
      };
         break;
      case 5: {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": pivots have different sign when factorizing supposedly definite matrix; " << info[25]
                      << " sign changes detected\n";
      };
         break;
      case 8: {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": inf norm of computed solution was zero\n";
      };
         break;
      case 10: {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": insufficient space in fact: " << lfact << "\n";
         rpessimism *= 2;

         lfact = static_cast<int>(rpessimism * lfact);
         fact.resize(lfact);
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << " resetting to " << lfact << "\n";

         error = true;
      };
         break;
      case 11: {
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << "WARNING MA57 " << name << ": insufficient factorization space in ifact: " << lifact << "\n";
         ipessimism *= 2;

         lifact = static_cast<int>(ipessimism * lifact);
         ifact.resize(lifact);
         if (gOoqpPrintLevel >= ooqp_print_level_warnings)
            std::cout << " resetting to " << lifact << "\n";

         error = true;
      };
         break;
      default : {
         assert(error_flag == 0);
      };
         break;
   }

   return error;
}

// TODO same as the one in MA27 - move somewhere else, some common MA_Solver thing maybe..
void Ma57Solver::getIndices(std::vector<int>& irowM, std::vector<int>& jcolM) const {
   const int* krowM = mat_storage->krowM;
   for (int i = 0; i < mat_storage->n; i++) {
      if (mat_storage->fortranIndexed()) {
         assert(krowM[i] - 1 >= 0);
         for (int k = krowM[i] - 1; k < krowM[i + 1] - 1; k++)
            irowM[k] = i + 1;
      }
      else
         for (int k = krowM[i]; k < krowM[i + 1]; k++)
            irowM[k] = i + 1;
   }

   for (int k = 0; k < nnz; k++) {
      if (!mat_storage->fortranIndexed())
         jcolM[k] = mat_storage->jcolM[k] + 1;
      else
         jcolM[k] = mat_storage->jcolM[k];
   }
}

std::tuple<unsigned int, unsigned int, unsigned int> Ma57Solver::get_inertia() const {
   assert(false && "TODO: Implement");
   return {0, 0, 0};
}
