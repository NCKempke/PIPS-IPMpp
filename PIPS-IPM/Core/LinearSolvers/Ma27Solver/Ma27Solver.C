/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#include "Ma27Solver.h"

extern int print_level;

Ma27Solver::Ma27Solver(const SparseSymmetricMatrix& sgm, std::string name_) : mat_storage(&sgm.getStorage()), n{mat_storage->n},
      nnz{sgm.numberOfNonZeros()}, n_threads{PIPSgetnOMPthreads()}, name(std::move(name_))//, scaler( new Mc30Scaler() )
{
   assert(n_threads >= 1);

   iter.resize(n * n_threads);
   iter_best.resize(n * n_threads);
   resid.resize(n * n_threads);

   if (n_threads > 1)
      info.resize(20 * n_threads);
   init();
}

void Ma27Solver::init() {
   assert(mat_storage->n == mat_storage->m);
   assert(n > 0);

   FNAME(ma27id)(icntl.data(), cntl.data());

   setThresholdPivoting(default_pivoting);
   setFratio(default_fratio);
   setSmallPivot(default_pivtol);

   icntl[0] = 0;
   icntl[1] = 0;
}

void Ma27Solver::firstCall() {
   /* convert to MA27 format */
   irowM.resize(nnz);
   jcolM.resize(nnz);
   getIndices(irowM, jcolM);

   /* set working arrays for iterative refinement */
   w_ma60.resize(3 * n);
   iw_ma60.resize(2 * n);

   /* set working stuff for factorization routine */
   liw = static_cast<int>(ipessimism * (2 * nnz + 3 * n + 1));
   iw = new int[liw];

   /* later iw1 will need to be of size nsteps * nsolvers with nsteps <= n -> set it to the max once */
   iw1 = new int[std::max(2 * n, n_threads * n)];
   ikeep = new int[3 * n];

   int iflag = 0; // set to 1 if ikeep contains pivot order
   double ops;

   /* do ordering */
   bool done{};
   int tries = 0;
   do {
      FNAME(ma27ad)(&n, &nnz, irowM.data(), jcolM.data(), iw, &liw, ikeep, iw1, &nsteps, &iflag, icntl.data(), cntl.data(), info.data(), &ops);
      done = !checkErrorsAndReact();
      ++tries;
   } while (!done && tries < max_tries);

   if (!done && tries > max_tries) {
      std::cerr << "ERROR MA27 " << name << ": could not get ordering of matrix after max " << max_tries << " tries\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   delete[] iw;

   la = static_cast<int>(rpessimism * minimumRealWorkspace() * 10.0);
   fact.resize(la);

   // set iw and in prep for calls to ma27bd and ma27cd
   liw = static_cast<int>(ipessimism * minimumIntWorkspace());
   iw = new int[liw];
}

void Ma27Solver::diagonalChanged(int /* idiag */, int /* extent */) {
   this->matrixChanged();
}

void Ma27Solver::matrixChanged() {
   if (fact.empty())
      this->firstCall();

   bool done;
   int tries = 0;
   do {
      // copy M to fact
      copyMatrixElements(fact, la);
      if (scaler)
         scaler->scaleMatrixTripletFormat(n, nnz, fact.data(), irowM.data(), jcolM.data(), true);

      FNAME(ma27bd)(&n, &nnz, irowM.data(), jcolM.data(), fact.data(), &la, iw, &liw, ikeep, &nsteps, &maxfrt, iw1, icntl.data(), cntl.data(),
            info.data());

      done = !checkErrorsAndReact();
      tries++;
   } while (!done && tries < max_tries);

   if (!done && tries > max_tries) {
      std::cout << "ERROR MA27 " << name << ": could not get factorization of matrix after max " << max_tries << " tries" << "\n";
      MPI_Abort(MPI_COMM_WORLD, -1);
   }

   assert(0 < nsteps && nsteps < n);
   assert(0 < maxfrt && nsteps < n);

   delete[] ww;

   ww = new double[n_threads * maxfrt];
   assert(ww);
}

void Ma27Solver::solve(int nrhss, double* rhss, int*) {
   const double tmp = max_n_iter_refinement;
   max_n_iter_refinement = 2;

#pragma omp parallel for schedule(dynamic, 1) num_threads(n_threads)
   for (int i = 0; i < nrhss; i++) {
      SimpleVector<double> v(rhss + i * n, n);
      solve(v);
      //solveIterRef(v); slower for some reason
   }
   max_n_iter_refinement = tmp;
}

void Ma27Solver::solve(Vector<double>& rhs_in) {
   auto& rhs = dynamic_cast<SimpleVector<double>&>(rhs_in);

#ifndef NDEBUG
   for (int i = 0; i < rhs.length(); ++i)
      if (std::fabs(rhs[i]) > 1e50) {
         std::cout << "Big entry in right hand side vector...\n";
         break;
      }
#endif
//   SimpleVector<double>* rhs_cpy = dynamic_cast<SimpleVector<double>*>(rhs_in.cloneFull());

//   /* sparsify rhs */
//   for( int i = 0; i < rhs.length(); ++i )
//      if( std::fabs(rhs[i]) < 1e-15 )
//         rhs[i] = 0.0;
   const int my_id = omp_get_thread_num();

   // define structures to save rhs and store residuals
   SimpleVector<double> iter_loc = SimpleVector<double>(iter.data() + my_id * n, n);
   iter_loc.setToZero();
   SimpleVector<double> best_iter_loc = SimpleVector<double>(iter_best.data() + my_id * n, n);
   best_iter_loc.setToZero();
   SimpleVector<double> residual_loc = SimpleVector<double>(resid.data() + my_id * n, n);
   residual_loc.copyFrom(rhs);

   double best_resid = std::numeric_limits<double>::infinity();


   const double rhsnorm = rhs.two_norm();

   bool done = false;
   int n_iter_ref = 0;

   double rnorm = -1.0;
   double res_last = -1.0;

   assert(my_id < n_threads);

   double* ww_loc = ww + my_id * maxfrt;
   int* iw1_loc = iw1 + my_id * n;
   int* info_loc = info.data() + my_id * 20;

   /* iterative refinement loop */
   while (!done && n_iter_ref < max_n_iter_refinement) {
      assert(maxfrt > 0);
      assert(nsteps > 0);
      assert(ww);
      assert(iw1_loc);

      /* solve DAD y = D*residual */
      if (scaler)
         scaler->scaleVector(residual_loc);
      FNAME(ma27cd)(&n, fact.data(), &la, iw, &liw, ww_loc, &maxfrt, residual_loc.elements(), iw1_loc, &nsteps, icntl.data(), info_loc);
      /* Dy = x */
      if (scaler)
         scaler->scaleVector(residual_loc);

      iter_loc.axpy(1.0, residual_loc);

      residual_loc.copyFrom(rhs);
      /* calculate residual and possibly new rhs */
      mat_storage->multSym(1.0, residual_loc.elements(), -1.0, iter_loc.elements());

      /* res = res - A * drhs where A * drhs_out = drhs_in */
      rnorm = residual_loc.two_norm();

      if (PIPSisEQ(rnorm, res_last, 1e-7) || rnorm > 100 * best_resid)
         done = true;
      else
         res_last = rnorm;

      if (rnorm < best_resid) {
         best_resid = rnorm;
         best_iter_loc.copyFrom(iter_loc);
      }

      if (rnorm < precision * (1.0 + rhsnorm))
         done = true;

      ++n_iter_ref;

      if (done && rnorm >= precision * (1.0 + rhsnorm)) {
         std::cout << "bad resnorm " << rnorm << " >= precision " << precision << " * 1+rhsnorm " << 1 + rhsnorm << " = "
                   << precision * (1.0 + rhsnorm) << "\n";
         if (thresholdPivoting() >= threshold_pivoting_max) {
            if (print_level >= ooqp_print_level_warnings) {
               std::cout << "WARNING MA27 " << name
                         << ": threshold_pivoting parameter is already at its max and iterative refinement steps are exceeded with unsifficient precision"
                         << "\n";
               std::cout << " did not converge but still keeping the iterate" << "\n";
            }
         }
         else {
            setThresholdPivoting(std::min(thresholdPivoting() * threshold_pivoting_factor, threshold_pivoting_max));

            if (print_level >= ooqp_print_level_warnings)
               std::cout << "STATUS MA27 " << name << ": Setting ThresholdPivoting parameter to " << thresholdPivoting()
                         << " refactorization suggested\n";
//
//            done = false;
//            n_iter_ref = 0;
//
//            residual->copyFrom(rhs);
//            iter->setToZero();
//            matrixChanged();
         }

      }
   }

   rnorm = best_resid;

   if (rnorm >= precision * (1.0 + rhsnorm)) {
//      std::cout << "ERROR " << rnorm/(1.0 + rhsnorm) << " > " << precision << " (after " << n_iter_ref << " iter refs)\n";
//      rhs_cpy->write_to_streamAll(std::cout);
//      std::cout << "ERROR " << rnorm << " vs " << precision * (1.0 + rhsnorm) << " required \n";
//      best_iter->write_to_streamAll(std::cout);
//      mat->write_to_streamDense(std::cout);
//      assert( false );
//      std::cout << "Writing K of local schur complement computation...\n";
//      std::ofstream myfile("../test.prb");
//
//      myfile << "n: " << n << "\n";
//      myfile << "nnz: " << nnz << "\n";
//
//      myfile << "ia: ";
//      for( int i = 0; i <= n; i++ )
//         myfile << mat->getStorage().krowM[i] << ", ";
//      myfile << "\n";
//
//      myfile << "ja: ";
//      for( int i = 0; i < nnz; i++ )
//         myfile << mat->getStorage().jcolM[i] << ", ";
//      myfile << "\n";
//
//      myfile << "a: ";
//      for( int i = 0; i < nnz; i++ )
//         myfile << mat->getStorage().M[i] << ", ";
//      myfile << "\n";
//
//      myfile.close();
//
//      std::cout << "Writing rhs from local schur complement computation...\n";
//      myfile.open("../test.rhs");
//
//      std::cout << "sizerhs " << size_t(1) * size_t(n) <<  "\n";
//
//      myfile << "nrhs: " << 1 << "\n";
//
//      myfile << "rhs: ";
//      for( int i = 0; i < n; i++ )
//         myfile << (*rhs_cpy)[i] << ", ";
//      myfile << "\n";
//
//      myfile.close();
//
//      assert(false);
   }

   rhs.copyFrom(best_iter_loc);

#ifndef NDEBUG
   auto pos = std::find_if(rhs.elements(), rhs.elements() + rhs.length(), [](double el) { return std::fabs(el) > 1e50; });
   if (pos != rhs.elements() + rhs.length()) {
      std::cout << *pos << "\n";
      assert(false && "Big entry in solution vector... ");

   }
#endif

   /* sparsify rhs */
   //std::transform( rhs.elements(), rhs.elements() + n, rhs.elements(), []( double el ){
   //   return std::fabs(el) < 1e-16 ? 0 : el;
   //});
}

void Ma27Solver::copyMatrixElements(std::vector<double>& afact, int lafact) const {
   assert(lafact >= nnz);
   const double* M = mat_storage->M;
   std::copy(M, M + nnz, afact.begin());

   if (lafact > nnz)
      std::fill(afact.begin() + nnz, afact.begin() + (lafact - nnz), 0.0);
}

// TODO same as the one in MA57 - move somewhere else, some common MA_Solver thing maybe..
void Ma27Solver::getIndices(std::vector<int>& irowM, std::vector<int>& jcolM) const {
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

Ma27Solver::~Ma27Solver() {
   delete scaler;
   freeWorkingArrays();
}

void Ma27Solver::freeWorkingArrays() {
   irowM.clear();
   jcolM.clear();
   fact.clear();

   delete[] ikeep;
   delete[] iw;
   delete[] iw1;
   delete[] ww;

   ikeep = iw = iw1 = nullptr;
   ww = nullptr;
}

bool Ma27Solver::checkErrorsAndReact() {
   bool error = false;

   const int error_flag = info[0];
   const int error_info = info[1];

   switch (error_flag) {
      case 0 :
         break;
      case -1 : {
         std::cout << "ERROR MA27 " << name << ": N out of range or < -1: " << n << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -2 : {
         std::cout << "ERROR MA27 " << name << ": NNZ out of range or < -1 : " << nnz << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case -3 : {
         if (print_level >= ooqp_print_level_warnings)
            std::cout << "WARNING MA27 " << name << ": insufficient space in iw: " << liw << " suggest reset to " << error_info << "\n";
         ipessimism *= 1.1;

         assert(iw);
         delete[] iw;

         liw = std::max(static_cast<int>(ipessimism * error_info), static_cast<int>(ipessimism * liw));
         iw = new int[liw];
         if (print_level >= ooqp_print_level_warnings)
            std::cout << " resetting to " << liw << "\n";

         error = true;
      };
         break;
      case -4 : {
         if (print_level >= ooqp_print_level_warnings)
            std::cout << "WARNING MA27 " << name << ": insufficient factorization space: " << la << "\n";;
         rpessimism *= 1.1;

         la = std::max(static_cast<int>(rpessimism * error_info), static_cast<int>(rpessimism * la));
         fact.resize(la);

         this->copyMatrixElements(fact, la);
         if (print_level >= ooqp_print_level_warnings)
            std::cout << " resetting to " << la << "\n";

         error = true;
      }
         break;
      case -5: {
         if (print_level >= ooqp_print_level_warnings)
            std::cout << "WARNING MA27 " << name << ": matrix apparently numerically singular, detected at stage " << error_info << "\n";

         if (getSmallPivot() <= threshold_pivtol) {
            std::cout << " cannot decrease pivtol anymore -- accepting factorization anyway" << "\n";
            assert(getSmallPivot() == threshold_pivtol);
         }
         else {
            const double curr_pivtol = getSmallPivot();
            const double new_pivtol = std::max(threshold_pivtol, curr_pivtol * threshold_pivtol_factor);
            std::cout << " decreasing pivtol from " << curr_pivtol << " to " << new_pivtol << "\n";

            setSmallPivot(new_pivtol);
         }
      };
         break;
      case -6: {
         if (print_level >= ooqp_print_level_warnings)
            std::cout << "WARNING MA27 " << name << ": change of sign of pivots detected at stage " << error_info << "\n";
      };
         break;
      case -7: {
         std::cerr << "ERROR MA27 " << name << ": value of NSTEPS out of range " << nsteps << " (should not happen..) " << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      };
         break;
      case 1 : {
         if (print_level >= ooqp_print_level_warnings)
            std::cout << "WARNING MA27 " << name << ": detected " << error_info << " entries out of range in irowM and jcolM; ignored" << "\n";
      };
         break;
      case 2 : {
         std::cerr << "ERROR MA27 " << name << ": change of sign in pivots detected when matrix is supposedly definite" << "\n";
         MPI_Abort(MPI_COMM_WORLD, -1);
      }
         break;
      case 3: {
         if (print_level >= ooqp_print_level_warnings)
            std::cout << "WARNING MA27 " << name << ": rank deficient matrix detected; apparent rank is " << error_info << " != n : " << this->n
                      << "\n";
         error = false;

      };
         break;
      default : {
         assert(0 == error_flag);
      };
         break;
   }

   return error;
}

std::tuple<unsigned int, unsigned int, unsigned int> Ma27Solver::get_inertia() const {
   assert(false && "TODO: Implement");
   return {0, 0, 0};
}
