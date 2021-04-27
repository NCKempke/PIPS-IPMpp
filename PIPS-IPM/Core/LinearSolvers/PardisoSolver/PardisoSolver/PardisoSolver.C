/* PIPS-IPM                                                             
 * Authors: Cosmin G. Petra, Miles Lubin, Murat Mut
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#include <iostream>
#include <algorithm>

#include "PardisoSolver.h"
#include "SparseStorage.h"
#include "SparseSymMatrix.h"
#include "SimpleVector.h"
#include "DenseGenMatrix.h"
#include "pipsdef.h"
#include <cstdlib>

#include "mpi.h"

PardisoSolver::PardisoSolver(const SparseSymMatrix* sgm) : Msys{sgm}, n{static_cast<int>(sgm->size())}, nnz{sgm->numberOfNonZeros()} {
   krowM = new int[n + 1];
   jcolM = new int[nnz];
   M = new double[nnz];

   nvec = new double[n];
}


PardisoSolver::PardisoSolver(const DenseSymMatrix* m) : Mdsys{m}, n{static_cast<int>(Mdsys->size())} {
   assert(m->size() < std::numeric_limits<int>::max());
   nvec = new double[n];
}

bool PardisoSolver::iparmUnchanged() const {
   bool same = true;
   bool print = true;

   std::vector<int> iparm_compare(64, 0);
   setIparm(iparm_compare.data());

   static const std::vector<int> to_compare{1, 7, 10, 12};

   for (int i = 0; i < 64; ++i) {
      // if entry should be compared
      if (std::find(to_compare.begin(), to_compare.end(), i) != to_compare.end()) {
         if (iparm[i] != iparm_compare[i]) {
            if (print)
               std::cout << "ERROR - PardisoSolver: elements in iparm changed at " << i << ": " << iparm[i] << " != " << iparm_compare[i] << "(new)"
                         << std::endl;
            same = false;
         }
      }

   }
   return same;
}

void PardisoSolver::initSystem() {
   if (Msys) {
      assert(Msys->isLower);

      //get the matrix in upper triangular
      Msys->getStorageRef().transpose(krowM, jcolM, M);

      //save the indices for diagonal entries for a streamlined later update
      int* krowMsys = Msys->getStorageRef().krowM;
      int* jcolMsys = Msys->getStorageRef().jcolM;
      for (int r = 0; r < n; r++) {
         // Msys - find the index in jcol for the diagonal (r,r)
         int idxDiagMsys = -1;
         for (int idx = krowMsys[r]; idx < krowMsys[r + 1]; idx++)
            if (jcolMsys[idx] == r) {
               idxDiagMsys = idx;
               break;
            }
         assert(idxDiagMsys >= 0);
         //must have all diagonals entries

         // aug - find the index in jcol for the diagonal (r,r)
         int idxDiagAug = -1;
         for (int idx = krowM[r]; idx < krowM[r + 1]; idx++)
            if (jcolM[idx] == r) {
               idxDiagAug = idx;
               break;
            }
         assert(idxDiagAug >= 0);

         diagMap.insert(std::pair<int, int>(idxDiagMsys, idxDiagAug));
      }
   }
   else if (Mdsys) {
      // the input is a dense matrix
      // the dense matrix is also processed in matrixChanged everytime the method is called

      // the input is a dense matrix
      const DenseSymMatrix& Md = (*Mdsys);
      nnz = Md.getNumberOfNonZeros();

      delete[] krowM;
      delete[] jcolM;
      delete[] M;
      krowM = new int[n + 1];
      jcolM = new int[nnz];
      M = new double[nnz];

      int nzIt = 0;
      for (int i = 0; i < n; i++) {

         krowM[i] = nzIt;
         //cout << i << " " << krowM[i] << endl;

         jcolM[nzIt] = i; // the diagonal
         M[nzIt] = Md[i][i];
         assert(nzIt < nnz);
         nzIt++;

         for (int j = i + 1; j < n; j++) {
            if (Md[i][j] != 0.0) {
               jcolM[nzIt] = j;
               M[nzIt] = Md[i][j];
               assert(nzIt < nnz);
               nzIt++;
            }
         }
      }
      //cout << "PardisoSolver::first call nzit=" << nzIt << endl;
      assert(nzIt == nnz);
      krowM[n] = nzIt;

   }
   else {
      assert(false);
   }

   // need Fortran indexes
   // extra method?
   for (int i = 0; i < n + 1; i++)
      krowM[i] += 1;
   for (int i = 0; i < nnz; i++)
      jcolM[i] += 1;

}

void PardisoSolver::diagonalChanged(int /* idiag */, int /* extent */) {
   this->matrixChanged();
}

void PardisoSolver::matrixChanged() {
   if (first) {
      firstCall();
      first = false;
   }
   else {
      if (Msys) {
         //update diagonal entries in the PARDISO aug sys (if the input is sparse)
         double* eltsMsys = Msys->getStorageRef().M;
         std::map<int, int>::iterator it;
         for (it = diagMap.begin(); it != diagMap.end(); it++)
            M[it->second] = eltsMsys[it->first];
      }
   }

   if (Mdsys) {
      // the input is a dense matrix
      const DenseSymMatrix& Md = (*Mdsys);
      //double tm=MPI_Wtime();
      int nzIt = 0;
      for (int i = 0; i < n; i++) {

         krowM[i] = nzIt;

         jcolM[nzIt] = i; // the diagonal
         M[nzIt] = Md[i][i];
         nzIt++;

         for (int j = i + 1; j < n; j++) {
            if (Md[i][j] != 0.0) {
               jcolM[nzIt] = j;
               M[nzIt] = Md[i][j];
               nzIt++;
            }
         }
      }
      //cout << "PardisoSolver::matrix changed nnzit=" << nzIt << endl;
      krowM[n] = nzIt;
      assert(nzIt == nnz);
      // need Fortran indexes
      for (int i = 0; i < n + 1; i++)
         krowM[i] += 1;
      for (int i = 0; i < nnz; i++)
         jcolM[i] += 1;
      //cout << "Forming the matrix took:" << MPI_Wtime()-tm << endl;
   }


   phase = 12; // Analysis, numerical factorization
   nrhs = 1;

   iparm[30] = 0; // do not specify sparse rhs
   // if one wants to use the sparse rhs and partial solves according to
   // the PARDISO user guide the perm vector has to be present during all
   // phases of pardiso

   assert(iparmUnchanged());

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

   if (error != 0) {
      printf("PardisoSolver - ERROR during factorization: %d\n", error);
      MPI_Abort(MPI_COMM_WORLD, -1);
   }
}

void PardisoSolver::solve(OoqpVector& rhs_in) {
   SimpleVector<double>& rhs = dynamic_cast<SimpleVector<double>&>(rhs_in);
   double* sol_local = nvec;

   //int maxRefinSteps=(gLackOfAccuracy==0?3:6);

   /* same for mkl_pardiso and pardiso */
   phase = 33; // solve and iterative refinement
   nrhs = 1;
   iparm[30] = 0; // do not specify sparse rhs at this point

   assert(iparmUnchanged());

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr, &nrhs, iparm, &msglvl, rhs.elements(), sol_local, &error);

   if (error != 0) {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }

   rhs.copyFromArray(sol_local);
}

void PardisoSolver::solve(GenMatrix& rhs_in) {
   DenseGenMatrix& rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

   int nrows, ncols;
   rhs.getSize(ncols, nrows);
   assert(nrows == n);

   if (static_cast<int>(sol.size()) < nrows * ncols)
      sol.resize(nrows * ncols);

   phase = 33; // solve and iterative refinement
   nrhs = ncols;
   iparm[30] = 0;
   assert(iparmUnchanged());

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr, &nrhs, iparm, &msglvl, &rhs[0][0], sol.data(), &error);

   if (error != 0)
      printf("PardisoSolver - ERROR during solve: %d", error);

   std::copy(&rhs[0][0], &rhs[0][0] + nrows * ncols, sol.begin());
}

void PardisoSolver::solve(GenMatrix& rhs_in, int* colSparsity) {
   std::cout << "PardisoSolver - using sparse rhs but might lead to numerical troubles .. \n";
   DenseGenMatrix& rhs = dynamic_cast<DenseGenMatrix&>(rhs_in);

   int nrows, ncols;
   rhs.getSize(ncols, nrows);
   assert(nrows == n);

   if (static_cast<int>(sol.size()) < nrows * ncols)
      sol.resize(nrows * ncols);

   /* same for mkl_pardiso and pardiso */
   phase = 33; // solve and iterative refinement
   iparm[30] = 1;

   nrhs = ncols;

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, colSparsity, &nrhs, iparm, &msglvl, &rhs[0][0], sol.data(), &error);

   if (error != 0)
      printf("PardisoSolver - ERROR during solve: %d", error);

   std::copy(&rhs[0][0], &rhs[0][0] + nrows * ncols, sol.begin());
}

void PardisoSolver::solve(int nrhss, double* rhss, int* /*colSparsity*/ ) {
   assert(rhss);
   assert(nrhss >= 1);

   if (static_cast<int>(sol.size()) < nrhss * n)
      sol.resize(nrhss * n);

   if (static_cast<int>(rhss_nonzero.size()) < nrhss * n) {
      rhss_nonzero.resize(nrhss * n);
      map_rhs_nonzero_original.resize(nrhss);
   }

   int nrhss_nnz{0};
   for (int i = 0; i < nrhss; ++i) {
      if (std::any_of(rhss + i * n, rhss + (i + 1) * n, [](double d) { return !PIPSisZero(d); })) {
         std::copy(rhss + i * n, rhss + (i + 1) * n, rhss_nonzero.begin() + nrhss_nnz * n);
         map_rhs_nonzero_original[nrhss_nnz] = i;
         ++nrhss_nnz;
      }
   }

#ifndef NDEBUG
//   if( colSparsity )
//   {]
//      for( int nr = 0; nr < nrhss; nr++ )
//      {
//         for( int i = 0; i < n; i++ )
//         {
//            const int rhspos = nr * n + i;
//            if( rhss[rhspos] != 0.0 )
//               assert(colSparsity[i] == 1);
//            else if( nrhss == 1 ) // does not work with zeroes in matrix, e.g. callback example
//               assert(colSparsity[i] == 0);
//         }
//      }
//   }
   // PARDISO cannot deal well with all zero rhs
   for (int i = 0; i < nrhss_nnz; ++i)
      assert(std::any_of(rhss_nonzero.begin() + i * n, rhss_nonzero.begin() + (i + 1) * n, [](double d) {
         return d != 0.0;
      }));
#endif

   /* same for mkl_pardiso and pardiso */
   phase = 33; // solve and iterative refinement
   int nrhss_local = nrhss;

   assert(iparmUnchanged());

// see notes on [30] earlier - cannot be specified on the go
// seems to detoriate parformance a lot - triggers many BiCGStab steps
//   if( colSparsity )
//   {
//      iparm[30] = 1; //sparse rhs
//   }
//   else
   {
      iparm[30] = 0;
   }

   assert(pt);
   assert(M);
   assert(krowM);
   assert(jcolM);
   assert(rhss);

   pardisoCall(pt, &maxfct, &mnum, &mtype, &phase, &n, M, krowM, jcolM, nullptr, &nrhss_local, iparm, &msglvl, rhss_nonzero.data(), sol.data(),
         &error);

   if (error != 0) {
      printf("PardisoSolver - ERROR during solve: %d", error);
      exit(1);
   }

   for (int i = 0; i < nrhss_nnz; ++i)
      std::copy(sol.begin() + i * n, sol.begin() + (i + 1) * n, rhss + map_rhs_nonzero_original[i] * n);
}

PardisoSolver::~PardisoSolver() {
   delete[] jcolM;
   delete[] krowM;
   delete[] M;
   delete[] nvec;
}
