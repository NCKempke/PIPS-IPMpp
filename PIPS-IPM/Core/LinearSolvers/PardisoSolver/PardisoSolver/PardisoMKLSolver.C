/*
 * PardisoMKLSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#include "PardisoMKLSolver.h"

#include "pipsdef.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"

PardisoMKLSolver::PardisoMKLSolver(const SparseSymmetricMatrix* sgm) : PardisoSolver(sgm) {
#ifdef TIMING
   if( PIPS_MPIgetRank() == 0 )
     std::cout << "PardisoMKLSolver::PardisoMKLSolver (sparse input)\n";
#endif
}

PardisoMKLSolver::PardisoMKLSolver(const DenseSymmetricMatrix* m) : PardisoSolver(m) {
#ifdef TIMING
   if( PIPS_MPIgetRank() == 0 )
     std::cout << "PardisoMKLSolver::PardisoMKLSolver (sparse input)\n";
#endif
}

void PardisoMKLSolver::firstCall() {
   iparm[0] = 0;
   pardisoinit(pt, &mtype, iparm);

   setIparm(iparm);

   initSystem();
}


/*
 * iparm[30] has to be set depending on the circumstances!
 */
void PardisoMKLSolver::setIparm(int* iparm) const {

   iparm[1] = 2; // 2 is for metis, 0 for min degree

   /* From INTEL (instead of iparm[2] which is not defined there):
   *  You can control the parallel execution of the solver by explicitly setting the MKL_NUM_THREADS environment variable.
   *  If fewer OpenMP threads are available than specified, the execution may slow down instead of speeding up.
   *  If MKL_NUM_THREADS is not defined, then the solver uses all available processors.
   */
   iparm[7] = 0; // mkl runs into numerical problems when setting iparm[7] too high
   iparm[10] = 1; // default, scaling for IPM KKT used with either mtype=11/13 or mtype=-2/-4/6 and iparm[12]=1
   iparm[12] = 1; // 0 disable matching, 1 enable matching, no other settings

#ifndef NDEBUG
   // enable matrix checker - mkl pardiso has no chkmatrix method
   iparm[26] = 1;
#endif
}

void PardisoMKLSolver::pardisoCall(void* pt, const int* maxfct, const int* mnum, const int* mtype, const int* phase, int* n, double* M, int* krowM,
      int* jcolM, int* perm, int* nrhs, int* iparm, const int* msglvl, double* rhs, double* sol, int* error) {
   pardiso(pt, maxfct, mnum, mtype, phase, n, M, krowM, jcolM, perm, nrhs, iparm, msglvl, rhs, sol, error);
}

PardisoMKLSolver::~PardisoMKLSolver() {
   phase = -1; // release internal memory
   nrhs = 1;

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, nullptr, krowM, jcolM, nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error);

   if (error != 0)
      printf("PardisoMKLSolver - ERROR in pardiso release: %d", error);
}
