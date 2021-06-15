/*
 * PardisoMKLIndefSolver.C
 *
 *  Created on: 11.12.2020
 *      Author: Nils-Christian Kempke
 */
#include "PardisoMKLIndefSolver.h"

#include "pipschecks.h"
#include "SimpleVector.hpp"
#include "pipsdef.h"
#include "PIPSIPMppOptions.h"

#include "mkl_pardiso.h"
#include "mkl_types.h"

PardisoMKLIndefSolver::PardisoMKLIndefSolver(DenseSymmetricMatrix* dm, bool solve_in_parallel, MPI_Comm mpi_comm) : PardisoIndefSolver(dm,
      solve_in_parallel, mpi_comm) {
   initPardiso();
}


PardisoMKLIndefSolver::PardisoMKLIndefSolver(SparseSymmetricMatrix* sm, bool solve_in_parallel, MPI_Comm mpi_comm) : PardisoIndefSolver(sm,
      solve_in_parallel, mpi_comm) {
   initPardiso();
}

void PardisoMKLIndefSolver::setIparm(int* iparm) const {

   iparm[9] = 13; /* pivot perturbation 10^{-xxx} */

   iparm[7] = 0; /* mkl runs into numerical problems when setting iparm[7] too high */

   /* enable matrix checker (default disabled) - mkl pardiso does not have chkmatrix */
#ifndef NDEBUG
   iparm[26] = 1;
#endif

   if (highAccuracy) {
      iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
      iparm[12] = 1;// MKL does not have a =2 option here
   }
   else {
      iparm[10] = 0;
      iparm[12] = 0;
   }
}

void PardisoMKLIndefSolver::initPardiso() {
   pardisoinit(pt, &mtype, iparm);

   setIparm(iparm);
}

void
PardisoMKLIndefSolver::pardisoCall(void* pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n, double* M, int* krowM, int* jcolM, int* perm,
      int* nrhs, int* iparm, int* msglvl, double* rhs, double* sol, int* error) {
   pardiso(pt, maxfct, mnum, mtype, phase, n, M, krowM, jcolM, perm, nrhs, iparm, msglvl, rhs, sol, error);
}

void PardisoMKLIndefSolver::checkMatrix() {
   return;
}

void PardisoMKLIndefSolver::getIparm(int* iparm) const {
   setIparm(iparm);
}

PardisoMKLIndefSolver::~PardisoMKLIndefSolver() {
   phase = -1; /* Release internal memory. */
   int error;
   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);

}
