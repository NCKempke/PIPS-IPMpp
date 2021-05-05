/*
 * PardisoProjectIndefSolver.C
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#include "PardisoProjectIndefSolver.h"
#include "pipsdef.h"

extern "C" void pardisoinit(void*, int*, int*, int*, double*, int*);
extern "C" void pardiso(void*, int*, int*, int*, int*, int*, double*, int*, int*, int*, int*, int*, int*, double*, double*, int*, double*);

extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats(int*, int*, double*, int*, int*, int*, double*, int*);

PardisoProjectIndefSolver::PardisoProjectIndefSolver(SparseSymmetricMatrix* sgm, bool solve_in_parallel, MPI_Comm mpi_comm) : PardisoIndefSolver(sgm,
      solve_in_parallel, mpi_comm) {
   assert(!sgm->isLower);

   num_threads = PIPSgetnOMPthreads();
   solver = 0; /* sparse direct solver */
   initPardiso();
}

PardisoProjectIndefSolver::PardisoProjectIndefSolver(DenseSymMatrix* m, bool solve_in_parallel, MPI_Comm mpi_comm) : PardisoIndefSolver(m,
      solve_in_parallel, mpi_comm) {
   num_threads = PIPSgetnOMPthreads();
   solver = 0; /* sparse direct solver */

   initPardiso();
}

void PardisoProjectIndefSolver::setIparm(int* iparm) const {
   iparm[9] = 13; /* pivot perturbation 10^{-xxx} */

   assert(num_threads >= 1 && pivotPerturbationExp >= 1);
   iparm[9] = pivotPerturbationExp;
   iparm[2] = num_threads;
   iparm[18] = 0; /* don't compute GFLOPS */
   iparm[7] = nIterativeRefins; /* max number of iterative refinement steps. */

   if (highAccuracy) {
      iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
      iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
   }
   else {
      iparm[10] = 1;
      iparm[12] = 0;
   }

   if (factorizationTwoLevel)
      iparm[23] = 1; // parallel Numerical Factorization (0=used in the last years, 1=two-level scheduling)
   else
      iparm[23] = 0;

   if (parallelForwardBackward)
      iparm[24] = 1; // parallelization for the forward and backward solve. 0=sequential, 1=parallel solve.
   else
      iparm[24] = 0;
}

void PardisoProjectIndefSolver::initPardiso() {
   int error = 0;
   solver = 0; /* use sparse direct solver */

//   #pragma omp critical
   pardisoinit(pt, &mtype, &solver, iparm, dparm, &error);

   if (error != 0) {
      if (error == -10)
         printf("PARDISO: No license file found \n");
      if (error == -11)
         printf("PARDISO: License is expired \n");
      if (error == -12)
         printf("PARDISO: Wrong username or hostname \n");

      PIPS_MPIabortIf(true, "Error in pardisoinit");
   }

   setIparm(iparm);
}

void PardisoProjectIndefSolver::getIparm(int* iparm) const {
   setIparm(iparm);
}

void
PardisoProjectIndefSolver::pardisoCall(void* pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n, double* M, int* krowM, int* jcolM, int* perm,
      int* nrhs, int* iparm, int* msglvl, double* rhs, double* sol, int* error) {
   pardiso(pt, maxfct, mnum, mtype, phase, n, M, krowM, jcolM, perm, nrhs, iparm, msglvl, rhs, sol, error, dparm);
}


void PardisoProjectIndefSolver::checkMatrix() {
#if !defined(NDEBUG) && defined(CHECK_PARDISO)
   pardiso_chkmatrix(&mtype, &n, a, ia, ja, &error);
   if( error != 0 )
   {
      printf("\nERROR in consistency of matrix: %d", error);
      MPI_Abort(MPI_COMM_WORLD);
   }
#endif
}

PardisoProjectIndefSolver::~PardisoProjectIndefSolver() {
   phase = -1; /* Release internal memory. */
   int error;
   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, &ddum, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);

   if (error != 0) {
      printf("PardisoIndefSolver - ERROR in pardiso release: %d", error);
   }

}
