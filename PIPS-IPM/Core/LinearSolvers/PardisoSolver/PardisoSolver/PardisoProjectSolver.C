/*
 * PardisoProjectSolver.C
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#include "PardisoProjectSolver.h"

extern "C" void pardisoinit(void*, const int*, int*, int*, double*, int*);
extern "C" void
pardiso(void*, const int*, const int*, const int*, const int*, int*, double*, int*, int*, int*, int*, int*, const int*, double*, double*, int*,
      double*);

extern "C" void pardiso_chkmatrix(int*, int*, double*, int*, int*, int*);
extern "C" void pardiso_chkvec(int*, int*, double*, int*);
extern "C" void pardiso_printstats(int*, int*, double*, int*, int*, int*, double*, int*);

PardisoProjectSolver::PardisoProjectSolver(const SparseSymmetricMatrix& sgm) : PardisoSolver(sgm) {
#ifdef TIMING
   if( PIPS_MPIgetRank() == 0 )
      std::cout << "PardisoProjectSolver::PardisoProjectSolver (sparse input)\n";
#endif

   num_threads = PIPSgetnOMPthreads();
   solver = 0; /* sparse direct solver */
}

PardisoProjectSolver::PardisoProjectSolver(const DenseSymmetricMatrix& m) : PardisoSolver(m) {
#ifdef TIMING
   if( myRank == PIPS_MPIgetRank() )
     std::cout << "PardisoProjectSolver created (dense input)\n";
#endif

   num_threads = PIPSgetnOMPthreads();
   solver = 0; /* sparse direct solver */
}

void PardisoProjectSolver::firstCall() {
   iparm[0] = 0;

   int error = 0;

   // the licence file read seems to be critical.. for PARDISO 6.0 and 6.2 - not 7.0/7.2 anymore?
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

   initSystem();
}


/*
 * iparm[30] has to be set depending on the circumstances!
 */
void PardisoProjectSolver::setIparm(int* iparm) const {
   iparm[1] = 2; // 2 is for metis, 0 for min degree
   iparm[2] = num_threads;

   iparm[7] = 2; // max number of iterative refinements
   iparm[10] = 1; // scaling for IPM KKT; used with IPARM(13)=1 or 2
   iparm[12] = 2; // improved accuracy for IPM KKT; used with IPARM(11)=1;
   // if needed, use 2 for advanced matchings and higher accuracy.
   iparm[18] = 0; /* don't compute GFLOPS */
}

void
PardisoProjectSolver::pardisoCall(void* pt, const int* maxfct, const int* mnum, const int* mtype, const int* phase, int* n, double* M, int* krowM,
      int* jcolM, int* perm, int* nrhs, int* iparm, const int* msglvl, double* rhs, double* sol, int* error) {
   pardiso(pt, maxfct, mnum, mtype, phase, n, M, krowM, jcolM, perm, nrhs, iparm, msglvl, rhs, sol, error, dparm);
}

PardisoProjectSolver::~PardisoProjectSolver() {
   phase = -1; // release internal memory
   nrhs = 1;

   pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, nullptr, krowM, jcolM, nullptr, &nrhs, iparm, &msglvl, nullptr, nullptr, &error, dparm);

   if (error != 0) {
      printf("PardisoSolver - ERROR in pardiso release: %d", error);
   }
}
