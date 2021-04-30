/*
 * PardisoProjectIndefSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#ifndef PARDISO_PROJECT_INDEF_LINSYS_H
#define PARDISO_PROJECT_INDEF_LINSYS_H

#include "PardisoIndefSolver.h"

class PardisoProjectIndefSolver : public PardisoIndefSolver {
public:
   PardisoProjectIndefSolver(SparseSymMatrix* sgm, bool solve_in_parallel, MPI_Comm mpi_comm);

   PardisoProjectIndefSolver(DenseSymMatrix* m, bool solve_in_parallel, MPI_Comm mpi_comm);

   ~PardisoProjectIndefSolver() override;

protected:
   void
   pardisoCall(void* pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n, double* M, int* krowM, int* jcolM,
      int* perm, int* nrhs, int* iparm,
      int* msglvl, double* rhs, double* sol, int* error) override;

   void checkMatrix() override;

   void getIparm(int* iparm) const override;

private:
   void initPardiso();

   void setIparm(int* iparm) const;

   double dparm[64];

   int solver{-1};
   int num_threads{1};

};

#endif /* PARDISO_PROJECT_INDEF_LINSYS_H */
