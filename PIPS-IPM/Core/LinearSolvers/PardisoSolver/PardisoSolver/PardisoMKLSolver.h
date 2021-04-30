/*
 * PardisoMKLSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#ifndef PARDISO_MKL_LINSYS_H
#define PARDISO_MKL_LINSYS_H

#include "PardisoSolver.h"
#include "pipsport.h"

class PardisoMKLSolver : public PardisoSolver {

public:
   void firstCall() override;

   PardisoMKLSolver(const SparseSymMatrix* sgm);
   PardisoMKLSolver(const DenseSymMatrix* m);


protected:
   void setIparm(int* iparm) const override;

   void
   pardisoCall(void* pt, const int* maxfct, const int* mnum, const int* mtype, const int* phase, int* n, double* M, int* krowM, int* jcolM, int* perm,
         int* nrhs, int* iparm, const int* msglvl, double* rhs, double* sol, int* error) override;
   ~PardisoMKLSolver() override;
};

#endif /* PARDISO_MKL_LINSYS_H */
