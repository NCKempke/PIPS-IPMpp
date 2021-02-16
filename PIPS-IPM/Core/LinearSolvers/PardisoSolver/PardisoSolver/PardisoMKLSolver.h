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

   PardisoMKLSolver( SparseSymMatrix * sgm );
   PardisoMKLSolver( DenseSymMatrix* m );


protected:
   void setIparm(int* iparm) const override;

   void pardisoCall(void *pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n, double* M, int* krowM, int* jcolM,
         int* perm, int* nrhs, int* iparm, int* msglvl, double* rhs, double* sol, int* error) override;
  ~PardisoMKLSolver() override;
};

#endif /* PARDISO_MKL_LINSYS_H */
