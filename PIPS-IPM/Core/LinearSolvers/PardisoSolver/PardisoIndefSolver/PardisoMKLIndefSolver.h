/*
 * PardisoMKLIndefSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: Nils-Christian Kempke
 */

#ifndef PARDISO_MKL_INDEF_LINSYS_H
#define PARDISO_MKL_INDEF_LINSYS_H

#include "PardisoIndefSolver.h"

class PardisoMKLIndefSolver : public PardisoIndefSolver
{
public:
   PardisoMKLIndefSolver( SparseSymMatrix * sgm, bool solve_in_parallel );
   PardisoMKLIndefSolver( DenseSymMatrix* m, bool solve_in_parallel );

 protected:
  void pardisoCall(void *pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n, double* M, int* krowM, int* jcolM,
        int* perm, int* nrhs, int* iparm, int* msglvl, double* rhs, double* sol, int* error) override;
  void checkMatrix() override;
  void getIparm( int* iparm ) const override;

  ~PardisoMKLIndefSolver() override;
 private:
  void initPardiso();
  void setIparm(int* iparm) const;
};

#endif /* PARDISO_MKL_INDEF_LINSYS_H */
