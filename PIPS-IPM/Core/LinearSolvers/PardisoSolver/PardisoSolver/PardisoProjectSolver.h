/*
 * PardisoProjectSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#ifndef PARDISO_PROJECT_LINSYS_H
#define PARDISO_PROJECT_LINSYS_H

#include "PardisoSolver.h"

#include <map>


class PardisoProjectSolver : public PardisoSolver
{
public:
  void firstCall() override;

  PardisoProjectSolver( const SparseSymMatrix * sgm );
  PardisoProjectSolver( const DenseSymMatrix* m);

 protected:
  void setIparm(int* iparm) const override;

  void pardisoCall(void *pt, const int* maxfct, const int* mnum, const int* mtype, const int* phase, int* n, double* M, int* krowM, int* jcolM,
        int* perm, int* nrhs, int* iparm, const int* msglvl, double* rhs, double* sol, int* error) override;

  ~PardisoProjectSolver() override;
 private:
  double dparm[64];

  int solver;
  int num_threads;

};

#endif /* PARDISO_PROJECT_LINSYS_H */
