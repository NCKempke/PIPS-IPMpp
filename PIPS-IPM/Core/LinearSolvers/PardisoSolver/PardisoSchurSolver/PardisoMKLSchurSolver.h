/*
 * PardisoMKLSchurSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: Nils-Christian Kempke
 */

#ifndef PARDISO_MKL_SCHUR_LINSYS_H
#define PARDISO_MKL_SCHUR_LINSYS_H

#include "PardisoSchurSolver.h"

class PardisoMKLSchurSolver : public PardisoSchurSolver
{
public:
   PardisoMKLSchurSolver( const SparseSymMatrix * sgm );
   void solve( OoqpVector& rhs ) override;
   using DoubleLinearSolver::solve;

 protected:
  ~PardisoMKLSchurSolver() override;

  void computeSC(int nSCO, /*const*/SparseGenMatrix& R,
     /*const*/SparseGenMatrix& A,
     /*const*/SparseGenMatrix& C,
     /*const*/SparseGenMatrix& F,
     /*const*/SparseGenMatrix& G, int*& rowptrSC, int*& colidxSC, double*& eltsSC) override;

  void setIparm(int* iparm) const override;
  void initPardiso() override;
};

#endif /* PARDISO_MKL_SCHUR_LINSYS_H */
