/*
 * PardisoMKLSchurSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: Nils-Christian Kempke
 */

#ifndef PARDISO_MKL_SCHUR_LINSYS_H
#define PARDISO_MKL_SCHUR_LINSYS_H

#include "PardisoSchurSolver.h"

class PardisoMKLSchurSolver : public PardisoSchurSolver {
public:
   PardisoMKLSchurSolver(const SparseSymmetricMatrix* sgm);

   void solve(Vector<double>& rhs) override;

   using DoubleLinearSolver::solve;

protected:
   ~PardisoMKLSchurSolver() override;

   void computeSC(int nSCO, const SparseMatrix& R, const SparseMatrix& A, const SparseMatrix& C, const SparseMatrix& F,
      const SparseMatrix& G, int*& rowptrSC, int*& colidxSC, double*& eltsSC) override;

   void setIparm(int* iparm) const override;

   void initPardiso() override;
};

#endif /* PARDISO_MKL_SCHUR_LINSYS_H */
