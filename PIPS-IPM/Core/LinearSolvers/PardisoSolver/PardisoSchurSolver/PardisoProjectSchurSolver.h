/*
 * PardisoProjectSchurSolver.h
 *
 *  Created on: 10.12.2020
 *      Author: bzfkempk
 */

#ifndef PARDISO_PROJECT_SCHUR_LINSYS_H
#define PARDISO_PROJECT_SCHUR_LINSYS_H

#include "PardisoSchurSolver.h"

class PardisoProjectSchurSolver : public PardisoSchurSolver {
public:
   explicit PardisoProjectSchurSolver(const SparseSymmetricMatrix& sgm);

   void solve(Vector<double>& rhs_in) override;
   using DoubleLinearSolver::solve;

protected:
   ~PardisoProjectSchurSolver() override;

   void computeSC(int nSCO, const SparseMatrix& R, const SparseMatrix& A, const SparseMatrix& C, const SparseMatrix& F,
         const SparseMatrix& G, int*& rowptrSC, int*& colidxSC, double*& eltsSC) override;
   void initPardiso() override;
   void setIparm(int* iparm) const override;

   double dparm[64]{};

   int solver{0};
   int num_threads{1};

};

#endif /* PARDISO_PROJECT_SCHUR_LINSYS_H */
