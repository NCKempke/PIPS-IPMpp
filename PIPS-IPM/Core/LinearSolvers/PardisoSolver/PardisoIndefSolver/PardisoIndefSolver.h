/*
 * PardisoIndefSolver.h
 *
 *  Created on: 21.03.2018
 *      Author: Daniel Rehfeldt
 */

#ifndef _PARDISOINDEFSOLVER_H_
#define _PARDISOINDEFSOLVER_H_

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "DenseSymMatrix.h"
#include "SparseSymMatrix.h"
#include "DenseStorageHandle.h"
#include "pipsport.h"

class PardisoIndefSolver : public DoubleLinearSolver
{
   public:
      DenseStorageHandle mStorage;
      SparseStorageHandle mStorageSparse;
   protected:

      constexpr static double precondDiagDomBound = 0.0001;
      constexpr static int pivotPerturbationExpDefault = 8;
      constexpr static int nIterativeRefinsDefault = 8;
      constexpr static bool highAccuracyDefault = true;
      constexpr static bool useSparseRhsDefault = true;
      constexpr static bool parallelForwardBackwardDefault = true;
      constexpr static bool factorizationTwoLevelDefault = true;


      double* x{}; /* solution vector */

      int mtype{-2};


      int n; /* size of the matrix */

      int nrhs{1}; /* Number of right hand sides. */

      void *pt[64];  /* Internal solver memory pointer pt */

      /* Pardiso control parameters. */
      int iparm[64];

      int maxfct{1};
      int mnum{1};
      int phase{11};
      int msglvl{0};
      int* ia{};
      int* ja{};
      double* a{};
      int idum{-1};
      double ddum{-1.0};
      int pivotPerturbationExp; // 10^-exp
      int nIterativeRefins;
      bool highAccuracy{highAccuracyDefault};
      bool parallelForwardBackward {parallelForwardBackwardDefault};
      bool factorizationTwoLevel{factorizationTwoLevelDefault};
      bool useSparseRhs{useSparseRhsDefault};
      bool deleteCSRpointers{false};
      bool solve_in_parallel;


      virtual void pardisoCall(void *pt, int* maxfct, int* mnum, int* mtype, int* phase, int* n, double* M, int* krowM, int* jcolM,
            int* perm, int* nrhs, int* iparm, int* msglvl, double* rhs, double* sol, int* error) = 0;
      virtual void checkMatrix() = 0;
      virtual void getIparm( int* iparm ) const = 0;
   public:
      PardisoIndefSolver(DenseSymMatrix * storage, bool solve_in_parallel);
      PardisoIndefSolver(SparseSymMatrix * storage, bool solve_in_parallel);
      void diagonalChanged(int idiag, int extent) override;
      void matrixChanged() override;
      void matrixRebuild( DoubleMatrix& matrixNew ) override;

      using DoubleLinearSolver::solve;
      void solve ( OoqpVector& vec ) override;
      void solve ( GenMatrix& vec ) override;
      void solveSynchronized( OoqpVector& vec ) override;
      ~PardisoIndefSolver() override;

   private:
      void initPardiso();
      void factorizeFromSparse();
      void factorizeFromSparse(SparseSymMatrix& matrix_fortran);
      void factorizeFromDense();
      void factorize();

      bool iparmUnchanged();

};

#endif /* _PARDISOINDEFSOLVER_H_ */
