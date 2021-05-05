/* PIPS-IPM
 * Authors: Cosmin G. Petra, Miles Lubin
 * (C) 2012 Argonne National Laboratory, see documentation for copyright
 */
#ifndef PARDISO_SCHUR_SOLVER
#define PARDISO_SCHUR_SOLVER

#include "DoubleLinearSolver.h"
#include "SparseSymmetricMatrix.h"
#include "SparseMatrix.h"
#include "DenseSymMatrix.h"
#include "Vector.hpp"
#include "SmartPointer.h"
#include "SparseStorage.h"
#include "pipsport.h"
#include <map>

/** implements the linear solver class using the Pardiso SC solver
 */

class PardisoSchurSolver : public DoubleLinearSolver {

   constexpr static int symbFactorIntervalDefault = 8;
   constexpr static int pivotPerturbationExpDefault = 8;
   constexpr static int nIterativeRefinsDefault = 8;
   constexpr static bool parallelForwardBackwardDefault = true;
   constexpr static bool factorizationTwoLevelDefault = true;


public:
   virtual void firstCall(); //first factorization call
   void
   firstSolveCall(SparseMatrix& R, SparseMatrix& A, SparseMatrix& C, SparseMatrix& F, SparseMatrix& G, int nSC0); //first solve call

   /** sets mStorage to refer to the argument sgm */
   PardisoSchurSolver(const SparseSymmetricMatrix* sgm);

   void diagonalChanged(int idiag, int extent) override;
   void matrixChanged() override;

   using DoubleLinearSolver::solve;
   void solve(Vector<double>& rhs) override = 0;
   void solve(GeneralMatrix&) override { assert(false && "Function not supported. Use PardisoSolver for this functionality."); };

   bool reports_inertia() const override { return false; };
   std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const override { return {0, 0, 0}; };

   /** Functions specific to the Schur approach. The last argument is the Schur first
    * stage matrix that will be updated.
    * 1. schur_solver( rhs, SC)
    *  - this is the generic function
    *
    * 2. schur_solve(R,A,B, SC)
    *  - avoids forming the matrix rhs [R' A' B']'
    *  - assumes rhs does not change
    */
   virtual void schur_solve(/*const*/ SparseMatrix& R,
         /*const*/ SparseMatrix& A,
         /*const*/ SparseMatrix& C,
         /*const*/ SparseMatrix& F,
         /*const*/ SparseMatrix& G, DenseSymMatrix& SC);

   virtual void schur_solve_sparse(/*const*/ SparseMatrix& R,
         /*const*/ SparseMatrix& A,
         /*const*/ SparseMatrix& C,
         /*const*/ SparseMatrix& F,
         /*const*/ SparseMatrix& G, SparseSymmetricMatrix& SC);

protected:

   const SparseSymmetricMatrix* Msys{}; // this is the (1,1) block in the augmented system
   bool first{true};
   bool firstSolve{true};
   void* pt[64];
   int iparm[64];
   bool useSparseRhs{true};
   int symbFactorInterval;

   int* shrinked2orgSC{};

   int pivotPerturbationExp; // 10^-exp
   int nIterativeRefins;
   bool parallelForwardBackward;
   bool factorizationTwoLevel;

   /* pardiso params */
   int maxfct{1};
   int mnum{1};
   int phase{0};
   int msglvl{0};
   int mtype{-2};
   int nrhs{1};

   /** dimension of the PARDISO augmented system */
   int n{-1};
   /** dimension of the Schur complement (# of rhs) */
   int nSC{-1};
   /** number of nonzeros in the PARDISO augmented matrix */
   int nnz{-1};
   /** storage for the upper triangular (in row-major format) */
   int* rowptrAug{};
   int* colidxAug{};
   double* eltsAug{};
   /** mapping from from the diagonals of the PIPS linear systems to
       the diagonal elements of the (1,1) block  in the augmented system */
   std::map<int, int> diagMap;

   //temporary vector of size n
   double* nvec{};
   double* nvec2{};
   int nvec_size{-1}; // to be save

   virtual void initPardiso() = 0;
   virtual void setIparm(int* iparm) const = 0;
   bool iparmUnchanged();

   virtual void computeSC(int nSCO,
         /*const*/SparseMatrix& R,
         /*const*/SparseMatrix& A,
         /*const*/SparseMatrix& C,
         /*const*/SparseMatrix& F,
         /*const*/SparseMatrix& G, int*& rowptrSC, int*& colidxSC, double*& eltsSC) = 0;

   ~PardisoSchurSolver() override;
};

#endif
