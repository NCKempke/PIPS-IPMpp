/*
 * MumpsSolverBase.h
 *
 *  Created on: 25.06.2019
 *      Author: bzfrehfe
 */

#ifndef PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERBASE_H_
#define PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERBASE_H_

#include "mpi.h"

#include "DoubleLinearSolver.h"
#include "SparseSymmetricMatrix.h"
#include "Vector.hpp"
#include "pipsport.h"
#include "dmumps_c.h"

enum MumpsVerbosity { verb_mute, verb_standard, verb_high };

#define ICNTL(I) icntl[(I)-1] // macro s.t. indices match documentation
#define INFOG(I) infog[(I)-1]
#define RINFOG(I) rinfog[(I)-1]

/** implements linear solver class that uses the MUMPS solver
 */

class MumpsSolverBase : public DoubleLinearSolver {

public:
   explicit MumpsSolverBase(const SparseSymmetricMatrix* sgm);
   MumpsSolverBase(MPI_Comm mpiCommPips_c, MPI_Comm mpiCommMumps_c, const SparseSymmetricMatrix* sgm);

   ~MumpsSolverBase() override;

   void diagonalChanged(int idiag, int extent) override;
   void matrixChanged() override = 0;

   using DoubleLinearSolver::solve;
   void solve(Vector<double>& rhs) override;

   static constexpr MumpsVerbosity defaultVerbosity = verb_mute;
//  static constexpr MumpsVerbosity defaultVerbosity = verb_standard;

   static constexpr unsigned defaultMaxNiterRefinments = 5;
   static constexpr int maxNreallocs = 5;

   bool reports_inertia() const override { return false; };
   std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const override {
      assert(false && "TODO : implement");
      return {0, 0, 0};
   };

protected:

   void setUpMpiData(MPI_Comm mpiCommPips_c, MPI_Comm mpiCommMumps_c);
   void setUpMumps();
   void processMumpsResultAnalysis(double starttime);
   void processMumpsResultFactor(double starttime);
   void processMumpsResultSolve(double starttime);

   void solve(double* vec);

   long long n{0};
   MumpsVerbosity verbosity;
   unsigned maxNiterRefinments{0};
   MPI_Comm mpiCommPips, mpiCommMumps;

   int rankMumps{-1};
   int rankPips{-1};

   DMUMPS_STRUC_C* mumps{};
   const SparseSymmetricMatrix* Msys{};

   // matrix pointer in triplet format
   int* tripletIrn{};
   int* tripletJcn{};
   double* tripletA{};
};


#endif /* PIPS_IPM_CORE_LINEARSOLVERS_MUMPSSOLVER_MUMPSSOLVERBASE_H_ */
