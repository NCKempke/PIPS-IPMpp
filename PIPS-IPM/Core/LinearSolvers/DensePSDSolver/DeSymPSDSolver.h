/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP
 *
 * C. Petra Added multiple Lsolve with multiple right hand sides      */

#ifndef DESYMPSDSOLVER_H
#define DESYMPSDSOLVER_H

#include "DoubleLinearSolver.h"
#include "Vector.hpp"
#include "DenseSymmetricMatrix.h"
#include "DenseMatrix.h"

/** A linear solver for dense, symmetric positive-definite systems.
 *  @ingroup DenseLinearAlgebra
 *  @ingroup LinearSolvers
 */
class DeSymPSDSolver : public DoubleLinearSolver {
protected:
   std::shared_ptr<DenseStorage> mStorage;
public:
   DeSymPSDSolver(const DenseSymmetricMatrix* dsm);
   void diagonalChanged(int idiag, int extent) override;
   void matrixChanged() override;

   using DoubleLinearSolver::solve;
   void solve(Vector<double>& x) override;

   //specialized method that uses BLAS-3 function DTRSM for the triagular solve.
   using DoubleLinearSolver::Lsolve;
   void Lsolve(GeneralMatrix& mat) override;

   ~DeSymPSDSolver() override = default;


   bool reports_inertia() const override { return true; };
   std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const override {
      assert(false && "TODO : implement");
      return {0, 0, 0};
   };

};

#endif
