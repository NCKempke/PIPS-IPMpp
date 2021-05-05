/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DESYMINDEFSOLVER_H
#define DESYMINDEFSOLVER_H

#include "DoubleLinearSolver.h"
#include "DenseSymMatrixHandle.h"
#include "SparseSymMatrix.h"
#include "DenseStorage.h"
#include "pipsport.h"

#include <vector>
#include <memory>

/** A linear solver for dense, symmetric indefinite systems
 * @ingroup DenseLinearAlgebra
 * @ingroup LinearSolvers
 */
class DeSymIndefSolver : public DoubleLinearSolver {
public:
   explicit DeSymIndefSolver(const DenseSymMatrix* storage);
   explicit DeSymIndefSolver(const SparseSymMatrix* storage);

   void diagonalChanged(int idiag, int extent) override;
   void matrixChanged() override;

   using DoubleLinearSolver::solve;
   void solve(Vector<double>& vec) override;
   void solve(GenMatrix& vec) override;

   ~DeSymIndefSolver() override = default;

   [[nodiscard]] bool reports_inertia() const override { return true; };
   [[nodiscard]] std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const override;

protected:

   void calculate_inertia_from_factorization() const;

   /* in PIPS symmetric matrices will be lower diagonal matrices which makes them upper diagonal in fortran access */
   const char fortranUplo = 'U';

   std::shared_ptr<DenseStorage> mStorage;
   std::vector<double> work;
   std::vector<int> ipiv;

   const SparseSymMatrix* sparseMat{};

   mutable int positive_eigenvalues{-1};
   mutable int negative_eigenvalues{-1};
   mutable int zero_eigenvalues{0};
};

#endif
