/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DOUBLELINEARSOLVER_H
#define DOUBLELINEARSOLVER_H

#include "Vector.hpp"
#include "SmartPointer.h"
#include "AbstractMatrix.h"
#include "pipsport.h"

/**
 * @defgroup LinearSolvers
 */

/**
 * Implements the main solver for linear systems that arise in
 * primal-dual interior-point methods.
 * @ingroup LinearSolvers
 * @ingroup AbstractLinearAlgebra
 */
class DoubleLinearSolver {
public:
   /** called if the diagonal elements of the matrix have
    *  changed. Triggers a refactorization of the matrix, if necessary.
    *
    *  @param idiag index of the first diagonal element that
    *               changed
    *  @param extent the number of diagonal element that changed.  */
   virtual void diagonalChanged(int idiag, int extent) = 0;

   /** called if some elements of the matrix have changed.  Triggers a
    *  refactorization of the matrix, if necessary.  */
   virtual void matrixChanged() = 0;

   /** called if new matrix (but same dimension) is to be used. Triggers factorization  */
   virtual void matrixRebuild(AbstractMatrix& /*matrixNew*/ ) { assert(0 && "Not implemented"); }

   /** solves a linear system.
    *
    * @param x on entry the right hand side of the system to be solved.
    *           On exit, the solution.  */
   virtual void solve(Vector<double>& x) = 0;

   /** does this solver report inertia of the factorized matrix back */
   [[nodiscard]] virtual bool reports_inertia() const = 0;

   /** get inertia of last factorized system */
   [[nodiscard]] virtual std::tuple<unsigned int, unsigned int, unsigned int> get_inertia() const = 0;

   /* override if necessary */
   virtual void solveSynchronized(Vector<double>& x) { solve(x); };

   // solve with multiple RHS
   virtual void solve(GeneralMatrix& /*rhs*/ ) { assert(0 && "Not implemented"); }

   // solve with multiple RHS and column sparsity array (can be nullptr)
   virtual void solve(int /*nrhss*/, double* /*rhss*/, int* /*colSparsity*/ ) { assert(0 && "Not implemented"); }

   // TODO: remove and only use solve
   void Lsolve(Vector<double>& /*x*/ ) { assert(false && "is always empty.. "); }
   virtual void Lsolve(GeneralMatrix& /*mat*/ ) { assert(0 && "Not implemented"); }
   void Dsolve(Vector<double>& x) { solve(x); }
   void Ltsolve(Vector<double>& /*x*/ ) { assert(false && "is always empty.. "); }

   /** Destructor  */
   virtual ~DoubleLinearSolver() = default;
protected:
   DoubleLinearSolver() = default;
};

class SymmetricLinearScaler {
public:
   /** compute scaling for a symmetric indefinite matrix and scale it */
   virtual void scaleMatrixTripletFormat(int n, int nnz, double* M, const int* rowM, const int* colM, bool fortran_indexed) = 0;

   virtual void scaleMatrixCSRFormat(int n, int nnz, double* M, const int* krowM, const int* jcolM, bool fortran_indexed) = 0;

   virtual const double* getScaling() const = 0;

   /* unscale a vector */
   virtual void scaleVector(Vector<double>& vec_in) const = 0;

   /* scale a vector */
   virtual void unscaleVector(Vector<double>& vec_in) const = 0;

   virtual ~SymmetricLinearScaler() = default;
};

/**
 * The abstract  interface to a mat-first operation required by
 * the iterative solvers.
 * @ingroup LinearSolvers
 * @ingroup AbstractLinearAlgebra
 */

class MatTimesVec {
public:
   /** y = beta * y + alpha * A * x */
   virtual void doIt(double beta, Vector<double>& y, double alpha, Vector<double>& x) = 0;
   virtual ~MatTimesVec() = default;
};

/**
 * The abstract interface for a linear iterative solver for linear systems
 * that arise in primal-dual interior-point methods for QP.
 * @ingroup LinearSolvers
 * @ingroup AbstractLinearAlgebra
 */

class DoubleIterativeLinearSolver : public DoubleLinearSolver {
public:
   DoubleIterativeLinearSolver(MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2 = nullptr);

   void diagonalChanged(int, int) override {};
   void matrixChanged() override {};

   ~DoubleIterativeLinearSolver() override = default;
protected:
   DoubleIterativeLinearSolver() = default;
   /** MatVec operation involving system matrix*/
   MatTimesVec* A{};

   /** MatVec ops for left and right preconditioner */
   MatTimesVec* ML{}, * MR{};

   /** Actual mat-first operations */
   void applyA(double beta, Vector<double>& res, double alpha, Vector<double>& x);
   void applyM1(double beta, Vector<double>& res, double alpha, Vector<double>& x);
   void applyM2(double beta, Vector<double>& res, double alpha, Vector<double>& x);
};

#endif
