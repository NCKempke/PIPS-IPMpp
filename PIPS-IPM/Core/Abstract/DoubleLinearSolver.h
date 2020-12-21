/* OOQP                                                               *
 * Authors: E. Michael Gertz, Stephen J. Wright                       *
 * (C) 2001 University of Chicago. See Copyright Notification in OOQP */

#ifndef DOUBLELINEARSOLVER_H
#define DOUBLELINEARSOLVER_H

#include "OoqpVectorHandle.h"
#include "DoubleMatrix.h"
#include "pipsport.h"

/**
 * @defgroup LinearSolvers
 */

/**
 * Implements the main solver for linear systems that arise in
 * primal-dual interior-point methods for QP.
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
  virtual void diagonalChanged( int idiag, int extent ) = 0;

  /** called if some elements of the matrix have changed.  Triggers a
   *  refactorization of the matrix, if necessary.  */
  virtual void matrixChanged() = 0;

  /** called if new matrix (but same dimension) is to be used. Triggers factorization  */
  virtual void matrixRebuild( DoubleMatrix& /*matrixNew*/ ) { assert(0 && "Not implemented"); }

  /** solves a linear system.
   *
   * @param x on entry the right hand side of the system to be solved.
   *           On exit, the solution.  */
   virtual void solve ( OoqpVector& x ) = 0;

   /* override if necessary */
   virtual void solveSynchronized( OoqpVector& x ) { solve(x); };

   // solve with multiple RHS
   virtual void solve ( GenMatrix& /*rhs*/ ) { assert(0 && "Not implemented"); }

   // solve with multiple RHS and column sparsity array (can be nullptr)
   virtual void solve ( int /*nrhss*/, double* /*rhss*/, int* /*colSparsity*/ ) { assert(0 && "Not implemented"); }

   // TODO: remove and only use solve
   void Lsolve( OoqpVector& /*x*/ ) { assert(false && "is always empty.. "); }
   virtual void Lsolve( GenMatrix& /*mat*/ ) { assert(0 && "Not implemented"); }
   void Dsolve( OoqpVector& x ) { solve(x);}
   void Ltsolve( OoqpVector& /*x*/ ){ assert(false && "is always empty.. "); }

  /** Destructor  */
  virtual ~DoubleLinearSolver() = default;
};

class SymmetricLinearScaler
{
   public:
      /** compute scaling for a symmetric indefinite matrix and scale it */
      virtual void scaleMatrixTripletFormat( int n, int nnz, double* M, const int* rowM, const int* colM, bool fortran_indexed ) = 0;

      virtual void scaleMatrixCSRFormat( int n, int nnz, double* M, const int* krowM, const int* jcolM, bool fortran_indexed ) = 0;

      virtual const double* getScaling() const = 0;

      /* unscale a vector */
      virtual void scaleVector( OoqpVector& vec_in ) const = 0;

      /* scale a vector */
      virtual void unscaleVector( OoqpVector& vec_in ) const = 0;

      virtual ~SymmetricLinearScaler() {};
};

/**
 * The abstract  interface to a mat-vec operation required by
 * the iterative solvers.
 * @ingroup LinearSolvers
 * @ingroup AbstractLinearAlgebra
 */

class MatTimesVec {
 public:
  /** y = beta * y + alpha * A * x */
  virtual void doIt(double beta, OoqpVector& y, double alpha, OoqpVector& x) = 0;
  virtual ~MatTimesVec();
};

/**
 * The abstract interface for a linear iterative solver for linear systems
 * that arise in primal-dual interior-point methods for QP.
 * @ingroup LinearSolvers
 * @ingroup AbstractLinearAlgebra
 */

class DoubleIterativeLinearSolver : public DoubleLinearSolver {
 public:
  DoubleIterativeLinearSolver( MatTimesVec* A, MatTimesVec* M1, MatTimesVec* M2=nullptr );

  void diagonalChanged( int, int ) override {};
  void matrixChanged() override {};

  virtual ~DoubleIterativeLinearSolver() = default;
 protected:
  DoubleIterativeLinearSolver() = default;
  /** MatVec operation involving system matrix*/
  MatTimesVec* A{};

  /** MatVec ops for left and right preconditioner */
  MatTimesVec *ML{}, *MR{};

  /** Actual mat-vec operations */
  void applyA (double beta, OoqpVector& res, double alpha, OoqpVector& x);
  void applyM1(double beta, OoqpVector& res, double alpha, OoqpVector& x);
  void applyM2(double beta, OoqpVector& res, double alpha, OoqpVector& x);
};

/**
 * An implementation of the abstract class @MatTimesVec that performs a mat-vec
 * with both the matrix and vector being on the same processor.
 *
 * It can use OOQP matrix and the implementation is based on SimpleVector class.
 */
class StoredMatTimesVec : public MatTimesVec {
 public:
  StoredMatTimesVec(DoubleMatrix* mat);
  virtual ~StoredMatTimesVec() {};

  void doIt(double beta, OoqpVector& y, double alpha, OoqpVector& x);
 protected:
  DoubleMatrix* mMat;
};
/**
 * An implementation of the abstract class @MatTimesVec that performs a
 * mat transpose-vec with both the matrix and vector being on the same processor.
 *
 * It can use OOQP matrix and the implementation is based on SimpleVector class.
 */
class StoredMatTransTimesVec : public MatTimesVec {
 public:
  StoredMatTransTimesVec(DoubleMatrix* mat);
  virtual ~StoredMatTransTimesVec() {};

  void doIt(double beta, OoqpVector& y, double alpha, OoqpVector& x);
 protected:
  DoubleMatrix* mMat;
};
#endif
